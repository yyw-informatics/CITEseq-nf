/*
 * Phase A: QC Foundation Subworkflow
 *
 * DAG (default, single branch):
 *   cellranger -> decontx_array -> qc_array -> scdblfinder_array -> select_branch
 *     -> [step02_dsb, step03_clr || step04_rna_norm]
 *     -> step05a_thresholdr
 *     -> [step05b_kde || step05c_fitness]
 *
 * DAG (--run_both_branches):
 *   cellranger -> [soupx_array, decontx_array] -> qc_array -> scdblfinder_array -> select_branch
 *     -> (same downstream as above)
 *
 * Fan-out pattern: every library is processed in parallel per step,
 * then aggregated before the next step begins. The library count comes
 * from --sample_metadata_json (or --n_libraries), never from a literal.
 *
 * By default only the selected branch (params.branch, default "decontx")
 * runs decontamination. Set --run_both_branches to run both SoupX and
 * DecontX for comparison — QC/scDblFinder still run on the selected
 * branch only.
 */

include { step00b_soupx_task;       step00b_soupx_aggregate       } from '../modules/step00b_soupx'
include { step00c_decontx_task;     step00c_decontx_aggregate     } from '../modules/step00c_decontx'
include { step01_qc_task;           step01_qc_aggregate           } from '../modules/step01_qc'
include { step01b_scdblfinder_task; step01b_scdblfinder_aggregate } from '../modules/step01b_scdblfinder'
include { select_branch }                                           from '../modules/select_branch'
include { step02_dsb }                                              from '../modules/step02_dsb'
include { step03_clr }                                              from '../modules/step03_clr'
include { step04_rna_norm }                                         from '../modules/step04_rna_norm'
include { step05a_thresholdr }                                      from '../modules/step05a_thresholdr'
include { step05b_kde }                                             from '../modules/step05b_kde'
include { step05c_fitness }                                         from '../modules/step05c_fitness'

/*
 * Number of libraries to fan out over. Derived from the sample metadata
 * JSON (one key per library, the same source scripts/config.py uses for
 * SAMPLE_DIRS), so the array steps match the study actually being run.
 * --n_libraries overrides it for pilot runs on a subset.
 */
def resolve_n_libraries() {
    if (params.n_libraries != null) {
        def n
        try {
            n = params.n_libraries as int
        }
        catch (Exception e) {
            error "--n_libraries must be a whole number, got: ${params.n_libraries}"
        }
        if (n < 1) {
            error "--n_libraries must be 1 or greater, got: ${params.n_libraries}"
        }
        return n
    }
    if (!params.sample_metadata_json) {
        error """
        Cannot determine how many libraries to process.
          Pass --sample_metadata_json /path/to/sample_metadata.json
          (a JSON object with one key per library, e.g. {"Lib01": {...}, ...}),
          or pass --n_libraries N to set the count directly.
        """.stripIndent()
    }
    def meta_file = file(params.sample_metadata_json)
    if (!meta_file.exists()) {
        error "--sample_metadata_json not found: ${params.sample_metadata_json}"
    }
    def meta = new groovy.json.JsonSlurper().parse(meta_file)
    if (!(meta instanceof Map) || meta.isEmpty()) {
        error "--sample_metadata_json must be a non-empty JSON object keyed by library id: ${params.sample_metadata_json}"
    }
    return meta.size()
}

workflow PHASE_A_QC {
    take:
    branch          // val: which decontamination branch to use ("decontx" or "soupx")

    main:
    def selected = params.branch ?: 'decontx'

    // Task indices are 0-based and shared by every array step below.
    def n_libraries = resolve_n_libraries()
    def library_ids = (0..<n_libraries).toList()
    log.info "Phase A: fanning out over ${n_libraries} librar${n_libraries == 1 ? 'y' : 'ies'}"

    // ── Decontamination: run selected branch (or both if requested) ─
    if (selected == 'soupx' || params.run_both_branches) {
        soupx_ids = Channel.fromList(library_ids)
        step00b_soupx_task(soupx_ids)
        step00b_soupx_aggregate(
            step00b_soupx_task.out.done.collect()
        )
    }

    if (selected == 'decontx' || params.run_both_branches) {
        decontx_ids = Channel.fromList(library_ids)
        step00c_decontx_task(decontx_ids)
        step00c_decontx_aggregate(
            step00c_decontx_task.out.done.collect()
        )
    }

    // ── Gate: wait for decontamination to complete ──────────────────
    if (params.run_both_branches) {
        decontam_done = step00b_soupx_aggregate.out.done
            .combine(step00c_decontx_aggregate.out.done)
    } else if (selected == 'soupx') {
        decontam_done = step00b_soupx_aggregate.out.done
    } else {
        decontam_done = step00c_decontx_aggregate.out.done
    }

    // ── Step 01: QC + Scrublet (one task per library, gated) ────────
    // flatMap emits the task IDs only after decontam_done fires
    qc_task_ids = decontam_done.flatMap { library_ids }
    step01_qc_task(qc_task_ids, branch)
    step01_qc_aggregate(
        step01_qc_task.out.done.collect(),
        branch
    )

    // ── Step 01b: scDblFinder (one task per library, gated) ─────────
    scdblfinder_task_ids = step01_qc_aggregate.out.done
        .flatMap { library_ids }
    step01b_scdblfinder_task(scdblfinder_task_ids, branch)
    step01b_scdblfinder_aggregate(
        step01b_scdblfinder_task.out.done.collect(),
        branch
    )

    // ── Branch selection: promote to canonical paths ────────────────
    select_branch(branch, step01b_scdblfinder_aggregate.out.done)

    // ── Path channels from canonical files on real disk ─────────────
    // select_branch.out.done is a process output — auto-forks in DSL2.
    // Each .map() derivation creates an independent gated path channel.
    def pd = "${params.project_root}/data/processed"

    // ── Steps 02, 03, 04: Parallel normalization block ──────────────
    // Step 02 (DSB) reads CellRanger h5 from disk — just needs gate
    step02_dsb(select_branch.out.done)

    // Step 03 (CLR) stages adata_qc + adata_adt_raw
    step03_clr(
        select_branch.out.done.map { file("${pd}/adata_qc.h5ad") },
        select_branch.out.done.map { file("${pd}/adata_adt_raw.h5ad") }
    )

    // Step 04 (RNA norm) stages adata_qc — parallel with step03
    step04_rna_norm(
        select_branch.out.done.map { file("${pd}/adata_qc.h5ad") }
    )

    // ── Step 05a: ThresholdR gating (after norm block) ──────────────
    // Collect all norm outputs as a gate — step05a starts only after
    // steps 02, 03, and 04 have all completed.
    norm_gate = step02_dsb.out.dsb
        .mix(step03_clr.out.clr, step04_rna_norm.out.normalized)
        .collect()

    step05a_thresholdr(
        norm_gate.map { file("${pd}/adt_counts_for_thresholdr.csv") }
    )

    // ── Steps 05b, 05c: Parallel diagnostics (after step05a) ───────
    // step05a.out.thresholds is a process output — auto-forks for
    // both step05b and step05c. The thresholds dependency transitively
    // gates these on the full norm block.

    // step05b needs CLR-normalized data (ADT_CD7 etc. in .obs).
    // step03_clr.out.clr is adata_qc_clr.h5ad; stage it as adata_qc.h5ad
    // so the script reads it via ADATA_QC = $PWD/adata_qc.h5ad.
    step05b_kde(
        step03_clr.out.clr,
        step05a_thresholdr.out.thresholds,
        step05a_thresholdr.out.thresholds.map { file("${pd}/adt_counts_for_thresholdr.csv") }
    )

    // step05c also needs CLR columns for its RNA concordance checks.
    step05c_fitness(
        step03_clr.out.clr,
        step04_rna_norm.out.normalized,
        step05a_thresholdr.out.thresholds
    )

    emit:
    clr        = step03_clr.out.clr
    normalized = step04_rna_norm.out.normalized
    hvg        = step04_rna_norm.out.hvg
    dsb        = step02_dsb.out.dsb
    thresholds = step05a_thresholdr.out.thresholds
    done       = step05b_kde.out.valleys
                     .mix(step05c_fitness.out.metrics)
                     .collect()
}

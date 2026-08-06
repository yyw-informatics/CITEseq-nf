#!/usr/bin/env nextflow

/*
 * citeseq-nf — Generic CITE-seq preprocessing + embedding pipeline
 *
 * This Nextflow pipeline covers the generic, method-driven portion of a
 * CITE-seq analysis: ambient-RNA correction, QC + doublet detection,
 * ADT + RNA normalization, marker gating, two-embedding integration
 * (totalVI + scVI), and embedding comparison. Downstream classification,
 * hypothesis testing, flow validation, and manuscript figures remain as
 * plain scripts in `scripts/python/` and `scripts/R/` — they are tightly
 * coupled to this study's biology and run rarely.
 *
 * Generic pipeline DAG (default workflow):
 *   cellranger → decontx/soupx → QC → scDblFinder → select_branch
 *     → [DSB ‖ CLR ‖ RNA_norm+HVG]
 *     → ThresholdR → [KDE ‖ fitness]
 *     → [totalVI ‖ scVI] → embedding_comparison → adata_qc_embedded
 *
 * Entry points:
 *   (default)          Full generic pipeline (Phase A QC + Phase B embed)
 *   from_cellranger    Alias for default — same as (default)
 *   from_qc            Skip decontam/QC/doublet — start from adata_qc.h5ad
 *   from_clr           Checkpoint staging for adata_qc_clr.h5ad
 *   clr_only           Phase 0 POC: step03_clr only
 *   from_hvg           Start from adata_qc_hvg.h5ad → totalVI + scVI + embed
 *   from_totalvi_scvi  Start from existing totalVI + scVI outputs → embed only
 */

nextflow.enable.dsl = 2

// ── Required-param validation ─────────────────────────────────────────
// Applied at script-parse time so every entry point and -preview hits it.
if (!params.project_root) {
    error """
    Required param missing: --project_root
      Point this at the calling study's project root. citeseq-nf exports
      CITESEQ_PROJECT_ROOT to every process, and scripts/config.py reads it
      to resolve paths.
    Example:
      ./run.sh run ./main.nf --project_root /path/to/your/study --outdir /path/to/your/study
    """.stripIndent().trim()
}
if (!params.outdir) {
    error """
    Required param missing: --outdir
      Where citeseq-nf should publish outputs (data/processed, results,
      models, figures, logs). Usually equal to --project_root.
    """.stripIndent().trim()
}
if (!params.conda_env) {
    error """
    Required param missing: --conda_env
      Path to the CITE-seq conda/mamba environment. Used by withLabel
      env_cite_seq and as the fallback for --python and --rscript
      (\${conda_env}/bin/python, \${conda_env}/bin/Rscript).
    Example:
      --conda_env /path/to/your/envs/cite_seq
    """.stripIndent().trim()
}
if (!params.panel) {
    error """
    Required param missing: --panel
      Path to an ADT panel YAML defining your antibody panel
      (cellranger_names, rename_map, functional, thresholdr,
      thresholdr_reference, plus optional fgprob_override / kde_validate /
      rna_concordance). See panels/example_panel.yaml for the schema.
    Example:
      --panel ${workflow.projectDir}/panels/example_panel.yaml
    """.stripIndent().trim()
}
if (!params.data_root) {
    error """
    Required param missing: --data_root
      Directory holding the CellRanger output, one subdirectory per library.
      Each library directory must contain raw_feature_bc_matrix.h5 and
      filtered_feature_bc_matrix.h5, plus graph-based clusters either as
      analysis/clustering/ or as analysis.tar.gz.
    Example:
      --data_root /path/to/your/study/data/cellranger
    """.stripIndent().trim()
}

include { PHASE_A_QC }    from './subworkflows/phase_a_qc'
include { PHASE_B_EMBED } from './subworkflows/phase_b_embed'
include { step02_dsb }         from './modules/step02_dsb'
include { step03_clr }         from './modules/step03_clr'
include { step04_rna_norm }    from './modules/step04_rna_norm'
include { step05a_thresholdr } from './modules/step05a_thresholdr'
include { step05b_kde }        from './modules/step05b_kde'
include { step05c_fitness }    from './modules/step05c_fitness'
include { step08_embedding }   from './modules/step08_embedding'

// ── Default workflow: full generic pipeline (Phase A QC → Phase B embed) ──
workflow {
    branch_ch = Channel.value(params.branch ?: 'decontx')
    PHASE_A_QC(branch_ch)
    PHASE_B_EMBED(PHASE_A_QC.out.hvg, PHASE_A_QC.out.clr)
}

// ── Named alias for the default: full generic pipeline from CellRanger ────
workflow from_cellranger {
    branch_ch = Channel.value(params.branch ?: 'decontx')
    PHASE_A_QC(branch_ch)
    PHASE_B_EMBED(PHASE_A_QC.out.hvg, PHASE_A_QC.out.clr)
}

// ── Checkpoint entry: skip decontam/QC/doublet, start from QC'd data ──
workflow from_qc {
    if (!params.adata_qc || !params.adata_adt_raw || !params.adt_counts_for_thresholdr) {
        error """
        Required params missing. Usage:
          nextflow run main.nf -entry from_qc \\
            --adata_qc /path/to/adata_qc.h5ad \\
            --adata_adt_raw /path/to/adata_adt_raw.h5ad \\
            --adt_counts_for_thresholdr /path/to/adt_counts_for_thresholdr.csv
        """.stripIndent().trim()
    }

    def pd = "${params.project_root}/data/processed"

    adata_qc_ch     = Channel.fromPath(params.adata_qc, checkIfExists: true)
    adata_adt_raw_ch = Channel.fromPath(params.adata_adt_raw, checkIfExists: true)
    adt_counts_ch   = Channel.fromPath(params.adt_counts_for_thresholdr, checkIfExists: true)

    // Parallel normalization block (steps 02, 03, 04)
    step02_dsb(adt_counts_ch.map { true })     // gate signal — reads CellRanger from disk
    step03_clr(adata_qc_ch, adata_adt_raw_ch)
    step04_rna_norm(adata_qc_ch)

    // Gate step05a after norm block
    norm_gate = step02_dsb.out.dsb
        .mix(step03_clr.out.clr, step04_rna_norm.out.normalized)
        .collect()

    step05a_thresholdr(
        norm_gate.map { file(params.adt_counts_for_thresholdr) }
    )

    // Parallel diagnostics (steps 05b, 05c)
    // Both need CLR-normalized data (ADT_<marker> columns in .obs)
    step05b_kde(
        step03_clr.out.clr,
        step05a_thresholdr.out.thresholds,
        step05a_thresholdr.out.thresholds.map { file(params.adt_counts_for_thresholdr) }
    )

    step05c_fitness(
        step03_clr.out.clr,
        step04_rna_norm.out.normalized,
        step05a_thresholdr.out.thresholds
    )

    // Chain into Phase B embedding
    PHASE_B_EMBED(step04_rna_norm.out.hvg, step03_clr.out.clr)
}

// ── Checkpoint entry: from HVG selection → totalVI + scVI + embedding ─
workflow from_hvg {
    def pd = "${params.project_root}/data/processed"
    def hvg = params.adata_qc_hvg ?: "${pd}/adata_qc_hvg.h5ad"
    def clr = params.adata_qc_clr ?: "${pd}/adata_qc_clr.h5ad"

    adata_qc_hvg_ch = Channel.fromPath(hvg, checkIfExists: true)
    adata_qc_clr_ch = Channel.fromPath(clr, checkIfExists: true)

    // Run steps 06 + 07 + 08 as a single subworkflow
    PHASE_B_EMBED(adata_qc_hvg_ch, adata_qc_clr_ch)
}

// ── Checkpoint entry: resume from trained embeddings, run only step 08 ─
workflow from_totalvi_scvi {
    def pd = "${params.project_root}/data/processed"
    def clr      = params.adata_qc_clr      ?: "${pd}/adata_qc_clr.h5ad"
    def totalvi  = params.adata_totalvi     ?: "${pd}/adata_totalvi.h5ad"
    def scvi     = params.adata_scvi        ?: "${pd}/adata_scvi.h5ad"

    adata_qc_clr_ch     = Channel.fromPath(clr,     checkIfExists: true)
    adata_totalvi_ch = Channel.fromPath(totalvi, checkIfExists: true)
    adata_scvi_ch    = Channel.fromPath(scvi,    checkIfExists: true)

    step08_embedding(adata_qc_clr_ch, adata_totalvi_ch, adata_scvi_ch)
}

// ── Checkpoint entry: CLR normalization only ────────────────────────
workflow clr_only {
    if (!params.adata_qc || !params.adata_adt_raw) {
        error """
        Required params missing. Usage:
          nextflow run main.nf -entry clr_only \\
            --adata_qc /path/to/adata_qc.h5ad \\
            --adata_adt_raw /path/to/adata_adt_raw.h5ad
        """.stripIndent().trim()
    }

    adata_qc_ch      = Channel.fromPath(params.adata_qc, checkIfExists: true)
    adata_adt_raw_ch  = Channel.fromPath(params.adata_adt_raw, checkIfExists: true)

    step03_clr(adata_qc_ch, adata_adt_raw_ch)

    step03_clr.out.clr
        .subscribe { file ->
            file.copyTo("${params.outdir}/data/processed/${file.name}")
        }
}

// ── Checkpoint entry: resume from CLR ───────────────────────────────
workflow from_clr {
    if (!params.adata_qc_clr) {
        error """
        Required param missing. Usage:
          nextflow run main.nf -entry from_clr \\
            --adata_qc_clr /path/to/adata_qc_clr.h5ad
        """.stripIndent().trim()
    }

    Channel.fromPath(params.adata_qc_clr, checkIfExists: true)
        .set { clr_ch }

    clr_ch.view { "Checkpoint entry: staged ${it.name} from ${it}" }
}

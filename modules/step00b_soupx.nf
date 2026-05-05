/*
 * Step 00b: SoupX Ambient RNA Decontamination
 *
 * Contract:
 *   Input:  CellRanger outputs (read from DATA_ROOT on disk)
 *   Output: soupx_parallel/<sample>/filtered_feature_bc_matrix.h5 (x16)
 *           results/phase1_qc/soupx/soupx_diag_<sample>.csv (x16)
 *   Aggregate: results/phase1_qc/soupx/soupx_summary_parallel.csv + barplot
 *
 * Array pattern: 16 per-library tasks -> 1 aggregate.
 * Scripts read from DATA_ROOT and write to PROCESSED_DIR/RESULTS_DIR
 * on the real filesystem (not work-dir staged), because CellRanger
 * inputs are too large/complex to stage via channels.
 *
 * Channel-pure: no publishDir. Publishing is done at workflow level.
 */

process step00b_soupx_task {
    label 'cpu_large'
    label 'env_cite_seq'

    input:
    val(task_id)

    output:
    val(task_id), emit: done

    script:
    // Script CLI uses 1-based --sample-index; channel emits 0-based task_id
    def sample_index = task_id + 1
    """
    ${params.python} ${workflow.projectDir}/scripts/python/step00b_soupx_decontamination.py \
        --sample-index ${sample_index}
    """
}

process step00b_soupx_aggregate {
    label 'cpu_small'
    label 'env_cite_seq'

    input:
    val(ready)

    output:
    val(true), emit: done

    script:
    """
    ${params.python} ${workflow.projectDir}/scripts/python/step00b_soupx_decontamination.py \
        --aggregate
    """
}

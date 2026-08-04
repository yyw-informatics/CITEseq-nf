/*
 * Step 00c: DecontX Ambient RNA Decontamination
 *
 * Contract:
 *   Input:  CellRanger outputs (read from DATA_ROOT on disk)
 *   Output: decontx/<sample>/filtered_feature_bc_matrix.h5 (one per library)
 *           step00c_temp/decontx_<sample>.csv (one per library)
 *   Aggregate: results/phase1_qc/decontx/decontx_summary.csv + barplot
 *
 * Array pattern: one task per library -> 1 aggregate.
 * Scripts read from DATA_ROOT and write to PROCESSED_DIR/RESULTS_DIR
 * on the real filesystem.
 *
 * Channel-pure: no publishDir. Publishing is done at workflow level.
 */

process step00c_decontx_task {
    label 'cpu_medium'
    label 'env_cite_seq'

    input:
    val(task_id)

    output:
    val(task_id), emit: done

    script:
    """
    ${params.rscript} ${workflow.projectDir}/scripts/R/step00c_decontx.R \
        --task ${task_id}
    """
}

process step00c_decontx_aggregate {
    label 'cpu_small'
    label 'env_cite_seq'

    input:
    val(ready)

    output:
    val(true), emit: done

    script:
    """
    ${params.rscript} ${workflow.projectDir}/scripts/R/step00c_decontx.R \
        --aggregate
    """
}

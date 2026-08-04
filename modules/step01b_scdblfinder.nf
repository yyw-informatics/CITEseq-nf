/*
 * Step 01b: scDblFinder Doublet Detection
 *
 * Contract:
 *   Input:  Decontaminated h5 + QC keep_barcodes from step 01
 *   Output per-library (--task):
 *     step01b_temp/<sample>_doublets_{branch}.csv
 *   Aggregate output:
 *     Updates adata_qc_{branch}.h5ad with scDblFinder annotations
 *     results/phase1_qc/doublet_comparison_{branch}.csv
 *
 * Array pattern: one task per library -> 1 aggregate.
 * The branch parameter selects soupx or decontx inputs.
 *
 * Channel-pure: no publishDir. Publishing is done at workflow level.
 */

process step01b_scdblfinder_task {
    label 'cpu_small'
    label 'env_cite_seq'

    input:
    val(task_id)
    val(branch)

    output:
    val(task_id), emit: done

    script:
    def branch_arg = branch ? "--branch ${branch}" : ''
    """
    ${params.rscript} ${workflow.projectDir}/scripts/R/step01b_scdblfinder.R \
        --task ${task_id} ${branch_arg}
    """
}

process step01b_scdblfinder_aggregate {
    label 'cpu_large'
    label 'env_cite_seq'

    input:
    val(ready)
    val(branch)

    output:
    val(true), emit: done

    script:
    def branch_arg = branch ? "--branch ${branch}" : ''
    """
    ${params.rscript} ${workflow.projectDir}/scripts/R/step01b_scdblfinder.R \
        --aggregate ${branch_arg}
    """
}

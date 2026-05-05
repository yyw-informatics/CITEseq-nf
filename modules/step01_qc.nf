/*
 * Step 01: Load Libraries, QC Filtering, and Scrublet
 *
 * Contract:
 *   Input:  Decontaminated h5 files from step 00b (soupx) or 00c (decontx)
 *   Output per-library (--task):
 *     step01_temp/<sample>_{rna,adt}_{branch}.h5ad
 *     step01_temp/<sample>_qc_{branch}.csv
 *     step01_temp/<sample>_keep_barcodes_{branch}.csv
 *   Aggregate output:
 *     adata_qc_{branch}.h5ad, adata_adt_raw_{branch}.h5ad
 *     adt_counts_for_thresholdr_{branch}.csv, qc_summary_{branch}.csv
 *
 * Array pattern: 16 per-library tasks -> 1 aggregate.
 * The branch parameter selects which decontamination input to use.
 *
 * Channel-pure: no publishDir. Publishing is done at workflow level.
 */

process step01_qc_task {
    label 'cpu_medium'
    label 'env_cite_seq'

    input:
    val(task_id)
    val(branch)

    output:
    val(task_id), emit: done

    script:
    def branch_arg = branch ? "--branch ${branch}" : ''
    """
    ${params.python} ${workflow.projectDir}/scripts/python/step01_qc.py \
        --task ${task_id} ${branch_arg}
    """
}

process step01_qc_aggregate {
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
    ${params.python} ${workflow.projectDir}/scripts/python/step01_qc.py \
        --aggregate ${branch_arg}
    """
}

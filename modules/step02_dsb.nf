/*
 * Step 02: DSB ADT Normalization
 *
 * Contract:
 *   Input:  CellRanger h5 files (read from fixed disk paths via config.R)
 *   Output: dsb_normalized_adt.csv
 *   Applies dsb Step I ambient correction to 4 functional ADT markers.
 *
 * Note: This step reads CellRanger filtered/raw matrices directly from
 * DATA_ROOT (10X_cloud/outs/), not from prior step outputs. The val(ready)
 * gate ensures upstream steps (select_branch) have completed.
 *
 * publishDir copies output to data/processed/ for downstream steps.
 */

process step02_dsb {
    label 'cpu_small'
    label 'env_cite_seq'

    publishDir "${params.outdir}/data/processed", mode: 'copy', overwrite: true

    input:
    val(ready)

    output:
    path('dsb_normalized_adt.csv', emit: dsb)

    script:
    """
    export SAGE_PROCESSED_DIR=\$PWD
    ${params.rscript} ${workflow.projectDir}/scripts/R/step02_dsb_normalization.R
    """
}

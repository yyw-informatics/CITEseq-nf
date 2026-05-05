/*
 * Step 04: RNA Normalization + HVG Selection
 *
 * Contract:
 *   Input:  adata_qc.h5ad (step01 aggregate)
 *   Output: adata_normalized.h5ad (log-normalized RNA)
 *           adata_qc_hvg.h5ad (raw counts with HVG flags for totalVI/scVI)
 *   Adds:   log-normalized .X, highly_variable flags, forced literature markers
 *
 * Independent of step03 (CLR) — both read adata_qc.h5ad and can run in parallel.
 *
 * publishDir copies output to data/processed/ for downstream steps.
 */

process step04_rna_norm {
    label 'cpu_large'
    label 'env_cite_seq'

    publishDir "${params.outdir}/data/processed", mode: 'copy', overwrite: true

    input:
    path(adata_qc, name: 'adata_qc.h5ad')

    output:
    path('adata_normalized.h5ad', emit: normalized)
    path('adata_qc_hvg.h5ad',    emit: hvg)

    script:
    """
    export SAGE_PROCESSED_DIR=\$PWD
    ${params.python} ${workflow.projectDir}/scripts/python/step04_rna_normalization.py
    """
}

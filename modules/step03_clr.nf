/*
 * Step 03: CLR Normalization and Protein Expression Matrix
 *
 * Contract:
 *   Input:  adata_qc.h5ad (step01), adata_adt_raw.h5ad (step01)
 *   Output: adata_qc_clr.h5ad
 *   Adds:   obsm['protein_expression'], uns['protein_names'],
 *           per-marker ADT_<marker> CLR columns in .obs (one per
 *           panel functional marker)
 *
 * publishDir copies output to data/processed/ for downstream steps.
 */

process step03_clr {
    label 'cpu_medium'
    label 'env_cite_seq'

    publishDir "${params.outdir}/data/processed", mode: 'copy', overwrite: true

    input:
    path(adata_qc,      name: 'adata_qc.h5ad')
    path(adata_adt_raw,  name: 'adata_adt_raw.h5ad')

    output:
    path('adata_qc_clr.h5ad', emit: clr)

    script:
    """
    export CITESEQ_PROCESSED_DIR=\$PWD
    ${params.python} ${workflow.projectDir}/scripts/python/step03_clr_protein.py
    """
}

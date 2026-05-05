/*
 * Step 07: scVI RNA-Only Model
 *
 * Contract:
 *   Input:  adata_qc_hvg.h5ad (step04 — HVG-flagged raw counts)
 *   Output: adata_scvi.h5ad, models/scvi_model/ directory, UMAP PNG
 *   Adds:   obsm['X_scVI'], leiden_scvi, UMAP
 *
 * Stochastic: outputs vary between runs. Validation uses invariants
 * (ELBO ±2%, cell count ±0.5%, latent dim match), not exact diff.
 *
 * Lighter than totalVI — no protein branch, smaller architecture.
 * CPU is fine (~20-40 min with early stopping).
 */

process step07_scvi {
    label 'gpu_scvi'
    label 'env_cite_seq'

    publishDir "${params.outdir}/data/processed", mode: 'copy', overwrite: true,
        pattern: '*.h5ad'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true,
        pattern: 'models/**'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true,
        pattern: 'figures/**'

    input:
    path(adata_qc_hvg, name: 'adata_qc_hvg.h5ad')

    output:
    path('adata_scvi.h5ad',                emit: scvi)
    path('models/scvi_model',               emit: model)
    path('figures/phase2_scvi_*.png',       emit: figures)

    script:
    """
    mkdir -p models results/phase2_embedding figures
    export SAGE_PROCESSED_DIR=\$PWD
    export SAGE_MODELS_DIR=\$PWD/models
    export SAGE_RESULTS_DIR=\$PWD/results
    export SAGE_FIGURES_DIR=\$PWD/figures
    ${params.python} ${workflow.projectDir}/scripts/python/step07_scvi.py
    """
}

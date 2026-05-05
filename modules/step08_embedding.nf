/*
 * Step 08: Embedding Comparison and Selection
 *
 * Contract:
 *   Input:   adata_qc_clr.h5ad       (step03 — CLR normalized, protein_expression)
 *            adata_totalvi.h5ad   (step06 — X_totalVI + denoised ADT + fg_prob)
 *            adata_scvi.h5ad      (step07 — X_scVI RNA-only latent)
 *   Output:  adata_qc_embedded.h5ad  (merged embeddings + neighbors + UMAP)
 *            embedding_comparison.csv (silhouette scores + primary embedding decision)
 *   Adds:    obsm['X_totalVI'], obsm['X_scVI'], obsm['X_umap'], obsm['X_umap_totalvi'],
 *            obs['ADT_*_denoised/_fg_prob'], obs['leiden_totalvi'],
 *            uns['primary_embedding']
 *
 * The script reads the three files itself and merges internally — no pre-run
 * graft needed. The Nextflow wrapper only has to stage all three with their
 * canonical names (step06 pattern: declared inputs, not merged upstream).
 *
 * Stochastic: silhouette scores vary between runs (subsample + stochastic
 * embeddings from steps 06/07). Validation uses invariants (±2% silhouette,
 * ±0.5% cell count, obsm key parity) per §19.7.
 *
 * Memory: CLR h5ad (~7.5 GB) drives peak — totalVI + scVI are small latent
 * files. cpu_large (32 GB) gives headroom for neighbors + UMAP on 167K cells.
 */

process step08_embedding {
    label 'cpu_large'
    label 'env_cite_seq'

    publishDir "${params.outdir}/data/processed", mode: 'copy', overwrite: true,
        pattern: '*.h5ad'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true,
        pattern: 'results/**'

    input:
    path(adata_qc_clr,     name: 'adata_qc_clr.h5ad')
    path(adata_totalvi, name: 'adata_totalvi.h5ad')
    path(adata_scvi,    name: 'adata_scvi.h5ad')

    output:
    path('adata_qc_embedded.h5ad',                         emit: embedded)
    path('results/phase2_embedding/embedding_comparison.csv', emit: comparison)

    script:
    // Scoped here (not in the global beforeScript) so that toggling or
    // sweeping --biological_separation_marker only invalidates step08's
    // cache, leaving Phase A / Phase B (steps 00-07) -resume-able.
    def bio_marker_export = params.biological_separation_marker ? "export SAGE_BIOLOGICAL_SEPARATION_MARKER=${params.biological_separation_marker}" : ""
    """
    mkdir -p results/phase2_embedding figures models
    export SAGE_PROCESSED_DIR=\$PWD
    export SAGE_MODELS_DIR=\$PWD/models
    export SAGE_RESULTS_DIR=\$PWD/results
    export SAGE_FIGURES_DIR=\$PWD/figures
    ${bio_marker_export}
    ${params.python} ${workflow.projectDir}/scripts/python/step08_embedding_comparison.py
    """
}

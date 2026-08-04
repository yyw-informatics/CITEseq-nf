/*
 * Step 06: totalVI Joint RNA+ADT Model
 *
 * Contract:
 *   Input:  adata_qc_hvg.h5ad (step04 — HVG-flagged raw counts)
 *           adata_qc_clr.h5ad (step03 — provides protein_expression in obsm)
 *   Output: adata_totalvi.h5ad, models/totalvi_model/ directory,
 *           training history CSVs, UMAP PNGs
 *   Adds:   obsm['X_totalVI'], obsm['protein_denoised'], obsm['protein_fg_prob'],
 *           per-marker _denoised and _fg_prob obs columns, leiden_totalvi, UMAP
 *
 * The script expects protein_expression already in the HVG adata. In the
 * original pipeline this happened via in-place writes; after the refactor,
 * protein_expression lives only in adata_qc_clr.h5ad. The merge step below
 * grafts it onto the HVG file before running the script.
 *
 * Stochastic: outputs vary between runs. Validation uses invariants
 * (ELBO ±2%, cell count ±0.5%, latent dim match), not exact diff.
 *
 * GPU recommended (~35 s/epoch) but works on CPU (~51 s/epoch).
 */

process step06_totalvi {
    label 'gpu_totalvi'
    label 'env_cite_seq'

    publishDir "${params.outdir}/data/processed", mode: 'copy', overwrite: true,
        pattern: '*.h5ad'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true,
        pattern: 'models/**'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true,
        pattern: 'results/**'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true,
        pattern: 'figures/**'

    input:
    path(adata_qc_hvg, name: 'adata_qc_hvg.h5ad')
    path(adata_qc_clr, name: 'adata_qc_clr.h5ad')

    output:
    path('adata_totalvi.h5ad',                            emit: totalvi)
    path('models/totalvi_model',                           emit: model)
    path('results/phase2_embedding/totalvi_*.csv',         emit: history)
    path('figures/phase2_totalvi_*.png',                   emit: figures)

    script:
    """
    # Merge protein_expression from CLR into HVG file.
    # After the in-place write refactor, these live in separate h5ad files.
    ${params.python} -c "
import anndata as ad, os, gc
hvg = ad.read_h5ad('adata_qc_hvg.h5ad')
clr = ad.read_h5ad('adata_qc_clr.h5ad', backed='r')
hvg.obsm['protein_expression'] = clr.obsm['protein_expression'][:]
clr.file.close()
del clr; gc.collect()
os.unlink('adata_qc_hvg.h5ad')
hvg.write('adata_qc_hvg.h5ad')
del hvg; gc.collect()
print('Merged protein_expression from CLR into HVG')
"

    mkdir -p models results/phase2_embedding figures
    export CITESEQ_PROCESSED_DIR=\$PWD
    export CITESEQ_MODELS_DIR=\$PWD/models
    export CITESEQ_RESULTS_DIR=\$PWD/results
    export CITESEQ_FIGURES_DIR=\$PWD/figures
    ${params.python} ${workflow.projectDir}/scripts/python/step06_totalvi.py
    """
}

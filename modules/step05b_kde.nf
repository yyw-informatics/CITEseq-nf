/*
 * Step 05b: KDE Gating Validation
 *
 * Contract:
 *   Input:  adata_qc_clr.h5ad (step03, staged as adata_qc.h5ad)
 *           thresholdr_thresholds.csv (step05a)
 *           adt_counts_for_thresholdr.csv (step01, raw ADT counts for
 *           the panel's KDE_MARKERS)
 *   Output: kde_valleys.csv, kde_stain_quality.csv (PROCESSED_DIR)
 *           kde_vs_thresholdr_comparison.csv (RESULTS_DIR)
 *           diagnostic PNGs (FIGURES_DIR)
 *
 * Diagnostic step — validates ThresholdR gating via independent KDE method.
 * Can run in parallel with step05c (fitness).
 *
 * Channel-pure: no publishDir. Publishing is done at workflow level.
 */

process step05b_kde {
    label 'cpu_small'
    label 'env_cite_seq'

    input:
    path(adata_qc,   name: 'adata_qc.h5ad')
    path(thresholds, name: 'thresholdr_thresholds.csv')
    path(adt_counts, name: 'adt_counts_for_thresholdr.csv')

    output:
    path('kde_valleys.csv',       emit: valleys)
    path('kde_stain_quality.csv', emit: stain_quality)
    path('results/phase1_qc/kde_vs_thresholdr_comparison.csv', emit: comparison)
    path('figures/phase1_kde_*.png',                           emit: figures)
    path('figures/phase1_kde_vs_thresholdr_grid.png',          emit: grid)

    script:
    """
    export CITESEQ_PROCESSED_DIR=\$PWD
    export CITESEQ_FIGURES_DIR=\$PWD/figures
    export CITESEQ_RESULTS_DIR=\$PWD/results
    mkdir -p figures results/phase1_qc
    ${params.python} ${workflow.projectDir}/scripts/python/step05b_kde_gating_validation.py
    """
}

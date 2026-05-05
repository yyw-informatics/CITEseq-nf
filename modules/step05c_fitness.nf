/*
 * Step 05c: ThresholdR Fitness Assessment
 *
 * Contract:
 *   Input:  adata_qc_clr.h5ad (step03, staged as adata_qc.h5ad — .obs only)
 *           adata_normalized.h5ad (step04, backed mode — subsetted to concordance genes)
 *           thresholdr_thresholds.csv (step05a)
 *   Output: thresholdr_fitness_metrics.csv (per-marker fitness verdicts)
 *           fitness summary PNG
 *
 * Diagnostic step — assesses ThresholdR gating quality via overlap, BIC,
 * stain quality, and RNA concordance metrics.
 * Can run in parallel with step05b (KDE).
 *
 * Channel-pure: no publishDir. Publishing is done at workflow level.
 */

process step05c_fitness {
    label 'cpu_medium'
    label 'env_cite_seq'

    input:
    path(adata_qc,   name: 'adata_qc.h5ad')
    path(adata_norm, name: 'adata_normalized.h5ad')
    path(thresholds, name: 'thresholdr_thresholds.csv')

    output:
    path('results/phase1_qc/thresholdr_fitness_metrics.csv',  emit: metrics)
    path('figures/phase1_thresholdr_fitness_summary.png',      emit: figure)

    script:
    """
    export SAGE_PROCESSED_DIR=\$PWD
    export SAGE_FIGURES_DIR=\$PWD/figures
    export SAGE_RESULTS_DIR=\$PWD/results
    mkdir -p figures results/phase1_qc
    ${params.python} ${workflow.projectDir}/scripts/python/step05c_thresholdr_fitness.py
    """
}

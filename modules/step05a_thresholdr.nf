/*
 * Step 05a: ThresholdR GMM Gating
 *
 * Contract:
 *   Input:  adt_counts_for_thresholdr.csv (raw ADT counts from step01)
 *   Output: thresholdr_thresholds.csv (per-marker GMM thresholds)
 *   Also writes diagnostic PNGs to FIGURES_DIR and PDFs to RESULTS_DIR
 *   (both on real disk — not redirected to work dir).
 *
 * Reads only raw ADT counts; does NOT depend on DSB (step02) output.
 * Gated after the parallel normalization block for DAG ordering.
 *
 * publishDir copies output to data/processed/ for downstream steps.
 */

process step05a_thresholdr {
    label 'cpu_small'
    label 'env_cite_seq'

    publishDir "${params.outdir}/data/processed", mode: 'copy', overwrite: true

    input:
    path(adt_counts, name: 'adt_counts_for_thresholdr.csv')

    output:
    path('thresholdr_thresholds.csv', emit: thresholds)

    script:
    """
    export CITESEQ_PROCESSED_DIR=\$PWD
    ${params.rscript} ${workflow.projectDir}/scripts/R/step05a_thresholdr_gating.R
    """
}

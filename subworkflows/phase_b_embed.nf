/*
 * Phase B: Embedding Subworkflow
 *
 * DAG:
 *   adata_qc_hvg ─┐
 *                 ├─→ step06_totalvi ─┐
 *   adata_qc_clr ─┤                   │
 *                 ├───────────────────┼─→ step08_embedding → adata_qc_embedded
 *                 ├─→ step07_scvi   ──┘
 *                 │
 *                 └─ (CLR is consumed by both step06 merge + step08)
 *
 * step06 needs protein_expression merged into the HVG file (see step06 module).
 * step07 is RNA-only, reads HVG only.
 * step08 reads all three files directly (CLR + totalVI output + scVI output).
 *
 * Entry surface: called from main.nf's default workflow (after Phase A emits
 * hvg + clr channels) and from `-entry from_hvg` (user supplies HVG + CLR
 * on disk).
 */

include { step06_totalvi }    from '../modules/step06_totalvi'
include { step07_scvi }       from '../modules/step07_scvi'
include { step08_embedding }  from '../modules/step08_embedding'

workflow PHASE_B_EMBED {
    take:
    adata_qc_hvg    // path channel: adata_qc_hvg.h5ad (from step04 HVG selection)
    adata_qc_clr    // path channel: adata_qc_clr.h5ad (from step03 CLR normalization)

    main:
    // Train both embeddings in parallel
    step06_totalvi(adata_qc_hvg, adata_qc_clr)
    step07_scvi(adata_qc_hvg)

    // Compare and select primary embedding (all three files needed)
    step08_embedding(
        adata_qc_clr,
        step06_totalvi.out.totalvi,
        step07_scvi.out.scvi,
    )

    emit:
    totalvi    = step06_totalvi.out.totalvi
    scvi       = step07_scvi.out.scvi
    embedded   = step08_embedding.out.embedded
    comparison = step08_embedding.out.comparison
}

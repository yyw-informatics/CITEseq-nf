#!/usr/bin/env python
"""Step 8: Embedding Comparison and Selection
================================================

Objective
---------
Head-to-head comparison of totalVI (joint RNA+protein, Step 6) and scVI
(RNA-only, Step 7) latent spaces using silhouette scores. Selects the
primary embedding for all downstream neighbor graphs, UMAP, and
clustering. When the panel carries few functional ADT markers, totalVI's
protein contribution may be negligible — this step detects that
empirically by computing the biological silhouette on
BIOLOGICAL_SEPARATION_MARKER.

Pipeline Position
-----------------
Phase:    2.2 (Tier 1)
Depends:  Steps 6 (totalVI), 7 (scVI)
Feeds:    Downstream ADT gating, RNA scoring, classification — any step
          that reads the neighbor graph / UMAP
Compute:  CPU, 8-16 GB, 5-15 min

Inputs
------
- data/processed/adata_qc_clr.h5ad     — From Step 3: QC + CLR-normalized ADT (h5ad)
- data/processed/adata_totalvi.h5ad     — totalVI latent + denoised ADT + fg_prob (h5ad)
- data/processed/adata_scvi.h5ad        — scVI latent, RNA-only (h5ad)

Outputs
-------
- results/phase2_embedding/embedding_comparison.csv  — Silhouette scores + decision (CSV)
- data/processed/adata_qc_embedded.h5ad              — Merged embeddings, neighbor graph, UMAP (h5ad)

Config Constants
----------------
- N_NEIGHBORS (15)                  — k for sc.pp.neighbors on selected embedding
- UMAP_MIN_DIST (0.3)               — min_dist for sc.tl.umap
- RANDOM_SEED (42)                   — Reproducible subsampling for silhouette computation
- ADT_FUNCTIONAL_MARKERS             — Panel-driven functional marker list
- BIOLOGICAL_SEPARATION_MARKER       — Canonical ADT name for the bio-silhouette
                                       proxy; None → skip bio-silhouette
- PROJECT_ROOT                       — Base path for all I/O
"""
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from config import *

import argparse
import scanpy as sc
import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score


def main(branch=None):
    paths = branched_paths(branch)
    suffix = f"_{branch}" if branch else ""
    if branch:
        print(f"=== Branch mode: '{branch}' ===")
    adata_rna = sc.read_h5ad(paths["ADATA_QC_CLR"])

    # Steps 6-7 write separate files so they can train in parallel on the cluster
    adata_totalvi = sc.read_h5ad(paths["ADATA_TOTALVI"])
    adata_scvi_file = sc.read_h5ad(paths["ADATA_SCVI"])

    # Merge totalVI outputs into master — unconditional to always use latest latent space
    adata_rna.obsm["X_totalVI"] = adata_totalvi.obsm["X_totalVI"]
    # Carry forward denoised ADT + fg_prob — Step 9 (gating) and Step 5b (KDE
    # validation, optional re-run) read these from master
    for col in adata_totalvi.obs.columns:
        if col.startswith("ADT_") and (col.endswith("_denoised") or col.endswith("_fg_prob")):
            adata_rna.obs[col] = adata_totalvi.obs[col].values
    if "leiden_totalvi" in adata_totalvi.obs:
        adata_rna.obs["leiden_totalvi"] = adata_totalvi.obs["leiden_totalvi"].values
    if "X_umap" in adata_totalvi.obsm:
        adata_rna.obsm["X_umap_totalvi"] = adata_totalvi.obsm["X_umap"]

    # Merge scVI — also unconditional to always use latest latent space
    adata_rna.obsm["X_scVI"] = adata_scvi_file.obsm["X_scVI"]

    latent_totalvi = adata_rna.obsm["X_totalVI"]
    latent_scvi = adata_rna.obsm["X_scVI"]
    assert latent_totalvi.shape[0] == latent_scvi.shape[0] == len(adata_rna), (
        f"Embedding row count mismatch: totalVI={latent_totalvi.shape[0]}, "
        f"scVI={latent_scvi.shape[0]}, adata={len(adata_rna)}"
    )

    # ── 8a: Batch-mixing silhouette ──────────────────────────────
    np.random.seed(RANDOM_SEED)  # RANDOM_SEED=42 from config
    subsample_size = min(5000, len(adata_rna))  # cap for silhouette runtime
    idx = np.random.choice(len(adata_rna), subsample_size, replace=False)

    # Batch = sample_id (per library); lower score = better cross-sample mixing
    sil_batch_totalvi = silhouette_score(
        latent_totalvi[idx], adata_rna.obs["sample_id"].values[idx]
    )
    sil_batch_scvi = silhouette_score(
        latent_scvi[idx], adata_rna.obs["sample_id"].values[idx]
    )

    print(f"Batch silhouette (lower = better mixing):")
    print(f"  totalVI: {sil_batch_totalvi:.4f}")
    print(f"  scVI:    {sil_batch_scvi:.4f}")

    # ── 8b: Biological separation (+/- on BIOLOGICAL_SEPARATION_MARKER) ──
    # Rough split via totalVI fg_prob > 0.5 on the configured marker;
    # skipped entirely when the marker is None.
    if BIOLOGICAL_SEPARATION_MARKER:
        fg_col = f"ADT_{BIOLOGICAL_SEPARATION_MARKER}_fg_prob"
        if fg_col not in adata_rna.obs.columns:
            raise KeyError(
                f"BIOLOGICAL_SEPARATION_MARKER='{BIOLOGICAL_SEPARATION_MARKER}' "
                f"but adata has no column '{fg_col}'. Set "
                f"CITESEQ_BIOLOGICAL_SEPARATION_MARKER to a marker in "
                f"ADT_FUNCTIONAL_MARKERS or unset it."
            )
        bio_label = (adata_rna.obs[fg_col] > 0.5).astype(str).values

        if len(np.unique(bio_label[idx])) > 1:
            sil_bio_totalvi = silhouette_score(latent_totalvi[idx], bio_label[idx])
            sil_bio_scvi = silhouette_score(latent_scvi[idx], bio_label[idx])
            print(f"\nBiological silhouette ({BIOLOGICAL_SEPARATION_MARKER}+/-, higher = better separation):")
            print(f"  totalVI: {sil_bio_totalvi:.4f}")
            print(f"  scVI:    {sil_bio_scvi:.4f}")
        else:
            sil_bio_totalvi = 0
            sil_bio_scvi = 0
            print(f"WARNING: {BIOLOGICAL_SEPARATION_MARKER} fg_prob does not produce two groups for silhouette")
    else:
        sil_bio_totalvi = np.nan
        sil_bio_scvi = np.nan
        print("\nBiological silhouette: SKIPPED (BIOLOGICAL_SEPARATION_MARKER is None)")

    # ── 8c: Embedding selection decision ─────────────────────────────
    # Prefer totalVI when it wins biological separation; scVI otherwise.
    # With a small functional-marker set, protein's contribution may be negligible.
    totalvi_better_batch = sil_batch_totalvi <= sil_batch_scvi
    if np.isnan(sil_bio_totalvi) or np.isnan(sil_bio_scvi):
        totalvi_better_bio = False
    else:
        totalvi_better_bio = sil_bio_totalvi >= sil_bio_scvi

    if not BIOLOGICAL_SEPARATION_MARKER:
        # No biological proxy available: fall back to batch-mixing only.
        if totalvi_better_batch:
            primary_embedding = "X_totalVI"
            decision = "totalVI: better batch mixing (biological silhouette skipped)"
        else:
            primary_embedding = "X_scVI"
            decision = "scVI: better batch mixing (biological silhouette skipped)"
    elif totalvi_better_batch and totalvi_better_bio:
        primary_embedding = "X_totalVI"
        decision = "totalVI: better batch mixing AND better protein-marker separation"
    elif totalvi_better_bio:
        primary_embedding = "X_totalVI"
        decision = "totalVI: better protein-marker separation (batch mixing similar)"
    else:
        primary_embedding = "X_scVI"
        decision = "scVI: RNA-only embedding preferred (totalVI protein contribution negligible)"

    print(f"\nDECISION: Primary embedding = {primary_embedding}")
    print(f"  Rationale: {decision}")

    # Audit trail CSV — consumed by tests/test_step08 and decision graph
    bio_row_label = (
        f"bio_silhouette_{BIOLOGICAL_SEPARATION_MARKER}"
        if BIOLOGICAL_SEPARATION_MARKER else "bio_silhouette_skipped"
    )
    comparison = pd.DataFrame({
        "metric": ["batch_silhouette", bio_row_label, "decision"],
        "totalvi": [sil_batch_totalvi, sil_bio_totalvi, primary_embedding],
        "scvi": [sil_batch_scvi, sil_bio_scvi, ""],
    })
    os.makedirs(RESULTS_PHASE2_EMBEDDING, exist_ok=True)
    comparison.to_csv(os.path.join(RESULTS_PHASE2_EMBEDDING, f"embedding_comparison{suffix}.csv"),
                      index=False)

    # Neighbor graph + UMAP on winner — all downstream steps inherit this layout
    sc.pp.neighbors(adata_rna, use_rep=primary_embedding,
                    n_neighbors=N_NEIGHBORS)   # N_NEIGHBORS=15 from config
    sc.tl.umap(adata_rna, min_dist=UMAP_MIN_DIST)  # UMAP_MIN_DIST=0.3 from config
    # Tag for downstream consumers (Step 9+ read this to know which rep to use)
    adata_rna.uns["primary_embedding"] = primary_embedding
    adata_rna.write_h5ad(paths["ADATA_QC_EMBEDDED"])

    print(f"\nStep 8 complete: primary_embedding = {primary_embedding}")


def validate(branch=None):
    paths = branched_paths(branch)
    suffix = f"_{branch}" if branch else ""
    results = {}

    # ── 8v-a: Comparison CSV exists ──────────────────────────────────
    _path = os.path.join(RESULTS_PHASE2_EMBEDDING, f"embedding_comparison{suffix}.csv")
    _exists = os.path.exists(_path)
    results["comparison_csv_exists"] = {
        "status": "PASS" if _exists else "FAIL",
        "detail": f"{_path} {'exists' if _exists else 'MISSING'}",
        "expected": "file exists",
        "actual": str(_exists),
    }

    # ── 8v-b: Decision row present ───────────────────────────────────
    if _exists:
        import pandas as pd
        comp = pd.read_csv(_path)
        has_decision = "decision" in comp["metric"].values
        if has_decision:
            decision_val = comp.loc[comp["metric"] == "decision", "totalvi"].values[0]
            results["decision_documented"] = {
                "status": "PASS" if decision_val else "FAIL",
                "detail": f"Decision: {decision_val}",
                "expected": "primary_embedding decision documented",
                "actual": str(decision_val),
            }
        else:
            results["decision_documented"] = {
                "status": "FAIL",
                "detail": "No decision row in comparison CSV",
                "expected": "decision row present",
                "actual": "missing",
            }

    # ── 8v-c: primary_embedding persisted in h5ad ────────────────────
    _h5ad_path = paths["ADATA_QC_EMBEDDED"]
    if os.path.exists(_h5ad_path):
        import anndata
        adata_check = anndata.read_h5ad(_h5ad_path, backed="r")
        has_key = "primary_embedding" in adata_check.uns
        has_umap = "X_umap" in adata_check.obsm
        adata_check.file.close()
        results["primary_embedding_persisted"] = {
            "status": "PASS" if has_key else "FAIL",
            "detail": f"uns['primary_embedding'] {'found' if has_key else 'MISSING'} in adata_qc.h5ad",
            "expected": "primary_embedding key in uns",
            "actual": str(has_key),
        }
        results["umap_persisted"] = {
            "status": "PASS" if has_umap else "FAIL",
            "detail": f"X_umap {'found' if has_umap else 'MISSING'} in adata_qc.h5ad",
            "expected": "X_umap in obsm",
            "actual": str(has_umap),
        }
    else:
        results["primary_embedding_persisted"] = {
            "status": "FAIL",
            "detail": f"{_h5ad_path} not found",
            "expected": "adata_qc.h5ad exists with primary_embedding in uns",
            "actual": "file missing",
        }

    # ── 8v-d: UMAP structure (requires visual inspection) ───────────
    results["umap_structure"] = {
        "status": "WARN",
        "detail": "Manual check: UMAP of selected embedding should show biologically meaningful structure",
        "expected": "biologically meaningful UMAP structure",
        "actual": "check visualization",
    }

    # ── 8v-e: Print summary ─────────────────────────────────────────
    print(f"\n{'='*60}")
    print(f"Validation: Step 8 - Embedding Comparison and Selection")
    print(f"{'='*60}")
    for name, result in results.items():
        status = result["status"]
        icon = {"PASS": "[OK]", "FAIL": "[!!]", "WARN": "[??]"}[status]
        print(f"  {icon} {status}: {name}")
        if status != "PASS":
            print(f"       Expected: {result['expected']}")
            print(f"       Actual:   {result['actual']}")
    n_pass = sum(1 for r in results.values() if r["status"] == "PASS")
    n_fail = sum(1 for r in results.values() if r["status"] == "FAIL")
    n_warn = sum(1 for r in results.values() if r["status"] == "WARN")
    print(f"\nSummary: {n_pass} PASS, {n_warn} WARN, {n_fail} FAIL")
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Step 8: Embedding Comparison and Selection")
    parser.add_argument("--branch", type=str, default=None,
                        help="Branch name for parallel decontamination comparison")
    args = parser.parse_args()
    main(branch=args.branch)
    results = validate(branch=args.branch)
    if any(r["status"] == "FAIL" for r in results.values()):
        print("\nValidation FAILED. Troubleshooting:")
        print("Default to totalVI (protein contribution is valuable for classification even if embedding is equivalent).")
        sys.exit(1)

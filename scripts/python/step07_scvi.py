#!/usr/bin/env python
"""Step 7: scVI RNA-Only Model
=================================================================

Objective
---------
Train scVI on the QC'd HVG-subset RNA matrix to produce an RNA-only latent
(X_scVI) for direct comparison against the joint totalVI embedding.

Design Rationale
----------------
The training schedule is set so KL warmup completes before early stopping
fires; otherwise the latent space is barely regularised when training halts:

  1. max_epochs sized to match the KL warmup schedule
  2. n_epochs_kl_warmup passed explicitly via plan_kwargs
  3. Relaxed early stopping (patience/min_delta) — scvi-tools defaults
  4. Architecture: n_hidden=256, n_layers=2 (KB recommended for this scale)
  5. DataLoader num_workers=8 to eliminate the bottleneck warning

Pipeline Position
-----------------
Phase:    2.2 (Tier 1)
Depends:  Step 5
Feeds:    Step 8 (compare embeddings)
Compute:  GPU recommended, 8-16 GB RAM, ~20-40 min (GPU) with early stopping

Inputs
------
- data/processed/adata_qc_hvg.h5ad

Outputs
-------
- models/scvi_model/
- data/processed/adata_scvi.h5ad
- figures/phase2_scvi_umap_overview.png

Config Constants
----------------
All SCVI_* constants from config.py.
"""
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from config import *

import argparse
import ctypes
import gc

import scanpy as sc
import scvi


def _malloc_trim():
    """Release glibc free-list memory back to the OS."""
    try:
        ctypes.CDLL("libc.so.6").malloc_trim(0)
    except Exception:
        pass


def main(branch=None):
    model_dir = os.path.join(MODELS_DIR, "scvi_model")
    adata_path = os.path.join(PROCESSED_DIR, "adata_scvi.h5ad")
    if branch:
        suffix = f"_{branch}"
        model_dir = os.path.join(MODELS_DIR, f"scvi_model{suffix}")
        adata_path = os.path.join(PROCESSED_DIR, f"adata_scvi{suffix}.h5ad")
        print(f"=== Branch mode: '{branch}' ===")
    else:
        suffix = ""

    paths = branched_paths(branch)
    adata_rna = sc.read_h5ad(paths["ADATA_QC_HVG"])
    adata_scvi = adata_rna[:, adata_rna.var["highly_variable"]].copy()
    del adata_rna
    gc.collect()

    # -- Register with scvi-tools ----------------------------------------------
    scvi.model.SCVI.setup_anndata(adata_scvi, layer="counts", batch_key="sample_id")

    # -- Initialize model with upgraded architecture ---------------------------
    os.makedirs(model_dir, exist_ok=True)
    os.makedirs(os.path.dirname(adata_path), exist_ok=True)

    model_scvi = scvi.model.SCVI(
        adata_scvi,
        n_latent=SCVI_N_LATENT,
        n_hidden=SCVI_N_HIDDEN,
        n_layers=SCVI_N_LAYERS,
        gene_likelihood=SCVI_GENE_LIKELIHOOD,
    )

    # -- Train -----------------------------------------------------------------
    print(f"\n=== Training Config ===")
    print(f"  n_latent={SCVI_N_LATENT}, n_hidden={SCVI_N_HIDDEN}, n_layers={SCVI_N_LAYERS}")
    print(f"  max_epochs={SCVI_MAX_EPOCHS}, lr={SCVI_LR}")
    print(f"  n_epochs_kl_warmup={SCVI_KL_WARMUP}")
    print(f"  early_stopping_patience={SCVI_EARLY_STOPPING_PATIENCE}, min_delta={SCVI_EARLY_STOPPING_MIN_DELTA}")

    model_scvi.train(
        max_epochs=SCVI_MAX_EPOCHS,
        early_stopping=SCVI_EARLY_STOPPING,
        early_stopping_patience=SCVI_EARLY_STOPPING_PATIENCE,
        early_stopping_min_delta=SCVI_EARLY_STOPPING_MIN_DELTA,
        early_stopping_monitor=SCVI_EARLY_STOPPING_MONITOR,
        plan_kwargs={
            "lr": SCVI_LR,
            "n_epochs_kl_warmup": SCVI_KL_WARMUP,
        },
    )

    gc.collect()
    _malloc_trim()

    # -- Print training diagnostics --------------------------------------------
    n_trained = len(model_scvi.history["elbo_train"])
    print(f"\n=== Training Diagnostics ===")
    print(f"Epochs trained: {n_trained}/{SCVI_MAX_EPOCHS}")
    for key in sorted(model_scvi.history.keys()):
        vals = model_scvi.history[key]
        print(f"  {key}: {len(vals)} entries, last={vals.iloc[-1].values[0]:.4f}")
    if n_trained < SCVI_MAX_EPOCHS:
        print(f"Early stopping triggered at epoch {n_trained}/{SCVI_MAX_EPOCHS}")
    else:
        print(f"WARNING: Trained full {SCVI_MAX_EPOCHS} epochs -- early stopping did NOT trigger")

    # -- Extract latent and build graph ----------------------------------------
    latent_scvi = model_scvi.get_latent_representation()
    adata_scvi.obsm["X_scVI"] = latent_scvi

    gc.collect()
    _malloc_trim()

    sc.pp.neighbors(adata_scvi, use_rep="X_scVI", n_neighbors=N_NEIGHBORS)
    sc.tl.umap(adata_scvi, min_dist=UMAP_MIN_DIST)
    sc.tl.leiden(adata_scvi, resolution=LEIDEN_RESOLUTION, key_added="leiden_scvi")

    # -- Figures ---------------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(FIGURES_DIR, exist_ok=True)
    sc.pl.umap(adata_scvi, color=["leiden_scvi", "sample_id"],
               ncols=2, show=False)
    plt.savefig(os.path.join(FIGURES_DIR, f"phase2_scvi_umap_overview{suffix}.png"),
                dpi=150, bbox_inches="tight")
    plt.close("all")

    # -- Save ------------------------------------------------------------------
    model_scvi.save(model_dir, overwrite=True)
    adata_scvi.write_h5ad(adata_path)

    print(f"Step 7 complete: scVI trained, latent dim={latent_scvi.shape[1]}")
    print(f"  Saved to {adata_path} (Step 8 will merge into master adata)")


def validate(branch=None):
    paths = branched_paths(branch)
    model_dir = paths["SCVI_MODEL_DIR"]
    adata_path = paths["ADATA_SCVI"]
    results = {}

    _exists = os.path.isdir(model_dir)
    results["scvi_model_dir_exists"] = {
        "status": "PASS" if _exists else "FAIL",
        "detail": f"{model_dir} {'exists' if _exists else 'MISSING'}",
        "expected": "directory exists",
        "actual": str(_exists),
    }

    _exists = os.path.exists(adata_path)
    results["adata_scvi_exists"] = {
        "status": "PASS" if _exists else "FAIL",
        "detail": f"{adata_path} {'exists' if _exists else 'MISSING'}",
        "expected": "file exists",
        "actual": str(_exists),
    }

    if _exists:
        adata = sc.read_h5ad(adata_path)

        has_latent = "X_scVI" in adata.obsm
        if has_latent:
            latent_shape = adata.obsm["X_scVI"].shape
            latent_ok = latent_shape[1] == SCVI_N_LATENT
            results["latent_shape"] = {
                "status": "PASS" if latent_ok else "FAIL",
                "detail": f"X_scVI shape: {latent_shape}",
                "expected": f"N x {SCVI_N_LATENT}",
                "actual": f"{latent_shape[0]} x {latent_shape[1]}",
            }
        else:
            results["latent_shape"] = {
                "status": "FAIL",
                "detail": "X_scVI missing from obsm",
                "expected": f"N x {SCVI_N_LATENT}",
                "actual": "missing",
            }

    results["training_convergence"] = {
        "status": "WARN",
        "detail": "Check: kl_weight should reach 1.0; early stopping should fire after warmup completes",
        "expected": "converged training with full KL, batch-mixed UMAP",
        "actual": "check console output",
    }

    print(f"\n{'='*60}")
    print(f"Validation: Step 7 - scVI")
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
    parser = argparse.ArgumentParser(description="Step 7: scVI RNA-Only Model")
    parser.add_argument("--branch", type=str, default=None,
                        help="Branch name for parallel decontamination comparison")
    args = parser.parse_args()
    main(branch=args.branch)
    results = validate(branch=args.branch)
    if any(r["status"] == "FAIL" for r in results.values()):
        print("\nValidation FAILED. Troubleshooting:")
        print("Check kl_weight in history — it should reach 1.0 before early stopping fires.")
        print("If ELBO does not converge, try reducing lr to 2e-3.")
        sys.exit(1)

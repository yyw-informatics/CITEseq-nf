#!/usr/bin/env python
"""Step 6: totalVI Joint RNA+ADT Model
=============================================================================

Objective
---------
Train totalVI on the QC'd HVG-subset RNA matrix plus the ADT panel, producing
a joint latent (X_totalVI), denoised protein expression, and per-marker
foreground probabilities.

Design Rationale
----------------
The training schedule is tuned to avoid the reconstruction plateau that
otherwise stalls this panel:

  1. LR scheduler monitors reconstruction_loss_validation directly (not elbo)
  2. lr_threshold aligned with early_stopping_min_delta — noise can't reset
     the LR-reduction patience window
  3. lr_patience kept short for faster LR cuts when reconstruction stalls
  4. pro_recons_weight upweights the small-cardinality protein loss branch
     against the ~4K-gene RNA branch
  5. reduce_lr_on_plateau passed via plan_kwargs (correct API for
     AdversarialTrainingPlan)
  6. Protein obsm preserved as a DataFrame with named columns so the
     marker-to-column mapping is unambiguous downstream.

Pipeline Position
-----------------
Phase:    2.1 (Tier 1)
Depends:  Steps 4 (CLR/protein_expression), 5 (HVGs)
Feeds:    Steps 8 (compare), 9 (gating), 11 (classify)
Compute:  GPU recommended, 16-32 GB RAM, ~30-60 min (GPU)

Inputs
------
- data/processed/adata_qc_hvg.h5ad

Outputs
-------
- models/totalvi_model/
- data/processed/adata_totalvi.h5ad
- results/phase2_embedding/totalvi_training_history.csv
- results/phase2_embedding/totalvi_full_history.csv

Config Constants
----------------
All TOTALVI_* constants from config.py.
"""
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from config import *

import argparse
import ctypes
import gc
import resource
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvi


def _malloc_trim():
    """Release glibc free-list memory back to the OS."""
    try:
        ctypes.CDLL("libc.so.6").malloc_trim(0)
    except Exception:
        pass


def _rss_gb():
    """Current process RSS in GB."""
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6


def main(branch=None):
    model_dir = os.path.join(MODELS_DIR, "totalvi_model")
    adata_path = os.path.join(PROCESSED_DIR, "adata_totalvi.h5ad")
    if branch:
        suffix = f"_{branch}"
        model_dir = os.path.join(MODELS_DIR, f"totalvi_model{suffix}")
        adata_path = os.path.join(PROCESSED_DIR, f"adata_totalvi{suffix}.h5ad")
        print(f"=== Branch mode: '{branch}' ===")
    else:
        suffix = ""

    paths = branched_paths(branch)
    adata_rna = sc.read_h5ad(paths["ADATA_QC_HVG"])

    # -- Subset to HVGs -------------------------------------------------------
    adata_totalvi = adata_rna[:, adata_rna.var["highly_variable"]].copy()
    del adata_rna
    gc.collect()

    # -- Validate protein matrix -----------------------------------------------
    assert "protein_expression" in adata_totalvi.obsm, "Missing protein_expression in obsm"
    n_adt = len(ADT_FUNCTIONAL_MARKERS)
    assert adata_totalvi.obsm["protein_expression"].shape[1] == n_adt, \
        f"Expected {n_adt} ADT markers, got {adata_totalvi.obsm['protein_expression'].shape[1]}"

    # -- Preserve protein column names so downstream code can index by marker --
    prot = adata_totalvi.obsm["protein_expression"]
    if not isinstance(prot, pd.DataFrame):
        adata_totalvi.obsm["protein_expression"] = pd.DataFrame(
            prot,
            index=adata_totalvi.obs_names,
            columns=ADT_FUNCTIONAL_MARKERS,
        )
        print(f"Protein obsm converted to DataFrame with columns: {ADT_FUNCTIONAL_MARKERS}")
    else:
        print(f"Protein obsm already a DataFrame, columns: {list(prot.columns)}")

    # -- Register with scvi-tools ----------------------------------------------
    scvi.model.TOTALVI.setup_anndata(
        adata_totalvi,
        layer="counts",
        protein_expression_obsm_key="protein_expression",
        batch_key="sample_id",
    )

    # -- Initialize model ------------------------------------------------------
    os.makedirs(model_dir, exist_ok=True)
    os.makedirs(os.path.dirname(adata_path), exist_ok=True)
    os.makedirs(RESULTS_PHASE2_EMBEDDING, exist_ok=True)

    model_totalvi = scvi.model.TOTALVI(
        adata_totalvi,
        n_latent=TOTALVI_N_LATENT,
        gene_likelihood=TOTALVI_GENE_LIKELIHOOD,
        empirical_protein_background_prior=TOTALVI_EMPIRICAL_PRIOR,
    )

    # -- Train with scheduler targeting reconstruction plateau -----------------
    print(f"Pre-training RSS: {_rss_gb():.2f} GB")
    print(f"\n=== Training Config ===")
    print(f"  max_epochs={TOTALVI_MAX_EPOCHS}, lr={TOTALVI_LR}")
    print(f"  n_epochs_kl_warmup={TOTALVI_KL_WARMUP_EPOCHS}")
    print(f"  early_stopping_patience={TOTALVI_EARLY_STOPPING_PATIENCE}, min_delta={TOTALVI_EARLY_STOPPING_MIN_DELTA}")
    print(f"  lr_patience={TOTALVI_LR_PATIENCE}, lr_factor={TOTALVI_LR_FACTOR}, lr_min={TOTALVI_LR_MIN}")
    print(f"  lr_threshold={TOTALVI_LR_THRESHOLD}, lr_scheduler_metric={TOTALVI_LR_SCHEDULER_METRIC}")
    print(f"  pro_recons_weight={TOTALVI_PRO_RECONS_WEIGHT}")
    t0 = time.time()

    model_totalvi.train(
        max_epochs=TOTALVI_MAX_EPOCHS,
        batch_size=TOTALVI_BATCH_SIZE_TRAIN,
        lr=TOTALVI_LR,
        n_epochs_kl_warmup=TOTALVI_KL_WARMUP_EPOCHS,
        early_stopping=TOTALVI_EARLY_STOPPING,
        early_stopping_patience=TOTALVI_EARLY_STOPPING_PATIENCE,
        early_stopping_min_delta=TOTALVI_EARLY_STOPPING_MIN_DELTA,
        early_stopping_monitor=TOTALVI_EARLY_STOPPING_MONITOR,
        plan_kwargs={
            "reduce_lr_on_plateau": True,
            "lr_patience": TOTALVI_LR_PATIENCE,
            "lr_factor": TOTALVI_LR_FACTOR,
            "lr_min": TOTALVI_LR_MIN,
            "lr_threshold": TOTALVI_LR_THRESHOLD,
            "lr_scheduler_metric": TOTALVI_LR_SCHEDULER_METRIC,
            "pro_recons_weight": TOTALVI_PRO_RECONS_WEIGHT,
        },
    )

    train_time = time.time() - t0
    print(f"\n=== Training Diagnostics ===")
    print(f"Training wall time: {train_time/60:.1f} min")
    print(f"Post-training RSS: {_rss_gb():.2f} GB")
    print(f"History keys: {list(model_totalvi.history.keys())}")
    for key in sorted(model_totalvi.history.keys()):
        vals = model_totalvi.history[key]
        print(f"  {key}: {len(vals)} entries, last={vals.iloc[-1].values[0]:.4f}")
    print(f"is_trained: {model_totalvi.is_trained}")

    n_trained = len(model_totalvi.history["elbo_train"])
    if n_trained < TOTALVI_MAX_EPOCHS:
        print(f"Early stopping triggered at epoch {n_trained}/{TOTALVI_MAX_EPOCHS}")
    else:
        print(f"WARNING: Trained full {TOTALVI_MAX_EPOCHS} epochs -- early stopping did NOT trigger")

    # -- Save training history -------------------------------------------------
    train_elbo = model_totalvi.history["elbo_train"]
    train_elbo.to_csv(os.path.join(RESULTS_PHASE2_EMBEDDING, f"totalvi_training_history{suffix}.csv"))

    all_history = {}
    for key in model_totalvi.history.keys():
        all_history[key] = model_totalvi.history[key].values.flatten()
    max_len = max(len(v) for v in all_history.values())
    for key in all_history:
        if len(all_history[key]) < max_len:
            all_history[key] = np.pad(all_history[key], (0, max_len - len(all_history[key])),
                                       constant_values=np.nan)
    pd.DataFrame(all_history).to_csv(
        os.path.join(RESULTS_PHASE2_EMBEDDING, f"totalvi_full_history{suffix}.csv"), index=False
    )

    gc.collect()
    _malloc_trim()
    print(f"Post-gc RSS: {_rss_gb():.2f} GB")

    # -- Save model immediately ------------------------------------------------
    print(f"\n=== Save model (RSS: {_rss_gb():.2f} GB) ===")
    model_totalvi.save(model_dir, overwrite=True)
    print("Model saved.")

    # -- Extract latent representation -----------------------------------------
    print(f"\n=== get_latent_representation (RSS: {_rss_gb():.2f} GB) ===")
    latent_totalvi = model_totalvi.get_latent_representation()
    adata_totalvi.obsm["X_totalVI"] = latent_totalvi
    gc.collect()
    print(f"Done (RSS: {_rss_gb():.2f} GB)")

    # -- Denoised protein expression -------------------------------------------
    print(f"\n=== get_normalized_expression n_samples=5 (RSS: {_rss_gb():.2f} GB) ===")
    _gene_norm, protein_norm_totalvi = model_totalvi.get_normalized_expression(
        return_mean=True, n_samples=5, batch_size=TOTALVI_BATCH_SIZE_INFERENCE
    )
    del _gene_norm
    adata_totalvi.obsm["protein_denoised"] = protein_norm_totalvi.values

    returned_cols = list(protein_norm_totalvi.columns)
    print(f"Protein columns from totalVI: {returned_cols}")
    print(f"ADT_FUNCTIONAL_MARKERS (config): {ADT_FUNCTIONAL_MARKERS}")
    for marker in ADT_FUNCTIONAL_MARKERS:
        if marker in returned_cols:
            adata_totalvi.obs[f"ADT_{marker}_denoised"] = protein_norm_totalvi[marker].values
        else:
            idx = ADT_FUNCTIONAL_MARKERS.index(marker)
            adata_totalvi.obs[f"ADT_{marker}_denoised"] = protein_norm_totalvi.iloc[:, idx].values
            print(f"  WARNING: '{marker}' not in returned columns, used position {idx}")
    gc.collect()
    _malloc_trim()
    print(f"Done (RSS: {_rss_gb():.2f} GB)")

    # -- Protein foreground probabilities --------------------------------------
    print(f"\n=== get_protein_foreground_probability n_samples=5 (RSS: {_rss_gb():.2f} GB) ===")
    protein_fg_prob = model_totalvi.get_protein_foreground_probability(
        n_samples=5, return_mean=True, batch_size=TOTALVI_BATCH_SIZE_INFERENCE
    )
    adata_totalvi.obsm["protein_fg_prob"] = protein_fg_prob.values
    fg_cols = list(protein_fg_prob.columns)
    for marker in ADT_FUNCTIONAL_MARKERS:
        if marker in fg_cols:
            adata_totalvi.obs[f"ADT_{marker}_fg_prob"] = protein_fg_prob[marker].values
        else:
            idx = ADT_FUNCTIONAL_MARKERS.index(marker)
            adata_totalvi.obs[f"ADT_{marker}_fg_prob"] = protein_fg_prob.iloc[:, idx].values
            print(f"  WARNING: '{marker}' not in fg_prob columns, used position {idx}")

    del protein_fg_prob
    gc.collect()
    _malloc_trim()
    print(f"Done (RSS: {_rss_gb():.2f} GB)")

    # Free model before UMAP
    del model_totalvi
    gc.collect()
    _malloc_trim()
    print(f"After model cleanup (RSS: {_rss_gb():.2f} GB)")

    # -- UMAP + Leiden ---------------------------------------------------------
    print(f"\n=== UMAP + Leiden (RSS: {_rss_gb():.2f} GB) ===")
    sc.pp.neighbors(adata_totalvi, use_rep="X_totalVI", n_neighbors=N_NEIGHBORS)
    sc.tl.umap(adata_totalvi, min_dist=UMAP_MIN_DIST)
    sc.tl.leiden(adata_totalvi, resolution=LEIDEN_RESOLUTION, key_added="leiden_totalvi")

    # -- Save ------------------------------------------------------------------
    print(f"\n=== Save adata (RSS: {_rss_gb():.2f} GB) ===")
    adata_totalvi.write_h5ad(adata_path)

    sc.pl.umap(adata_totalvi, color=["leiden_totalvi", "sample_id"],
               show=False)
    plt.savefig(os.path.join(FIGURES_DIR, f"phase2_totalvi_umap_overview{suffix}.png"),
                dpi=150, bbox_inches="tight")
    plt.close()

    sc.pl.umap(
        adata_totalvi,
        color=[f"ADT_{m}_denoised" for m in ADT_FUNCTIONAL_MARKERS],
        show=False,
    )
    plt.savefig(os.path.join(FIGURES_DIR, f"phase2_totalvi_umap_adt_denoised{suffix}.png"),
                dpi=150, bbox_inches="tight")
    plt.close()

    print(f"\nStep 6 complete: totalVI trained, latent dim={latent_totalvi.shape[1]}")
    print(f"Final RSS: {_rss_gb():.2f} GB")
    print(f"Total wall time: {(time.time() - t0)/60:.1f} min")


def validate(branch=None):
    paths = branched_paths(branch)
    suffix = f"_{branch}" if branch else ""
    model_dir = paths["TOTALVI_MODEL_DIR"]
    adata_path = paths["ADATA_TOTALVI"]
    results = {}

    _exists = os.path.isdir(model_dir)
    results["totalvi_model_dir_exists"] = {
        "status": "PASS" if _exists else "FAIL",
        "detail": f"{model_dir} {'exists' if _exists else 'MISSING'}",
        "expected": "directory exists",
        "actual": str(_exists),
    }

    _exists = os.path.exists(adata_path)
    results["adata_totalvi_exists"] = {
        "status": "PASS" if _exists else "FAIL",
        "detail": f"{adata_path} {'exists' if _exists else 'MISSING'}",
        "expected": "file exists",
        "actual": str(_exists),
    }

    if _exists:
        adata = sc.read_h5ad(adata_path)

        has_latent = "X_totalVI" in adata.obsm
        if has_latent:
            latent_shape = adata.obsm["X_totalVI"].shape
            latent_ok = latent_shape[1] == TOTALVI_N_LATENT
            results["latent_shape"] = {
                "status": "PASS" if latent_ok else "FAIL",
                "detail": f"X_totalVI shape: {latent_shape}",
                "expected": f"N x {TOTALVI_N_LATENT}",
                "actual": f"{latent_shape[0]} x {latent_shape[1]}",
            }
        else:
            results["latent_shape"] = {
                "status": "FAIL",
                "detail": "X_totalVI missing from obsm",
                "expected": f"N x {TOTALVI_N_LATENT}",
                "actual": "missing",
            }

        n_adt = len(ADT_FUNCTIONAL_MARKERS)
        has_denoised = "protein_denoised" in adata.obsm
        if has_denoised:
            dn_shape = adata.obsm["protein_denoised"].shape
            dn_ok = dn_shape[1] == n_adt
            results["denoised_protein_shape"] = {
                "status": "PASS" if dn_ok else "FAIL",
                "detail": f"protein_denoised shape: {dn_shape}",
                "expected": f"N x {n_adt}",
                "actual": f"{dn_shape[0]} x {dn_shape[1]}",
            }
        else:
            results["denoised_protein_shape"] = {
                "status": "FAIL",
                "detail": "protein_denoised missing from obsm",
                "expected": f"N x {n_adt}",
                "actual": "missing",
            }

        has_fg = "protein_fg_prob" in adata.obsm
        if has_fg:
            fg_shape = adata.obsm["protein_fg_prob"].shape
            fg_ok = fg_shape[1] == n_adt
            results["fg_prob_shape"] = {
                "status": "PASS" if fg_ok else "FAIL",
                "detail": f"protein_fg_prob shape: {fg_shape}",
                "expected": f"N x {n_adt}",
                "actual": f"{fg_shape[0]} x {fg_shape[1]}",
            }
        else:
            results["fg_prob_shape"] = {
                "status": "FAIL",
                "detail": "protein_fg_prob missing from obsm",
                "expected": f"N x {n_adt}",
                "actual": "missing",
            }

    _hist_path = os.path.join(RESULTS_PHASE2_EMBEDDING, f"totalvi_training_history{suffix}.csv")
    _hist_exists = os.path.exists(_hist_path)
    results["training_history_exists"] = {
        "status": "PASS" if _hist_exists else "FAIL",
        "detail": f"{_hist_path} {'exists' if _hist_exists else 'MISSING'}",
        "expected": "file exists",
        "actual": str(_hist_exists),
    }

    results["umap_batch_mixing"] = {
        "status": "WARN",
        "detail": "Manual check: UMAP should show batch mixing (samples interleaved)",
        "expected": "samples interleaved in UMAP",
        "actual": "check PNG plots",
    }

    results["fg_prob_bimodal"] = {
        "status": "WARN",
        "detail": "Manual check: functional-marker fg_prob histograms should show bimodal distributions (near 0 and near 1) for well-stained markers.",
        "expected": "bimodal distribution",
        "actual": "check PNG plots",
    }

    print(f"\n{'='*60}")
    print(f"Validation: Step 6 - totalVI")
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
    parser = argparse.ArgumentParser(description="Step 6: totalVI Joint RNA+ADT Model")
    parser.add_argument("--branch", type=str, default=None,
                        help="Branch name for parallel decontamination comparison")
    args = parser.parse_args()
    main(branch=args.branch)
    results = validate(branch=args.branch)
    if any(r["status"] == "FAIL" for r in results.values()):
        print("\nValidation FAILED. Troubleshooting:")
        print("Check lr-Adam in full_history CSV -- it should drop from TOTALVI_LR during training.")
        print("If recon loss plateaus high, the marker-count ceiling may be fundamental for this panel.")
        print("Compare against scVI (Step 8) to determine if totalVI adds value.")
        sys.exit(1)

#!/usr/bin/env python
"""Step 3: CLR Normalization and Protein Expression Matrix
==========================================================

Objective
---------
Compute CLR-normalized ADT values for the functional markers and store raw
ADT counts in obsm['protein_expression'] for totalVI input (Step 6).

Pipeline Position
-----------------
Phase:    1.2 (Tier 1)
Depends:  Step 1
Feeds:    Steps 6 (totalVI), 9 (gating)
Compute:  CPU, <1 GB, <1 min

Inputs
------
- data/processed/adata_qc.h5ad       -- QC-passed RNA from Step 1
- data/processed/adata_adt_raw.h5ad   -- Raw ADT counts synced to QC barcodes

Outputs
-------
- data/processed/adata_qc_clr.h5ad   -- Distinct artifact: adata_qc + obsm['protein_expression'],
                                         uns['protein_names'], and ADT_* columns in .obs

Config Constants
----------------
- ADT_FUNCTIONAL_MARKERS  -- Panel-driven; non-functional markers are excluded upstream via panel YAML / config.
- ADT_CANONICAL_NAMES     -- Short names used throughout pipeline
- PROJECT_ROOT (path)     -- Top-level project directory
"""
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from config import *

import argparse
import scanpy as sc
import numpy as np
import scipy.sparse as sp


def main(branch=None):
    paths = branched_paths(branch)
    if branch:
        print(f"=== Branch mode: '{branch}' ===")
    adata_rna = sc.read_h5ad(paths["ADATA_QC"])
    adata_adt = sc.read_h5ad(paths["ADATA_ADT_RAW"])

    print(f"ADT markers found: {list(adata_adt.var_names)}")

    # Barcode alignment: RNA and ADT must match exactly
    assert (adata_rna.obs_names == adata_adt.obs_names).all(), (
        "Barcode mismatch between RNA and ADT anndata objects"
    )

    # Map CellRanger names to canonical names, preserving paired order
    functional_pairs = [(name, ADT_RENAME_MAP.get(name, name))
                        for name in adata_adt.var_names
                        if ADT_RENAME_MAP.get(name, name) in ADT_FUNCTIONAL_MARKERS]
    functional_adt_names = [p[0] for p in functional_pairs]
    canonical_names = [p[1] for p in functional_pairs]

    print(f"Functional ADT markers: {functional_adt_names} -> {canonical_names}")
    assert len(functional_adt_names) == len(ADT_FUNCTIONAL_MARKERS), (
        f"Expected {len(ADT_FUNCTIONAL_MARKERS)} functional markers, got {len(functional_adt_names)}"
    )

    adt_counts = adata_adt[:, functional_adt_names].X
    if sp.issparse(adt_counts):
        adt_counts = adt_counts.toarray()
    adt_counts = adt_counts.astype(np.float64)

    # CLR: center log-ratio per cell -- standard ADT normalization baseline
    log_adt = np.log1p(adt_counts)
    mean_log = np.mean(log_adt, axis=1, keepdims=True)  # = log(geometric mean of x+1)
    clr_adt = log_adt - mean_log

    for i, marker in enumerate(canonical_names):
        adata_rna.obs[f"ADT_{marker}"] = clr_adt[:, i]

    # Raw counts stored for totalVI (Step 6); it needs untransformed integers
    adata_rna.obsm["protein_expression"] = adt_counts
    adata_rna.uns["protein_names"] = canonical_names

    adata_rna.write_h5ad(paths["ADATA_QC_CLR"])

    print(f"Step 3 complete: {len(canonical_names)} ADT markers CLR-normalized and stored")
    for marker in canonical_names:
        vals = adata_rna.obs[f"ADT_{marker}"]
        print(f"  ADT_{marker}: mean={vals.mean():.3f}, std={vals.std():.3f}, "
              f"min={vals.min():.3f}, max={vals.max():.3f}")


def validate(branch=None):
    paths = branched_paths(branch)
    results = {}

    _path = paths["ADATA_QC_CLR"]
    _exists = os.path.exists(_path)
    results["adata_qc_clr_exists"] = {
        "status": "PASS" if _exists else "FAIL",
        "detail": f"{_path} {'exists' if _exists else 'MISSING'}",
        "expected": "file exists",
        "actual": str(_exists),
    }

    if _exists:
        import scanpy as sc
        adata = sc.read_h5ad(_path)
        has_protein = "protein_expression" in adata.obsm
        results["protein_expression_in_obsm"] = {
            "status": "PASS" if has_protein else "FAIL",
            "detail": f"protein_expression {'found' if has_protein else 'MISSING'} in obsm",
            "expected": "obsm['protein_expression'] exists",
            "actual": str(has_protein),
        }

        if has_protein:
            shape = adata.obsm["protein_expression"].shape
            n_expected = len(ADT_FUNCTIONAL_MARKERS)
            shape_ok = shape[1] == n_expected
            results["protein_expression_shape"] = {
                "status": "PASS" if shape_ok else "FAIL",
                "detail": f"protein_expression shape: {shape}",
                "expected": f"N x {n_expected}",
                "actual": f"{shape[0]} x {shape[1]}",
            }

            pe = adata.obsm["protein_expression"]
            nonneg = (pe >= 0).all()
            results["protein_counts_nonneg"] = {
                "status": "PASS" if nonneg else "FAIL",
                "detail": f"protein_expression non-negative: {nonneg}",
                "expected": "all values >= 0",
                "actual": f"min={pe.min():.3f}",
            }

        adt_cols = [c for c in adata.obs.columns if c.startswith("ADT_")]
        n_expected_cols = len(ADT_FUNCTIONAL_MARKERS)
        has_adt_cols = len(adt_cols) >= n_expected_cols
        results["adt_obs_columns"] = {
            "status": "PASS" if has_adt_cols else "FAIL",
            "detail": f"ADT columns in obs: {adt_cols}",
            "expected": f"{n_expected_cols} ADT columns in .obs",
            "actual": f"{len(adt_cols)} columns: {adt_cols}",
        }

        if has_adt_cols:
            clr_cols = np.column_stack([adata.obs[c].values for c in adt_cols[:n_expected_cols]])
            row_means = np.mean(clr_cols, axis=1)
            centered = np.allclose(row_means, 0, atol=1e-6)
            results["clr_zero_centered"] = {
                "status": "PASS" if centered else "FAIL",
                "detail": f"CLR row means: max abs = {np.max(np.abs(row_means)):.2e}",
                "expected": "CLR values zero-centered per cell (mean ~ 0)",
                "actual": f"max |row mean| = {np.max(np.abs(row_means)):.2e}",
            }

    if "clr_zero_centered" not in results:
        results["clr_zero_centered"] = {
            "status": "WARN",
            "detail": "Could not check CLR centering (missing data)",
            "expected": "CLR values zero-centered",
            "actual": "skipped",
        }

    print(f"\n{'='*60}")
    print(f"Validation: Step 3 - CLR Normalization and Protein Expression Matrix")
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
    parser = argparse.ArgumentParser(description="Step 3: CLR Normalization and Protein Expression Matrix")
    parser.add_argument("--branch", type=str, default=None,
                        help="Branch name for parallel decontamination comparison")
    args = parser.parse_args()
    main(branch=args.branch)
    results = validate(branch=args.branch)
    if any(r["status"] == "FAIL" for r in results.values()):
        print("\nValidation FAILED. Troubleshooting:")
        print("Check that every ADT_CELLRANGER_NAMES entry in the panel YAML matches a row of adata_adt_raw.var_names; ADT_RENAME_MAP must cover all of them.")
        sys.exit(1)

#!/usr/bin/env python
"""Step 4: RNA Normalization and HVG Selection
===============================================

Objective
---------
Log-normalize RNA counts, select N_HVG batch-aware HVGs (seurat_v3 flavor),
and force-include any literature marker genes supplied via FORCED_HVGS.
The forced-HVG mechanism exists because rare populations (a few percent of
the dataset) often have their lineage markers outranked by more abundant
programs; studies pre-declare those genes so downstream classification and
hypothesis testing can still rely on them.

Pipeline Position
-----------------
Phase:    1.3 (Tier 1)
Depends:  Step 1
Feeds:    Steps 6 (totalVI), 7 (scVI), downstream classification / scoring
Compute:  CPU, 8-16 GB RAM, 5-10 min

Inputs
------
- data/processed/adata_qc.h5ad         -- QC-passed RNA with raw counts in .X and .layers['counts']

Outputs
-------
- data/processed/adata_normalized.h5ad  -- Log-normalized RNA with HVG annotations
- data/processed/adata_qc_hvg.h5ad     -- Raw counts + .var['highly_variable'] for Steps 6/7

Config Constants
----------------
- N_HVG (4000)                          -- Number of top HVGs to select
- FORCED_HVGS (list, study-supplied)    -- Markers force-included as HVGs
- HVG_SENTINEL_GENES (list, optional)   -- Sentinel genes the HVG set must retain; empty → skip check
- ADATA_QC, ADATA_QC_HVG, ADATA_NORMALIZED -- I/O paths
"""
import sys
import os
import logging

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from config import (
    N_HVG, FORCED_HVGS, HVG_SENTINEL_GENES,
    ADATA_QC, ADATA_QC_HVG, ADATA_NORMALIZED,
    branched_paths,
)

import argparse
import scanpy as sc

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ── 4c: Force-include literature marker genes ────────────────
# Loaded at module level so validate() can also reference them.
forced_genes = FORCED_HVGS


def main(branch=None):
    paths = branched_paths(branch)
    if branch:
        log.info(f"=== Branch mode: '{branch}' ===")
    adata_rna = sc.read_h5ad(paths["ADATA_QC"])

    # ── 4a: Log-normalize for scanpy-based analyses ────────────
    adata_norm = adata_rna.copy()
    sc.pp.normalize_total(adata_norm, target_sum=1e4)
    sc.pp.log1p(adata_norm)

    # ── 4b: Batch-aware HVG selection ──────────────────────────
    # seurat_v3 requires raw counts; per-batch ranking avoids one
    # library dominating the HVG set (N_HVG = 4000 from config).
    sc.pp.highly_variable_genes(
        adata_norm,
        n_top_genes=N_HVG,
        flavor="seurat_v3",
        batch_key="sample_id",
        layer="counts",
    )

    n_forced = 0
    for gene in forced_genes:
        if gene in adata_norm.var_names:
            if not adata_norm.var.loc[gene, "highly_variable"]:
                adata_norm.var.loc[gene, "highly_variable"] = True
                n_forced += 1

    n_hvg = adata_norm.var["highly_variable"].sum()
    log.info(f"HVGs: {n_hvg} (including {n_forced} additionally forced literature markers)")

    # Propagate HVG flags back to raw-count object for Steps 6/7
    adata_rna.var["highly_variable"] = False
    hvg_genes = adata_norm.var_names[adata_norm.var["highly_variable"]]
    adata_rna.var.loc[adata_rna.var_names.isin(hvg_genes), "highly_variable"] = True

    # ── 4d: Save ──────────────────────────────────────────────
    # Write normalized AnnData (log-normalized .X with HVG annotations)
    adata_norm.write_h5ad(paths["ADATA_NORMALIZED"])

    # Write HVG-annotated raw-count object for Steps 6/7 (totalVI/scVI).
    # Saved alongside (not overwriting) adata_qc.h5ad so the Step 1
    # checkpoint remains intact and re-runs are safe.
    adata_rna.write_h5ad(paths["ADATA_QC_HVG"])
    log.info(f"HVG-annotated raw counts saved to {paths['ADATA_QC_HVG']}")

    genes_found = [g for g in forced_genes if g in adata_norm.var_names]
    genes_missing = [g for g in forced_genes if g not in adata_norm.var_names]
    log.info(f"Forced genes found: {len(genes_found)}/{len(forced_genes)}")
    if genes_missing:
        log.warning(f"Missing genes: {genes_missing}")

    log.info(f"Step 4 complete: {adata_rna.n_obs} cells, {n_hvg} HVGs")


def validate(branch=None):
    paths = branched_paths(branch)
    _adata_normalized = paths["ADATA_NORMALIZED"]
    results = {}

    _exists = os.path.exists(_adata_normalized)
    results["adata_normalized_exists"] = {
        "status": "PASS" if _exists else "FAIL",
        "detail": f"{_adata_normalized} {'exists' if _exists else 'MISSING'}",
        "expected": "file exists",
        "actual": str(_exists),
    }

    if _exists:
        adata = sc.read_h5ad(_adata_normalized, backed="r")
        n_hvg = adata.var["highly_variable"].sum()
        hvg_ok = n_hvg >= 4000
        results["hvg_count"] = {
            "status": "PASS" if hvg_ok else "WARN",
            "detail": f"{n_hvg} HVGs selected",
            "expected": "~4000+ HVGs",
            "actual": str(n_hvg),
        }

        # Sentinel genes — see HVG_SENTINEL_GENES in config
        if HVG_SENTINEL_GENES:
            hvg_names = set(adata.var_names[adata.var["highly_variable"]])
            markers_present = [g for g in HVG_SENTINEL_GENES if g in hvg_names]
            markers_missing = [g for g in HVG_SENTINEL_GENES if g not in hvg_names]
            all_present = len(markers_missing) == 0
            results["key_markers_in_hvg"] = {
                "status": "PASS" if all_present else "WARN",
                "detail": f"Sentinels in HVG: {markers_present}, missing: {markers_missing}",
                "expected": f"all of {HVG_SENTINEL_GENES} in HVG set",
                "actual": f"{len(markers_present)}/{len(HVG_SENTINEL_GENES)} present",
            }

        # Steps 6/7 (scvi-tools) require raw counts in this layer
        has_counts = "counts" in adata.layers
        results["raw_counts_layer"] = {
            "status": "PASS" if has_counts else "FAIL",
            "detail": f"layers['counts'] {'found' if has_counts else 'MISSING'}",
            "expected": "layers['counts'] exists with raw counts",
            "actual": str(has_counts),
        }
        adata.file.close()

    log.info(f"\n{'='*60}")
    log.info(f"Validation: Step 4 - RNA Normalization and HVG Selection")
    log.info(f"{'='*60}")
    for name, result in results.items():
        status = result["status"]
        icon = {"PASS": "[OK]", "FAIL": "[!!]", "WARN": "[??]"}[status]
        log.info(f"  {icon} {status}: {name}")
        if status != "PASS":
            log.info(f"       Expected: {result['expected']}")
            log.info(f"       Actual:   {result['actual']}")
    n_pass = sum(1 for r in results.values() if r["status"] == "PASS")
    n_fail = sum(1 for r in results.values() if r["status"] == "FAIL")
    n_warn = sum(1 for r in results.values() if r["status"] == "WARN")
    log.info(f"\nSummary: {n_pass} PASS, {n_warn} WARN, {n_fail} FAIL")
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Step 4: RNA Normalization and HVG Selection")
    parser.add_argument("--branch", type=str, default=None,
                        help="Branch name for parallel decontamination comparison")
    args = parser.parse_args()
    main(branch=args.branch)
    results = validate(branch=args.branch)
    if any(r["status"] == "FAIL" for r in results.values()):
        log.error("Validation FAILED. Troubleshooting:")
        log.error("If seurat_v3 fails (requires raw counts in `.layers['counts']`), fall back to `flavor='cell_ranger'`. If many forced genes are missing, check gene name format (HUGO vs aliases).")
        sys.exit(1)

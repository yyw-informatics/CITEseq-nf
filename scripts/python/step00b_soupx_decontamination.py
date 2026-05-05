#!/usr/bin/env python
"""Step 0b: SoupX Ambient RNA Decontamination
==============================================

Objective
---------
Remove ambient RNA contamination from each CellRanger library using SoupX.
Empty droplets from the raw (unfiltered) matrix estimate the "soup" profile;
CellRanger graph-based clusters provide the per-cell prior. Corrected
counts are written as CellRanger-format h5 files with RNA corrected and
ADT rows unchanged, so Step 1 can read them transparently.

Pipeline Position
-----------------
Phase:    0 (pre-QC)
Depends:  CellRanger outputs (raw + filtered h5, analysis.tar.gz)
Feeds:    Step 1 (QC)
Compute:  CPU, 32 GB RAM, ~2-4 min/library sequential

Inputs
------
- 10X_cloud/outs/<sample>/raw_feature_bc_matrix.h5        -- unfiltered droplets
- 10X_cloud/outs/<sample>/filtered_feature_bc_matrix.h5   -- filtered cells
- 10X_cloud/outs/<sample>/analysis.tar.gz                 -- CellRanger clusters

Outputs
-------
- data/processed/soupx/<sample>/filtered_feature_bc_matrix.h5  -- corrected h5 (x16)
- results/phase1_qc/soupx/soupx_summary.csv                   -- rho + diagnostics
- results/phase1_qc/soupx/soupx_rho_barplot.pdf               -- rho bar chart

Config Constants
----------------
- SAMPLE_METADATA, DATA_ROOT, PROJECT_ROOT
- SOUPX_DIR, SOUPX_RESULTS_DIR, SOUPX_METHOD, SOUPX_ROUND_TO_INT
"""
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from config import (
    SAMPLE_METADATA, DATA_ROOT, PROJECT_ROOT,
    SOUPX_DIR, SOUPX_PARALLEL_DIR, SOUPX_RESULTS_DIR,
    SOUPX_METHOD, SOUPX_ROUND_TO_INT,
    SOUPX_ANCHOR_GENES,
    ensure_dirs,
)

import gc
import io
import tarfile
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import h5py
import soupx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

def extract_clusters(sample_dir):
    """Extract CellRanger graph-based clusters from analysis.tar.gz.

    Returns a DataFrame with columns ['Barcode', 'Cluster'].
    """
    tar_path = os.path.join(sample_dir, "analysis.tar.gz")
    target = "clustering/gene_expression_graphclust/clusters.csv"
    with tarfile.open(tar_path, "r:gz") as tar:
        member = tar.getmember(target)
        f = tar.extractfile(member)
        clusters_df = pd.read_csv(io.BytesIO(f.read()))
    return clusters_df


def anchor_detection_rates(adata_rna):
    """Fraction of cells with non-zero counts for each SOUPX_ANCHOR_GENES entry.

    Returns an empty dict when SOUPX_ANCHOR_GENES is empty.
    """
    rates = {}
    for gene in SOUPX_ANCHOR_GENES:
        if gene in adata_rna.var_names:
            col = adata_rna[:, gene].X
            if sp.issparse(col):
                col = col.toarray()
            rates[gene] = float((col > 0).mean())
        else:
            rates[gene] = np.nan
    return rates


def write_corrected_h5(original_h5_path, output_h5_path, corrected_rna_csc,
                       n_gex_features):
    """Write a CellRanger-format h5 with corrected RNA + original ADT.

    Parameters
    ----------
    original_h5_path : str
        Path to the original filtered_feature_bc_matrix.h5.
    output_h5_path : str
        Path for the output corrected h5.
    corrected_rna_csc : scipy.sparse.csc_matrix
        Corrected RNA counts, shape (n_gex_features, n_cells), CSC format.
    n_gex_features : int
        Number of Gene Expression features (first N rows in the h5 matrix).
    """
    os.makedirs(os.path.dirname(output_h5_path), exist_ok=True)

    with h5py.File(original_h5_path, "r") as src, \
         h5py.File(output_h5_path, "w") as dst:
        # Copy the entire file structure first
        for key in src.keys():
            src.copy(key, dst)

        # Extract original ADT rows from the CSC matrix
        orig_data = np.array(dst["matrix/data"])
        orig_indices = np.array(dst["matrix/indices"])
        orig_indptr = np.array(dst["matrix/indptr"])
        orig_shape = np.array(dst["matrix/shape"])
        n_features, n_cells = orig_shape

        # Reconstruct original full matrix as CSC to extract ADT rows
        original_csc = sp.csc_matrix(
            (orig_data, orig_indices, orig_indptr),
            shape=(n_features, n_cells),
        )
        adt_csc = original_csc[n_gex_features:, :]  # last rows = ADT

        # Stack corrected RNA (top) + original ADT (bottom)
        combined = sp.vstack([corrected_rna_csc, adt_csc], format="csc")

        # Overwrite the sparse matrix datasets in-place
        del dst["matrix/data"]
        del dst["matrix/indices"]
        del dst["matrix/indptr"]
        del dst["matrix/shape"]

        dst.create_dataset("matrix/data", data=combined.data)
        dst.create_dataset("matrix/indices", data=combined.indices)
        dst.create_dataset("matrix/indptr", data=combined.indptr)
        dst.create_dataset("matrix/shape", data=np.array(combined.shape,
                                                          dtype=np.int32))


def process_library(sample_name, output_dir=None):
    """Run SoupX on a single library. Returns a dict of diagnostics.

    Parameters
    ----------
    sample_name : str
        Key from SAMPLE_METADATA.
    output_dir : str, optional
        Base directory for corrected h5 output. Defaults to SOUPX_DIR.
    """
    if output_dir is None:
        output_dir = SOUPX_DIR
    sample_dir = os.path.join(DATA_ROOT, sample_name)
    raw_h5 = os.path.join(sample_dir, "raw_feature_bc_matrix.h5")
    filtered_h5 = os.path.join(sample_dir, "filtered_feature_bc_matrix.h5")

    print(f"\n{'='*60}")
    print(f"Processing: {sample_name}")
    print(f"{'='*60}")

    # 1. Extract CellRanger clusters
    clusters_df = extract_clusters(sample_dir)
    print(f"  Clusters: {clusters_df['Cluster'].nunique()} groups, "
          f"{len(clusters_df)} barcodes")

    # 2. Load RNA-only matrices (gex_only=True excludes ADT)
    print("  Loading raw matrix...")
    adata_raw = sc.read_10x_h5(raw_h5, gex_only=True)
    adata_raw.var_names_make_unique()
    print(f"  Raw: {adata_raw.shape[0]} droplets x {adata_raw.shape[1]} genes")

    print("  Loading filtered matrix...")
    adata_filt = sc.read_10x_h5(filtered_h5, gex_only=True)
    adata_filt.var_names_make_unique()
    print(f"  Filtered: {adata_filt.shape[0]} cells x {adata_filt.shape[1]} genes")

    # Anchor-gene detection rates BEFORE correction
    rates_before = anchor_detection_rates(adata_filt)

    # 3. Align clusters to filtered barcodes
    clusters_df = clusters_df.set_index("Barcode")
    common_barcodes = adata_filt.obs_names.intersection(clusters_df.index)
    if len(common_barcodes) < adata_filt.n_obs:
        print(f"  WARNING: {adata_filt.n_obs - len(common_barcodes)} filtered "
              f"barcodes missing from clusters; subsetting")
        adata_filt = adata_filt[common_barcodes].copy()
    clusters_array = clusters_df.loc[adata_filt.obs_names, "Cluster"].values

    # 4. Ensure gene sets match between raw and filtered
    # SoupX requires identical gene ordering
    common_genes = adata_raw.var_names.intersection(adata_filt.var_names)
    if len(common_genes) < adata_filt.n_vars:
        print(f"  Subsetting to {len(common_genes)} common genes")
        adata_raw = adata_raw[:, common_genes].copy()
        adata_filt = adata_filt[:, common_genes].copy()

    # 5. Run SoupX
    # SoupX expects genes x cells/droplets in CSR format
    tod = adata_raw.X.T.tocsr()   # genes x droplets
    toc = adata_filt.X.T.tocsr()  # genes x cells
    del adata_raw; gc.collect()

    print("  Creating SoupChannel...")
    sc_obj = soupx.SoupChannel(tod=tod, toc=toc)
    del tod; gc.collect()

    print("  Setting clusters...")
    sc_obj.setClusters(clusters_array)

    print("  Estimating contamination fraction (autoEstCont)...")
    sc_obj = soupx.autoEstCont(sc_obj)
    rho = sc_obj.fit["rhoEst"]
    print(f"  Estimated rho = {rho:.4f}")

    print(f"  Adjusting counts (method={SOUPX_METHOD}, "
          f"roundToInt={SOUPX_ROUND_TO_INT})...")
    corrected = soupx.adjustCounts(
        sc_obj,
        method=SOUPX_METHOD,
        roundToInt=SOUPX_ROUND_TO_INT,
    )
    # corrected is CSR, genes x cells
    del sc_obj; gc.collect()

    # 6. Compute anchor-gene detection rates AFTER correction
    # Build a temporary AnnData from corrected (transpose to cells x genes)
    corrected_t = corrected.T.tocsr()  # cells x genes
    adata_corrected = sc.AnnData(
        X=corrected_t,
        obs=adata_filt.obs.copy(),
        var=adata_filt.var.copy(),
    )
    rates_after = anchor_detection_rates(adata_corrected)
    del adata_corrected, corrected_t; gc.collect()

    # 7. Write corrected h5 with RNA corrected + original ADT unchanged
    # Convert corrected (genes x cells, CSR) to CSC for h5 storage
    corrected_csc = corrected.tocsc()
    del corrected; gc.collect()

    # Determine n_gex_features from the original h5
    with h5py.File(filtered_h5, "r") as f:
        feature_types = [x.decode() for x in f["matrix/features/feature_type"][:]]
    n_gex = sum(1 for ft in feature_types if ft == "Gene Expression")

    # SoupX was run on gex_only subset — need to map back to full gene set
    # The filtered h5 has n_gex GEX features; SoupX corrected len(common_genes)
    # If they match, direct use; otherwise we need to reconstruct
    if corrected_csc.shape[0] == n_gex:
        corrected_rna_csc = corrected_csc
    else:
        # Build a full-size GEX matrix, filling uncorrected genes with original
        print(f"  Mapping {corrected_csc.shape[0]} corrected genes back to "
              f"{n_gex} GEX features...")
        adata_filt_full = sc.read_10x_h5(filtered_h5, gex_only=True)
        adata_filt_full.var_names_make_unique()
        full_csc = adata_filt_full.X.T.tocsc()  # genes x cells
        gene_idx = [list(adata_filt_full.var_names).index(g)
                    for g in adata_filt.var_names]
        full_lil = full_csc.tolil()
        corrected_lil = corrected_csc.tolil()
        for i_new, i_full in enumerate(gene_idx):
            full_lil[i_full] = corrected_lil[i_new]
        corrected_rna_csc = full_lil.tocsc()
        del adata_filt_full, full_csc, full_lil, corrected_lil; gc.collect()

    output_h5 = os.path.join(output_dir, sample_name,
                             "filtered_feature_bc_matrix.h5")
    write_corrected_h5(filtered_h5, output_h5, corrected_rna_csc, n_gex)
    print(f"  Wrote: {output_h5}")

    del corrected_csc, corrected_rna_csc, adata_filt; gc.collect()

    return {
        "sample": sample_name,
        "rho": rho,
        "n_cells": len(clusters_array),
        "n_clusters": int(clusters_df["Cluster"].nunique()),
        **{f"anchor_{g}_before": rates_before.get(g, np.nan)
           for g in SOUPX_ANCHOR_GENES},
        **{f"anchor_{g}_after": rates_after.get(g, np.nan)
           for g in SOUPX_ANCHOR_GENES},
    }


def _save_per_library_csv(diag):
    """Save a single-library diagnostics dict as a CSV row."""
    os.makedirs(SOUPX_RESULTS_DIR, exist_ok=True)
    out = os.path.join(SOUPX_RESULTS_DIR,
                       f"soupx_diag_{diag['sample']}.csv")
    pd.DataFrame([diag]).to_csv(out, index=False)
    print(f"  Per-library diagnostics: {out}")


def _generate_summary(parallel=False):
    """Aggregate per-library CSVs into summary CSV + rho barplot.

    Parameters
    ----------
    parallel : bool
        If True, write to ``soupx_summary_parallel.csv`` and
        ``soupx_rho_barplot_parallel.pdf`` to avoid collisions with
        a concurrent sequential run.
    """
    suffix = "_parallel" if parallel else ""
    csvs = sorted(
        [os.path.join(SOUPX_RESULTS_DIR, f)
         for f in os.listdir(SOUPX_RESULTS_DIR)
         if f.startswith("soupx_diag_") and f.endswith(".csv")]
    )
    if not csvs:
        print("ERROR: No per-library diagnostics CSVs found.")
        sys.exit(1)

    summary_df = pd.concat([pd.read_csv(c) for c in csvs], ignore_index=True)

    # Ensure consistent sample ordering (match SAMPLE_METADATA key order)
    sample_order = list(SAMPLE_METADATA.keys())
    summary_df["_sort"] = summary_df["sample"].map(
        {s: i for i, s in enumerate(sample_order)}
    )
    summary_df = summary_df.sort_values("_sort").drop(columns="_sort")

    summary_csv = os.path.join(SOUPX_RESULTS_DIR,
                               f"soupx_summary{suffix}.csv")
    summary_df.to_csv(summary_csv, index=False)
    print(f"Summary saved: {summary_csv} ({len(summary_df)} libraries)")

    # Rho bar plot
    fig, ax = plt.subplots(figsize=(12, 5))
    short_names = [s.split("_sample_")[0].replace("donor_", "")
                   for s in summary_df["sample"]]
    bars = ax.bar(range(len(summary_df)), summary_df["rho"], color="#4292c6")
    ax.set_xticks(range(len(summary_df)))
    ax.set_xticklabels(short_names, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Estimated rho (contamination fraction)")
    ax.set_title("SoupX Contamination Estimate per Library")
    for bar, rho in zip(bars, summary_df["rho"]):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.002,
                f"{rho:.3f}", ha="center", va="bottom", fontsize=7)
    plt.tight_layout()
    rho_plot = os.path.join(SOUPX_RESULTS_DIR,
                            f"soupx_rho_barplot{suffix}.pdf")
    plt.savefig(rho_plot, bbox_inches="tight")
    plt.close()
    print(f"Rho plot saved: {rho_plot}")

    # Anchor-gene summary
    if SOUPX_ANCHOR_GENES:
        print("\nAnchor Gene Detection Rates (before -> after):")
        for gene in SOUPX_ANCHOR_GENES:
            before_col = f"anchor_{gene}_before"
            after_col = f"anchor_{gene}_after"
            if before_col in summary_df.columns:
                mean_before = summary_df[before_col].mean()
                mean_after = summary_df[after_col].mean()
                print(f"  {gene}: {mean_before:.3f} -> {mean_after:.3f}")


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="Step 0b: SoupX ambient RNA decontamination")
    parser.add_argument("--sample-index", type=int, default=None,
                        help="1-based sample index (for SGE array jobs). "
                             "Processes a single library.")
    parser.add_argument("--aggregate", action="store_true",
                        help="Aggregate per-library CSVs into summary + plot, "
                             "then validate.")
    args = parser.parse_args()

    ensure_dirs()
    os.makedirs(SOUPX_RESULTS_DIR, exist_ok=True)

    if args.aggregate:
        # Aggregation mode: collect per-library CSVs, build summary, validate
        # Uses "_parallel" suffix to avoid overwriting sequential-mode outputs
        print("SoupX Aggregation Mode (parallel)")
        _generate_summary(parallel=True)
        print(f"\nAggregation complete.")
        return

    all_samples = list(SAMPLE_METADATA.keys())

    if args.sample_index is not None:
        # Single-library mode (SGE array job)
        idx = args.sample_index - 1  # convert 1-based to 0-based
        if idx < 0 or idx >= len(all_samples):
            print(f"ERROR: --sample-index {args.sample_index} out of range "
                  f"(1-{len(all_samples)})")
            sys.exit(1)
        sample_name = all_samples[idx]
        print(f"SoupX Single-Library Mode: index={args.sample_index}, "
              f"sample={sample_name}")
        print(f"  Method: {SOUPX_METHOD}, Round to int: {SOUPX_ROUND_TO_INT}")
        print(f"  Output dir: {SOUPX_PARALLEL_DIR}")
        diag = process_library(sample_name, output_dir=SOUPX_PARALLEL_DIR)
        _save_per_library_csv(diag)
        print(f"\nLibrary {args.sample_index}/{len(SAMPLE_METADATA)} complete.")
        return

    # Sequential mode (all libraries)
    print("SoupX Ambient RNA Decontamination (sequential mode)")
    print(f"  Method: {SOUPX_METHOD}")
    print(f"  Round to int: {SOUPX_ROUND_TO_INT}")
    print(f"  Libraries: {len(SAMPLE_METADATA)}")
    print(f"  Output dir: {SOUPX_DIR}")

    for sample_name in all_samples:
        diag = process_library(sample_name)
        _save_per_library_csv(diag)
        gc.collect()

    _generate_summary()
    print(f"\nStep 0b complete. {len(all_samples)} libraries processed.")


def validate(parallel=False):
    """Validate Step 0b outputs.

    Parameters
    ----------
    parallel : bool
        If True, check ``SOUPX_PARALLEL_DIR`` for h5 files and
        ``soupx_summary_parallel.csv`` for diagnostics.
    """
    h5_dir = SOUPX_PARALLEL_DIR if parallel else SOUPX_DIR
    suffix = "_parallel" if parallel else ""
    mode_label = "parallel" if parallel else "sequential"
    results = {}

    # Check: every library has a corrected h5 file
    n_expected = len(SAMPLE_METADATA)
    missing = []
    for sample_name in SAMPLE_METADATA:
        h5_path = os.path.join(h5_dir, sample_name,
                               "filtered_feature_bc_matrix.h5")
        if not os.path.exists(h5_path):
            missing.append(sample_name)
    n_found = n_expected - len(missing)
    results["corrected_h5_files"] = {
        "status": "PASS" if len(missing) == 0 else "FAIL",
        "detail": f"{n_found}/{n_expected} corrected h5 files in {h5_dir}",
        "expected": f"{n_expected} corrected h5 files",
        "actual": f"{n_found} found, {len(missing)} missing",
    }

    # Check: readable by scanpy, feature count matches original
    if len(missing) == 0:
        sample0 = list(SAMPLE_METADATA.keys())[0]
        h5_path = os.path.join(h5_dir, sample0,
                               "filtered_feature_bc_matrix.h5")
        try:
            adata = sc.read_10x_h5(h5_path, gex_only=False)
            orig_h5 = os.path.join(DATA_ROOT, sample0,
                                   "filtered_feature_bc_matrix.h5")
            adata_orig = sc.read_10x_h5(orig_h5, gex_only=False)

            features_match = adata.n_vars == adata_orig.n_vars
            results["feature_count"] = {
                "status": "PASS" if features_match else "FAIL",
                "detail": f"{sample0}: corrected={adata.n_vars}, "
                          f"original={adata_orig.n_vars}",
                "expected": "feature count matches original",
                "actual": f"corrected={adata.n_vars}, original={adata_orig.n_vars}",
            }

            cells_match = adata.n_obs == adata_orig.n_obs
            results["cell_count_match"] = {
                "status": "PASS" if cells_match else "FAIL",
                "detail": f"{sample0}: corrected={adata.n_obs}, "
                          f"original={adata_orig.n_obs}",
                "expected": "cell counts match original",
                "actual": f"corrected={adata.n_obs}, original={adata_orig.n_obs}",
            }
            del adata, adata_orig
        except Exception as e:
            results["feature_count"] = {
                "status": "FAIL",
                "detail": f"Failed to read corrected h5: {e}",
                "expected": "readable by scanpy",
                "actual": str(e),
            }

    # Check: summary CSV
    summary_name = f"soupx_summary{suffix}.csv"
    summary_path = os.path.join(SOUPX_RESULTS_DIR, summary_name)
    if os.path.exists(summary_path):
        summary_df = pd.read_csv(summary_path)
        n_rows = len(summary_df)
        results["summary_csv"] = {
            "status": "PASS" if n_rows == n_expected else "FAIL",
            "detail": f"{summary_name} has {n_rows} rows",
            "expected": f"{n_expected} rows",
            "actual": str(n_rows),
        }

        # Check rho range
        rho_min = summary_df["rho"].min()
        rho_max = summary_df["rho"].max()
        rho_mean = summary_df["rho"].mean()
        rho_ok = (0 <= rho_min) and (rho_max <= 1) and (0.01 <= rho_mean <= 0.5)
        results["rho_range"] = {
            "status": "PASS" if rho_ok else "WARN",
            "detail": f"rho: min={rho_min:.4f}, max={rho_max:.4f}, "
                      f"mean={rho_mean:.4f}",
            "expected": "all in [0,1], mean in [0.01, 0.5]",
            "actual": f"min={rho_min:.4f}, max={rho_max:.4f}, "
                      f"mean={rho_mean:.4f}",
        }

        # Check anchor-gene reduction (skipped when SOUPX_ANCHOR_GENES is empty)
        reductions = []
        for gene in SOUPX_ANCHOR_GENES:
            bc = f"anchor_{gene}_before"
            ac = f"anchor_{gene}_after"
            if bc in summary_df.columns and ac in summary_df.columns:
                mean_b = summary_df[bc].mean()
                mean_a = summary_df[ac].mean()
                if not np.isnan(mean_b) and mean_b > 0:
                    reductions.append(mean_a < mean_b)
        if reductions:
            all_reduced = all(reductions)
            results["anchor_gene_reduction"] = {
                "status": "PASS" if all_reduced else "WARN",
                "detail": f"{sum(reductions)}/{len(reductions)} anchor genes "
                          f"show reduced detection",
                "expected": "all anchor genes reduced",
                "actual": f"{sum(reductions)}/{len(reductions)} reduced",
            }
    else:
        results["summary_csv"] = {
            "status": "FAIL",
            "detail": f"{summary_name} not found",
            "expected": "file exists",
            "actual": "MISSING",
        }

    # Print results
    print(f"\n{'='*60}")
    print(f"Validation: Step 0b - SoupX ({mode_label} mode)")
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
    main()
    # Validate only in aggregate or sequential mode (not single-library)
    if "--aggregate" in sys.argv:
        results = validate(parallel=True)
    elif "--sample-index" not in sys.argv:
        results = validate(parallel=False)
    else:
        results = None
    if results and any(r["status"] == "FAIL" for r in results.values()):
        print("\nValidation FAILED.")
        sys.exit(1)

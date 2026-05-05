#!/usr/bin/env python
"""Step 1: Load CITE-seq Libraries and QC Filtering
===================================================

Objective
---------
Load all CellRanger filtered matrices, separate RNA and ADT modalities,
apply per-library adaptive QC thresholds (MAD-based), detect doublets with
Scrublet, and save raw count layers for downstream normalization.

Pipeline Position
-----------------
Phase:    1.1 (Tier 1)
Depends:  Step 0b (SoupX decontamination) — uses corrected h5 when available
Feeds:    Steps 2 (dsb), 3 (CLR), 4 (RNA norm), 5a (ThresholdR)
Compute:  CPU, 16-32 GB RAM, 10-20 min

Inputs
------
- 10X_cloud/outs/<sample>/filtered_feature_bc_matrix.h5 -- CellRanger output, one per library

Outputs
-------
- data/processed/adata_qc.h5ad                          -- QC-passed RNA with raw counts in .X and .layers['counts']
- data/processed/adata_adt_raw.h5ad                      -- Raw ADT counts synced to QC-passing barcodes
- data/processed/adt_counts_for_thresholdr.csv           -- ADT counts with canonical column names for R/ThresholdR
- results/phase1_qc/qc_summary.csv                      -- Per-library cell retention stats

Config Constants
----------------
- SAMPLE_METADATA (dict)                -- library_id → donor/age/timepoint map
- DATA_ROOT (path)                      -- CellRanger outs directory
- PROJECT_ROOT (path)                   -- Top-level project directory
- QC_MAD_N_MADS (5)                     -- MAD multiplier for adaptive filtering
- QC_MIN_GENES_ABSOLUTE (100)           -- Hard floor for gene count filter
- QC_MIN_UMI_ABSOLUTE (200)             -- Hard floor for UMI count filter
- QC_MAX_MT_PCT (20)                    -- Maximum mitochondrial fraction
- QC_MIN_CELLS_PER_GENE (3)            -- Minimum cells expressing a gene
- SCRUBLET_EXPECTED_DOUBLET_RATE (0.06) -- Prior doublet rate for Scrublet
- SAMPLE_GROUP_FIELD (optional)         -- SAMPLE_METADATA key for coloring/ordering
"""
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from config import (
    SAMPLE_METADATA, DATA_ROOT, PROJECT_ROOT, FIGURES_DIR, PROCESSED_DIR,
    RESULTS_PHASE1_QC, SOUPX_DIR, SOUPX_PARALLEL_DIR, DECONTX_DIR,
    QC_MAD_N_MADS, QC_MIN_GENES_ABSOLUTE, QC_MIN_UMI_ABSOLUTE,
    QC_MAX_MT_PCT, QC_MIN_CELLS_PER_GENE, SCRUBLET_EXPECTED_DOUBLET_RATE,
    SAMPLE_GROUP_FIELD,
    ADT_RENAME_MAP, branched_paths,
)

import argparse
import scanpy as sc
import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)


# ── 1a: Sample metadata dictionary ────────────────────────────
sample_metadata = SAMPLE_METADATA


def _load_one_library(sample_name, meta, branch=None):
    """Load a single library, split into RNA and ADT, return both.

    When branch is set, force the corresponding decontamination source
    instead of the default priority cascade.
    """
    if branch == "decontx":
        h5_path = os.path.join(DECONTX_DIR, sample_name, "filtered_feature_bc_matrix.h5")
        source = "DecontX-corrected (forced)"
    elif branch == "soupx":
        h5_path = os.path.join(SOUPX_PARALLEL_DIR, sample_name, "filtered_feature_bc_matrix.h5")
        source = "SoupX-corrected (forced)"
    else:
        # Default priority cascade
        decontx_h5 = os.path.join(DECONTX_DIR, sample_name, "filtered_feature_bc_matrix.h5")
        soupx_h5 = os.path.join(SOUPX_DIR, sample_name, "filtered_feature_bc_matrix.h5")
        soupx_par_h5 = os.path.join(SOUPX_PARALLEL_DIR, sample_name, "filtered_feature_bc_matrix.h5")
        if os.path.exists(decontx_h5):
            h5_path = decontx_h5
            source = "DecontX-corrected"
        elif os.path.exists(soupx_h5):
            h5_path = soupx_h5
            source = "SoupX-corrected"
        elif os.path.exists(soupx_par_h5):
            h5_path = soupx_par_h5
            source = "SoupX-corrected, parallel"
        else:
            h5_path = os.path.join(DATA_ROOT, sample_name, "filtered_feature_bc_matrix.h5")
            source = "WARNING: no decontaminated output found, using original"
    adata = sc.read_10x_h5(h5_path, gex_only=False)
    adata.obs["sample_id"] = sample_name
    for k, v in meta.items():
        adata.obs[k] = str(v) if v is not None else "none"
    adata.var_names_make_unique()
    is_adt = adata.var["feature_types"].isin(["Antibody Capture"])
    rna = adata[:, ~is_adt].copy()
    adt = adata[:, is_adt].copy()
    del adata
    return sample_name, source, rna, adt


def _scrublet_one_library(sample, sub):
    """Run Scrublet on one library, return (sample, predicted_doublet, doublet_score) or None on failure."""
    try:
        sc.pp.scrublet(sub, expected_doublet_rate=SCRUBLET_EXPECTED_DOUBLET_RATE)
        return sample, sub.obs["predicted_doublet"].values, sub.obs["doublet_score"].values
    except Exception as e:
        print(f"  WARNING: Scrublet failed for {sample}: {e}")
        return None


def _mad_filter(series, n_mads=QC_MAD_N_MADS):
    median = np.median(series)
    mad = np.median(np.abs(series - median))
    return median - n_mads * mad, median + n_mads * mad


# ── Per-library task (for SGE array parallelism) ────────────────
STEP01_TEMP_DIR = os.path.join(PROCESSED_DIR, "step01_temp")


def _ordered_samples_and_palette():
    """Return (ordered_sample_ids, {sample_id: color}).

    When SAMPLE_GROUP_FIELD is set to a key present in SAMPLE_METADATA
    entries, samples are sorted by that field's value (stable within group)
    and colored from matplotlib tab10 per unique group value. When unset,
    samples appear in SAMPLE_METADATA insertion order with a single color.
    """
    all_ids = list(sample_metadata.keys())
    if SAMPLE_GROUP_FIELD:
        def _key(sid):
            return (str(sample_metadata.get(sid, {}).get(SAMPLE_GROUP_FIELD, "")), sid)
        ordered = sorted(all_ids, key=_key)
        groups = []
        seen = set()
        for sid in ordered:
            g = str(sample_metadata[sid].get(SAMPLE_GROUP_FIELD, ""))
            if g not in seen:
                seen.add(g)
                groups.append(g)
        cmap = plt.get_cmap("tab10")
        group_colors = {g: cmap(i % 10) for i, g in enumerate(groups)}
        palette = {
            sid: group_colors[str(sample_metadata[sid].get(SAMPLE_GROUP_FIELD, ""))]
            for sid in ordered
        }
    else:
        ordered = list(all_ids)
        palette = {sid: "#4292c6" for sid in ordered}
    return ordered, palette


def _qc_figure_title():
    base = "Per-Library QC Metrics"
    if SAMPLE_GROUP_FIELD:
        return f"{base} (colored by {SAMPLE_GROUP_FIELD})"
    return base


def _plot_qc_violins(adata_rna):
    """Render the per-library QC violin triptych to FIGURES_DIR."""
    ordered, palette = _ordered_samples_and_palette()
    adata_rna.obs["sample_id"] = pd.Categorical(
        adata_rna.obs["sample_id"], categories=ordered, ordered=True
    )
    metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    titles = ["Genes Detected", "Total UMI Counts", "% Mitochondrial"]
    fig, axes = plt.subplots(3, 1, figsize=(14, 12))
    for ax, key, title in zip(axes, metrics, titles):
        sc.pl.violin(adata_rna, keys=[key], groupby="sample_id",
                     palette=palette, rotation=90, ax=ax, show=False)
        ax.set_title(title, fontsize=12, fontweight="bold")
        ax.set_xlabel("")
        ax.tick_params(axis="x", labelsize=8)
    fig.suptitle(_qc_figure_title(), fontsize=13, fontweight="bold", y=1.01)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "phase1_qc_violin_metrics.png"),
                dpi=150, bbox_inches="tight")
    plt.close()


def main_task(task_index, branch=None):
    """Process a single library: load, QC, Scrublet, save intermediate."""
    import gc
    sample_names = list(sample_metadata.keys())
    sample_name = sample_names[task_index]
    meta = sample_metadata[sample_name]
    print(f"=== Task {task_index}: {sample_name} ===")

    _, source, rna, adt = _load_one_library(sample_name, meta, branch)
    print(f"  Loaded ({source}): {rna.n_obs} cells, {rna.n_vars} genes")

    # QC metrics
    rna.var["mt"] = rna.var_names.str.startswith("MT-")
    rna.var["ribo"] = rna.var_names.str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(
        rna, qc_vars=["mt", "ribo"], percent_top=None, log1p=True, inplace=True
    )

    # MAD-based filtering
    lo_genes, hi_genes = _mad_filter(rna.obs["n_genes_by_counts"])
    lo_genes = max(lo_genes, QC_MIN_GENES_ABSOLUTE)
    lo_umi, hi_umi = _mad_filter(rna.obs["total_counts"])
    lo_umi = max(lo_umi, QC_MIN_UMI_ABSOLUTE)
    cell_mask = (
        (rna.obs["n_genes_by_counts"] >= lo_genes)
        & (rna.obs["n_genes_by_counts"] <= hi_genes)
        & (rna.obs["total_counts"] >= lo_umi)
        & (rna.obs["total_counts"] <= hi_umi)
        & (rna.obs["pct_counts_mt"] < QC_MAX_MT_PCT)
    )
    n_total = rna.n_obs
    rna = rna[cell_mask].copy()
    n_pass = rna.n_obs
    print(f"  QC: {n_pass}/{n_total} cells pass (genes>={lo_genes:.0f})")

    # Scrublet
    result = _scrublet_one_library(sample_name, rna)
    if result is not None:
        _, pred, scores = result
        rna.obs["predicted_doublet"] = pred
        rna.obs["doublet_score"] = scores
    else:
        rna.obs["predicted_doublet"] = False
        rna.obs["doublet_score"] = 0.0
    n_doublets = int(rna.obs["predicted_doublet"].sum())
    print(f"  Scrublet: {n_doublets} doublets ({n_doublets/n_pass*100:.1f}%)")

    # Sync ADT to QC-passing cells
    adt = adt[rna.obs_names].copy()

    # Save intermediates
    os.makedirs(STEP01_TEMP_DIR, exist_ok=True)
    suffix = f"_{branch}" if branch else ""
    rna.write_h5ad(os.path.join(STEP01_TEMP_DIR, f"{sample_name}_rna{suffix}.h5ad"))
    adt.write_h5ad(os.path.join(STEP01_TEMP_DIR, f"{sample_name}_adt{suffix}.h5ad"))
    qc_row = {"sample": sample_name, "pass": n_pass, "total": n_total,
              "pct": n_pass / n_total * 100, "min_genes": lo_genes,
              "doublets": n_doublets, "pass_post_doublet": n_pass - n_doublets}
    pd.DataFrame([qc_row]).to_csv(
        os.path.join(STEP01_TEMP_DIR, f"{sample_name}_qc{suffix}.csv"), index=False
    )
    # Export QC-passing barcodes for scDblFinder (step01b)
    pd.Series(rna.obs_names).to_csv(
        os.path.join(STEP01_TEMP_DIR, f"{sample_name}_keep_barcodes{suffix}.csv"),
        index=False, header=False
    )
    print(f"  Saved to {STEP01_TEMP_DIR}")


# ── Aggregation (runs after all tasks complete) ────────────────
def main_aggregate(branch=None):
    """Combine per-library intermediates into final outputs."""
    import gc
    paths = branched_paths(branch)
    suffix = f"_{branch}" if branch else ""

    print("=== Aggregating per-library results ===")
    rna_list = {}
    adt_list = {}
    qc_frames = []
    for sample_name in sample_metadata:
        rna_path = os.path.join(STEP01_TEMP_DIR, f"{sample_name}_rna{suffix}.h5ad")
        adt_path = os.path.join(STEP01_TEMP_DIR, f"{sample_name}_adt{suffix}.h5ad")
        qc_path = os.path.join(STEP01_TEMP_DIR, f"{sample_name}_qc{suffix}.csv")
        if not os.path.exists(rna_path):
            print(f"  MISSING: {rna_path} — run task for this library first")
            sys.exit(1)
        rna_list[sample_name] = sc.read_h5ad(rna_path)
        adt_list[sample_name] = sc.read_h5ad(adt_path)
        qc_frames.append(pd.read_csv(qc_path))
        print(f"  Loaded {sample_name}: {rna_list[sample_name].n_obs} cells")

    # Concatenate
    adata_rna = sc.concat(rna_list, label="sample_id", index_unique="-")
    del rna_list; gc.collect()
    adata_adt = sc.concat(adt_list, label="sample_id", index_unique="-")
    del adt_list; gc.collect()
    qc_summary = pd.concat(qc_frames, ignore_index=True)

    print(f"Total: {adata_rna.n_obs} cells, {adata_rna.n_vars} genes")
    print(f"ADT: {adata_adt.n_obs} cells, {adata_adt.n_vars} markers -> {list(adata_adt.var_names)}")

    # Remove doublets
    n_doublets = int(adata_rna.obs["predicted_doublet"].sum())
    print(f"Doublets detected: {n_doublets} ({n_doublets/len(adata_rna)*100:.1f}%)")
    adata_rna = adata_rna[~adata_rna.obs["predicted_doublet"]].copy()

    # Filter genes across all libraries
    sc.pp.filter_genes(adata_rna, min_cells=QC_MIN_CELLS_PER_GENE)

    # Sync ADT to final RNA barcodes
    adata_adt = adata_adt[adata_rna.obs_names].copy()

    # Save raw counts layer
    adata_rna.layers["counts"] = adata_rna.X.copy()

    # QC summary CSV
    os.makedirs(RESULTS_PHASE1_QC, exist_ok=True)
    qc_summary.to_csv(os.path.join(RESULTS_PHASE1_QC, "qc_summary.csv"), index=False)

    _plot_qc_violins(adata_rna)

    adata_rna.write_h5ad(paths["ADATA_QC"])
    adata_adt.write_h5ad(paths["ADATA_ADT_RAW"])

    # Export raw ADT counts for ThresholdR
    if not branch:
        adt_counts_dense = adata_adt.X
        if sp.issparse(adt_counts_dense):
            adt_counts_dense = adt_counts_dense.toarray()
        adt_export = pd.DataFrame(
            adt_counts_dense,
            index=adata_adt.obs_names,
            columns=adata_adt.var_names
        )
        adt_export.rename(columns=ADT_RENAME_MAP, inplace=True)
        adt_export.to_csv(os.path.join(PROCESSED_DIR, "adt_counts_for_thresholdr.csv"))
        print(f"Exported ADT counts for ThresholdR: {adt_export.shape}")
    else:
        print(f"Branch mode: skipping ThresholdR CSV export (ADT unchanged across branches)")

    print(f"\nStep 1 complete: {adata_rna.n_obs} cells, {adata_rna.n_vars} genes")
    print(f"ADT matrix: {adata_adt.n_obs} cells, {adata_adt.n_vars} markers")


# ── Original monolithic mode (kept for backward compatibility) ──
def main(branch=None):
    paths = branched_paths(branch)
    if branch:
        print(f"=== Branch mode: '{branch}' ===")
        print(f"  Output QC: {paths['ADATA_QC']}")
        print(f"  Output ADT: {paths['ADATA_ADT_RAW']}")

    # ── 1b: Load all GEX + ADT libraries (parallel) ─────────────
    import gc
    from concurrent.futures import ThreadPoolExecutor

    rna_list = {}
    adt_list = {}
    with ThreadPoolExecutor(max_workers=8) as pool:
        futures = {pool.submit(_load_one_library, name, meta, branch): name
                   for name, meta in sample_metadata.items()}
        for future in futures:
            sample_name, source, rna, adt = future.result()
            print(f"Loading {sample_name} ({source})...")
            rna_list[sample_name] = rna
            adt_list[sample_name] = adt
    gc.collect()

    # ── 1c: Concatenate RNA and ADT separately ──────────────────
    adata_rna = sc.concat(rna_list, label="sample_id", index_unique="-")
    del rna_list; gc.collect()
    adata_adt = sc.concat(adt_list, label="sample_id", index_unique="-")
    del adt_list; gc.collect()
    print(f"Loaded {adata_rna.n_obs} total cells across {len(sample_metadata)} libraries")

    print(f"RNA features: {adata_rna.n_vars}")
    print(f"ADT features: {adata_adt.n_vars} -> {list(adata_adt.var_names)}")

    # ── 1d: RNA quality control metrics ─────────────────────────
    adata_rna.var["mt"] = adata_rna.var_names.str.startswith("MT-")
    adata_rna.var["ribo"] = adata_rna.var_names.str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(
        adata_rna, qc_vars=["mt", "ribo"], percent_top=None, log1p=True, inplace=True
    )

    # ── 1e: Adaptive MAD-based filtering per library ────────────
    # Per-library thresholds accommodate the depth range across libraries.
    keep_cells = []
    qc_summary = []
    for sample in adata_rna.obs["sample_id"].unique():
        mask = adata_rna.obs["sample_id"] == sample
        sub = adata_rna[mask]
        lo_genes, hi_genes = _mad_filter(sub.obs["n_genes_by_counts"])
        lo_genes = max(lo_genes, QC_MIN_GENES_ABSOLUTE)
        lo_umi, hi_umi = _mad_filter(sub.obs["total_counts"])
        lo_umi = max(lo_umi, QC_MIN_UMI_ABSOLUTE)
        cell_mask = (
            (sub.obs["n_genes_by_counts"] >= lo_genes)
            & (sub.obs["n_genes_by_counts"] <= hi_genes)
            & (sub.obs["total_counts"] >= lo_umi)
            & (sub.obs["total_counts"] <= hi_umi)
            & (sub.obs["pct_counts_mt"] < QC_MAX_MT_PCT)
        )
        keep_cells.extend(sub.obs_names[cell_mask].tolist())
        n_pass = cell_mask.sum()
        n_total = mask.sum()
        qc_summary.append({"sample": sample, "pass": n_pass, "total": n_total,
                            "pct": n_pass / n_total * 100, "min_genes": lo_genes})
        print(f"  {sample}: {n_pass}/{n_total} cells pass QC (genes>={lo_genes:.0f})")

    adata_rna = adata_rna[keep_cells].copy()
    sc.pp.filter_genes(adata_rna, min_cells=QC_MIN_CELLS_PER_GENE)

    # ── 1f: Doublet detection (per library, parallel) ────────────
    # NOTE: Scrublet may fail on very small libraries; failures are
    # caught and skipped so those cells remain in the dataset.
    adata_rna.obs["predicted_doublet"] = False
    adata_rna.obs["doublet_score"] = 0.0
    scrublet_inputs = {}
    for sample in adata_rna.obs["sample_id"].unique():
        mask = adata_rna.obs["sample_id"] == sample
        scrublet_inputs[sample] = (mask, adata_rna[mask].copy())

    with ThreadPoolExecutor(max_workers=8) as pool:
        futures = {pool.submit(_scrublet_one_library, sample, sub): sample
                   for sample, (mask, sub) in scrublet_inputs.items()}
        for future in futures:
            result = future.result()
            if result is not None:
                sample, pred, scores = result
                mask = scrublet_inputs[sample][0]
                adata_rna.obs.loc[mask, "predicted_doublet"] = pred
                adata_rna.obs.loc[mask, "doublet_score"] = scores
    del scrublet_inputs

    n_doublets = adata_rna.obs["predicted_doublet"].sum()
    doublet_pct = n_doublets / len(adata_rna) * 100
    print(f"\nDoublets detected: {n_doublets} ({doublet_pct:.1f}%)")

    # Add per-library post-doublet counts to qc_summary
    doublet_counts = adata_rna.obs.groupby("sample_id")["predicted_doublet"].sum()
    qc_lookup = {row["sample"]: i for i, row in enumerate(qc_summary)}
    for sample, n_dbl in doublet_counts.items():
        idx = qc_lookup[sample]
        qc_summary[idx]["doublets"] = int(n_dbl)
        qc_summary[idx]["pass_post_doublet"] = qc_summary[idx]["pass"] - int(n_dbl)

    adata_rna = adata_rna[~adata_rna.obs["predicted_doublet"]].copy()

    # ── 1g: Sync ADT matrix to QC-passing cells ────────────────
    assert adata_rna.obs_names.isin(adata_adt.obs_names).all(), (
        "Barcode mismatch: some RNA barcodes missing from ADT matrix"
    )
    adata_adt = adata_adt[adata_rna.obs_names].copy()

    # ── 1h: Save raw counts layer ──────────────────────────────
    # Steps 5/6/7 require raw integer counts in layers['counts'].
    adata_rna.layers["counts"] = adata_rna.X.copy()

    # ── 1i: Save QC results and figures ─────────────────────────
    os.makedirs(RESULTS_PHASE1_QC, exist_ok=True)
    pd.DataFrame(qc_summary).to_csv(
        os.path.join(RESULTS_PHASE1_QC, "qc_summary.csv"), index=False
    )

    _plot_qc_violins(adata_rna)

    adata_rna.write_h5ad(paths["ADATA_QC"])
    adata_adt.write_h5ad(paths["ADATA_ADT_RAW"])

    # ── 1j: Export raw ADT counts for ThresholdR ───────────────
    # Step 5a (R) reads this CSV; columns renamed from CellRanger
    # naming convention to canonical short names via ADT_RENAME_MAP.
    # Only export when running without branch (shared output, ADT unchanged).
    if not branch:
        adt_counts_dense = adata_adt.X
        if sp.issparse(adt_counts_dense):
            adt_counts_dense = adt_counts_dense.toarray()
        adt_export = pd.DataFrame(
            adt_counts_dense,
            index=adata_adt.obs_names,
            columns=adata_adt.var_names
        )
        adt_export.rename(columns=ADT_RENAME_MAP, inplace=True)
        adt_export.to_csv(os.path.join(PROCESSED_DIR, "adt_counts_for_thresholdr.csv"))
        print(f"Exported ADT counts for ThresholdR: {adt_export.shape}")
    else:
        print(f"Branch mode: skipping ThresholdR CSV export (ADT unchanged across branches)")

    print(f"\nStep 1 complete: {adata_rna.n_obs} cells, {adata_rna.n_vars} genes")
    print(f"ADT matrix: {adata_adt.n_obs} cells, {adata_adt.n_vars} markers")


def validate(branch=None):
    # In branch mode, main() writes branched filenames (adata_qc_<branch>.h5ad);
    # the unbranched canonical name is produced later by the select_branch
    # module. Use branched_paths() so validate checks the actual outputs of
    # this invocation rather than what a downstream step will eventually
    # promote.
    paths = branched_paths(branch)
    results = {}

    _path = paths["ADATA_QC"]
    _exists = os.path.exists(_path)
    results["adata_qc_exists"] = {
        "status": "PASS" if _exists else "FAIL",
        "detail": f"{_path} {'exists' if _exists else 'MISSING'}",
        "expected": "file exists",
        "actual": str(_exists),
    }

    _path = paths["ADATA_ADT_RAW"]
    _exists = os.path.exists(_path)
    results["adata_adt_raw_exists"] = {
        "status": "PASS" if _exists else "FAIL",
        "detail": f"{_path} {'exists' if _exists else 'MISSING'}",
        "expected": "file exists",
        "actual": str(_exists),
    }

    # adt_counts_for_thresholdr.csv is only written in non-branch mode (the ADT
    # matrix is identical across decontamination branches, so this CSV is
    # shared and emitted once by the unbranched run).
    if not branch:
        _path = os.path.join(PROCESSED_DIR, "adt_counts_for_thresholdr.csv")
        _exists = os.path.exists(_path)
        results["adt_counts_csv_exists"] = {
            "status": "PASS" if _exists else "FAIL",
            "detail": f"{_path} {'exists' if _exists else 'MISSING'}",
            "expected": "file exists",
            "actual": str(_exists),
        }

    _qc_path = paths["ADATA_QC"]
    if os.path.exists(_qc_path):
        import scanpy as sc
        adata = sc.read_h5ad(_qc_path)
        results["adata_qc_shape"] = {
            "status": "PASS",
            "detail": f"{adata.n_obs} cells x {adata.n_vars} genes",
            "expected": "informational",
            "actual": f"{adata.n_obs} cells x {adata.n_vars} genes",
        }
    else:
        results["adata_qc_shape"] = {
            "status": "FAIL",
            "detail": "Cannot check shape: adata_qc.h5ad missing",
            "expected": "file exists",
            "actual": "N/A",
        }

    _path = os.path.join(RESULTS_PHASE1_QC, "qc_summary.csv")
    _exists = os.path.exists(_path)
    results["qc_summary_exists"] = {
        "status": "PASS" if _exists else "FAIL",
        "detail": f"{_path} {'exists' if _exists else 'MISSING'}",
        "expected": "file exists",
        "actual": str(_exists),
    }

    # WARNING: <50% retention in any library may indicate overly aggressive filtering
    if os.path.exists(_path):
        import pandas as pd
        qc_df = pd.read_csv(_path)
        n_samples = len(qc_df)
        expected_samples = len(sample_metadata)
        low_retention = qc_df[qc_df["pct"] < 50]
        samples_ok = n_samples == expected_samples and len(low_retention) == 0
        results["sample_retention"] = {
            "status": "PASS" if samples_ok else "WARN",
            "detail": f"{n_samples} samples, {len(low_retention)} with <50% retention",
            "expected": f"{expected_samples} samples, all >50% retention",
            "actual": f"{n_samples} samples, {len(low_retention)} low retention",
        }

        # Doublet rate check from saved qc_summary
        if "doublets" in qc_df.columns and "pass" in qc_df.columns:
            total_doublets = qc_df["doublets"].sum()
            total_pre_doublet = qc_df["pass"].sum()
            doublet_rate = total_doublets / total_pre_doublet * 100 if total_pre_doublet > 0 else 0
            rate_ok = 3 <= doublet_rate <= 8
            results["doublet_rate_check"] = {
                "status": "PASS" if rate_ok else "WARN",
                "detail": f"{doublet_rate:.1f}% overall doublet rate ({int(total_doublets)} doublets)",
                "expected": "3-8% doublet rate",
                "actual": f"{doublet_rate:.1f}%",
            }
        else:
            results["doublet_rate_check"] = {
                "status": "WARN",
                "detail": "qc_summary.csv missing doublet columns (re-run step 1)",
                "expected": "3-8% doublet rate",
                "actual": "N/A",
            }
    else:
        results["doublet_rate_check"] = {
            "status": "WARN",
            "detail": "Cannot check doublet rate: qc_summary.csv missing",
            "expected": "3-8% doublet rate",
            "actual": "N/A",
        }

    print(f"\n{'='*60}")
    print(f"Validation: Step 1 - Load CITE-seq Libraries and QC Filtering")
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
    parser = argparse.ArgumentParser(description="Step 1: Load CITE-seq Libraries and QC Filtering")
    parser.add_argument("--branch", type=str, default=None,
                        help="Branch name for parallel decontamination comparison (e.g. 'decontx', 'soupx')")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--task", type=int, default=None, metavar="N",
                      help="Process library N (0-based) only; save to step01_temp/. For SGE array jobs.")
    mode.add_argument("--aggregate", action="store_true",
                      help="Combine per-library results from step01_temp/ into final outputs.")
    args = parser.parse_args()

    if args.task is not None:
        main_task(args.task, branch=args.branch)
    elif args.aggregate:
        main_aggregate(branch=args.branch)
        results = validate(branch=args.branch)
        if any(r["status"] == "FAIL" for r in results.values()):
            print("\nValidation FAILED. Troubleshooting:")
            print("If a library loses >80% of cells, lower the `lo_genes` minimum (QC_MIN_GENES_ABSOLUTE). If Scrublet fails on any sample, that sample's doublets are skipped but its cells remain in the dataset.")
            sys.exit(1)
    else:
        main(branch=args.branch)
        results = validate(branch=args.branch)
        if any(r["status"] == "FAIL" for r in results.values()):
            print("\nValidation FAILED. Troubleshooting:")
            print("If a library loses >80% of cells, lower the `lo_genes` minimum (QC_MIN_GENES_ABSOLUTE). If Scrublet fails on any sample, that sample's doublets are skipped but its cells remain in the dataset.")
            sys.exit(1)

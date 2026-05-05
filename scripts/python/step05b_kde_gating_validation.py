#!/usr/bin/env python
"""Step 5b: KDE-Based ADT Gating Validation
============================================

Objective
---------
Validate ThresholdR gates using non-parametric KDE peak/valley detection
and ADTnorm-style stain quality scoring (pure Python / scipy). Provides a
complementary, distribution-shape-aware check that does not assume Gaussian
mixture components -- like holding a second ruler alongside ThresholdR's GMM.
Output is diagnostic only; it does not replace ThresholdR or fg_prob gates.

Pipeline Position
-----------------
Phase:    1.2b (Tier 1)
Depends:  Step 5a (ThresholdR thresholds), Step 3 (CLR values in adata_qc),
          Step 1 (QC-filtered cells)
Feeds:    Step 5c (ThresholdR fitness metrics), downstream threshold-override
          decisions
Compute:  CPU, ~4 GB, ~2 min

Inputs
------
- data/processed/adata_qc.h5ad               -- QC-filtered AnnData with CLR
                                                 values in .obs and sample_id
                                                 metadata (h5ad)
- data/processed/thresholdr_thresholds.csv    -- Per-marker bi_thr from Step 5a
                                                 (CSV)
- data/processed/adt_counts_for_thresholdr.csv -- Raw ADT counts for KDE-only
                                                  markers (standalone log-
                                                  transform path)

Outputs
-------
- data/processed/kde_valleys.csv              -- Per-marker, per-sample KDE
                                                 valley positions (CSV)
- data/processed/kde_stain_quality.csv        -- Stain quality scores for every
                                                 marker in KDE_MARKERS (CSV)
- results/phase1_qc/kde_vs_thresholdr_comparison.csv
                                              -- Threshold positions, Jaccard
                                                 similarities, reclassification
                                                 counts for ADT_THRESHOLDR_MARKERS
                                                 (CSV)
- figures/phase1_kde_density_{marker}.png     -- Per-marker density with KDE
                                                 landmarks and ThresholdR
                                                 overlay (one PNG per marker)
- figures/phase1_kde_vs_thresholdr_grid.png   -- Side-by-side bar + Jaccard
                                                 heatmap comparison (PNG)

Config Constants
----------------
- ADATA_QC                   -- Path to QC-filtered AnnData
- THRESHOLDR_THRESHOLDS      -- Path to ThresholdR threshold CSV
- ADT_COUNTS_FOR_THRESHOLDR  -- Path to raw ADT counts CSV (includes KDE-only
                                markers that bypass joint CLR)
- PROCESSED_DIR              -- Base directory for intermediate outputs
- RESULTS_PHASE1_QC          -- Phase 1 QC results directory
- FIGURES_DIR                -- Figures output directory
- ADT_THRESHOLDR_MARKERS     -- Markers with ThresholdR GMM gates
- KDE_MARKERS                -- Markers inspected via KDE (superset of
                                ADT_THRESHOLDR_MARKERS; may include
                                negative-control markers)
"""

import sys
import os
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from config import (
    ADATA_QC,
    THRESHOLDR_THRESHOLDS,
    ADT_COUNTS_FOR_THRESHOLDR,
    PROCESSED_DIR,
    RESULTS_PHASE1_QC,
    FIGURES_DIR,
    ADT_THRESHOLDR_MARKERS,
    KDE_MARKERS,
)

warnings.filterwarnings("ignore", category=FutureWarning)

# Fine grid for smooth density estimation; 2000 points is cheap at ~161K cells
KDE_GRID_POINTS = 2000


# ── 5a: Load data ────────────────────────────────────────────────────────

def load_data():
    """Load adata_qc, ThresholdR thresholds, and raw ADT counts."""
    print("Loading adata_qc.h5ad (obs only, backed mode)...")
    # QC-filtered AnnData with CLR values in .obs (from Steps 1 + 3)
    # Only .obs is needed (CLR columns + sample_id); backed mode avoids
    # loading the full ~7 GB expression matrix into memory.
    adata_qc = sc.read_h5ad(ADATA_QC, backed="r")

    print("Loading ThresholdR thresholds...")
    # GMM-derived bi_thr per marker from Step 5a
    thr_df = pd.read_csv(THRESHOLDR_THRESHOLDS, index_col="marker")

    print("Loading raw ADT counts (for markers outside ADT_THRESHOLDR_MARKERS)...")
    # Includes any markers in KDE_MARKERS that bypass joint CLR normalization.
    raw_adt = pd.read_csv(ADT_COUNTS_FOR_THRESHOLDR, index_col=0)

    # Doublet removal and other post-export filtering may have reduced
    # adata_qc relative to the raw ADT export, so intersect barcodes.
    common = raw_adt.index.intersection(adata_qc.obs_names)
    if len(common) < adata_qc.n_obs:
        print(f"  WARNING: only {len(common)}/{adata_qc.n_obs} adata_qc barcodes found in raw_adt")
    raw_adt = raw_adt.reindex(adata_qc.obs_names).dropna()
    print(f"  raw_adt aligned to {len(raw_adt)} cells (matching adata_qc barcodes)")

    return adata_qc, thr_df, raw_adt


# ── 5b: KDE peak/valley detection ─────────────────────────────────────

def detect_peaks_valleys(values, bw_method="scott", _retries_left=2):
    """
    Detect peaks and valleys in a distribution using KDE.

    If fewer than 2 peaks are found, retries with halved bandwidth (up to
    ``_retries_left`` times). This handles over-smoothing from Scott's rule
    on skewed or heavy-tailed distributions.

    Returns dict with keys: x_grid, kde_y, peaks, valleys, peak_heights,
    valley_heights, bandwidth.
    """
    kde = gaussian_kde(values, bw_method=bw_method)
    x_min, x_max = values.min(), values.max()
    margin = (x_max - x_min) * 0.1
    x_grid = np.linspace(x_min - margin, x_max + margin, KDE_GRID_POINTS)
    kde_y = kde(x_grid)

    # Prominence thresholds tuned to avoid spurious peaks from noise
    peak_indices, _ = find_peaks(kde_y, height=0, prominence=0.001)

    # Valleys = local minima (peaks of negated density)
    valley_indices, _ = find_peaks(-kde_y, prominence=0.0001)

    # Scott's rule over-smooths skewed, heavy-right-tail distributions;
    # halving bandwidth recovers bimodality. Skip for large N where Scott is
    # already well-calibrated and bandwidth reduction amplifies noise.
    if len(peak_indices) < 2 and _retries_left > 0 and len(values) < 10000:
        return detect_peaks_valleys(
            values, bw_method=kde.factor * 0.5, _retries_left=_retries_left - 1
        )

    return {
        "x_grid": x_grid,
        "kde_y": kde_y,
        "peaks": x_grid[peak_indices],
        "peak_heights": kde_y[peak_indices],
        "valleys": x_grid[valley_indices],
        "valley_heights": kde_y[valley_indices],
        "bandwidth": kde.factor,
    }


def get_primary_valley(result, values=None):
    """
    Get the primary valley (gating threshold) — the deepest valley that best
    separates negative and positive populations.

    Parameters
    ----------
    result : dict from detect_peaks_valleys
    values : np.ndarray, optional
        Raw CLR values, used to compute SD for shoulder detection minimum distance.

    Strategy:
    - For 2+ peaks: find the valley with the greatest "contrast" — defined as
      the minimum height of its flanking peaks minus the valley height. This
      favors deep valleys between well-separated populations over shallow
      notches within a single population.
    - For 1 peak: shoulder detection via slope inflection point on the right
      tail (where the density transitions from concave to convex).
    """
    peaks = result["peaks"]
    peak_heights = result["peak_heights"]
    valleys = result["valleys"]
    valley_heights = result["valley_heights"]

    if len(peaks) < 2:
        # Unimodal: no valley between peaks, so use shoulder detection --
        # find the inflection point on the right tail (concave-to-convex
        # transition). Requires minimum distance of 1 SD from the peak to
        # avoid gating inside the main negative population.
        kde_y = result["kde_y"]
        x_grid = result["x_grid"]
        if len(peaks) == 1:
            peak_x = peaks[0]
            data_sd = np.std(values) if values is not None else 1.0
            min_distance = max(data_sd, 0.5)

            mask = x_grid > (peak_x + min_distance)
            if mask.any():
                d2 = np.gradient(np.gradient(kde_y[mask], x_grid[mask]), x_grid[mask])
                sign_changes = np.where(np.diff(np.sign(d2)))[0]
                if len(sign_changes) > 0:
                    idx = sign_changes[0]
                    return x_grid[mask][idx], kde_y[mask][idx]
                # Fallback: steepest descent beyond min_distance
                slope_after = np.gradient(kde_y[mask], x_grid[mask])
                min_slope_idx = np.argmin(slope_after)
                return x_grid[mask][min_slope_idx], kde_y[mask][min_slope_idx]
        return np.nan, np.nan

    # Bimodal+: pick the valley with greatest "contrast" (min flanking peak
    # height minus valley height). Favors deep separations over shallow
    # notches within a single population.
    best_valley_x = np.nan
    best_valley_h = np.nan
    best_contrast = -np.inf

    for vi, (vx, vh) in enumerate(zip(valleys, valley_heights)):
        left_peaks = peaks[peaks < vx]
        right_peaks = peaks[peaks > vx]
        if len(left_peaks) == 0 or len(right_peaks) == 0:
            continue

        left_idx = np.argmin(np.abs(peaks - left_peaks[-1]))
        right_idx = np.argmin(np.abs(peaks - right_peaks[0]))
        flanking_min_height = min(peak_heights[left_idx], peak_heights[right_idx])

        contrast = flanking_min_height - vh
        if contrast > best_contrast:
            best_contrast = contrast
            best_valley_x = vx
            best_valley_h = vh

    if np.isnan(best_valley_x):
        return (peaks[0] + peaks[-1]) / 2, np.nan

    return best_valley_x, best_valley_h


# ── 5c: Stain quality score ───────────────────────────────────────────

def compute_stain_quality(values, result):
    """
    ADTnorm-style stain quality score.

    For bimodal (2+ peaks):
      (pos_mode - neg_mode) / (SD_neg + SD_pos)
        × (pos_height - valley_height + 1)
        × (AUC_pos + 1)

    For unimodal (1 peak):
      (valley - peak_mode) / SD_all
        × (0 - valley_height + 1)
        × (AUC_right_tail + 1)

    WARNING: bimodal and unimodal scores are on different scales (bimodal
    typically 2-10, unimodal typically 0.01-0.5). Do not compare scores
    across peak-count categories directly. For cross-marker ranking, the
    number of peaks is the primary signal; the score is secondary.
    """
    peaks = result["peaks"]
    peak_heights = result["peak_heights"]
    valley_x, valley_height = get_primary_valley(result, values)

    if len(peaks) >= 2:
        order = np.argsort(peaks)
        sorted_peaks = peaks[order]
        sorted_heights = peak_heights[order]

        neg_mode = sorted_peaks[0]
        pos_mode = sorted_peaks[-1]
        pos_height = sorted_heights[-1]

        # Split cells at the valley to estimate per-population SD
        if not np.isnan(valley_x):
            neg_cells = values[values < valley_x]
            pos_cells = values[values >= valley_x]
        else:
            mid = (neg_mode + pos_mode) / 2
            neg_cells = values[values < mid]
            pos_cells = values[values >= mid]

        sd_neg = np.std(neg_cells) if len(neg_cells) > 1 else 1.0
        sd_pos = np.std(pos_cells) if len(pos_cells) > 1 else 1.0
        auc_pos = len(pos_cells) / len(values) if len(values) > 0 else 0.0

        if (sd_neg + sd_pos) == 0:
            return 0.0

        valley_h = valley_height if not np.isnan(valley_height) else 0.0
        score = ((pos_mode - neg_mode) / (sd_neg + sd_pos)
                 * (pos_height - valley_h + 1)
                 * (auc_pos + 1))
        return score

    elif len(peaks) == 1:
        # Unimodal formula: score will be low, reflecting poor separation
        peak_mode = peaks[0]
        sd_all = np.std(values) if len(values) > 1 else 1.0
        if sd_all == 0 or np.isnan(valley_x):
            return 0.0

        right_tail = values[values > valley_x]
        auc_right = len(right_tail) / len(values) if len(values) > 0 else 0.0
        valley_h = valley_height if not np.isnan(valley_height) else 0.0

        score = ((valley_x - peak_mode) / sd_all
                 * (0 - valley_h + 1)
                 * (auc_right + 1))
        return score

    return 0.0


# ── 5d: Method comparison ─────────────────────────────────────────────

def jaccard(set_a, set_b):
    """Jaccard similarity between two boolean arrays."""
    intersection = (set_a & set_b).sum()
    union = (set_a | set_b).sum()
    return intersection / union if union > 0 else 0.0


def compare_methods(adata_qc, thr_df, valleys_df):
    """
    Compare positive calls across 3 methods: KDE valley, ThresholdR, fg_prob.

    Runs only on ADT_THRESHOLDR_MARKERS; KDE-only markers have no ThresholdR
    or fg_prob data to compare against.
    """
    rows = []
    for marker in ADT_THRESHOLDR_MARKERS:
        clr_col = f"ADT_{marker}"
        clr_values = adata_qc.obs[clr_col].values.astype(float)

        kde_valley = valleys_df.loc[
            valleys_df["marker"] == marker, "valley_aggregated"
        ].values[0]

        # GMM-derived threshold from Step 5a
        thr_val = thr_df.loc[marker, "bi_thr"]

        # totalVI fg_prob only available after Step 6; N/A in early pipeline runs
        fg_col = f"ADT_{marker}_fg_prob"
        has_fg = fg_col in adata_qc.obs.columns
        if has_fg:
            fg_values = adata_qc.obs[fg_col].values.astype(float)
            fg_pos = fg_values > 0.5
        else:
            fg_pos = np.zeros(len(clr_values), dtype=bool)

        kde_pos = clr_values > kde_valley if not np.isnan(kde_valley) else np.zeros(len(clr_values), dtype=bool)
        thr_pos = clr_values > thr_val

        n_cells = len(clr_values)
        rows.append({
            "marker": marker,
            "kde_valley": round(kde_valley, 4) if not np.isnan(kde_valley) else np.nan,
            "thresholdr_bi_thr": round(thr_val, 4),
            "n_kde_pos": int(kde_pos.sum()),
            "n_thr_pos": int(thr_pos.sum()),
            "n_fg_pos": int(fg_pos.sum()) if has_fg else np.nan,
            "pct_kde_pos": round(100 * kde_pos.sum() / n_cells, 2),
            "pct_thr_pos": round(100 * thr_pos.sum() / n_cells, 2),
            "pct_fg_pos": round(100 * fg_pos.sum() / n_cells, 2) if has_fg else np.nan,
            "jaccard_kde_thr": round(jaccard(kde_pos, thr_pos), 4),
            "jaccard_kde_fg": round(jaccard(kde_pos, fg_pos), 4) if has_fg else np.nan,
            "jaccard_thr_fg": round(jaccard(thr_pos, fg_pos), 4) if has_fg else np.nan,
            "reclassified_kde_vs_thr": int((kde_pos != thr_pos).sum()),
        })

    return pd.DataFrame(rows)


# ── 5e: Visualization ────────────────────────────────────────────────

def plot_per_marker_density(marker, values, result, valley_x, thr_val, sample_ids, output_path):
    """Per-marker density plot with KDE, peaks, valleys, and ThresholdR threshold."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left panel: aggregated KDE with peak/valley landmarks + threshold lines
    ax = axes[0]
    ax.fill_between(result["x_grid"], result["kde_y"], alpha=0.3, color="steelblue")
    ax.plot(result["x_grid"], result["kde_y"], color="steelblue", linewidth=2, label="KDE (all cells)")

    for px, ph in zip(result["peaks"], result["peak_heights"]):
        ax.plot(px, ph, "v", color="green", markersize=10, zorder=5)
    for vx in result["valleys"]:
        ax.axvline(vx, color="gray", linestyle=":", alpha=0.5)

    if not np.isnan(valley_x):
        ax.axvline(valley_x, color="blue", linestyle="--", linewidth=2, label=f"KDE valley={valley_x:.2f}")
    # ThresholdR line absent for markers outside ADT_THRESHOLDR_MARKERS
    if not np.isnan(thr_val):
        ax.axvline(thr_val, color="red", linestyle="--", linewidth=2, label=f"ThresholdR={thr_val:.2f}")

    ax.set_title(f"{marker} — All Cells (n={len(values):,})")
    ax.set_xlabel("CLR")
    ax.set_ylabel("Density")
    ax.legend(fontsize=8)

    # Right panel: per-sample KDE overlays reveal batch-specific threshold shifts
    ax = axes[1]
    unique_samples = np.unique(sample_ids)
    colors = plt.cm.tab20(np.linspace(0, 1, len(unique_samples)))
    for i, sid in enumerate(unique_samples):
        mask = sample_ids == sid
        sample_vals = values[mask]
        if len(sample_vals) < 50:
            continue
        try:
            kde_s = gaussian_kde(sample_vals)
            x_grid = result["x_grid"]
            ax.plot(x_grid, kde_s(x_grid), color=colors[i], alpha=0.5, linewidth=1)
        except Exception:
            pass

    if not np.isnan(valley_x):
        ax.axvline(valley_x, color="blue", linestyle="--", linewidth=2)
    if not np.isnan(thr_val):
        ax.axvline(thr_val, color="red", linestyle="--", linewidth=2)
    ax.set_title(f"{marker} — Per Sample (n={len(unique_samples)})")
    ax.set_xlabel("CLR")
    ax.set_ylabel("Density")

    fig.suptitle(f"Step 5b: KDE Gating Validation — {marker}", fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_comparison_grid(comparison_df, output_path):
    """Comparison grid: bar chart of positive percentages + Jaccard heatmap."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    markers = comparison_df["marker"].values

    # Left panel: positive call rate comparison across 3 gating methods
    ax = axes[0]
    x = np.arange(len(markers))
    w = 0.25
    ax.bar(x - w, comparison_df["pct_kde_pos"], w, label="KDE", color="#2196F3")
    ax.bar(x, comparison_df["pct_thr_pos"], w, label="ThresholdR", color="#FF5722")
    fg_pcts = comparison_df["pct_fg_pos"].fillna(0)
    has_fg = comparison_df["pct_fg_pos"].notna().any()
    ax.bar(x + w, fg_pcts, w, label="fg_prob>0.5" + (" (N/A)" if not has_fg else ""), color="#4CAF50", alpha=0.3 if not has_fg else 1.0)
    ax.set_xticks(x)
    ax.set_xticklabels(markers)
    ax.set_ylabel("% Positive")
    ax.set_title("Positive Call Rates by Method")
    ax.legend()

    # Right panel: Jaccard similarity shows inter-method agreement per marker
    ax = axes[1]
    jaccard_data = np.zeros((len(markers), 3))
    for i, (_, row) in enumerate(comparison_df.iterrows()):
        jaccard_data[i] = [
            row["jaccard_kde_thr"],
            row["jaccard_kde_fg"] if not pd.isna(row["jaccard_kde_fg"]) else 0.0,
            row["jaccard_thr_fg"] if not pd.isna(row["jaccard_thr_fg"]) else 0.0,
        ]

    im = ax.imshow(jaccard_data, cmap="YlGn", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["KDE vs ThresholdR", "KDE vs fg_prob", "ThresholdR vs fg_prob"], fontsize=9)
    ax.set_yticks(range(len(markers)))
    ax.set_yticklabels(markers)
    ax.set_title("Jaccard Similarity")

    for i in range(len(markers)):
        for j in range(3):
            ax.text(j, i, f"{jaccard_data[i, j]:.2f}", ha="center", va="center", fontsize=10)

    fig.colorbar(im, ax=ax, shrink=0.8)
    fig.suptitle("Step 5b: KDE vs ThresholdR vs fg_prob Comparison", fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ── 5f: Main ─────────────────────────────────────────────────────────

def main():
    print("=" * 60)
    print("Step 5b: KDE-Based ADT Gating Validation")
    print("=" * 60)

    adata_qc, thr_df, raw_adt = load_data()
    sample_ids = adata_qc.obs["sample_id"].values

    # KDE_MARKERS may extend ADT_THRESHOLDR_MARKERS with non-functional
    # controls (e.g. negative/control antibodies listed in the panel YAML).
    valleys_rows = []
    quality_rows = []

    for marker in KDE_MARKERS:
        print(f"\n--- {marker} ---")

        if marker not in ADT_THRESHOLDR_MARKERS:
            # Markers outside ADT_THRESHOLDR_MARKERS are not part of the joint
            # CLR normalization. Use a standalone mean-centered log1p for KDE
            # shape analysis; stain-quality magnitudes are not comparable
            # to functional markers.
            raw_vals = raw_adt[marker].values.astype(float)
            log_vals = np.log1p(raw_vals)
            clr_values = log_vals - log_vals.mean()
            # raw_adt may have fewer cells than adata_qc (dropna after reindex),
            # so use matching sample_ids for per-sample KDE
            marker_sample_ids = adata_qc.obs.loc[raw_adt.index, "sample_id"].values
            print(f"  {marker}: standalone log-transform from raw counts, n={len(clr_values)}")
        else:
            clr_col = f"ADT_{marker}"
            clr_values = adata_qc.obs[clr_col].values.astype(float)
            marker_sample_ids = sample_ids

        result = detect_peaks_valleys(clr_values)
        valley_x, valley_h = get_primary_valley(result, clr_values)
        print(f"  Peaks: {result['peaks']}")
        print(f"  Primary valley: {valley_x:.4f}" if not np.isnan(valley_x) else "  Primary valley: N/A")

        # ADTnorm-style stain quality (see KB: ADTnorm theory.md)
        sq = compute_stain_quality(clr_values, result)
        print(f"  Stain quality: {sq:.4f}")

        # Per-sample KDE detects batch-specific threshold shifts per library
        per_sample_valleys = []
        for sid in np.unique(marker_sample_ids):
            mask = marker_sample_ids == sid
            sample_vals = clr_values[mask]
            if len(sample_vals) < 50:
                continue
            try:
                s_result = detect_peaks_valleys(sample_vals)
                s_valley, _ = get_primary_valley(s_result, sample_vals)
                per_sample_valleys.append({
                    "marker": marker,
                    "sample_id": sid,
                    "valley": round(s_valley, 4) if not np.isnan(s_valley) else np.nan,
                    "n_peaks": len(s_result["peaks"]),
                    "n_cells": int(mask.sum()),
                })
            except Exception as e:
                print(f"    WARNING: KDE failed for {marker}/{sid} ({mask.sum()} cells): {e}")
                per_sample_valleys.append({
                    "marker": marker,
                    "sample_id": sid,
                    "valley": np.nan,
                    "n_peaks": 0,
                    "n_cells": int(mask.sum()),
                })

        valleys_rows.extend(per_sample_valleys)

        # Median across per-sample valleys; falls back to all-cell valley
        sample_valleys = [v["valley"] for v in per_sample_valleys if not np.isnan(v.get("valley", np.nan))]
        agg_valley = np.median(sample_valleys) if sample_valleys else valley_x

        quality_rows.append({
            "marker": marker,
            "n_peaks": len(result["peaks"]),
            "stain_quality_score": round(sq, 4),
            "valley_all_cells": round(valley_x, 4) if not np.isnan(valley_x) else np.nan,
            "valley_median_per_sample": round(agg_valley, 4) if not np.isnan(agg_valley) else np.nan,
            "valley_sd_per_sample": round(np.std(sample_valleys), 4) if len(sample_valleys) > 1 else np.nan,
        })

        thr_val = thr_df.loc[marker, "bi_thr"] if marker in thr_df.index else np.nan
        fig_path = os.path.join(FIGURES_DIR, f"phase1_kde_density_{marker}.png")
        plot_per_marker_density(marker, clr_values, result, valley_x, thr_val, marker_sample_ids, fig_path)
        print(f"  Saved {fig_path}")

    valleys_df = pd.DataFrame(valleys_rows)
    valleys_path = os.path.join(PROCESSED_DIR, "kde_valleys.csv")
    valleys_df.to_csv(valleys_path, index=False)
    print(f"\nSaved per-sample valleys to {valleys_path}")

    quality_df = pd.DataFrame(quality_rows)
    quality_path = os.path.join(PROCESSED_DIR, "kde_stain_quality.csv")
    quality_df.to_csv(quality_path, index=False)
    print(f"Saved stain quality scores to {quality_path}")

    # Comparison restricted to ADT_THRESHOLDR_MARKERS; KDE-only markers have
    # no ThresholdR/fg_prob to compare. Comparison CSV therefore has fewer rows
    # than kde_stain_quality.csv when KDE_MARKERS is a strict superset.
    # All-cell valley preferred over per-sample median: more stable at the
    # full-cohort cell count; per-sample valleys have high variance for
    # unimodal markers.
    agg_valleys = []
    for _, row in quality_df.iterrows():
        agg_valleys.append({
            "marker": row["marker"],
            "valley_aggregated": row["valley_all_cells"],
        })
    agg_df = pd.DataFrame(agg_valleys)

    comparison_df = compare_methods(adata_qc, thr_df, agg_df)
    comparison_path = os.path.join(RESULTS_PHASE1_QC, "kde_vs_thresholdr_comparison.csv")
    os.makedirs(RESULTS_PHASE1_QC, exist_ok=True)
    comparison_df.to_csv(comparison_path, index=False)
    print(f"Saved comparison to {comparison_path}")

    grid_path = os.path.join(FIGURES_DIR, "phase1_kde_vs_thresholdr_grid.png")
    plot_comparison_grid(comparison_df, grid_path)
    print(f"Saved comparison grid to {grid_path}")

    print("\n" + "=" * 60)
    print("STAIN QUALITY SUMMARY")
    print("=" * 60)
    for _, row in quality_df.iterrows():
        print(f"  {row['marker']:8s}  peaks={row['n_peaks']}  "
              f"stain_quality={row['stain_quality_score']:.3f}  "
              f"valley={row['valley_all_cells']:.3f}" if not pd.isna(row['valley_all_cells'])
              else f"  {row['marker']:8s}  peaks={row['n_peaks']}  "
              f"stain_quality={row['stain_quality_score']:.3f}  valley=N/A")

    # Non-functional markers (KDE_MARKERS \ ADT_THRESHOLDR_MARKERS) are
    # expected to score lowest, since they are intentionally retained as
    # negative controls in the panel YAML.
    # WARNING: multimodal scores (noise fragmentation) are not directly
    # comparable to unimodal scores.
    non_functional = set(KDE_MARKERS) - set(ADT_THRESHOLDR_MARKERS)
    functional_df = quality_df[~quality_df["marker"].isin(non_functional)]
    for nf_marker in non_functional:
        nf_row = quality_df[quality_df["marker"] == nf_marker]
        if not nf_row.empty and not functional_df.empty:
            nf_sq = nf_row.iloc[0]["stain_quality_score"]
            other_min = functional_df["stain_quality_score"].min()
            if nf_sq < other_min:
                print(f"\n  ✓ {nf_marker} has lowest stain quality (expected — non-functional antibody)")
            else:
                print(f"\n  ⚠ {nf_marker} stain quality ({nf_sq:.3f}) is NOT lowest — "
                      f"functional min is {other_min:.3f}")

    print("\n" + "=" * 60)
    print("COMPARISON SUMMARY")
    print("=" * 60)
    for _, row in comparison_df.iterrows():
        fg_str = f"fg_prob={row['pct_fg_pos']:5.1f}%" if not pd.isna(row['pct_fg_pos']) else "fg_prob=  N/A"
        jfg_str = f"J(KDE,FG)={row['jaccard_kde_fg']:.3f}" if not pd.isna(row['jaccard_kde_fg']) else "J(KDE,FG)=  N/A"
        print(f"  {row['marker']:8s}  KDE={row['pct_kde_pos']:5.1f}%  "
              f"ThresholdR={row['pct_thr_pos']:5.1f}%  "
              f"{fg_str}  "
              f"J(KDE,THR)={row['jaccard_kde_thr']:.3f}  "
              f"{jfg_str}")


if __name__ == "__main__":
    main()

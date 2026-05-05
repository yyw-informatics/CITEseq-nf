#!/usr/bin/env python
"""Step 5c: ThresholdR Fitness Metrics
=====================================

Objective
---------
Quantify how well each ThresholdR GMM fit separates negative from positive
cells, producing a per-marker PASS/WARN/FAIL verdict.  Think of it as a
"report card" for each ADT gate: overlap measures how much the two Gaussian
bells bleed into each other, stain quality captures antibody signal strength,
and RNA concordance checks whether the ADT gate agrees with the transcriptome.
Markers that FAIL here are candidates for totalVI fg_prob override.

Pipeline Position
-----------------
Phase:    1.2b (Tier 1)
Depends:  Step 5a (ThresholdR thresholds), Step 4 (normalized RNA for concordance)
Feeds:    Step 9 (informs fg_prob override decisions), documentation / QC audit
Compute:  CPU, 32 GB h_vmem (two 7 GB h5ad via backed mode), <2 min

Inputs
------
- data/processed/thresholdr_thresholds.csv  — Per-marker GMM thresholds from R ThresholdR (CSV)
- data/processed/adata_qc.h5ad              — CLR values in .obs (h5ad, backed read)
- data/processed/adata_normalized.h5ad      — Log-normalized RNA for concordance genes (h5ad, gene-subset)

Outputs
-------
- results/phase1_qc/thresholdr_fitness_metrics.csv  — Per-marker fitness metrics + verdict (CSV)
- figures/phase1_thresholdr_fitness_summary.png      — GMM overlay plots + metrics table (PNG)

Config Constants
----------------
- ADT_THRESHOLDR_MARKERS       — Panel-driven list of markers to evaluate
- ADT_RNA_CONCORDANCE          (dict) — RNA gene sets per marker
- THRESHOLDR_OVERLAP_THRESHOLD (0.15)  — PASS ceiling for overlap coefficient
- THRESHOLDR_OVERLAP_FAIL      (0.30)  — FAIL floor for overlap coefficient
- THRESHOLDR_STAIN_QUALITY_THRESHOLD (2.0) — PASS floor for stain quality score
- THRESHOLDR_STAIN_QUALITY_FAIL      (1.0) — FAIL ceiling for stain quality score
- THRESHOLDR_RNA_CONCORDANCE_THRESHOLD (1.5) — PASS floor for RNA enrichment ratio
- THRESHOLDR_RNA_CONCORDANCE_FAIL      (1.2) — FAIL ceiling for RNA enrichment ratio
- THRESHOLDR_BIC_MARGIN_FAIL   (0.0)  — BIC margin below this means k=1 fits better
"""

import sys
import os
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.stats import norm, gaussian_kde
from sklearn.mixture import GaussianMixture

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from config import (
    ADATA_QC,
    ADATA_NORMALIZED,
    THRESHOLDR_THRESHOLDS,
    RESULTS_PHASE1_QC,
    FIGURES_DIR,
    ADT_THRESHOLDR_MARKERS,
    ADT_RNA_CONCORDANCE,
    THRESHOLDR_OVERLAP_THRESHOLD,
    THRESHOLDR_STAIN_QUALITY_THRESHOLD,
    THRESHOLDR_RNA_CONCORDANCE_THRESHOLD,
    THRESHOLDR_OVERLAP_FAIL,
    THRESHOLDR_STAIN_QUALITY_FAIL,
    THRESHOLDR_RNA_CONCORDANCE_FAIL,
    THRESHOLDR_BIC_MARGIN_FAIL,
)

warnings.filterwarnings("ignore", category=FutureWarning)


# ── 5c-a: Load data ───────────────────────────────────────────────────

def load_data():
    """Load ThresholdR thresholds, adata_qc (.obs only), and adata_normalized.

    adata_qc is opened in backed mode and only .obs is extracted — the full
    X/layers are never loaded, keeping peak memory low and avoiding any
    HDF5 corruption in layers from concurrent writes by other steps.

    adata_normalized is also backed-loaded and subsetted to the small set
    of RNA concordance genes before reading into memory.
    """
    print("Loading ThresholdR thresholds...")
    thr_df = pd.read_csv(THRESHOLDR_THRESHOLDS, index_col="marker")

    print("Loading adata_qc.h5ad (obs only, backed mode)...")
    adata_qc = sc.read_h5ad(ADATA_QC, backed="r")

    # Collect the small set of RNA genes needed for concordance
    rna_genes_needed = set()
    for info in ADT_RNA_CONCORDANCE.values():
        rna_genes_needed.update(info.get("positive", []))
    rna_genes_needed = sorted(rna_genes_needed)
    print(f"Loading adata_normalized.h5ad (subsetting to {len(rna_genes_needed)} concordance genes)...")

    adata_norm_full = sc.read_h5ad(ADATA_NORMALIZED, backed="r")
    available = [g for g in rna_genes_needed if g in adata_norm_full.var_names]
    adata_norm = adata_norm_full[:, available].to_memory()
    adata_norm_full.file.close()

    return thr_df, adata_qc, adata_norm


# ── 5c-b: Fit GMM and compute metrics ────────────────────────────────

def compute_overlap_coefficient(means, stds, weights):
    """
    Compute overlap coefficient between two Gaussian components.
    overlap = ∫ min(w1·φ1(x), w2·φ2(x)) dx, estimated numerically.
    """
    mu1, mu2 = means
    s1, s2 = stds
    w1, w2 = weights
    x_range = np.linspace(
        min(mu1 - 4 * s1, mu2 - 4 * s2),
        max(mu1 + 4 * s1, mu2 + 4 * s2),
        2000,
    )
    pdf1 = w1 * norm.pdf(x_range, mu1, s1)
    pdf2 = w2 * norm.pdf(x_range, mu2, s2)
    overlap = np.trapezoid(np.minimum(pdf1, pdf2), x_range)
    return overlap


def compute_stain_quality_score(clr_values, means, stds, weights, valley):
    """
    ADTnorm-style stain quality score for bimodal distribution:
    (pos_mode - neg_mode) / (SD_neg + SD_pos)
      × (pos_height - valley_height + 1)
      × (AUC_pos + 1)
    """
    # Identify negative (lower mean) and positive (higher mean) components
    idx_neg, idx_pos = (0, 1) if means[0] < means[1] else (1, 0)
    neg_mode = means[idx_neg]
    pos_mode = means[idx_pos]
    sd_neg = stds[idx_neg]
    sd_pos = stds[idx_pos]

    # Separation term
    separation = (pos_mode - neg_mode) / (sd_neg + sd_pos)

    # Heights at modes and valley (from KDE for robustness)
    try:
        kde = gaussian_kde(clr_values)
        pos_height = kde(pos_mode)[0]
        valley_height = kde(valley)[0] if valley is not None else 0.0
    except Exception:
        pos_height = 1.0
        valley_height = 0.0

    # AUC of positive component
    auc_pos = weights[idx_pos]

    score = separation * (pos_height - valley_height + 1) * (auc_pos + 1)
    return score


def compute_rna_concordance(adata_norm, adata_qc, marker, threshold):
    """
    RNA concordance ratio: mean(RNA_positive | ADT+) / mean(RNA_positive | ADT-).
    Uses positive RNA genes from ADT_RNA_CONCORDANCE config.
    """
    if marker not in ADT_RNA_CONCORDANCE:
        return np.nan

    rna_genes = ADT_RNA_CONCORDANCE[marker].get("positive", [])
    if not rna_genes:
        return np.nan

    # Filter to genes present in the normalized data
    available_genes = [g for g in rna_genes if g in adata_norm.var_names]
    if not available_genes:
        return np.nan

    clr_col = f"ADT_{marker}"
    clr_values = adata_qc.obs[clr_col].values
    adt_pos = clr_values > threshold
    adt_neg = ~adt_pos

    if adt_pos.sum() < 10 or adt_neg.sum() < 10:
        return np.nan

    # Per-gene mean expression in ADT+ vs ADT- cells, then average across genes.
    # This prevents highly-expressed genes (e.g., NKG7) from dominating the ratio.
    rna_matrix = adata_norm[:, available_genes].X
    if hasattr(rna_matrix, "toarray"):
        rna_matrix = rna_matrix.toarray()

    # Shape: (n_cells, n_genes) → per-gene means → average across genes
    per_gene_mean_pos = rna_matrix[adt_pos].mean(axis=0)  # (n_genes,)
    per_gene_mean_neg = rna_matrix[adt_neg].mean(axis=0)  # (n_genes,)

    # Per-gene ratio, then median (robust to outlier genes)
    with np.errstate(divide="ignore", invalid="ignore"):
        per_gene_ratio = per_gene_mean_pos / per_gene_mean_neg
    valid = np.isfinite(per_gene_ratio)
    if not valid.any():
        return np.inf

    return float(np.median(per_gene_ratio[valid]))


def compute_fitness_metrics(thr_df, adata_qc, adata_norm):
    """Compute all fitness metrics for each ThresholdR marker.

    Returns (metrics_df, gmm_params) where gmm_params is a dict of
    {marker: (means, stds, weights)} for reuse in plotting.
    """
    results = []
    gmm_params = {}

    for marker in ADT_THRESHOLDR_MARKERS:
        print(f"\n--- {marker} ---")
        clr_col = f"ADT_{marker}"
        clr_values = adata_qc.obs[clr_col].values.astype(float)

        # Re-fit 2-component GMM in Python (R ThresholdR uses mclust; we mirror with sklearn)
        gmm = GaussianMixture(n_components=2, random_state=42)
        gmm.fit(clr_values.reshape(-1, 1))
        means = gmm.means_.flatten()
        stds = np.sqrt(gmm.covariances_.flatten())
        weights = gmm.weights_

        # Sort by mean (negative component first)
        order = np.argsort(means)
        means = means[order]
        stds = stds[order]
        weights = weights[order]

        gmm_params[marker] = (means.copy(), stds.copy(), weights.copy())

        print(f"  GMM: μ=[{means[0]:.3f}, {means[1]:.3f}], "
              f"σ=[{stds[0]:.3f}, {stds[1]:.3f}], "
              f"λ=[{weights[0]:.3f}, {weights[1]:.3f}]")

        # BIC margin: (BIC_k1 - BIC_k2) / |BIC_k1|  — positive means k=2 is better.
        # Negative margin → unimodal fits better → ThresholdR bimodal gate is unsupported.
        gmm_k1 = GaussianMixture(n_components=1, random_state=42)
        gmm_k1.fit(clr_values.reshape(-1, 1))
        bic_k1 = gmm_k1.bic(clr_values.reshape(-1, 1))
        bic_k2 = gmm.bic(clr_values.reshape(-1, 1))
        bic_margin = (bic_k1 - bic_k2) / abs(bic_k1)
        print(f"  BIC margin: {bic_margin:.4f}")

        # Overlap coefficient
        overlap = compute_overlap_coefficient(means, stds, weights)
        print(f"  Overlap: {overlap:.4f}")

        # Threshold-to-valley distance — how far ThresholdR's bi_thr (mu_neg + 3*sigma_neg)
        # deviates from Python GMM's expected threshold, in units of sigma_neg.
        thr_val = thr_df.loc[marker, "bi_thr"]
        expected_thr = means[0] + 3 * stds[0]
        thr_valley_dist = abs(expected_thr - thr_val) / stds[0] if stds[0] > 0 else np.nan
        print(f"  Thr-valley dist: {thr_valley_dist:.4f} "
              f"(expected={expected_thr:.3f}, actual={thr_val:.3f})")

        # Min component fraction
        min_comp_frac = min(weights)
        print(f"  Min component frac: {min_comp_frac:.4f}")

        # ADTnorm-style stain quality — uses bi_cut (GMM intersection) as valley
        valley = thr_df.loc[marker, "bi_cut"]
        stain_quality = compute_stain_quality_score(
            clr_values, means, stds, weights, valley
        )
        print(f"  Stain quality: {stain_quality:.4f}")

        # RNA concordance
        rna_concordance = compute_rna_concordance(
            adata_norm, adata_qc, marker, thr_val
        )
        print(f"  RNA concordance: {rna_concordance:.4f}" if not np.isnan(rna_concordance) else "  RNA concordance: N/A")

        # Fitness verdict
        verdict = compute_verdict(overlap, stain_quality, rna_concordance, bic_margin)
        print(f"  Verdict: {verdict}")

        results.append({
            "marker": marker,
            "k_selected": int(thr_df.loc[marker, "k_selected"]),
            "bic_margin": round(bic_margin, 4),
            "overlap_coeff": round(overlap, 4),
            "stain_quality_score": round(stain_quality, 4),
            "min_component_frac": round(min_comp_frac, 4),
            "thr_valley_dist": round(thr_valley_dist, 4),
            "rna_concordance_ratio": round(rna_concordance, 4) if not np.isnan(rna_concordance) else np.nan,
            "fitness_verdict": verdict,
        })

    return pd.DataFrame(results), gmm_params


def compute_verdict(overlap, stain_quality, rna_concordance, bic_margin):
    """
    Fitness verdict logic:
      FAIL: bic_margin < 0 (bimodality unsupported) OR overlap > 0.30
            OR stain_quality < 1.0 OR rna_concordance < 1.2
      PASS: overlap < threshold AND stain_quality > threshold AND rna_concordance > threshold
      WARN: otherwise
    """
    # FAIL conditions (any one triggers)
    if bic_margin < THRESHOLDR_BIC_MARGIN_FAIL:
        return "FAIL"
    if overlap > THRESHOLDR_OVERLAP_FAIL:
        return "FAIL"
    if stain_quality < THRESHOLDR_STAIN_QUALITY_FAIL:
        return "FAIL"
    if not np.isnan(rna_concordance) and rna_concordance < THRESHOLDR_RNA_CONCORDANCE_FAIL:
        return "FAIL"

    # PASS conditions (all must hold)
    if (overlap < THRESHOLDR_OVERLAP_THRESHOLD
            and stain_quality > THRESHOLDR_STAIN_QUALITY_THRESHOLD
            and (np.isnan(rna_concordance) or rna_concordance > THRESHOLDR_RNA_CONCORDANCE_THRESHOLD)):
        return "PASS"

    return "WARN"


# ── 5c-c: Visualization ──────────────────────────────────────────────

def plot_fitness_summary(metrics_df, adata_qc, thr_df, gmm_params):
    """Generate fitness summary figure: metrics heatmap + per-marker GMM plots."""
    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(2, 4, figure=fig, hspace=0.35, wspace=0.3)

    # Top row: per-marker CLR distribution with GMM overlay
    for i, marker in enumerate(ADT_THRESHOLDR_MARKERS):
        ax = fig.add_subplot(gs[0, i])
        clr_values = adata_qc.obs[f"ADT_{marker}"].values.astype(float)

        # Histogram
        ax.hist(clr_values, bins=100, density=True, alpha=0.5, color="gray", label="CLR")

        # Reuse GMM params from metrics computation (avoid redundant re-fit)
        means, stds, weights = gmm_params[marker]

        x_range = np.linspace(clr_values.min(), clr_values.max(), 500)
        for j, color in enumerate(["#2196F3", "#FF5722"]):
            ax.plot(x_range, weights[j] * norm.pdf(x_range, means[j], stds[j]),
                    color=color, linewidth=2, label=f"Comp {j+1}")

        # ThresholdR threshold line
        thr_val = thr_df.loc[marker, "bi_thr"]
        ax.axvline(thr_val, color="red", linestyle="--", linewidth=1.5, label=f"bi_thr={thr_val:.2f}")

        row = metrics_df[metrics_df["marker"] == marker].iloc[0]
        verdict = row["fitness_verdict"]
        verdict_colors = {"PASS": "green", "WARN": "orange", "FAIL": "red"}
        ax.set_title(f"{marker} [{verdict}]", fontsize=12, fontweight="bold",
                     color=verdict_colors.get(verdict, "black"))
        ax.set_xlabel("CLR")
        ax.set_ylabel("Density")
        if i == 0:
            ax.legend(fontsize=7)

    # Bottom row: metrics table as heatmap
    ax_table = fig.add_subplot(gs[1, :])
    ax_table.axis("off")

    metric_cols = ["bic_margin", "overlap_coeff", "stain_quality_score",
                   "min_component_frac", "thr_valley_dist", "rna_concordance_ratio"]
    display_names = ["BIC margin", "Overlap", "Stain quality", "Min comp frac",
                     "Thr-valley dist", "RNA concordance"]

    table_data = []
    cell_colors = []
    for _, row in metrics_df.iterrows():
        row_data = [row["marker"]]
        row_colors = ["white"]
        for col in metric_cols:
            val = row[col]
            row_data.append(f"{val:.3f}" if not pd.isna(val) else "N/A")
            row_colors.append("white")
        row_data.append(row["fitness_verdict"])
        verdict_bg = {"PASS": "#C8E6C9", "WARN": "#FFF9C4", "FAIL": "#FFCDD2"}
        row_colors.append(verdict_bg.get(row["fitness_verdict"], "white"))
        table_data.append(row_data)
        cell_colors.append(row_colors)

    col_labels = ["Marker"] + display_names + ["Verdict"]
    table = ax_table.table(
        cellText=table_data,
        colLabels=col_labels,
        cellColours=cell_colors,
        loc="center",
        cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.0, 1.5)
    # Bold header
    for j in range(len(col_labels)):
        table[0, j].set_text_props(fontweight="bold")

    fig.suptitle("Step 5c: ThresholdR Fitness Metrics", fontsize=14, fontweight="bold", y=0.98)
    return fig


# ── 5c-d: Main ──────────────────────────────────────────────────────

def main():
    print("=" * 60)
    print("Step 5c: ThresholdR Fitness Metrics")
    print("=" * 60)

    thr_df, adata_qc, adata_norm = load_data()

    # Compute metrics
    metrics_df, gmm_params = compute_fitness_metrics(thr_df, adata_qc, adata_norm)

    # Save CSV
    os.makedirs(RESULTS_PHASE1_QC, exist_ok=True)
    output_csv = os.path.join(RESULTS_PHASE1_QC, "thresholdr_fitness_metrics.csv")
    metrics_df.to_csv(output_csv, index=False)
    print(f"\nSaved fitness metrics to {output_csv}")

    # Plot (reuses GMM params — no redundant re-fitting)
    fig = plot_fitness_summary(metrics_df, adata_qc, thr_df, gmm_params)
    output_fig = os.path.join(FIGURES_DIR, "phase1_thresholdr_fitness_summary.png")
    fig.savefig(output_fig, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved fitness summary figure to {output_fig}")

    # Close backed adata_qc
    if hasattr(adata_qc, "file") and adata_qc.file is not None:
        adata_qc.file.close()

    # Summary
    print("\n" + "=" * 60)
    print("FITNESS SUMMARY")
    print("=" * 60)
    for _, row in metrics_df.iterrows():
        rna_str = f"{row['rna_concordance_ratio']:.3f}" if not pd.isna(row['rna_concordance_ratio']) else "N/A"
        print(f"  {row['marker']:8s}  {row['fitness_verdict']:4s}  "
              f"overlap={row['overlap_coeff']:.3f}  "
              f"stain_q={row['stain_quality_score']:.3f}  "
              f"rna_conc={rna_str}")


if __name__ == "__main__":
    main()

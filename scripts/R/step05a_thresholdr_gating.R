#!/usr/bin/env Rscript
# ================================================================
# Step 5a: ThresholdR Binary Gating
# ================================================================
#
# Objective
# ---------
# Fit Gaussian mixture models (GMM) to CLR-normalized ADT
# distributions for the panel's THRESHOLDR_MARKERS. Determine binary
# positive/negative thresholds via the noise mean + 3*sigma rule.
#
# Pipeline Position
# -----------------
# Phase:    1.2 (Tier 1)
# Depends:  Step 1 (QC)
# Feeds:    Step 9 (ADT gating integration)
# Compute:  CPU, <1 GB, <2 min
#
# Inputs
# ------
# - data/processed/adt_counts_for_thresholdr.csv — Raw ADT counts from Step 1
#
# Outputs
# -------
# - data/processed/thresholdr_thresholds.csv        — Per-marker threshold table
# - results/phase1_qc/thresholdr_<marker>_fit.pdf   — Individual GMM fit per marker (intermediate)
# - figures/phase1_thresholdr_fit_<marker>.png      — Individual GMM fit PNG (final)
# - figures/phase1_thresholdr_fit_grid.png          — Combined grid PNG (final)
#
# Config Constants
# ----------------
# - THRESHOLDR_MARKERS — Markers entered into ThresholdR GMM gating
#     (panel YAML).
# - THRESHOLDR_REFERENCE_MARKERS — Subset of THRESHOLDR_MARKERS used to
#     drive the BIC threshold (panel YAML); should be markers known to
#     be cleanly bimodal in the data.
# - THRESHOLDR_MIN_COMPONENT_FRAC (0.05) — Warn if smallest GMM
#     component falls below this fraction (fit instability risk).
# - RANDOM_SEED (42) — Reproducible GMM fits
# ================================================================

# ── 5a.1: Setup ───────────────────────────────────────────────
this_script <- commandArgs(trailingOnly = FALSE)
this_script <- this_script[grep("--file=", this_script)]
this_script <- sub("--file=", "", this_script)
source(file.path(dirname(this_script), "..", "config.R"))

ensure_dirs()

library(ThresholdR)
library(ggplot2)

# ── 5a.2: Load and subset ADT data ───────────────────────────
adt_counts <- read.csv(ADT_COUNTS_FOR_THRESHOLDR, row.names = 1)

# From config.R: panel-defined ThresholdR markers
functional_markers <- THRESHOLDR_MARKERS
adt_data <- adt_counts[, intersect(functional_markers, colnames(adt_counts))]

# ── 5a.3: CLR normalisation ──────────────────────────────────
# ThresholdR expects CLR-transformed input: log(x+1) - row mean
log_adt <- log1p(adt_data)
geom_mean <- rowMeans(log_adt)
clr_adt <- log_adt - geom_mean

# ── 5a.4: ThresholdR GMM workflow ──────────────────────────────
bic_values <- calculateBIC(data = as.data.frame(clr_adt), k = 3)

# Reference markers: panel-defined subset known to be cleanly bimodal;
# drives the adaptive BIC threshold for all THRESHOLDR_MARKERS.
ref_markers <- intersect(THRESHOLDR_REFERENCE_MARKERS, colnames(clr_adt))
if (length(ref_markers) == 0) ref_markers <- colnames(clr_adt)[1]
k_list <- get_k_list(bic_values, markers = ref_markers)

fittings <- fit_k(data = as.data.frame(clr_adt), k_list = k_list, seed = RANDOM_SEED)
thresholds <- get_thresholds(fittings, k_list)

# ── 5a.5: Component-size QC ──────────────────────────────────
# Warn if the smallest GMM component < THRESHOLDR_MIN_COMPONENT_FRAC
# (fit instability risk, especially for markers gating rare populations).
for (marker in names(k_list)) {
    k <- k_list[[marker]][1]
    if (k >= 2) {
        fit <- fittings[[marker]][[1]]  # Index into BIC-selected fit (nested list)
        if (!is.null(fit$lambda)) {
            min_frac <- min(fit$lambda)
            if (min_frac < THRESHOLDR_MIN_COMPONENT_FRAC) {
                cat(sprintf("WARN: %s smallest GMM component = %.1f%% of cells (< %.0f%% threshold). Fit may be unstable.\n",
                            marker, min_frac * 100, THRESHOLDR_MIN_COMPONENT_FRAC * 100))
            } else {
                cat(sprintf("  OK: %s smallest GMM component = %.1f%% of cells.\n",
                            marker, min_frac * 100))
            }
        }
    }
}

# ── 5a.6: Visualise and export ────────────────────────────────
# publish_fittings_thr writes individual per-marker PDFs named
# {marker}_{Bimodal|Trimodal}.pdf into the output directory.
qc_dir <- file.path(PROJECT_ROOT, "results/phase1_qc")
publish_fittings_thr(fittings, k_list, thresholds,
                      path = qc_dir,
                      seed = RANDOM_SEED)

# ── 5a.6-zoom: Overwrite with x-axis-trimmed versions ──────────
# publish_fittings_thr() auto-scales x to max(data), so outlier
# cells create long tails that compress the GMM curves. Overwrite with
# coord_cartesian()-trimmed versions that zoom to the GMM structure
# region without altering the density estimate.
for (marker in names(k_list)) {
    for (j in seq_along(k_list[[marker]])) {
        k <- k_list[[marker]][j]
        if (k < 2) next
        fit <- fittings[[marker]][[j]]
        x <- fit$x;  mu <- fit$mu;  sig <- fit$sigma;  lam <- fit$lambda

        # X-limit: rightmost component mean + 4*sigma, at least data 99th pctl
        max_i <- which.max(mu)
        x_upper <- max(mu[max_i] + 4 * sig[max_i], quantile(x, 0.99))
        x_lower <- min(x)
        pad <- (x_upper - x_lower) * 0.02
        xlim <- c(x_lower - pad, x_upper + pad)

        # Shared theme matching ThresholdR's publish_fittings_thr() style
        base_theme <- theme(
            legend.text = element_text(size = 16),
            legend.title = element_text(size = 18),
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 20),
            plot.title = element_text(size = 24),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

        if (k == 2) {
            df <- data.frame(x = x,
                d1 = dnorm(x, min(mu), sig[which.min(mu)]) * lam[which.min(mu)],
                d2 = dnorm(x, max(mu), sig[which.max(mu)]) * lam[which.max(mu)])
            p <- ggplot(df, aes(x = x)) +
                geom_density(aes(y = after_stat(density)), adjust = 1) +
                ggtitle(paste("Bimodal Fitting for", marker)) +
                geom_line(aes(y = d1, color = "C1")) +
                geom_line(aes(y = d2, color = "C2")) +
                labs(x = marker, y = "Density", colour = "Components") +
                scale_color_manual(values = c("C1" = "blue", "C2" = "red")) +
                geom_vline(xintercept = thresholds[marker, "bi_thr"],
                           linetype = "dashed", linewidth = 0.5) +
                annotate("label", x = thresholds[marker, "bi_thr"], y = Inf, vjust = 2,
                         label = paste("M1+3SD=", round(thresholds[marker, "bi_thr"], 2))) +
                geom_vline(xintercept = thresholds[marker, "bi_cut"],
                           linetype = "dashed", linewidth = 0.5) +
                annotate("label", x = thresholds[marker, "bi_cut"],
                         y = max(density(x)$y) * 0.5, vjust = 2,
                         label = paste("cut=", round(thresholds[marker, "bi_cut"], 2))) +
                base_theme + coord_cartesian(xlim = xlim)
            fit_label <- "Bimodal"

        } else {
            mid_i <- order(mu)[2]  # middle component by rank (safe with tied means)
            df <- data.frame(x = x,
                d1 = dnorm(x, min(mu), sig[which.min(mu)]) * lam[which.min(mu)],
                d2 = dnorm(x, mu[mid_i], sig[mid_i]) * lam[mid_i],
                d3 = dnorm(x, max(mu), sig[which.max(mu)]) * lam[which.max(mu)])
            p <- ggplot(df, aes(x = x)) +
                geom_density(aes(y = after_stat(density)), adjust = 1) +
                ggtitle(paste("Trimodal Fitting for", marker)) +
                geom_line(aes(y = d1, color = "C1")) +
                geom_line(aes(y = d2, color = "C2")) +
                geom_line(aes(y = d3, color = "C3")) +
                labs(x = marker, y = "Density", colour = "Components") +
                scale_color_manual(values = c("C1" = "blue", "C2" = "red", "C3" = "purple")) +
                geom_vline(xintercept = thresholds[marker, "tri_thr1"],
                           linetype = "dashed", linewidth = 0.5) +
                annotate("label", x = thresholds[marker, "tri_thr1"], y = Inf, vjust = 2,
                         label = paste("thr1=", round(thresholds[marker, "tri_thr1"], 2))) +
                geom_vline(xintercept = thresholds[marker, "tri_cut1"],
                           linetype = "dashed", linewidth = 0.5) +
                annotate("label", x = thresholds[marker, "tri_cut1"],
                         y = max(density(x)$y) * 0.8, vjust = 2,
                         label = paste("cut1=", round(thresholds[marker, "tri_cut1"], 2))) +
                geom_vline(xintercept = thresholds[marker, "tri_thr2"],
                           linetype = "dashed", linewidth = 0.5) +
                annotate("label", x = thresholds[marker, "tri_thr2"],
                         y = max(density(x)$y) * 0.6, vjust = 2,
                         label = paste("thr2=", round(thresholds[marker, "tri_thr2"], 2))) +
                geom_vline(xintercept = thresholds[marker, "tri_cut2"],
                           linetype = "dashed", linewidth = 0.5) +
                annotate("label", x = thresholds[marker, "tri_cut2"],
                         y = max(density(x)$y) * 0.4, vjust = 2,
                         label = paste("cut2=", round(thresholds[marker, "tri_cut2"], 2))) +
                base_theme + coord_cartesian(xlim = xlim)
            fit_label <- "Trimodal"
        }

        pdf(file.path(qc_dir, paste0(marker, "_", fit_label, ".pdf")), width = 7, height = 5)
        print(p)
        dev.off()
        cat(sprintf("  Zoomed figure: %s_%s.pdf\n", marker, fit_label))
    }
}

# Remove non-selected fit PDFs so only the BIC-chosen model survives.
# k_list[[marker]][1] is the BIC-selected k for each marker.
for (marker in names(k_list)) {
    selected_k <- k_list[[marker]][1]
    selected_label <- if (selected_k == 2) "Bimodal" else "Trimodal"
    for (alt_label in c("Bimodal", "Trimodal")) {
        if (alt_label != selected_label) {
            alt_pdf <- file.path(qc_dir, paste0(marker, "_", alt_label, ".pdf"))
            if (file.exists(alt_pdf)) {
                file.remove(alt_pdf)
                cat(sprintf("  Removed non-selected: %s_%s.pdf\n", marker, alt_label))
            }
        }
    }
}

# Rename surviving PDFs to thresholdr_<marker>_fit.pdf for clarity.
raw_pdfs <- list.files(qc_dir, pattern = "_(Bimodal|Trimodal)\\.pdf$", full.names = TRUE)
for (raw_pdf in raw_pdfs) {
    marker_name <- sub("_(Bimodal|Trimodal)\\.pdf$", "", basename(raw_pdf))
    dest <- file.path(qc_dir, paste0("thresholdr_", marker_name, "_fit.pdf"))
    file.copy(raw_pdf, dest, overwrite = TRUE)
    cat(sprintf("  Individual figure: %s\n", dest))
}

# Tile per-marker PDFs into one grid via Python helper (R magick unavailable)
combine_script <- file.path(PROJECT_ROOT, "scripts/python/combine_pdf_grid.py")
if (file.exists(combine_script)) {
    py_bin <- PYTHON_BIN
    # Save combined grid as PNG directly to FIGURES_DIR
    rc <- system2(py_bin, args = c(combine_script,
        file.path(PROJECT_ROOT, "results/phase1_qc"),
        shQuote("_(Bimodal|Trimodal)\\.pdf$"),
        file.path("..", "..", "figures", "phase1_thresholdr_fit_grid.png"),
        "--cols", "2"),
        stdout = TRUE, stderr = TRUE)
    cat(paste(rc, collapse = "\n"), "\n")
} else {
    cat("WARN: combine_pdf_grid.py not found, skipping combined figure.\n")
}

# Clean up raw ThresholdR-named files (the renamed copies are the keepers)
for (raw_pdf in raw_pdfs) {
    file.remove(raw_pdf)
}

# ── 5a.7: Convert PDFs to PNG in figures/ directory ────────────
# ThresholdR can only output PDF; convert to PNG for centralized figures.
for (marker in rownames(thresholds)) {
    pdf_src <- file.path(qc_dir, paste0("thresholdr_", marker, "_fit.pdf"))
    png_dest <- file.path(FIGURES_DIR, paste0("phase1_thresholdr_fit_", tolower(marker), ".png"))
    if (file.exists(pdf_src)) {
        gs_rc <- system2("gs", args = c("-dBATCH", "-dNOPAUSE", "-dQUIET",
                               "-sDEVICE=png16m", "-r150",
                               paste0("-sOutputFile=", png_dest), pdf_src))
        if (gs_rc == 0 && file.exists(png_dest)) {
            cat(sprintf("  PNG figure: %s\n", png_dest))
            file.remove(pdf_src)
        } else {
            cat(sprintf("WARN: Ghostscript failed for %s (exit %d). Keeping PDF.\n", marker, gs_rc))
        }
    } else {
        cat(sprintf("WARN: Expected PDF missing: %s\n", pdf_src))
    }
}

k_sel <- sapply(k_list, function(x) x[1])
thr_df <- data.frame(
    marker = rownames(thresholds),
    bi_thr = thresholds$bi_thr,    # noise_mean + 3*sigma (bimodal gate)
    bi_cut = thresholds$bi_cut,    # intersection of two Gaussians (bimodal alternative)
    k_selected = k_sel
)
# For trimodal fits, include tri_thr1 (noise boundary) and tri_thr2 (mid boundary)
if ("tri_thr1" %in% colnames(thresholds)) thr_df$tri_thr1 <- thresholds$tri_thr1
if ("tri_thr2" %in% colnames(thresholds)) thr_df$tri_thr2 <- thresholds$tri_thr2
# threshold_used: the recommended gate value regardless of k.
# Bimodal → bi_thr (noise mean + 3*sigma).
# Trimodal → tri_thr1 (first noise boundary, analogous to bi_thr).
# Unimodal → NA (Step 9 falls back to fg_prob).
thr_df$threshold_used <- ifelse(
    k_sel == 2, thresholds$bi_thr,
    ifelse(k_sel >= 3 & !is.null(thresholds$tri_thr1), thresholds$tri_thr1, NA)
)
write.csv(thr_df,
          THRESHOLDR_THRESHOLDS,
          row.names = FALSE)

cat("ThresholdR gating complete.\n")
print(thr_df)

# ── 5a.8: Validation ──────────────────────────────────────────
validate <- function() {
    results <- list()

    path <- THRESHOLDR_THRESHOLDS
    exists <- file.exists(path)
    results[["thresholdr_thresholds_csv_exists"]] <- list(
        status = if (exists) "PASS" else "FAIL",
        detail = sprintf("File %s: %s", if (exists) "exists" else "MISSING", path),
        expected = "file exists",
        actual = as.character(exists)
    )

    # Individual per-marker PNGs (final figures in centralized directory)
    for (m in THRESHOLDR_MARKERS) {
        png_file <- file.path(FIGURES_DIR, paste0("phase1_thresholdr_fit_", tolower(m), ".png"))
        png_exists <- file.exists(png_file)
        results[[paste0("thresholdr_", m, "_png")]] <- list(
            status = if (png_exists) "PASS" else "WARN",
            detail = sprintf("File %s: %s", if (png_exists) "exists" else "MISSING", png_file),
            expected = "file exists",
            actual = as.character(png_exists)
        )
    }

    # Combined grid PNG (final figure in centralized directory)
    grid_png <- file.path(FIGURES_DIR, "phase1_thresholdr_fit_grid.png")
    grid_exists <- file.exists(grid_png)
    results[["thresholdr_all_fits_png_exists"]] <- list(
        status = if (grid_exists) "PASS" else "WARN",
        detail = sprintf("File %s: %s", if (grid_exists) "exists" else "MISSING", grid_png),
        expected = "file exists",
        actual = as.character(grid_exists)
    )

    # Figures must be at least as recent as the thresholds CSV.
    # Allow 5s tolerance: PNGs are written before the CSV in the same run,
    # so sub-second or ~1s differences are expected and not staleness.
    staleness_tol <- 5  # seconds
    if (exists) {
        csv_mtime <- file.mtime(path)
        for (m in THRESHOLDR_MARKERS) {
            png_file <- file.path(FIGURES_DIR, paste0("phase1_thresholdr_fit_", tolower(m), ".png"))
            if (file.exists(png_file)) {
                png_mtime <- file.mtime(png_file)
                is_current <- as.numeric(difftime(csv_mtime, png_mtime, units = "secs")) <= staleness_tol
                results[[paste0("thresholdr_", m, "_png_current")]] <- list(
                    status = if (is_current) "PASS" else "FAIL",
                    detail = sprintf("%s mtime %s vs CSV mtime %s (tol %ds)",
                                     basename(png_file), format(png_mtime), format(csv_mtime), staleness_tol),
                    expected = "PNG within 5s of CSV mtime",
                    actual = sprintf("PNG %s, CSV %s", format(png_mtime), format(csv_mtime))
                )
            }
        }
        grid_png <- file.path(FIGURES_DIR, "phase1_thresholdr_fit_grid.png")
        if (file.exists(grid_png)) {
            grid_mtime <- file.mtime(grid_png)
            is_current <- as.numeric(difftime(csv_mtime, grid_mtime, units = "secs")) <= staleness_tol
            results[["thresholdr_grid_png_current"]] <- list(
                status = if (is_current) "PASS" else "FAIL",
                detail = sprintf("Grid mtime %s vs CSV mtime %s (tol %ds)",
                                 format(grid_mtime), format(csv_mtime), staleness_tol),
                expected = "grid PNG within 5s of CSV mtime",
                actual = sprintf("grid %s, CSV %s", format(grid_mtime), format(csv_mtime))
            )
        }
    }

    # Expect one row per THRESHOLDR_MARKERS entry (4 markers)
    if (exists) {
        thr_data <- read.csv(path)
        n_rows <- nrow(thr_data)
        n_expected <- length(THRESHOLDR_MARKERS)
        row_ok <- n_rows == n_expected
        results[["thresholdr_row_count"]] <- list(
            status = if (row_ok) "PASS" else "FAIL",
            detail = sprintf("Row count: %d", n_rows),
            expected = sprintf("%d rows (one per ThresholdR marker)", n_expected),
            actual = as.character(n_rows)
        )

        # Columns needed by Step 9
        required_cols <- c("bi_thr", "bi_cut", "threshold_used", "k_selected")
        has_cols <- all(required_cols %in% colnames(thr_data))
        results[["thresholdr_columns"]] <- list(
            status = if (has_cols) "PASS" else "FAIL",
            detail = sprintf("Columns: %s", paste(colnames(thr_data), collapse = ", ")),
            expected = paste(required_cols, collapse = ", "),
            actual = paste(colnames(thr_data), collapse = ", ")
        )

        # threshold_used must be non-NA for every marker with k >= 2
        if (has_cols) {
            na_thresholds <- thr_data$marker[is.na(thr_data$threshold_used) & thr_data$k_selected >= 2]
            results[["thresholdr_threshold_used_complete"]] <- list(
                status = if (length(na_thresholds) == 0) "PASS" else "FAIL",
                detail = sprintf("threshold_used NA for k>=2 markers: %s",
                                 if (length(na_thresholds) == 0) "none" else paste(na_thresholds, collapse = ", ")),
                expected = "no NA threshold_used when k >= 2",
                actual = if (length(na_thresholds) == 0) "all populated" else paste(na_thresholds, collapse = ", ")
            )
        }
    }

    cat("\n", strrep("=", 60), "\n", sep = "")
    cat(sprintf("Validation: Step %s - %s\n", "5a", "ThresholdR Binary Gating"))
    cat(strrep("=", 60), "\n", sep = "")
    for (name in names(results)) {
        r <- results[[name]]
        icon <- switch(r$status, PASS = "[OK]", FAIL = "[!!]", WARN = "[??]")
        cat(sprintf("  %s %s: %s\n", icon, r$status, name))
        if (r$status != "PASS") {
            cat(sprintf("       Expected: %s\n", r$expected))
            cat(sprintf("       Actual:   %s\n", r$actual))
        }
    }
    n_pass <- sum(sapply(results, function(r) r$status == "PASS"))
    n_fail <- sum(sapply(results, function(r) r$status == "FAIL"))
    n_warn <- sum(sapply(results, function(r) r$status == "WARN"))
    cat(sprintf("\nSummary: %d PASS, %d WARN, %d FAIL\n", n_pass, n_warn, n_fail))
    invisible(results)
}

validation_results <- validate()
if (any(sapply(validation_results, function(r) r$status == "FAIL"))) {
    cat("\nValidation FAILED. Troubleshooting:\n")
    cat("If ThresholdR produces unimodal fits for any THRESHOLDR_REFERENCE_MARKERS (expected bimodal), fall back to totalVI foreground probability for that marker (panel YAML: fgprob_override). If GMM fitting fails, check that CLR normalization produces non-degenerate distributions.\n")
    quit(status = 1)
}

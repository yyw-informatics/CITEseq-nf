#!/usr/bin/env Rscript
# ================================================================
# Step 0c: DecontX Ambient RNA Decontamination
# ================================================================
#
# Objective
# ---------
# Apply Bayesian ambient RNA removal (DecontX) to the RNA modality
# of each library.  DecontX models each cell's counts as a mixture
# of native + contamination multinomials and estimates per-cell
# contamination fractions via variational EM.  This is an
# alternative/complement to SoupX (Step 0b) — both address ambient
# transcript leakage into cell-containing droplets.
#
# Pipeline Position
# -----------------
# Phase:    0 (pre-QC)
# Depends:  CellRanger outputs (raw + filtered h5, analysis.tar.gz)
# Feeds:    Step 1 (QC) — reads corrected h5 from DECONTX_DIR
# Compute:  CPU, 8-16 GB per library, ~2 min per library
#
# Usage
# -----
# Per-library:  Rscript step00c_decontx.R --task N        (0-based index)
# Aggregate:    Rscript step00c_decontx.R --aggregate
# Sequential:   Rscript step00c_decontx.R                 (all 16, legacy)
#
# Inputs
# ------
# - 10X_cloud/outs/<sample>/filtered_feature_bc_matrix.h5  — cell calls
# - 10X_cloud/outs/<sample>/raw_feature_bc_matrix.h5       — all droplets
# - 10X_cloud/outs/<sample>/analysis.tar.gz                — CellRanger clusters
#
# Outputs
# -------
# Per-library (--task):
# - data/processed/decontx/<sample>/filtered_feature_bc_matrix.h5
# - data/processed/step00c_temp/decontx_<sample>.csv
#
# Aggregate (--aggregate):
# - results/phase1_qc/decontx/decontx_summary.csv
# - results/phase1_qc/decontx/decontx_contamination_barplot.pdf
#
# Config Constants
# ----------------
# - DECONTX_DIR, DECONTX_RESULTS_DIR    — output directories
# - DECONTX_MAX_ITER (500)              — EM max iterations
# - DECONTX_CONVERGENCE (0.001)         — log-likelihood threshold
# - DECONTX_DELTA (c(10,10))            — Beta prior hyperparameters
#
# Note: pDC anchor detection-rate diagnostic was moved to
# explorations/pdc_anchor_contamination_check.py (post-classification,
# where cell-type labels make it informative).
# ================================================================

# ── 0a: Setup ─────────────────────────────────────────────────
this_script <- commandArgs(trailingOnly = FALSE)
this_script <- this_script[grep("--file=", this_script)]
this_script <- sub("--file=", "", this_script)
source(file.path(dirname(this_script), "..", "config.R"))

ensure_dirs()

suppressPackageStartupMessages({
    library(decontX)
    library(Matrix)
    library(rhdf5)
})

sample_dirs <- SAMPLE_DIRS

TEMP_DIR <- file.path(PROCESSED_DIR, "step00c_temp")


# ── CLI argument parsing ─────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
    parsed <- list(task = NULL, aggregate = FALSE)
    i <- 1
    while (i <= length(args)) {
        if (args[i] == "--task") {
            i <- i + 1
            parsed$task <- as.integer(args[i])
        } else if (args[i] == "--aggregate") {
            parsed$aggregate <- TRUE
        }
        i <- i + 1
    }
    parsed
}

opts <- parse_args(args)


# ── Helper: extract CellRanger clusters from analysis.tar.gz ──
extract_clusters <- function(sample_dir) {
    tar_path <- file.path(DATA_ROOT, sample_dir, "analysis.tar.gz")
    target <- "clustering/gene_expression_graphclust/clusters.csv"
    tmp <- tempdir()
    untar(tar_path, files = target, exdir = tmp)
    clusters <- read.csv(file.path(tmp, target), stringsAsFactors = FALSE)
    unlink(file.path(tmp, "clustering"), recursive = TRUE)
    clusters
}


# ── Helper: write corrected h5 in CellRanger format ──────────
# Copies the original h5 structure, replaces only GEX data rows
# with decontaminated counts, leaves ADT rows untouched.
write_corrected_h5 <- function(original_h5, output_h5, corrected_gex_dgc,
                                n_gex_features) {
    dir.create(dirname(output_h5), recursive = TRUE, showWarnings = FALSE)
    if (file.exists(output_h5)) file.remove(output_h5)

    # Copy entire h5 structure
    file.copy(original_h5, output_h5)

    # Read original sparse matrix components
    orig_data    <- as.numeric(h5read(output_h5, "matrix/data"))
    orig_indices <- as.integer(h5read(output_h5, "matrix/indices"))
    orig_indptr  <- as.integer(h5read(output_h5, "matrix/indptr"))
    orig_shape   <- as.integer(h5read(output_h5, "matrix/shape"))
    n_features <- orig_shape[1]
    n_cells    <- orig_shape[2]

    # Reconstruct original as dgCMatrix (features x cells, CSC)
    original_csc <- sparseMatrix(
        i = orig_indices + 1L,
        p = orig_indptr,
        x = orig_data,
        dims = c(n_features, n_cells),
        repr = "C"
    )

    # Extract ADT rows (everything after GEX)
    adt_csc <- original_csc[(n_gex_features + 1):n_features, , drop = FALSE]

    # Stack: corrected GEX (top) + original ADT (bottom)
    combined <- rbind(corrected_gex_dgc, adt_csc)

    # Force CSC and ensure integer-compatible values
    combined <- as(combined, "CsparseMatrix")

    # Overwrite sparse matrix datasets
    h5delete(output_h5, "matrix/data")
    h5delete(output_h5, "matrix/indices")
    h5delete(output_h5, "matrix/indptr")
    h5delete(output_h5, "matrix/shape")

    h5write(combined@x,                 output_h5, "matrix/data")
    h5write(combined@i,                 output_h5, "matrix/indices")
    h5write(combined@p,                 output_h5, "matrix/indptr")
    h5write(as.integer(combined@Dim),   output_h5, "matrix/shape")

    H5close()
    invisible(NULL)
}


# ── Process a single library ─────────────────────────────────
process_library <- function(sdir) {
    cat(sprintf("\n%s\nProcessing: %s\n%s\n",
                strrep("=", 60), sdir, strrep("=", 60)))

    # 1. Load filtered matrix (RNA + ADT)
    filtered_h5 <- file.path(DATA_ROOT, sdir, "filtered_feature_bc_matrix.h5")
    filtered_all <- Seurat::Read10X_h5(filtered_h5)

    # Split RNA vs ADT
    if (is.list(filtered_all)) {
        filtered_rna <- filtered_all[["Gene Expression"]]
        n_gex <- nrow(filtered_rna)
    } else {
        # Fallback: identify features by type from h5 metadata
        feat_type <- h5read(filtered_h5, "matrix/features/feature_type")
        gex_mask <- feat_type == "Gene Expression"
        n_gex <- sum(gex_mask)
        filtered_rna <- filtered_all[gex_mask, , drop = FALSE]
    }
    n_cells <- ncol(filtered_rna)
    gene_names <- rownames(filtered_rna)
    cat(sprintf("  Filtered RNA: %d genes x %d cells\n", nrow(filtered_rna), n_cells))

    # 2. Load raw matrix for empty droplet background (RNA only)
    raw_h5 <- file.path(DATA_ROOT, sdir, "raw_feature_bc_matrix.h5")
    raw_all <- Seurat::Read10X_h5(raw_h5)
    if (is.list(raw_all)) {
        raw_rna <- raw_all[["Gene Expression"]]
    } else {
        feat_type_raw <- h5read(raw_h5, "matrix/features/feature_type")
        raw_rna <- raw_all[feat_type_raw == "Gene Expression", , drop = FALSE]
    }
    cat(sprintf("  Raw RNA: %d genes x %d droplets\n", nrow(raw_rna), ncol(raw_rna)))

    # Identify empty droplets: in raw but not filtered, 1 < UMI < 100
    cell_bcs <- colnames(filtered_rna)
    empty_bcs <- setdiff(colnames(raw_rna), cell_bcs)
    empty_rna <- raw_rna[, empty_bcs, drop = FALSE]
    empty_umi <- colSums(empty_rna)
    empty_rna <- empty_rna[, empty_umi > 1 & empty_umi < 100, drop = FALSE]
    cat(sprintf("  Empty droplets (1-100 UMI): %d\n", ncol(empty_rna)))

    rm(raw_all, raw_rna); gc(verbose = FALSE)

    # 3. Extract CellRanger graph-based clusters
    clusters_df <- extract_clusters(sdir)
    # Align to filtered barcodes
    clusters_vec <- clusters_df$Cluster[match(cell_bcs, clusters_df$Barcode)]
    # Fill unmatched barcodes with a catch-all cluster
    if (any(is.na(clusters_vec))) {
        max_cl <- max(clusters_vec, na.rm = TRUE)
        clusters_vec[is.na(clusters_vec)] <- max_cl + 1L
        cat(sprintf("  WARNING: %d barcodes unmatched to clusters, assigned to cluster %d\n",
                    sum(is.na(clusters_vec)), max_cl + 1L))
    }
    cat(sprintf("  Clusters: %d groups\n", length(unique(clusters_vec))))

    # 4. Run DecontX
    cat("  Running DecontX (this may take a few minutes)...\n")
    t0 <- proc.time()

    result <- decontX(
        x            = filtered_rna,
        z            = clusters_vec,
        background   = if (ncol(empty_rna) >= 100) empty_rna else NULL,
        maxIter      = DECONTX_MAX_ITER,
        delta        = DECONTX_DELTA,
        convergence  = DECONTX_CONVERGENCE,
        seed         = RANDOM_SEED,
        verbose      = TRUE
    )

    elapsed <- (proc.time() - t0)[["elapsed"]]
    cat(sprintf("  DecontX complete in %.1f s\n", elapsed))

    # Extract decontaminated counts
    decontam_counts <- result$decontXcounts
    contamination   <- result$contamination  # per-cell contamination fraction

    cat(sprintf("  Contamination: median=%.4f, mean=%.4f, range=[%.4f, %.4f]\n",
                median(contamination), mean(contamination),
                min(contamination), max(contamination)))

    # 5. Round to integer and write corrected h5
    decontam_int <- round(decontam_counts)
    decontam_int <- as(decontam_int, "CsparseMatrix")

    output_h5 <- file.path(DECONTX_DIR, sdir, "filtered_feature_bc_matrix.h5")
    write_corrected_h5(filtered_h5, output_h5, decontam_int, n_gex)
    cat(sprintf("  Wrote: %s\n", output_h5))

    # 6. Build diagnostics row
    summary_row <- data.frame(
        sample_id           = sdir,
        n_cells             = n_cells,
        n_clusters          = length(unique(clusters_vec)),
        n_empty_drops       = ncol(empty_rna),
        contamination_median = median(contamination),
        contamination_mean   = mean(contamination),
        contamination_min    = min(contamination),
        contamination_max    = max(contamination),
        elapsed_sec          = elapsed,
        stringsAsFactors     = FALSE
    )

    # Memory cleanup
    rm(filtered_all, filtered_rna, empty_rna, result, decontam_counts,
       decontam_int, contamination, clusters_df, clusters_vec)
    gc(verbose = FALSE)

    summary_row
}


# ── Save per-library CSV to temp dir ─────────────────────────
save_per_library_csv <- function(summary_row) {
    dir.create(TEMP_DIR, recursive = TRUE, showWarnings = FALSE)
    csv_path <- file.path(TEMP_DIR, paste0("decontx_", summary_row$sample_id, ".csv"))
    write.csv(summary_row, csv_path, row.names = FALSE)
    cat(sprintf("  Per-library CSV: %s\n", csv_path))
}


# ── Aggregate: collect per-library CSVs, build summary + plot ─
run_aggregate <- function() {
    cat("DecontX Aggregation Mode\n")
    dir.create(DECONTX_RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

    # Collect per-library CSVs
    csv_files <- list.files(TEMP_DIR, pattern = "^decontx_.*\\.csv$",
                            full.names = TRUE)
    if (length(csv_files) == 0) {
        stop("No per-library CSVs found in ", TEMP_DIR,
             "\n  Run --task for each library first.")
    }
    cat(sprintf("  Found %d per-library CSVs\n", length(csv_files)))

    summary_df <- do.call(rbind, lapply(csv_files, read.csv))
    # Sort to match SAMPLE_DIRS order
    summary_df <- summary_df[match(SAMPLE_DIRS, summary_df$sample_id), ]

    summary_csv <- file.path(DECONTX_RESULTS_DIR, "decontx_summary.csv")
    write.csv(summary_df, summary_csv, row.names = FALSE)
    cat(sprintf("  Summary written: %s\n", summary_csv))

    # Contamination barplot
    pdf_path <- file.path(DECONTX_RESULTS_DIR, "decontx_contamination_barplot.pdf")
    pdf(pdf_path, width = 12, height = 5)
    par(mar = c(10, 4, 3, 1))
    bp <- barplot(summary_df$contamination_median,
                  names.arg = sub("_sample_.*", "", summary_df$sample_id),
                  las = 2, col = "steelblue",
                  main = "DecontX: Median Contamination per Library",
                  ylab = "Median contamination fraction")
    # Add mean as red points
    points(bp, summary_df$contamination_mean, pch = 4, col = "red", cex = 1.2)
    legend("topright", legend = c("Median", "Mean"), pch = c(15, 4),
           col = c("steelblue", "red"), bty = "n")
    dev.off()
    cat(sprintf("  Barplot written: %s\n", pdf_path))

    # Validate
    validation_results <- validate()
    if (any(sapply(validation_results, function(r) r$status == "FAIL"))) {
        cat("\nValidation FAILED. Check outputs and re-run.\n")
        quit(status = 1)
    }

    cat(sprintf("\nAggregation complete. %d libraries summarized.\n",
                nrow(summary_df)))
}


# ── 0d: Validation ───────────────────────────────────────────
validate <- function() {
    results <- list()

    summary_csv <- file.path(DECONTX_RESULTS_DIR, "decontx_summary.csv")
    pdf_path <- file.path(DECONTX_RESULTS_DIR, "decontx_contamination_barplot.pdf")

    # Check 1: summary CSV exists
    csv_exists <- file.exists(summary_csv)
    results[["summary_csv_exists"]] <- list(
        status = if (csv_exists) "PASS" else "FAIL",
        detail = sprintf("%s %s", summary_csv, if (csv_exists) "exists" else "MISSING"),
        expected = "file exists",
        actual = as.character(csv_exists)
    )

    # Check 2: 16 corrected h5 files
    expected_h5 <- file.path(DECONTX_DIR, SAMPLE_DIRS,
                              "filtered_feature_bc_matrix.h5")
    h5_exist <- file.exists(expected_h5)
    n_h5 <- sum(h5_exist)
    results[["corrected_h5_count"]] <- list(
        status = if (n_h5 == 16) "PASS" else if (n_h5 > 0) "WARN" else "FAIL",
        detail = sprintf("%d/16 corrected h5 files exist", n_h5),
        expected = "16",
        actual = as.character(n_h5)
    )

    # Check 3: summary has 16 rows
    if (csv_exists) {
        summ <- read.csv(summary_csv)
        n_rows <- nrow(summ)
        results[["summary_row_count"]] <- list(
            status = if (n_rows == 16) "PASS" else "WARN",
            detail = sprintf("Summary has %d rows", n_rows),
            expected = "16",
            actual = as.character(n_rows)
        )

        # Check 4: contamination in plausible range
        med_range_ok <- all(summ$contamination_median >= 0 &
                            summ$contamination_median <= 1)
        mean_plausible <- mean(summ$contamination_median) > 0.01 &
                          mean(summ$contamination_median) < 0.5
        results[["contamination_range"]] <- list(
            status = if (med_range_ok && mean_plausible) "PASS" else "WARN",
            detail = sprintf("Mean of medians: %.4f", mean(summ$contamination_median)),
            expected = "all in [0,1], mean of medians in (0.01, 0.5)",
            actual = sprintf("range [%.4f, %.4f], mean %.4f",
                             min(summ$contamination_median),
                             max(summ$contamination_median),
                             mean(summ$contamination_median))
        )

    }

    # Check 5: barplot PDF exists
    pdf_exists <- file.exists(pdf_path)
    results[["barplot_exists"]] <- list(
        status = if (pdf_exists) "PASS" else "WARN",
        detail = sprintf("%s %s", pdf_path, if (pdf_exists) "exists" else "MISSING"),
        expected = "file exists",
        actual = as.character(pdf_exists)
    )

    # Print summary
    cat(sprintf("\n%s\n", strrep("=", 60)))
    cat("Validation: Step 0c - DecontX Ambient RNA Decontamination\n")
    cat(sprintf("%s\n", strrep("=", 60)))
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


# ── Main dispatch ─────────────────────────────────────────────
if (!is.null(opts$task)) {
    # Single-library mode (SGE array job)
    task_idx <- opts$task  # 0-based
    if (task_idx < 0 || task_idx >= length(sample_dirs)) {
        cat(sprintf("ERROR: --task %d out of range (0-%d)\n",
                    task_idx, length(sample_dirs) - 1))
        quit(status = 1)
    }
    sdir <- sample_dirs[task_idx + 1]  # R is 1-based
    cat(sprintf("DecontX Single-Library Mode: task=%d, sample=%s\n",
                task_idx, sdir))
    summary_row <- process_library(sdir)
    save_per_library_csv(summary_row)
    cat(sprintf("\nLibrary %d/16 complete.\n", task_idx + 1))

} else if (opts$aggregate) {
    # Aggregation mode
    run_aggregate()

} else {
    # Sequential mode (all 16 libraries, legacy)
    cat("DecontX Ambient RNA Decontamination (sequential mode)\n")
    cat(sprintf("  Libraries: %d\n", length(sample_dirs)))

    all_summaries <- list()
    for (sdir in sample_dirs) {
        summary_row <- process_library(sdir)
        all_summaries[[sdir]] <- summary_row
    }

    # Save summary and plot directly (no temp CSVs needed)
    dir.create(DECONTX_RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
    summary_df <- do.call(rbind, all_summaries)
    summary_csv <- file.path(DECONTX_RESULTS_DIR, "decontx_summary.csv")
    write.csv(summary_df, summary_csv, row.names = FALSE)
    cat(sprintf("\nSummary written: %s\n", summary_csv))

    pdf_path <- file.path(DECONTX_RESULTS_DIR, "decontx_contamination_barplot.pdf")
    pdf(pdf_path, width = 12, height = 5)
    par(mar = c(10, 4, 3, 1))
    bp <- barplot(summary_df$contamination_median,
                  names.arg = sub("_sample_.*", "", summary_df$sample_id),
                  las = 2, col = "steelblue",
                  main = "DecontX: Median Contamination per Library",
                  ylab = "Median contamination fraction")
    points(bp, summary_df$contamination_mean, pch = 4, col = "red", cex = 1.2)
    legend("topright", legend = c("Median", "Mean"), pch = c(15, 4),
           col = c("steelblue", "red"), bty = "n")
    dev.off()
    cat(sprintf("Barplot written: %s\n", pdf_path))

    validation_results <- validate()
    if (any(sapply(validation_results, function(r) r$status == "FAIL"))) {
        cat("\nValidation FAILED. Check outputs and re-run.\n")
        quit(status = 1)
    }
}

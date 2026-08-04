#!/usr/bin/env Rscript
# ================================================================
# Step 2: dsb ADT Normalization (Step I only)
# ================================================================
#
# Objective
# ---------
# Apply dsb Step I ambient correction to the panel's functional ADT
# markers using empty droplets from raw matrices.  Step II (isotype-
# control denoising) is skipped — assumes no isotype controls in panel.
#
# Pipeline Position
# -----------------
# Phase:    1.2 (Tier 1)
# Depends:  Step 1 (QC)
# Feeds:    Diagnostic comparison only; Step 3 applies CLR independently
# Compute:  CPU, 4-8 GB, 5-10 min
#
# Inputs
# ------
# - <DATA_ROOT>/<sample>/filtered_feature_bc_matrix.h5 — CellRanger cell calls
# - <DATA_ROOT>/<sample>/raw_feature_bc_matrix.h5       — All barcodes (for empties)
#
# Outputs
# -------
# - data/processed/dsb_normalized_adt.csv — dsb-normalized ADT values per cell
#
# Config Constants
# ----------------
# - PROJECT_ROOT              — Base directory for all I/O
# - DATA_ROOT                 — CellRanger outs root
# - FUNCTIONAL_MARKERS — ADT panel's functional subset (panel YAML)
# - ADT_RENAME — CellRanger feature name -> canonical short name (panel YAML)
# ================================================================

# ── 1a: Setup ─────────────────────────────────────────────────
this_script <- commandArgs(trailingOnly = FALSE)
this_script <- this_script[grep("--file=", this_script)]
this_script <- sub("--file=", "", this_script)
source(file.path(dirname(this_script), "..", "config.R"))

ensure_dirs()

library(dsb)
library(Matrix)

# ── 1b: Per-sample dsb normalization loop ─────────────────────
sample_dirs <- SAMPLE_DIRS

functional_markers <- FUNCTIONAL_MARKERS
# Panel-defined CellRanger -> canonical rename happens below before
# marker intersection (see ADT_RENAME in config.R).

all_dsb_results <- list()

for (sdir in sample_dirs) {
    cat(sprintf("Processing %s...\n", sdir))

    filtered_path <- file.path(DATA_ROOT, sdir, "filtered_feature_bc_matrix.h5")
    filtered <- Seurat::Read10X_h5(filtered_path)

    # CellRanger CITE-seq returns a named list; extract ADT modality
    if (is.list(filtered)) {
        cell_adt <- as.matrix(filtered[["Antibody Capture"]])
    } else {
        adt_idx <- grep("TotalSeq", rownames(filtered), ignore.case = TRUE)
        cell_adt <- as.matrix(filtered[adt_idx, ])
    }

    # Raw matrix provides empty-droplet background for dsb ambient estimation
    raw_path <- file.path(DATA_ROOT, sdir, "raw_feature_bc_matrix.h5")
    raw <- Seurat::Read10X_h5(raw_path)
    if (is.list(raw)) {
        raw_adt <- as.matrix(raw[["Antibody Capture"]])
    } else {
        adt_idx <- grep("TotalSeq", rownames(raw), ignore.case = TRUE)
        raw_adt <- as.matrix(raw[adt_idx, ])
    }

    # Empty drops: in raw but not filtered, 1 < UMI < 100 (dsb convention)
    cell_bcs <- colnames(cell_adt)
    empty_bcs <- setdiff(colnames(raw_adt), cell_bcs)
    empty_adt <- raw_adt[, empty_bcs]
    empty_umi <- colSums(empty_adt)
    empty_adt <- empty_adt[, empty_umi > 1 & empty_umi < 100]

    # Apply panel-defined CellRanger -> canonical name remap (ADT_RENAME from config.R).
    .remap <- function(rn) {
        idx <- match(rn, names(ADT_RENAME))
        rn[!is.na(idx)] <- unname(ADT_RENAME[idx[!is.na(idx)]])
        rn
    }
    rownames(cell_adt)  <- .remap(rownames(cell_adt))
    rownames(empty_adt) <- .remap(rownames(empty_adt))

    cat(sprintf("  Cells: %d, Empty droplets: %d\n", ncol(cell_adt), ncol(empty_adt)))

    if (ncol(empty_adt) < 100) {
        cat("  WARNING: Too few empty droplets, skipping dsb for this sample\n")
        next
    }

    avail_markers <- intersect(rownames(cell_adt), functional_markers)
    if (length(avail_markers) == 0) {
        # Fallback: partial matching against any panel marker name handles
        # residual naming mismatches (built from ADT_RENAME keys and values).
        .panel_pattern <- paste(c(names(ADT_RENAME), unname(ADT_RENAME)),
                                collapse = "|")
        avail_markers <- rownames(cell_adt)[grep(.panel_pattern,
                                                   rownames(cell_adt))]
    }
    # Ensure markers exist in both cell and empty matrices
    avail_markers <- intersect(avail_markers, rownames(empty_adt))

    if (length(avail_markers) == 0) {
        cat("  WARNING: No matching markers found, skipping\n")
        next
    }

    # Step I only: ambient correction via empty-drop background subtraction.
    # Step II (isotype-control denoising) skipped.
    dsb_norm <- DSBNormalizeProtein(
        cell_protein_matrix = cell_adt[avail_markers, , drop = FALSE],
        empty_drop_matrix = empty_adt[avail_markers, , drop = FALSE],
        denoise.counts = FALSE,
        use.isotype.control = FALSE
    )

    dsb_df <- as.data.frame(t(dsb_norm))
    dsb_df$barcode <- rownames(dsb_df)
    dsb_df$sample_id <- sdir
    all_dsb_results[[sdir]] <- dsb_df
}

# ── 1c: Combine and export ────────────────────────────────────
if (length(all_dsb_results) == 0) {
    cat("\nERROR: No samples produced dsb results. Exiting.\n")
    quit(status = 1)
}
dsb_combined <- do.call(rbind, all_dsb_results)
write.csv(dsb_combined, DSB_NORMALIZED_ADT, row.names = FALSE)

cat(sprintf("\ndsb normalization complete: %d cells, %d markers\n",
            nrow(dsb_combined), ncol(dsb_combined) - 2))

# ── 1d: Validation ────────────────────────────────────────────
validate <- function() {
    results <- list()

    path <- DSB_NORMALIZED_ADT
    exists <- file.exists(path)
    results[["dsb_normalized_adt_csv_exists"]] <- list(
        status = if (exists) "PASS" else "FAIL",
        detail = sprintf("File %s: %s", if (exists) "exists" else "MISSING", path),
        expected = "file exists",
        actual = as.character(exists)
    )

    # Shape check: one row per normalized cell, one column per functional
    # panel marker plus barcode and sample_id.  Cohort size is study-
    # specific, so only require a non-empty table.
    n_expected_cols <- length(FUNCTIONAL_MARKERS) + 2  # + barcode + sample_id
    if (exists) {
        dsb_data <- read.csv(path)
        n_rows <- nrow(dsb_data)
        n_cols <- ncol(dsb_data)
        row_ok <- n_rows > 0
        results[["dsb_row_count"]] <- list(
            status = if (row_ok) "PASS" else "WARN",
            detail = sprintf("Row count: %d", n_rows),
            expected = ">0 rows (cells)",
            actual = as.character(n_rows)
        )
        col_ok <- n_cols == n_expected_cols
        results[["dsb_col_count"]] <- list(
            status = if (col_ok) "PASS" else "WARN",
            detail = sprintf("Column count: %d (includes barcode + sample_id)", n_cols),
            expected = sprintf("%d columns (%d panel markers + barcode + sample_id)",
                               n_expected_cols, length(FUNCTIONAL_MARKERS)),
            actual = as.character(n_cols)
        )
    }

    cat("\n", strrep("=", 60), "\n", sep = "")
    cat(sprintf("Validation: Step %d - %s\n", 2, "dsb ADT Normalization"))
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
    cat("If empty droplets are insufficient for any sample, fall back to CLR normalization (Step 3) for that sample. If dsb produces all-zero output for a marker, check marker name matching between Cell Ranger output and the functional_markers vector.\n")
    quit(status = 1)
}

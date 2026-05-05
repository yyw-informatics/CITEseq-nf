#!/usr/bin/env Rscript
# Step 1b: scDblFinder Doublet Detection
# =======================================
#
# Objective
# ---------
# Run scDblFinder on each QC-passed library to detect doublets using a
# gradient-boosted classifier with cluster-aware artificial doublet generation.
# Replaces Scrublet as the primary doublet caller (Scrublet retained as fallback).
#
# Pipeline Position
# -----------------
# Phase:    1.1b (Tier 1)
# Depends:  Step 01 (QC) -- uses keep_barcodes CSVs + decontaminated h5
# Feeds:    Step 02 (DSB normalization)
# Compute:  CPU, ~8 GB per library (array job), ~32 GB aggregate
#
# Usage
# -----
# Per-library:  Rscript step01b_scdblfinder.R --task N [--branch decontx|soupx]
# Aggregate:    Rscript step01b_scdblfinder.R --aggregate [--branch decontx|soupx]

# ── Source config ────────────────────────────────────────────────
this_script <- commandArgs(trailingOnly = FALSE)
this_script <- this_script[grep("--file=", this_script)]
this_script <- sub("--file=", "", this_script)
source(file.path(dirname(this_script), "..", "config.R"))
ensure_dirs()

suppressPackageStartupMessages({
  library(rhdf5)
  library(Matrix)
  library(SingleCellExperiment)
  library(scDblFinder)
})

# ── CLI argument parsing ────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  parsed <- list(task = NULL, aggregate = FALSE, branch = "")
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--task") {
      i <- i + 1
      parsed$task <- as.integer(args[i])
    } else if (args[i] == "--aggregate") {
      parsed$aggregate <- TRUE
    } else if (args[i] == "--branch") {
      i <- i + 1
      parsed$branch <- args[i]
    }
    i <- i + 1
  }
  parsed
}

opts <- parse_args(args)

# ── Resolve branch-specific paths ───────────────────────────────
get_h5_dir <- function(sample_name, branch) {
  if (branch == "soupx") {
    file.path(PROCESSED_DIR, "soupx_parallel", sample_name)
  } else {
    # Default and "decontx" both read from decontx directory
    file.path(PROCESSED_DIR, "decontx", sample_name)
  }
}

branch_suffix <- function(branch) {
  if (!is.null(branch) && branch != "") paste0("_", branch) else ""
}

# ── Read 10x h5 as SingleCellExperiment (RNA only) ──────────────
read_10x_h5_as_sce <- function(h5_path) {
  barcodes    <- h5read(h5_path, "/matrix/barcodes")
  gene_names  <- h5read(h5_path, "/matrix/features/name")
  gene_ids    <- h5read(h5_path, "/matrix/features/id")
  feat_types  <- h5read(h5_path, "/matrix/features/feature_type")
  data_vals   <- h5read(h5_path, "/matrix/data")
  indices     <- h5read(h5_path, "/matrix/indices")
  indptr      <- h5read(h5_path, "/matrix/indptr")
  mat_shape   <- h5read(h5_path, "/matrix/shape")

  # Build sparse matrix (CSC format: features x barcodes)
  mat <- sparseMatrix(
    i = indices + 1L,
    p = indptr,
    x = as.numeric(data_vals),
    dims = mat_shape,
    dimnames = list(gene_names, barcodes)
  )

  # Subset to Gene Expression only (exclude Antibody Capture)
  rna_idx <- which(feat_types == "Gene Expression")
  mat <- mat[rna_idx, , drop = FALSE]
  gene_names <- gene_names[rna_idx]
  gene_ids   <- gene_ids[rna_idx]

  # Handle duplicate gene names
  if (any(duplicated(gene_names))) {
    gene_names <- make.unique(gene_names)
    rownames(mat) <- gene_names
  }

  sce <- SingleCellExperiment(
    assays = list(counts = mat),
    rowData = data.frame(ID = gene_ids, Symbol = gene_names, stringsAsFactors = FALSE),
    colData = data.frame(Barcode = barcodes, row.names = barcodes, stringsAsFactors = FALSE)
  )
  sce
}

# ── Per-library task mode ───────────────────────────────────────
run_task <- function(task_index, branch) {
  suffix <- branch_suffix(branch)
  sample_name <- SAMPLE_DIRS[task_index + 1]  # 0-based to 1-based
  cat(sprintf("=== Task %d: %s (branch=%s) ===\n", task_index, sample_name, branch))

  # Resolve input paths
  h5_dir <- get_h5_dir(sample_name, branch)
  h5_path <- file.path(h5_dir, "filtered_feature_bc_matrix.h5")
  # Try branch-suffixed keep_barcodes first, fall back to unsuffixed
  # (step01 may have been run without --branch)
  keep_path <- file.path(STEP01_TEMP_DIR,
                         paste0(sample_name, "_keep_barcodes", suffix, ".csv"))
  if (!file.exists(keep_path)) {
    keep_path <- file.path(STEP01_TEMP_DIR,
                           paste0(sample_name, "_keep_barcodes.csv"))
  }

  if (!file.exists(h5_path)) stop("h5 not found: ", h5_path)
  if (!file.exists(keep_path)) stop("keep_barcodes not found: ", keep_path)

  # Read keep barcodes
  keep_barcodes <- readLines(keep_path)
  keep_barcodes <- trimws(keep_barcodes)
  keep_barcodes <- keep_barcodes[nchar(keep_barcodes) > 0]
  cat(sprintf("  Keep barcodes: %d\n", length(keep_barcodes)))

  # Load h5 as SCE
  cat(sprintf("  Loading: %s\n", h5_path))
  sce <- read_10x_h5_as_sce(h5_path)
  cat(sprintf("  Loaded: %d cells, %d genes\n", ncol(sce), nrow(sce)))

  # Subset to QC-passing barcodes
  shared <- intersect(colnames(sce), keep_barcodes)
  if (length(shared) < length(keep_barcodes)) {
    cat(sprintf("  WARNING: %d/%d keep_barcodes found in h5 (%.1f%%)\n",
                length(shared), length(keep_barcodes),
                100 * length(shared) / length(keep_barcodes)))
  }
  sce <- sce[, shared]
  cat(sprintf("  After subsetting: %d cells\n", ncol(sce)))

  # Run scDblFinder
  set.seed(RANDOM_SEED)
  cat("  Running scDblFinder...\n")
  sce <- scDblFinder(
    sce,
    clusters = SCDBLFINDER_CLUSTERS,
    dbr = SCDBLFINDER_DBR,
    verbose = TRUE
  )

  # Extract results
  results <- data.frame(
    barcode = colnames(sce),
    scDblFinder.score = sce$scDblFinder.score,
    scDblFinder.class = sce$scDblFinder.class,
    stringsAsFactors = FALSE
  )

  n_doublets <- sum(results$scDblFinder.class == "doublet")
  cat(sprintf("  scDblFinder: %d doublets (%.1f%%)\n",
              n_doublets, 100 * n_doublets / nrow(results)))

  # Write output
  dir.create(SCDBLFINDER_TEMP_DIR, recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(SCDBLFINDER_TEMP_DIR,
                        paste0(sample_name, "_doublets", suffix, ".csv"))
  write.csv(results, out_path, row.names = FALSE)
  cat(sprintf("  Saved: %s\n", out_path))
}

# ── Aggregate mode ──────────────────────────────────────────────
run_aggregate <- function(branch) {
  suffix <- branch_suffix(branch)
  cat(sprintf("=== Aggregating scDblFinder results (branch=%s) ===\n", branch))

  # Load all 16 doublet CSVs
  all_doublets <- list()
  for (sample_name in SAMPLE_DIRS) {
    csv_path <- file.path(SCDBLFINDER_TEMP_DIR,
                          paste0(sample_name, "_doublets", suffix, ".csv"))
    if (!file.exists(csv_path)) {
      stop("Missing doublet CSV: ", csv_path, "\n  Run --task for this library first.")
    }
    df <- read.csv(csv_path, stringsAsFactors = FALSE)
    df$sample_id <- sample_name
    all_doublets[[sample_name]] <- df
    n_dbl <- sum(df$scDblFinder.class == "doublet")
    cat(sprintf("  %s: %d cells, %d doublets (%.1f%%)\n",
                sample_name, nrow(df), n_dbl, 100 * n_dbl / nrow(df)))
  }
  doublet_df <- do.call(rbind, all_doublets)
  rownames(doublet_df) <- NULL
  cat(sprintf("  Total: %d cells, %d doublets (%.1f%%)\n",
              nrow(doublet_df), sum(doublet_df$scDblFinder.class == "doublet"),
              100 * sum(doublet_df$scDblFinder.class == "doublet") / nrow(doublet_df)))

  # Load adata_qc.h5ad via Python helper (anndata -> CSV bridge)
  adata_qc_path <- file.path(PROCESSED_DIR, paste0("adata_qc", suffix, ".h5ad"))
  adata_adt_path <- file.path(PROCESSED_DIR, paste0("adata_adt_raw", suffix, ".h5ad"))

  if (!file.exists(adata_qc_path)) stop("adata_qc not found: ", adata_qc_path)
  if (!file.exists(adata_adt_path)) stop("adata_adt_raw not found: ", adata_adt_path)

  # Use Python to update the h5ad files
  # Build script via paste to avoid sprintf/f-string escaping issues
  all_csv_path_py <- file.path(SCDBLFINDER_TEMP_DIR, paste0("all_doublets", suffix, ".csv"))
  comp_csv_name <- paste0("doublet_comparison", suffix, ".csv")

  python_lines <- c(
    'import sys, os',
    paste0('sys.path.insert(0, os.path.join("', PROJECT_ROOT, '", "scripts"))'),
    'import scanpy as sc',
    'import pandas as pd',
    'import scipy.sparse as sp',
    '',
    '# Load doublet results',
    paste0('doublet_df = pd.read_csv("', all_csv_path_py, '")'),
    '# Reconstruct full barcode: raw_barcode + "-" + sample_id (matches sc.concat index_unique="-")',
    'doublet_df["barcode"] = doublet_df["barcode"] + "-" + doublet_df["sample_id"]',
    'doublet_df = doublet_df.set_index("barcode")',
    '',
    '# Load adata_qc',
    paste0('adata = sc.read_h5ad("', adata_qc_path, '")'),
    'print(f"Loaded adata_qc: {adata.n_obs} cells")',
    '',
    '# Rename existing Scrublet columns',
    'if "predicted_doublet" in adata.obs.columns:',
    '    adata.obs["scrublet_doublet"] = adata.obs["predicted_doublet"]',
    'if "doublet_score" in adata.obs.columns:',
    '    adata.obs["scrublet_score"] = adata.obs["doublet_score"]',
    '',
    '# Match doublet results to adata barcodes',
    'shared = adata.obs_names.intersection(doublet_df.index)',
    'print(f"Matched {len(shared)}/{adata.n_obs} barcodes to scDblFinder results")',
    '',
    '# Add scDblFinder columns as primary doublet columns',
    'adata.obs["scDblFinder_score"] = 0.0',
    'adata.obs["scDblFinder_class"] = "singlet"',
    'adata.obs.loc[shared, "scDblFinder_score"] = doublet_df.loc[shared, "scDblFinder.score"].values',
    'adata.obs.loc[shared, "scDblFinder_class"] = doublet_df.loc[shared, "scDblFinder.class"].values',
    '',
    '# Set primary doublet columns from scDblFinder',
    'adata.obs["predicted_doublet"] = adata.obs["scDblFinder_class"] == "doublet"',
    'adata.obs["doublet_score"] = adata.obs["scDblFinder_score"]',
    '',
    '# Count pre-removal',
    'n_scrublet = int(adata.obs["scrublet_doublet"].sum()) if "scrublet_doublet" in adata.obs.columns else 0',
    'n_scdblfinder = int(adata.obs["predicted_doublet"].sum())',
    'print(f"Scrublet doublets: {n_scrublet}")',
    'print(f"scDblFinder doublets: {n_scdblfinder}")',
    '',
    '# Build comparison CSV before removal',
    'comparison_rows = []',
    'for sample in adata.obs["sample_id"].unique():',
    '    mask = adata.obs["sample_id"] == sample',
    '    n_cells = int(mask.sum())',
    '    n_scrub = int(adata.obs.loc[mask, "scrublet_doublet"].sum()) if "scrublet_doublet" in adata.obs.columns else 0',
    '    n_scdb = int(adata.obs.loc[mask, "predicted_doublet"].sum())',
    '    comparison_rows.append({',
    '        "sample": sample, "n_cells": n_cells,',
    '        "scrublet_doublets": n_scrub, "scrublet_pct": round(100*n_scrub/n_cells, 2),',
    '        "scDblFinder_doublets": n_scdb, "scDblFinder_pct": round(100*n_scdb/n_cells, 2),',
    '    })',
    'comp_df = pd.DataFrame(comparison_rows)',
    paste0('comp_path = os.path.join("', RESULTS_PHASE1_QC, '", "', comp_csv_name, '")'),
    'comp_df.to_csv(comp_path, index=False)',
    'print(f"Comparison CSV: {comp_path}")',
    'print(comp_df.to_string(index=False))',
    '',
    '# Remove scDblFinder-flagged doublets',
    'adata = adata[~adata.obs["predicted_doublet"]].copy()',
    'print(f"After doublet removal: {adata.n_obs} cells")',
    '',
    '# Save updated adata_qc',
    paste0('adata.write_h5ad("', adata_qc_path, '")'),
    paste0('print("Saved: ', adata_qc_path, '")'),
    '',
    '# Sync ADT to final barcodes',
    paste0('adata_adt = sc.read_h5ad("', adata_adt_path, '")'),
    'adata_adt = adata_adt[adata.obs_names].copy()',
    paste0('adata_adt.write_h5ad("', adata_adt_path, '")'),
    'print(f"Saved synced ADT: {adata_adt.n_obs} cells")'
  )
  python_script <- paste(python_lines, collapse = "\n")

  # First save concatenated doublet CSV for the Python helper
  all_csv_path <- file.path(SCDBLFINDER_TEMP_DIR, paste0("all_doublets", suffix, ".csv"))
  write.csv(doublet_df, all_csv_path, row.names = FALSE)
  cat(sprintf("  Wrote combined CSV: %s\n", all_csv_path))

  # Write Python helper script and execute
  helper_path <- file.path(SCDBLFINDER_TEMP_DIR, paste0("aggregate_helper", suffix, ".py"))
  writeLines(python_script, helper_path)
  cat("  Running Python aggregate helper...\n")
  exit_code <- system2(PYTHON_BIN, helper_path, stdout = "", stderr = "")
  if (exit_code != 0) stop("Python aggregate helper failed with exit code ", exit_code)

  cat("=== Aggregation complete ===\n")
}

# ── Main dispatch ───────────────────────────────────────────────
if (!is.null(opts$task)) {
  run_task(opts$task, opts$branch)
} else if (opts$aggregate) {
  run_aggregate(opts$branch)
} else {
  cat("Usage:\n")
  cat("  Per-library:  Rscript step01b_scdblfinder.R --task N [--branch decontx|soupx]\n")
  cat("  Aggregate:    Rscript step01b_scdblfinder.R --aggregate [--branch decontx|soupx]\n")
  quit(status = 1)
}

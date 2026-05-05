# citeseq-nf — generic CITE-seq preprocessing + embedding pipeline config (R).
#
# Scope: generic + panel-specific only. Study overlays (sample metadata,
# forced HVG gene lists, Tier-2 classification thresholds) live in the
# calling repo, NOT here.
#
# Environment-variable contract — see scripts/config.py for full reference.
#   SAGE_PROJECT_ROOT           REQUIRED. Errors if unset.
#   SAGE_PROCESSED_DIR          Override for {PROJECT_ROOT}/data/processed
#   SAGE_MODELS_DIR             Override for {PROJECT_ROOT}/models
#   SAGE_RESULTS_DIR            Override for {PROJECT_ROOT}/results
#   SAGE_FIGURES_DIR            Override for {PROJECT_ROOT}/figures
#   SAGE_EXTERNAL_DIR           Override for {PROJECT_ROOT}/data/external
#   SAGE_DATA_ROOT              Override for {PROJECT_ROOT}/10X_cloud/outs
#   SAGE_SAMPLE_METADATA_JSON   JSON file path; maps library_id -> metadata
#   SAGE_FORCED_HVGS_FILE       Text file path; one gene symbol per line
#   SAGE_PANEL_YAML             REQUIRED. Panel YAML path (see
#                               panels/example_panel.yaml for the schema).

# ═════════════════════════════════════════════════════════════════════
# Project root — REQUIRED from env
# ═════════════════════════════════════════════════════════════════════
PROJECT_ROOT <- Sys.getenv("SAGE_PROJECT_ROOT", unset = NA_character_)
if (is.na(PROJECT_ROOT) || PROJECT_ROOT == "") {
    stop(paste(
        "citeseq-nf requires SAGE_PROJECT_ROOT to be set. Nextflow exports",
        "it automatically via beforeScript; for direct script invocation,",
        "set SAGE_PROJECT_ROOT to the calling study's project root."
    ))
}

# ═════════════════════════════════════════════════════════════════════
# Base paths (env-overridable, defaults derived from PROJECT_ROOT)
# ═════════════════════════════════════════════════════════════════════
DATA_ROOT     <- Sys.getenv("SAGE_DATA_ROOT",     file.path(PROJECT_ROOT, "10X_cloud/outs"))
PROCESSED_DIR <- Sys.getenv("SAGE_PROCESSED_DIR", file.path(PROJECT_ROOT, "data/processed"))
EXTERNAL_DIR  <- Sys.getenv("SAGE_EXTERNAL_DIR",  file.path(PROJECT_ROOT, "data/external"))
MODELS_DIR    <- Sys.getenv("SAGE_MODELS_DIR",    file.path(PROJECT_ROOT, "models"))
RESULTS_DIR   <- Sys.getenv("SAGE_RESULTS_DIR",   file.path(PROJECT_ROOT, "results"))
FIGURES_DIR   <- Sys.getenv("SAGE_FIGURES_DIR",   file.path(PROJECT_ROOT, "figures"))

RESULTS_PHASE1_QC        <- file.path(RESULTS_DIR, "phase1_qc")
RESULTS_PHASE2_EMBEDDING <- file.path(RESULTS_DIR, "phase2_embedding")

# ═════════════════════════════════════════════════════════════════════
# Sample metadata + sample dirs — loaded from env-pointed JSON file
# ═════════════════════════════════════════════════════════════════════
.sample_meta_path <- Sys.getenv("SAGE_SAMPLE_METADATA_JSON", unset = "")
if (nzchar(.sample_meta_path) && file.exists(.sample_meta_path)) {
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
        stop("jsonlite is required to load SAGE_SAMPLE_METADATA_JSON")
    }
    SAMPLE_METADATA <- jsonlite::fromJSON(
        .sample_meta_path, simplifyVector = FALSE
    )
} else {
    SAMPLE_METADATA <- list()
}
SAMPLE_DIRS <- names(SAMPLE_METADATA)
ALL_SAMPLES <- names(SAMPLE_METADATA)

# ═════════════════════════════════════════════════════════════════════
# DecontX ambient RNA decontamination (Step 0c)
# ═════════════════════════════════════════════════════════════════════
DECONTX_DIR         <- file.path(PROCESSED_DIR, "decontx")
DECONTX_RESULTS_DIR <- file.path(RESULTS_PHASE1_QC, "decontx")
DECONTX_MAX_ITER    <- 500L
DECONTX_CONVERGENCE <- 0.001
DECONTX_DELTA       <- c(10, 10)

# ═════════════════════════════════════════════════════════════════════
# Canonical artifact names (Python ↔ R bridge files)
# ═════════════════════════════════════════════════════════════════════
ADT_COUNTS_FOR_THRESHOLDR <- file.path(PROCESSED_DIR, "adt_counts_for_thresholdr.csv")
DSB_NORMALIZED_ADT        <- file.path(PROCESSED_DIR, "dsb_normalized_adt.csv")
THRESHOLDR_THRESHOLDS     <- file.path(PROCESSED_DIR, "thresholdr_thresholds.csv")

# ═════════════════════════════════════════════════════════════════════
# scDblFinder doublet detection (Step 1b)
# ═════════════════════════════════════════════════════════════════════
SCDBLFINDER_ENABLED  <- TRUE
SCDBLFINDER_CLUSTERS <- TRUE
SCDBLFINDER_DBR      <- NULL   # auto: n^2/10^5
SCDBLFINDER_TEMP_DIR <- file.path(PROCESSED_DIR, "step01b_temp")
STEP01_TEMP_DIR      <- file.path(PROCESSED_DIR, "step01_temp")

# ═════════════════════════════════════════════════════════════════════
# HVG — forced-include gene list (loaded from env-pointed text file)
# ═════════════════════════════════════════════════════════════════════
.forced_hvgs_path <- Sys.getenv("SAGE_FORCED_HVGS_FILE", unset = "")
if (nzchar(.forced_hvgs_path) && file.exists(.forced_hvgs_path)) {
    .lines <- readLines(.forced_hvgs_path, warn = FALSE)
    .lines <- trimws(.lines)
    FORCED_HVGS <- .lines[nzchar(.lines) & !startsWith(.lines, "#")]
} else {
    FORCED_HVGS <- character(0)
}

# ═════════════════════════════════════════════════════════════════════
# ADT panel — loaded from required YAML (SAGE_PANEL_YAML)
# ═════════════════════════════════════════════════════════════════════
# SAGE_PANEL_YAML must point at a panel YAML file (see
# panels/example_panel.yaml for the schema). Schema-validated at load.
.panel_yaml_path <- Sys.getenv("SAGE_PANEL_YAML", unset = "")
if (!nzchar(.panel_yaml_path)) {
    stop(paste(
        "citeseq-nf requires SAGE_PANEL_YAML to be set.",
        "Pass --panel /path/to/your_panel.yaml to nextflow, or export",
        "SAGE_PANEL_YAML when invoking R scripts directly.",
        "See panels/example_panel.yaml for the schema."
    ))
}
if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("R package 'yaml' is required to load SAGE_PANEL_YAML.")
}
.panel <- yaml::read_yaml(.panel_yaml_path)
for (.k in c("cellranger_names", "rename_map", "functional",
             "thresholdr", "thresholdr_reference")) {
    if (is.null(.panel[[.k]])) {
        stop(sprintf("Panel YAML %s missing required key: %s",
                     .panel_yaml_path, .k))
    }
}
if (length(.panel$thresholdr_reference) == 0L) {
    stop(sprintf(
        "Panel YAML %s: thresholdr_reference must be non-empty.",
        .panel_yaml_path))
}

ADT_RENAME                   <- unlist(.panel$rename_map)
FUNCTIONAL_MARKERS           <- unlist(.panel$functional)
THRESHOLDR_MARKERS           <- unlist(.panel$thresholdr)
THRESHOLDR_REFERENCE_MARKERS <- unlist(.panel$thresholdr_reference)

# ThresholdR GMM stability warning threshold (not panel-specific)
THRESHOLDR_MIN_COMPONENT_FRAC <- 0.05

# ═════════════════════════════════════════════════════════════════════
# Python binary (for R scripts that spawn Python subprocesses)
# Set by nextflow.config beforeScript via params.python.
# ═════════════════════════════════════════════════════════════════════
PYTHON_BIN <- Sys.getenv("SAGE_PYTHON_BIN", unset = NA_character_)
if (is.na(PYTHON_BIN) || PYTHON_BIN == "") {
    stop("SAGE_PYTHON_BIN env var not set. Pass --python (or --conda_env) to nextflow, or export SAGE_PYTHON_BIN before running R directly.")
}

# ═════════════════════════════════════════════════════════════════════
# General infrastructure
# ═════════════════════════════════════════════════════════════════════
RANDOM_SEED <- 42

# ═════════════════════════════════════════════════════════════════════
# Utility
# ═════════════════════════════════════════════════════════════════════
ensure_dirs <- function() {
    dirs <- c(
        PROCESSED_DIR, EXTERNAL_DIR, MODELS_DIR,
        RESULTS_DIR, FIGURES_DIR,
        RESULTS_PHASE1_QC, RESULTS_PHASE2_EMBEDDING,
        DECONTX_DIR, DECONTX_RESULTS_DIR,
        SCDBLFINDER_TEMP_DIR, STEP01_TEMP_DIR
    )
    for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

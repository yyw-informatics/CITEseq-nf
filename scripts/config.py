#!/usr/bin/env python
"""
citeseq-nf — generic CITE-seq preprocessing + embedding pipeline config.

Scope
-----
This file ships with citeseq-nf and covers only the GENERIC surface: paths,
canonical artifact names, and decontamination / QC / HVG / embedding
hyperparameters. The antibody panel is provided by the caller via a panel
YAML (--panel / CITESEQ_PANEL_YAML); see panels/example_panel.yaml for the
schema.

Study overlays (sample metadata, forced HVG gene lists, classification
thresholds, hypothesis gene sets, manuscript figure constants, etc.) live
in the calling repo, NOT here.

Environment-variable contract
-----------------------------
Nextflow's beforeScript (nextflow.config) exports these on every process;
direct-script callers set them by hand.

  CITESEQ_PROJECT_ROOT          REQUIRED. Caller study's project root.
                             Errors if unset.

  CITESEQ_DATA_ROOT             REQUIRED. CellRanger input root: one
                             subdirectory per library, each holding that
                             library's raw_feature_bc_matrix.h5 and
                             filtered_feature_bc_matrix.h5. Errors if unset.

  CITESEQ_PROCESSED_DIR         Override for {PROJECT_ROOT}/data/processed
  CITESEQ_MODELS_DIR            Override for {PROJECT_ROOT}/models
  CITESEQ_RESULTS_DIR           Override for {PROJECT_ROOT}/results
  CITESEQ_FIGURES_DIR           Override for {PROJECT_ROOT}/figures
  CITESEQ_EXTERNAL_DIR          Override for {PROJECT_ROOT}/data/external

  CITESEQ_SAMPLE_METADATA_JSON  Path to JSON file. Maps library_id -> dict of
                             {donor, age, age_group, timepoint, replicate}.
                             If unset or missing → SAMPLE_METADATA = {}.

  CITESEQ_FORCED_HVGS_FILE      Path to text file. One gene symbol per line;
                             blank lines and '#' comments ignored.
                             If unset or missing → FORCED_HVGS = [].

  CITESEQ_QC_MAD_N_MADS         Override for QC_MAD_N_MADS (default 5).

  CITESEQ_PANEL_YAML            REQUIRED. Path to an ADT panel YAML (see
                             panels/example_panel.yaml for the schema).
                             Defines the ADT_* constants used by every
                             downstream step. Errors if unset.

  CITESEQ_SOUPX_ANCHOR_GENES    Comma-separated list of anchor genes for the
                             SoupX post-correction validation step. Empty
                             default → validation step is skipped.

  CITESEQ_HVG_SENTINEL_GENES    Comma-separated list of sentinel genes the
                             forced-HVG injection should retain. Empty
                             default → sentinel check is skipped.

  CITESEQ_BIOLOGICAL_SEPARATION_MARKER
                             Canonical ADT name used for the step08
                             biological-silhouette metric. Unset → block
                             is skipped.

  CITESEQ_SAMPLE_GROUP_FIELD    Key in SAMPLE_METADATA entries used to order
                             and color samples in step01 QC figures (e.g.,
                             "age_group"). Unset → insertion order,
                             single color.
"""
import json
import os

import scanpy as sc
import yaml

# ═════════════════════════════════════════════════════════════════════
# Project root — REQUIRED from env; citeseq-nf does not embed a default
# ═════════════════════════════════════════════════════════════════════
PROJECT_ROOT = os.environ.get("CITESEQ_PROJECT_ROOT")
if PROJECT_ROOT is None:
    raise RuntimeError(
        "citeseq-nf requires CITESEQ_PROJECT_ROOT to be set. "
        "Nextflow exports it automatically via beforeScript; for direct "
        "script invocation, set CITESEQ_PROJECT_ROOT to the calling study's "
        "project root (the dir holding data/processed/, models/, etc.)."
    )

# ═════════════════════════════════════════════════════════════════════
# Input data root — REQUIRED from env; citeseq-nf does not embed a default
# ═════════════════════════════════════════════════════════════════════
DATA_ROOT = os.environ.get("CITESEQ_DATA_ROOT")
if not DATA_ROOT:
    raise RuntimeError(
        "citeseq-nf requires CITESEQ_DATA_ROOT to be set. "
        "Pass --data_root /path/to/your/cellranger_outs to nextflow, or "
        "export CITESEQ_DATA_ROOT when invoking scripts directly. It must "
        "point at a directory holding one subdirectory per library, each "
        "containing that library's CellRanger raw_feature_bc_matrix.h5 and "
        "filtered_feature_bc_matrix.h5."
    )

# ═════════════════════════════════════════════════════════════════════
# Base paths (env-overridable, defaults derived from PROJECT_ROOT)
# ═════════════════════════════════════════════════════════════════════
PROCESSED_DIR = os.environ.get("CITESEQ_PROCESSED_DIR", os.path.join(PROJECT_ROOT, "data/processed"))
EXTERNAL_DIR  = os.environ.get("CITESEQ_EXTERNAL_DIR",  os.path.join(PROJECT_ROOT, "data/external"))
MODELS_DIR    = os.environ.get("CITESEQ_MODELS_DIR",    os.path.join(PROJECT_ROOT, "models"))
RESULTS_DIR   = os.environ.get("CITESEQ_RESULTS_DIR",   os.path.join(PROJECT_ROOT, "results"))
FIGURES_DIR   = os.environ.get("CITESEQ_FIGURES_DIR",   os.path.join(PROJECT_ROOT, "figures"))

TOTALVI_MODEL_DIR = os.path.join(MODELS_DIR, "totalvi_model")
SCVI_MODEL_DIR    = os.path.join(MODELS_DIR, "scvi_model")

# Generic phase-result subdirs (Tier 1 only; Tier 2+ subdirs stay in study overlay)
RESULTS_PHASE1_QC        = os.path.join(RESULTS_DIR, "phase1_qc")
RESULTS_PHASE2_EMBEDDING = os.path.join(RESULTS_DIR, "phase2_embedding")

# ═════════════════════════════════════════════════════════════════════
# Canonical artifact names (generic — filenames are the handoff contract)
# ═════════════════════════════════════════════════════════════════════
ADATA_QC           = os.path.join(PROCESSED_DIR, "adata_qc.h5ad")
ADATA_ADT_RAW      = os.path.join(PROCESSED_DIR, "adata_adt_raw.h5ad")
ADATA_QC_HVG       = os.path.join(PROCESSED_DIR, "adata_qc_hvg.h5ad")
ADATA_NORMALIZED   = os.path.join(PROCESSED_DIR, "adata_normalized.h5ad")
ADATA_QC_CLR       = os.path.join(PROCESSED_DIR, "adata_qc_clr.h5ad")
ADATA_QC_EMBEDDED  = os.path.join(PROCESSED_DIR, "adata_qc_embedded.h5ad")  # ← terminal artifact
ADATA_TOTALVI      = os.path.join(PROCESSED_DIR, "adata_totalvi.h5ad")
ADATA_SCVI         = os.path.join(PROCESSED_DIR, "adata_scvi.h5ad")

# Python ↔ R bridge files (generic)
ADT_COUNTS_FOR_THRESHOLDR = os.path.join(PROCESSED_DIR, "adt_counts_for_thresholdr.csv")
DSB_NORMALIZED_ADT        = os.path.join(PROCESSED_DIR, "dsb_normalized_adt.csv")
THRESHOLDR_THRESHOLDS     = os.path.join(PROCESSED_DIR, "thresholdr_thresholds.csv")

# ═════════════════════════════════════════════════════════════════════
# Sample metadata — loaded from env-pointed JSON file (empty default)
# ═════════════════════════════════════════════════════════════════════
# Format:
#   {
#     "Lib01": {"donor": "Donor_A", "age": 25,
#         "age_group": "young", "timepoint": "d0", "replicate": 1},
#     ...
#   }
_sample_meta_path = os.environ.get("CITESEQ_SAMPLE_METADATA_JSON")
if _sample_meta_path and os.path.exists(_sample_meta_path):
    with open(_sample_meta_path) as _f:
        SAMPLE_METADATA = json.load(_f)
else:
    SAMPLE_METADATA = {}

# CITESEQ_N_LIBRARIES (--n_libraries) runs a pilot on the first N libraries.
# Truncate the metadata itself, not just the Nextflow fan-out, so that every
# step agrees on which libraries exist. Steps that iterate SAMPLE_DIRS (DSB)
# or SAMPLE_METADATA (the QC and doublet aggregates) would otherwise still
# expect all of them and pull in stale files from an earlier full run.
_n_libraries = os.environ.get("CITESEQ_N_LIBRARIES", "").strip()
if _n_libraries and SAMPLE_METADATA:
    try:
        _n = int(_n_libraries)
    except ValueError:
        raise RuntimeError(
            f"CITESEQ_N_LIBRARIES must be a whole number, got: {_n_libraries!r}"
        )
    if _n < 1:
        raise RuntimeError(
            f"CITESEQ_N_LIBRARIES must be 1 or greater, got: {_n}"
        )
    if _n > len(SAMPLE_METADATA):
        raise RuntimeError(
            f"CITESEQ_N_LIBRARIES is {_n}, but the sample metadata lists only "
            f"{len(SAMPLE_METADATA)} librar"
            f"{'y' if len(SAMPLE_METADATA) == 1 else 'ies'}. "
            "Lower --n_libraries, or point --sample_metadata_json at a file "
            "with more entries."
        )
    SAMPLE_METADATA = dict(list(SAMPLE_METADATA.items())[:_n])

ALL_SAMPLES = list(SAMPLE_METADATA.keys())
SAMPLE_DIRS = list(SAMPLE_METADATA.keys())


def _env_csv_list(var):
    """Parse env var as comma-separated list. Empty/unset → []."""
    raw = os.environ.get(var, "")
    return [x.strip() for x in raw.split(",") if x.strip()]


# Optional sample-grouping field for step01 QC figures. When set to a key
# present in SAMPLE_METADATA entries, step01 sorts samples by that field's
# value and colors them from matplotlib tab10.
SAMPLE_GROUP_FIELD = os.environ.get("CITESEQ_SAMPLE_GROUP_FIELD") or None

# ═════════════════════════════════════════════════════════════════════
# SoupX ambient RNA decontamination (Step 0b)
# ═════════════════════════════════════════════════════════════════════
SOUPX_DIR           = os.path.join(PROCESSED_DIR, "soupx")
SOUPX_PARALLEL_DIR  = os.path.join(PROCESSED_DIR, "soupx_parallel")
SOUPX_RESULTS_DIR   = os.path.join(RESULTS_PHASE1_QC, "soupx")
SOUPX_METHOD        = "subtraction"
SOUPX_ROUND_TO_INT  = True

# Anchor genes for the post-SoupX validation pass (genes expected to be
# exclusively ambient in the input data; study-specific). Empty → skip
# validation.
SOUPX_ANCHOR_GENES = _env_csv_list("CITESEQ_SOUPX_ANCHOR_GENES")

# ═════════════════════════════════════════════════════════════════════
# DecontX ambient RNA decontamination (Step 0c)
# ═════════════════════════════════════════════════════════════════════
DECONTX_DIR         = os.path.join(PROCESSED_DIR, "decontx")
DECONTX_RESULTS_DIR = os.path.join(RESULTS_PHASE1_QC, "decontx")
DECONTX_MAX_ITER    = 500
DECONTX_CONVERGENCE = 0.001
DECONTX_DELTA       = [10, 10]

# ═════════════════════════════════════════════════════════════════════
# QC + doublet detection (Steps 1, 1b)
# ═════════════════════════════════════════════════════════════════════
QC_MAD_N_MADS              = int(os.environ.get("CITESEQ_QC_MAD_N_MADS", "5"))
QC_MIN_GENES_ABSOLUTE      = 100
QC_MIN_UMI_ABSOLUTE        = 200
QC_MAX_MT_PCT              = 20
QC_MIN_CELLS_PER_GENE      = 3
SCRUBLET_EXPECTED_DOUBLET_RATE = 0.06

SCDBLFINDER_ENABLED  = True
SCDBLFINDER_CLUSTERS = True
SCDBLFINDER_DBR      = None  # auto: n^2/10^5
SCDBLFINDER_TEMP_DIR = os.path.join(PROCESSED_DIR, "step01b_temp")

# ═════════════════════════════════════════════════════════════════════
# HVG selection (Step 4)
# ═════════════════════════════════════════════════════════════════════
N_HVG = 4000  # Seurat-v3 flavor, batch_key="sample_id"

# FORCED_HVGS — loaded from env-pointed text file (empty default).
# Format: one gene symbol per line; blank lines and '#'-comments ignored.
# Study overlay provides its own gene list; citeseq-nf is study-neutral.
_forced_hvgs_path = os.environ.get("CITESEQ_FORCED_HVGS_FILE")
if _forced_hvgs_path and os.path.exists(_forced_hvgs_path):
    with open(_forced_hvgs_path) as _f:
        FORCED_HVGS = [
            line.strip() for line in _f
            if line.strip() and not line.lstrip().startswith("#")
        ]
else:
    FORCED_HVGS = []

# Sentinel genes for the forced-HVG injection check. Empty → skip check.
HVG_SENTINEL_GENES = _env_csv_list("CITESEQ_HVG_SENTINEL_GENES")

# ═════════════════════════════════════════════════════════════════════
# ThresholdR + KDE (ADT gating — Step 5a/5b/5c)
# ═════════════════════════════════════════════════════════════════════
# Stability
THRESHOLDR_MIN_COMPONENT_FRAC       = 0.05
# Fitness PASS thresholds (all must hold)
THRESHOLDR_OVERLAP_THRESHOLD        = 0.15
THRESHOLDR_STAIN_QUALITY_THRESHOLD  = 2.0
THRESHOLDR_RNA_CONCORDANCE_THRESHOLD = 1.5
# Fitness FAIL triggers (any one triggers FAIL)
THRESHOLDR_OVERLAP_FAIL             = 0.30
THRESHOLDR_STAIN_QUALITY_FAIL       = 1.0
THRESHOLDR_RNA_CONCORDANCE_FAIL     = 1.2
THRESHOLDR_BIC_MARGIN_FAIL          = 0.0
# totalVI fg_prob fallback when ThresholdR unimodal or overridden
ADT_FG_PROB_FALLBACK                = 0.5

# ═════════════════════════════════════════════════════════════════════
# ADT panel — loaded from required YAML (CITESEQ_PANEL_YAML)
# ═════════════════════════════════════════════════════════════════════
# CITESEQ_PANEL_YAML must point at a panel YAML file (see
# panels/example_panel.yaml for the schema). The ADT_* constants below
# are loaded from it and schema-validated at import time.
def _load_panel_yaml(path):
    """Load and validate a panel YAML. Returns the parsed dict."""
    with open(path) as f:
        panel = yaml.safe_load(f)
    if not isinstance(panel, dict):
        raise ValueError(f"Panel YAML {path} did not parse to a mapping.")

    required = ["cellranger_names", "rename_map", "functional",
                "thresholdr", "thresholdr_reference"]
    missing = [k for k in required if k not in panel]
    if missing:
        raise ValueError(
            f"Panel YAML {path} missing required keys: {missing}"
        )

    cellranger_names = list(panel["cellranger_names"])
    rename_map       = dict(panel["rename_map"])
    functional       = list(panel["functional"])
    thresholdr       = list(panel["thresholdr"])
    thresholdr_ref   = list(panel["thresholdr_reference"])
    fgprob_override  = list(panel.get("fgprob_override", []))
    kde_validate     = list(panel.get("kde_validate", []))
    rna_concordance  = dict(panel.get("rna_concordance", {}))

    # rename_map keys ⊆ cellranger_names
    stray_keys = [k for k in rename_map if k not in cellranger_names]
    if stray_keys:
        raise ValueError(
            f"Panel YAML {path}: rename_map keys not in cellranger_names: "
            f"{stray_keys}"
        )

    canonical = set(rename_map.values())

    # functional/thresholdr/thresholdr_ref/fgprob_override/kde_validate
    # must appear as a value in rename_map (i.e., be canonical names).
    for section, values in [
        ("functional", functional),
        ("thresholdr", thresholdr),
        ("thresholdr_reference", thresholdr_ref),
        ("fgprob_override", fgprob_override),
        ("kde_validate", kde_validate),
    ]:
        stray = [v for v in values if v not in canonical]
        if stray:
            raise ValueError(
                f"Panel YAML {path}: {section} contains markers not in "
                f"rename_map values: {stray}"
            )

    # thresholdr_reference ⊆ thresholdr, and non-empty
    if not thresholdr_ref:
        raise ValueError(
            f"Panel YAML {path}: thresholdr_reference must be non-empty."
        )
    stray_ref = [m for m in thresholdr_ref if m not in thresholdr]
    if stray_ref:
        raise ValueError(
            f"Panel YAML {path}: thresholdr_reference not subset of "
            f"thresholdr: {stray_ref}"
        )

    # rna_concordance keys ⊆ functional
    stray_rna = [k for k in rna_concordance if k not in functional]
    if stray_rna:
        raise ValueError(
            f"Panel YAML {path}: rna_concordance keys not in functional: "
            f"{stray_rna}"
        )

    return {
        "cellranger_names": cellranger_names,
        "rename_map": rename_map,
        "functional": functional,
        "thresholdr": thresholdr,
        "thresholdr_reference": thresholdr_ref,
        "fgprob_override": fgprob_override,
        "kde_validate": kde_validate,
        "rna_concordance": rna_concordance,
    }


_panel_yaml_path = os.environ.get("CITESEQ_PANEL_YAML")
if not _panel_yaml_path:
    raise RuntimeError(
        "citeseq-nf requires CITESEQ_PANEL_YAML to be set. "
        "Pass --panel /path/to/your_panel.yaml to nextflow, or export "
        "CITESEQ_PANEL_YAML when invoking scripts directly. "
        "See panels/example_panel.yaml for the schema."
    )
_panel = _load_panel_yaml(_panel_yaml_path)
ADT_CELLRANGER_NAMES             = _panel["cellranger_names"]
ADT_RENAME_MAP                   = _panel["rename_map"]
ADT_FUNCTIONAL_MARKERS           = _panel["functional"]
ADT_THRESHOLDR_MARKERS           = _panel["thresholdr"]
ADT_THRESHOLDR_REFERENCE_MARKERS = _panel["thresholdr_reference"]
ADT_FGPROB_OVERRIDE_MARKERS      = _panel["fgprob_override"]
KDE_MARKERS                      = _panel["kde_validate"]
ADT_RNA_CONCORDANCE              = _panel["rna_concordance"]
# Canonical short names downstream = functional markers in panel order.
ADT_CANONICAL_NAMES              = list(_panel["functional"])

# ═════════════════════════════════════════════════════════════════════
# Embedding — totalVI hyperparameters (Step 6 — step06_totalvi.py)
# ═════════════════════════════════════════════════════════════════════
# Targets recon plateau diagnosed in earlier iterations: LR scheduler on
# reconstruction_loss_validation, lr_threshold aligned with early-stopping
# min_delta, faster LR cuts, and pro_recons_weight upweight for the 4-protein
# loss branch vs ~4K genes.
TOTALVI_N_LATENT              = 15
TOTALVI_GENE_LIKELIHOOD       = "nb"
TOTALVI_MAX_EPOCHS            = 100
TOTALVI_EARLY_STOPPING        = True
TOTALVI_EARLY_STOPPING_PATIENCE = 16
TOTALVI_EARLY_STOPPING_MIN_DELTA = 1.0
TOTALVI_EARLY_STOPPING_MONITOR   = "elbo_validation"
TOTALVI_LR                    = 2e-3
TOTALVI_EMPIRICAL_PRIOR       = False
TOTALVI_BATCH_SIZE_TRAIN      = 128
TOTALVI_BATCH_SIZE_INFERENCE  = 64
TOTALVI_KL_WARMUP_EPOCHS      = 30
TOTALVI_LR_PATIENCE           = 4
TOTALVI_LR_FACTOR             = 0.5
TOTALVI_LR_MIN                = 1e-4
TOTALVI_LR_THRESHOLD          = 1.0
TOTALVI_LR_SCHEDULER_METRIC   = "reconstruction_loss_validation"
TOTALVI_PRO_RECONS_WEIGHT     = 2.0

# ═════════════════════════════════════════════════════════════════════
# Embedding — scVI hyperparameters (Step 7 — step07_scvi.py)
# ═════════════════════════════════════════════════════════════════════
# Wider architecture (n_hidden=256, n_layers=2) and full KL ramp aligned to
# max_epochs so the latent space is actually regularised before early stopping.
SCVI_N_LATENT              = 15
SCVI_N_HIDDEN              = 256
SCVI_N_LAYERS              = 2
SCVI_GENE_LIKELIHOOD       = "nb"
SCVI_MAX_EPOCHS            = 400
SCVI_EARLY_STOPPING        = True
SCVI_EARLY_STOPPING_PATIENCE = 30
SCVI_EARLY_STOPPING_MIN_DELTA = 0.0
SCVI_EARLY_STOPPING_MONITOR   = "elbo_validation"
SCVI_LR                    = 4e-3
SCVI_KL_WARMUP             = 400

# Neighbor graph / UMAP / Leiden (Step 8)
N_NEIGHBORS       = 15
UMAP_MIN_DIST     = 0.3
LEIDEN_RESOLUTION = 0.8

# Canonical ADT name used as a biological-separation proxy in step08's
# silhouette metric. When None, that block is skipped (batch-mixing
# silhouette still runs).
BIOLOGICAL_SEPARATION_MARKER = os.environ.get("CITESEQ_BIOLOGICAL_SEPARATION_MARKER") or None

# ═════════════════════════════════════════════════════════════════════
# General infrastructure
# ═════════════════════════════════════════════════════════════════════
RANDOM_SEED      = 42
SCANPY_VERBOSITY = 2

sc.settings.verbosity = SCANPY_VERBOSITY
sc.settings.figdir    = FIGURES_DIR


# ═════════════════════════════════════════════════════════════════════
# Utility functions
# ═════════════════════════════════════════════════════════════════════
def ensure_dirs():
    """Create all citeseq-nf output directories (idempotent)."""
    for d in [
        PROCESSED_DIR, EXTERNAL_DIR, MODELS_DIR, TOTALVI_MODEL_DIR,
        SCVI_MODEL_DIR, RESULTS_DIR, FIGURES_DIR,
        RESULTS_PHASE1_QC, RESULTS_PHASE2_EMBEDDING,
        SOUPX_DIR, SOUPX_PARALLEL_DIR, SOUPX_RESULTS_DIR,
        DECONTX_DIR, DECONTX_RESULTS_DIR,
        SCDBLFINDER_TEMP_DIR,
    ]:
        os.makedirs(d, exist_ok=True)


def get_sample_path(sample_id, matrix_type="filtered"):
    """Get path to a sample's feature-barcode matrix h5 file."""
    return os.path.join(DATA_ROOT, sample_id, f"{matrix_type}_feature_bc_matrix.h5")


def branched_paths(branch):
    """Return canonical or branch-suffixed artifact paths.

    The pipeline can run on SoupX- or DecontX-decontaminated input (or both
    side-by-side). When ``branch`` is set, intermediate artifacts get
    ``_{branch}`` spliced in before the file extension (or appended to
    directory names) so the two branches don't stomp each other. When
    ``branch`` is None, the canonical unsuffixed paths are returned.
    """
    canonical = {
        "ADATA_QC":          ADATA_QC,
        "ADATA_ADT_RAW":     ADATA_ADT_RAW,
        "ADATA_QC_CLR":      ADATA_QC_CLR,
        "ADATA_QC_HVG":      ADATA_QC_HVG,
        "ADATA_NORMALIZED":  ADATA_NORMALIZED,
        "ADATA_TOTALVI":     ADATA_TOTALVI,
        "ADATA_SCVI":        ADATA_SCVI,
        "ADATA_QC_EMBEDDED": ADATA_QC_EMBEDDED,
        "TOTALVI_MODEL_DIR": TOTALVI_MODEL_DIR,
        "SCVI_MODEL_DIR":    SCVI_MODEL_DIR,
    }
    if not branch:
        return canonical
    out = {}
    for key, path in canonical.items():
        root, ext = os.path.splitext(path)
        out[key] = f"{root}_{branch}{ext}"
    return out


if __name__ == "__main__":
    ensure_dirs()
    print(f"PROJECT_ROOT:         {PROJECT_ROOT}")
    print(f"DATA_ROOT:            {DATA_ROOT}")
    print(f"PROCESSED_DIR:        {PROCESSED_DIR}")
    print(f"ADATA_QC_EMBEDDED:    {ADATA_QC_EMBEDDED}")
    print(f"Samples loaded:       {len(SAMPLE_METADATA)}")
    print(f"Functional ADT:       {ADT_FUNCTIONAL_MARKERS}")
    print(f"FORCED_HVGS loaded:   {len(FORCED_HVGS)}")
    print(f"N_HVG:                {N_HVG}")
    print("Directories created.")

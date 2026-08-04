#!/usr/bin/env Rscript
# =====================================================================
# citeseq-nf - install the R packages that conda does not provide
# =====================================================================
#
# What this is for
# ----------------
# The pipeline runs on R 4.4.3 inside a conda environment. Conda covers
# most of the R stack (Seurat, Matrix, ggplot2, mclust, scDblFinder,
# SingleCellExperiment, ...). A handful of packages were added on top of
# the conda environment by hand. This script reinstalls exactly those, so
# the R side can be rebuilt without guessing.
#
# How to run it
# -------------
#   Rscript env/r_extra.R
#
# Run it from the repository root, with the project's conda environment
# already active. The script installs into the first writable library in
# .libPaths(), which is the conda environment's own R library when the
# environment is active.
#
# It is safe to run more than once. A package that is already installed
# at the pinned version is left alone.
#
# It stops with a non-zero exit status if anything fails to install or
# fails to load afterwards.
#
# READ THIS ABOUT ThresholdR
# --------------------------
# ThresholdR is not on CRAN or Bioconductor. The build that produced the
# published results was a LOCALLY PATCHED copy of upstream, and R kept no
# record of where it came from. A rebuilt source copy of that exact build
# is vendored in this repository at env/vendor/ThresholdR_0.1.0.tar.gz.
#
# By default this script installs ThresholdR from upstream GitHub, and
# then checks whether the result matches the vendored build. If it does
# not, it prints a loud warning telling you which numbers can change. To
# install the vendored copy instead, run:
#
#   THRESHOLDR_SOURCE=vendor Rscript env/r_extra.R
#
# See env/vendor/README.md for the full story.
# =====================================================================

options(warn = 1)
options(repos = c(CRAN = "https://cloud.r-project.org"))

BIOC_RELEASE   <- "3.20"          # Bioconductor release matching R 4.4.x
R_VERSION_USED <- "4.4.3"         # R the pipeline was validated on

# ---------------------------------------------------------------------
# Locate the repository root so every path below stays relative to it
# ---------------------------------------------------------------------
find_repo_root <- function() {
    env <- Sys.getenv("CITESEQ_NF_ROOT", "")
    if (nzchar(env)) return(normalizePath(env, mustWork = TRUE))
    args <- commandArgs(trailingOnly = FALSE)
    hit <- grep("^--file=", args, value = TRUE)
    if (length(hit)) {
        # this file lives at <repo root>/env/r_extra.R
        return(normalizePath(file.path(dirname(sub("^--file=", "", hit[1])), ".."),
                             mustWork = TRUE))
    }
    # sourced interactively: fall back to the working directory
    normalizePath(".", mustWork = TRUE)
}

REPO_ROOT     <- find_repo_root()
VENDOR_DIR    <- file.path(REPO_ROOT, "env", "vendor")
THRESHOLDR_TARBALL <- file.path(VENDOR_DIR, "ThresholdR_0.1.0.tar.gz")

say <- function(...) cat(sprintf(...), "\n", sep = "")
rule <- function() cat(strrep("-", 70), "\n")

say("citeseq-nf: installing the R packages conda does not provide")
rule()
say("R in use          : %s", getRversion())
say("Library target    : %s", .libPaths()[1])
say("Bioconductor pin  : %s", BIOC_RELEASE)
rule()

if (getRversion() < "4.4.0" || getRversion() >= "4.5.0") {
    warning(sprintf(paste0("This pipeline was validated on R %s. You are on R %s. ",
                           "Bioconductor %s is built for R 4.4.x, so some ",
                           "installs below may fail or give different versions."),
                    R_VERSION_USED, getRversion(), BIOC_RELEASE),
            immediate. = TRUE)
}

failures <- character(0)
record_failure <- function(pkg, e) {
    failures <<- c(failures, sprintf("%s: %s", pkg, conditionMessage(e)))
    message(sprintf("  FAILED  %s -- %s", pkg, conditionMessage(e)))
}

installed_version <- function(pkg) {
    v <- tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) NA_character_)
    v
}

# ---------------------------------------------------------------------
# 0. Bootstrap the two installer helpers
# ---------------------------------------------------------------------
say("[0/5] bootstrap installers")
for (helper in c("BiocManager", "remotes")) {
    if (!requireNamespace(helper, quietly = TRUE)) {
        say("  installing %s", helper)
        utils::install.packages(helper)
    }
    if (!requireNamespace(helper, quietly = TRUE))
        stop(sprintf("cannot install %s, and nothing else can proceed without it", helper))
    say("  %-14s %s", helper, installed_version(helper))
}

# ---------------------------------------------------------------------
# 1. Bioconductor packages
#    The version on each line is what the pipeline actually ran with.
#    A Bioconductor release channel keeps serving patch updates inside
#    the release, so asking for release 3.20 does NOT pin the exact
#    version. Expect small patch differences and treat them as normal.
#    scDblFinder is deliberately absent here: the recorded 1.23.4 is a
#    Bioconductor devel build and release 3.20 serves 1.20.2, so this
#    channel can never produce it. env/environment.yml pins
#    bioconductor-scdblfinder=1.23.4 through conda instead.
# ---------------------------------------------------------------------
bioc_pkgs <- c(
    "Rhdf5lib"             = "1.28.0",   # HDF5 C library used by rhdf5
    "rhdf5filters"         = "1.18.1",   # compression filters for rhdf5
    "rhdf5"                = "2.50.2",   # reads the 10x .h5 matrices
    "SingleCellExperiment" = "1.28.0",   # container used by scDblFinder
    "celda"                = "1.22.1",   # hard dependency of decontX
    "decontX"              = "1.4.1"     # ambient RNA removal, step 0c
)

say("")
say("[1/5] Bioconductor %s packages", BIOC_RELEASE)
for (pkg in names(bioc_pkgs)) {
    want <- bioc_pkgs[[pkg]]
    have <- installed_version(pkg)
    if (!is.na(have) && have == want) {
        say("  %-22s %-10s already installed", pkg, have)
        next
    }
    say("  %-22s installing (want %s, have %s)", pkg, want,
        ifelse(is.na(have), "nothing", have))
    tryCatch(
        BiocManager::install(pkg, version = BIOC_RELEASE, update = FALSE,
                             ask = FALSE, checkBuilt = FALSE),
        error = function(e) record_failure(pkg, e))
    got <- installed_version(pkg)
    if (!is.na(got) && got != want)
        say(paste0("  note: %s installed as %s, the pipeline ran %s. ",
                   "Bioconductor serves patch updates inside a release, ",
                   "so a small difference here is expected."),
            pkg, got, want)
}

# ---------------------------------------------------------------------
# 2. CRAN packages, pinned to the exact recorded versions
#    mixtools, mclust and foreach are required by ThresholdR at run time.
#    ThresholdR calls require() for them from inside its functions, so a
#    missing one only shows up when step 5a is already running.
# ---------------------------------------------------------------------
cran_pkgs <- c(
    "mixtools" = "2.0.0.1",   # normalmixEM, the GMM fit behind step 5a
    "mclust"   = "6.1.2",     # BIC model selection in ThresholdR
    "foreach"  = "1.5.2",     # loop backend used by ThresholdR::calculateBIC
    "ggplot2"  = "4.0.2"      # step 5a fit plots
)

say("")
say("[2/5] CRAN packages (pinned versions)")
for (pkg in names(cran_pkgs)) {
    want <- cran_pkgs[[pkg]]
    have <- installed_version(pkg)
    if (!is.na(have) && have == want) {
        say("  %-22s %-10s already installed", pkg, have)
        next
    }
    say("  %-22s installing %s (have %s)", pkg, want,
        ifelse(is.na(have), "nothing", have))
    tryCatch(
        remotes::install_version(pkg, version = want, upgrade = "never",
                                 dependencies = TRUE),
        error = function(e) record_failure(pkg, e))
}

# ---------------------------------------------------------------------
# 3. dsb - CITE-seq ADT normalisation, step 2
#    The installed copy recorded RemoteType: github, MattPM/dsb at the
#    commit below, so that pin is exact. dsb 2.0.1 is also on CRAN and is
#    used as the fallback if GitHub is unreachable.
# ---------------------------------------------------------------------
DSB_VERSION <- "2.0.1"
DSB_REPO    <- "MattPM/dsb"
DSB_SHA     <- "aab21279fe8f96ab480cc27503716410ecca1310"

say("")
say("[3/5] dsb")
have <- installed_version("dsb")
if (!is.na(have) && have == DSB_VERSION) {
    say("  %-22s %-10s already installed", "dsb", have)
} else {
    say("  installing dsb from GitHub %s@%s", DSB_REPO, DSB_SHA)
    ok <- tryCatch({
        remotes::install_github(paste0(DSB_REPO, "@", DSB_SHA),
                                upgrade = "never", dependencies = TRUE)
        TRUE
    }, error = function(e) {
        message("  GitHub install of dsb failed: ", conditionMessage(e))
        FALSE
    })
    if (!ok) {
        say("  falling back to CRAN dsb %s", DSB_VERSION)
        tryCatch(
            remotes::install_version("dsb", version = DSB_VERSION,
                                     upgrade = "never", dependencies = TRUE),
            error = function(e) record_failure("dsb", e))
    }
}

# ---------------------------------------------------------------------
# 4. ThresholdR - GMM gating of ADT markers, step 5a
#
#    Not on CRAN or Bioconductor. Upstream is MDMotlagh/ThresholdR (MIT).
#    The installed copy in the original environment carried no Repository,
#    RemoteUrl or RemoteSha field, so R recorded no origin for it.
#
#    Set THRESHOLDR_SOURCE=vendor to skip GitHub and install the vendored
#    source copy of the exact build the pipeline ran.
# ---------------------------------------------------------------------
THRESHOLDR_VERSION <- "0.1.0"
THRESHOLDR_REPO    <- "MDMotlagh/ThresholdR"
# Newest upstream commit on main that predates the recorded build date of
# the original install (2026-03-29). Confirmed to exist upstream. It is
# NOT confirmed to be the commit the original build came from, and it is
# known to differ from it (see env/vendor/README.md).
THRESHOLDR_REF     <- "9b639507d76e5348afda9ff633ca6c652bbd01b3"

say("")
say("[4/5] ThresholdR")
source_choice <- tolower(Sys.getenv("THRESHOLDR_SOURCE", "github"))
have <- installed_version("ThresholdR")

install_thresholdr_vendored <- function() {
    if (!file.exists(THRESHOLDR_TARBALL))
        stop(sprintf("vendored tarball not found at env/vendor/%s (looked under the repo root)",
                     basename(THRESHOLDR_TARBALL)))
    say("  installing vendored env/vendor/%s", basename(THRESHOLDR_TARBALL))
    utils::install.packages(THRESHOLDR_TARBALL, repos = NULL, type = "source")
    # install.packages() reports a failed source build as a WARNING, not an
    # error, so tryCatch(error=) around it never fires. Check the result
    # ourselves, otherwise a failed reinstall can leave an older ThresholdR
    # in place and this script would still exit 0.
    got <- installed_version("ThresholdR")
    if (is.na(got))
        stop("vendored ThresholdR install failed: package is not installed afterwards")
    if (got != THRESHOLDR_VERSION)
        stop(sprintf("vendored ThresholdR install failed: found %s, expected %s",
                     got, THRESHOLDR_VERSION))
    if (!(".intersect2mix" %in% ls(getNamespace("ThresholdR"), all.names = TRUE)))
        stop("vendored ThresholdR installed but .intersect2mix is missing, so this is not the recorded build")
    invisible(TRUE)
}

if (source_choice == "vendor") {
    # Explicit request: always reinstall from the vendored copy, even if
    # some other ThresholdR is already present.
    tryCatch(install_thresholdr_vendored(),
             error = function(e) record_failure("ThresholdR", e))
} else if (!is.na(have) && have == THRESHOLDR_VERSION) {
    say("  %-22s %-10s already installed", "ThresholdR", have)
} else {
    say("  installing from GitHub %s@%s", THRESHOLDR_REPO, THRESHOLDR_REF)
    ok <- tryCatch({
        remotes::install_github(paste0(THRESHOLDR_REPO, "@", THRESHOLDR_REF),
                                upgrade = "never", dependencies = TRUE)
        TRUE
    }, error = function(e) {
        message("  GitHub install of ThresholdR failed: ", conditionMessage(e))
        FALSE
    })
    if (!ok) {
        say("  falling back to the vendored source copy")
        tryCatch(install_thresholdr_vendored(),
                 error = function(e) record_failure("ThresholdR", e))
    }
}

# Tell the user plainly whether the ThresholdR they now have is the one
# the published results came from. The vendored build carries a private
# helper .intersect2mix(); upstream uses AdaptGauss::Intersect2Mixtures()
# instead, which is a different algorithm for the same quantity.
if (requireNamespace("ThresholdR", quietly = TRUE)) {
    is_vendored_build <- ".intersect2mix" %in%
        ls(asNamespace("ThresholdR"), all.names = TRUE)
    if (is_vendored_build) {
        say("  check: this is the build the pipeline ran (private .intersect2mix present)")
    } else {
        message(strrep("=", 70))
        message("WARNING: this ThresholdR is upstream, not the build the pipeline ran.")
        message("")
        message("  Upstream get_thresholds() computes the component cut point with")
        message("  AdaptGauss::Intersect2Mixtures(); the build used here computed it")
        message("  with a private uniroot-based helper. The two can disagree, so the")
        message("  bi_cut, tri_cut1 and tri_cut2 columns of the threshold table can")
        message("  change. The bi_thr / tri_thr columns (noise mean + 3*sigma), which")
        message("  are what step 5a gates on, are computed identically in both.")
        message("")
        message("  Upstream also needs the AdaptGauss package; install it with")
        message("  install.packages(\"AdaptGauss\") if get_thresholds() errors.")
        message("")
        message("  To reproduce the original numbers exactly, run:")
        message("    THRESHOLDR_SOURCE=vendor Rscript env/r_extra.R")
        message(strrep("=", 70))
    }
}

# ---------------------------------------------------------------------
# 5. Verify: every package must load, not just be present on disk
# ---------------------------------------------------------------------
say("")
say("[5/5] verifying that everything loads")

must_load <- c(names(bioc_pkgs), names(cran_pkgs), "dsb", "ThresholdR")
# These come from conda, not from this script, but the R steps stop dead
# without them, so check them here and fail early if the conda env is
# incomplete. scripts/config.R needs jsonlite and yaml; the step scripts
# need Matrix, Seurat and scDblFinder.
must_load <- c(must_load, "Matrix", "Seurat", "scDblFinder", "jsonlite", "yaml")

for (pkg in must_load) {
    ok <- suppressWarnings(suppressMessages(
        requireNamespace(pkg, quietly = TRUE)))
    if (ok) {
        say("  OK      %-22s %s", pkg, installed_version(pkg))
    } else {
        failures <- c(failures, sprintf("%s: installed but will not load", pkg))
        message(sprintf("  FAILED  %-22s will not load", pkg))
    }
}

# The five functions step 5a calls must actually exist.
if (requireNamespace("ThresholdR", quietly = TRUE)) {
    need_fun <- c("calculateBIC", "get_k_list", "fit_k",
                  "get_thresholds", "publish_fittings_thr")
    missing_fun <- need_fun[!need_fun %in% getNamespaceExports("ThresholdR")]
    if (length(missing_fun))
        failures <- c(failures,
                      sprintf("ThresholdR is missing exported functions: %s",
                              paste(missing_fun, collapse = ", ")))
}

rule()
if (length(failures)) {
    message("R environment setup FAILED. ", length(failures), " problem(s):")
    for (f in failures) message("  - ", f)
    message("")
    message("Nothing was rolled back. Fix the problems above and run this ",
            "script again; it skips whatever already installed cleanly.")
    quit(save = "no", status = 1)
}

say("SUCCESS: the R side of citeseq-nf is ready.")
say("")
say("Installed and loadable:")
say("  %s", paste(sort(must_load), collapse = ", "))
say("")
say("Recorded provenance for this environment is in env/r_session_info.txt.")
quit(save = "no", status = 0)

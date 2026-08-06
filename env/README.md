# Building the citeseq-nf environment

This folder holds everything you need to build the software environment for the
pipeline. You build it once. After that you only pass its path to Nextflow.

The environment has two parts:

1. `environment.yml`: almost everything, installed by conda.
2. `r_extra.R`: the R packages that conda cannot supply.

You must run both steps. Step 1 does the heavy work. Step 2 is short.

---

## What you need first

- A conda-style package manager. Any one of these works:
  `micromamba`, `mamba`, or `conda`. `micromamba` and `mamba` are much faster.
- About 10 GB of free disk space and a working internet connection. Most of
  that is the CUDA runtime that ships inside the pinned `torch` wheel.
- Linux on a 64-bit Intel or AMD machine (`linux-64`). The version pins were
  tested there.
- A GPU if you want to run step 06 (totalVI) and step 07 (scVI) at full speed.
  The CPU also works, but it is slow.

All commands below are run from the top folder of this repository.

---

## Step 1: create the environment

Pick the line that matches your tool.

```bash
# micromamba
micromamba create -f env/environment.yml

# mamba
mamba env create -f env/environment.yml

# conda
conda env create -f env/environment.yml
```

This creates an environment named `citeseq-nf`. It installs about 440 packages
and downloads roughly 370 MB. Expect 5 to 15 minutes.

Now turn the environment on:

```bash
micromamba activate citeseq-nf     # or: mamba activate / conda activate
```

---

## Step 2: add the extra R packages

One R package, `ThresholdR`, is not published on any conda channel. It lives on
GitHub only. The script `r_extra.R` fetches it for you and then checks that the
whole R side loads.

With the environment still active, run:

```bash
Rscript env/r_extra.R
```

The script also re-checks the R packages that step 1 already installed. It
skips anything that is already at the right version. On a fresh build five
packages always need work, because `environment.yml` cannot supply them at the
version the pipeline ran with:

| Package | `environment.yml` pins | The run used | Gap |
|---|---|---|---|
| `Rhdf5lib` | not pinned at all | 1.28.0 | conda pin missing |
| `rhdf5filters` | not pinned at all | 1.18.1 | conda pin missing |
| `rhdf5` | 2.50.0 | 2.50.2 | two patch releases |
| `celda` | 1.22.0 | 1.22.1 | one patch release |
| `decontX` | 1.4.0 | 1.4.1 | one patch release |

That is expected. Allow ten minutes or more, because `Rhdf5lib` compiles
the HDF5 C library from source.

`ThresholdR` needs `mixtools`, `mclust`, `foreach`, and `ggplot2`. Step 1
installs all four at the right versions, so the script leaves them alone.

Read `env/r_extra.R` for the ThresholdR options. There is a vendored copy in
`env/vendor/` if you need the exact build the published results came from.

---

## Step 3: check that it worked

Run these two checks. Both should print `OK` and no error.

Python side:

```bash
python -c "import scanpy, anndata, scvi, torch, soupx, skmisc, skimage, h5py, yaml; print('OK')"
```

R side:

```bash
Rscript -e 'suppressPackageStartupMessages({library(dsb); library(decontX); library(ThresholdR); library(scDblFinder); library(SingleCellExperiment); library(rhdf5); library(Matrix); library(ggplot2)}); cat("OK\n")'
```

If a line fails, the error names the missing package. Install that one package
and run the check again.

---

## Step 4: tell the pipeline where the environment is

Nextflow does not find the environment by itself. You pass the folder path with
`--conda_env`. Print the path like this:

```bash
python -c "import sys, pathlib; print(pathlib.Path(sys.prefix))"
```

Then use it when you start a run:

```bash
./run.sh run ./main.nf \
    --conda_env  <the path printed above> \
    --project_root <your project folder> \
    --outdir       <your output folder> \
    --panel        <your panel file>
```

The pipeline builds `--python` and `--rscript` from that path, so you do not
need to set them.

---

## Notes on the pinned versions

- **Versions are exact on purpose.** Almost every package is pinned to the
  version used for the validated pipeline run. The exceptions are the five
  R packages listed in step 2, which `r_extra.R` installs instead. This keeps
  results reproducible. If you loosen a pin, re-run the tests in `tests/`.
- **Scrublet is not a separate package.** The doublet finder Scrublet is built
  into `scanpy` 1.11 and is called as `sc.pp.scrublet`. Do not install a
  standalone `scrublet` package. It would conflict.
- **`torch` is pinned to 2.5.1.** That is the torch, Lightning and
  `scvi-tools` combination the pipeline was validated with. If you change it,
  re-run the tests in `tests/`. On Linux the pinned `torch==2.5.1` wheel is
  the CUDA 12.4 build, so it uses the GPU when one is present.
- **A few pip packages are used instead of conda ones.** `torch`,
  `scvi-tools`, `mudata`, `soupx-python`, `scikit-misc`, and `scikit-image`
  come from PyPI, because that is the combination that was tested.
- **No `prefix:` line.** The file has no hard-coded install location, so your
  package manager puts the environment in your own envs folder.
- **`AdaptGauss` is missing on purpose.** conda has no build of it for R 4.4.
  You only need it if `ThresholdR::get_thresholds()` raises an error. In that
  case run `Rscript -e 'install.packages("AdaptGauss")'`.
- **`rpy2` is installed but not imported.** No step script calls it today. It
  was part of the validated environment, so it stays pinned.

---

## Rebuilding after a change

If you edit `environment.yml`, delete the old environment first. Updating in
place often leaves stale packages behind.

```bash
micromamba env remove -n citeseq-nf
micromamba create -f env/environment.yml
micromamba activate citeseq-nf
Rscript env/r_extra.R
```

Before you build for real, you can ask the solver to check the file without
installing anything:

```bash
micromamba create --dry-run -n citeseq-nf-check -f env/environment.yml
```

The conda part must solve cleanly. The pip part is skipped in a dry run, and
the warning it prints about a missing prefix is expected.

---

## What is in the environment

Python is 3.10.19. Every version below is quoted from `env/environment.yml`.

Imported by name in `scripts/config.py` or `scripts/python/*.py`:

| Package | Version | Where it is used |
|---|---|---|
| `scanpy` | 1.11.5 | every step script (`import scanpy as sc`) |
| `anndata` | 0.11.4 | step 08 |
| `numpy` | 2.2.6 | most step scripts |
| `pandas` | 2.3.3 | most step scripts |
| `scipy` | 1.15.2 | sparse matrices, `gaussian_kde` (step 05b, 05c), `find_peaks` (step 05b) |
| `scikit-learn` | 1.7.2 | `GaussianMixture` (step 05c), `silhouette_score` (step 08) |
| `matplotlib-base` | 3.10.8 | all figures |
| `h5py` | 3.15.1 | step 00b reads `.h5` matrices directly |
| `pyyaml` | 6.0.3 | `scripts/config.py` reads the panel YAML |
| `scvi-tools` | 1.3.3 | step 06 totalVI, step 07 scVI (`import scvi`) |
| `soupx-python` | 0.3.0 | step 00b (`import soupx`). This is the **Python** SoupX port, not the R `SoupX` package. |

Never imported by name, but a step fails without them:

| Package | Version | Why it is needed |
|---|---|---|
| `torch` | 2.5.1 | the deep-learning backend `scvi-tools` runs on. Pinned below 2.6 so the Lightning version `scvi-tools` pulls in keeps working. |
| `mudata` | 0.3.3 | data container `scvi-tools` uses |
| `scikit-misc` | 0.5.2 | needed by `sc.pp.highly_variable_genes(flavor="seurat_v3")` in step 04 |
| `scikit-image` | 0.25.2 | threshold search inside `sc.pp.scrublet` in step 01 |
| `leidenalg` | 0.11.0 | `sc.tl.leiden` in steps 06 and 07 |
| `python-igraph` | 1.0.0 | graph backend for `leidenalg` |
| `umap-learn` | 0.5.11 | `sc.tl.umap` in steps 06, 07 and 08 |
| `pynndescent` | 0.5.13 | nearest-neighbour backend for `umap-learn` |
| `statsmodels` | 0.14.6 | required by `soupx-python` |
| `numba` | 0.63.1 | just-in-time compiler that scanpy and `umap-learn` use |
| `seaborn-base` | 0.13.2 | used by some scanpy plotting calls |

`env/environment.yml` pins a few more small helpers that scanpy and `h5py`
need: `natsort` 8.4.0, `tqdm` 4.67.3, `session-info2` 0.3 and `hdf5` 1.14.6.

`rpy2` 3.6.6 is also pinned. No step script imports it today. It stays because
it was part of the validated environment.

### R packages

R is 4.4.3, with Bioconductor release 3.20. Versions come from
`env/environment.yml` unless the table says otherwise.

| Package | Version | Loaded by |
|---|---|---|
| `decontX` | 1.4.1 (from `env/r_extra.R`) | step 00c, ambient-RNA removal |
| `rhdf5` | 2.50.2 (from `env/r_extra.R`) | steps 00c and 01b, reads 10x `.h5` matrices |
| `scDblFinder` | 1.23.4 | step 01b, doublet calling |
| `SingleCellExperiment` | 1.28.0 | step 01b container |
| `dsb` | 2.0.1 | step 02, ADT denoising |
| `ThresholdR` | 0.1.0 (vendored, see below) | step 05a, per-marker gating |
| `ggplot2` | 4.0.2 | step 05a fit plots |
| `Matrix` | 1.7-4 | steps 00c, 01b, 02 |
| `Seurat` | 5.4.0 | `Seurat::Read10X_h5` in steps 00c and 02 |
| `jsonlite` | 2.0.0 | `scripts/config.R` |
| `yaml` | 2.3.12 | `scripts/config.R` |

These packages are installed as dependencies. No script loads them.

| Package | Version | Needed by |
|---|---|---|
| `celda` | 1.22.1 (from `env/r_extra.R`) | `decontX` |
| `mclust` | 6.1.2 | called from inside `ThresholdR` |
| `mixtools` | 2.0.0.1 | called from inside `ThresholdR` |
| `foreach` | 1.5.2 | called from inside `ThresholdR` |
| `limma` | 3.62.1 | `dsb` |
| `magrittr` | 2.0.4 | `dsb` |
| `remotes` | 2.5.0 | `env/r_extra.R` |

Two notes on the R side:

- **`ThresholdR` is on no conda channel.** It lives on GitHub only.
  `env/r_extra.R` installs it. `env/vendor/ThresholdR_0.1.0.tar.gz` is a
  rebuild of the build the reference run of 2026-04-28 used, not the original
  source. No original source tarball was archived, so the package was rebuilt
  by deparsing the installed package database. The executed code was checked
  function by function and is identical. The source comments and the code
  formatting are not the original ones. Read `env/vendor/README.md` before you
  choose between the two sources. This is the most fragile dependency in the
  whole stack.
- **`env/r_extra.R` reinstalls five R packages** that `env/environment.yml`
  does not supply at the recorded version. This is expected.

  | Package | `env/environment.yml` pins | Recorded version | Gap |
  |---|---|---|---|
  | `rhdf5` | 2.50.0 | 2.50.2 | two patch releases |
  | `celda` | 1.22.0 | 1.22.1 | one patch release |
  | `decontX` | 1.4.0 | 1.4.1 | one patch release |
  | `Rhdf5lib` | not pinned at all | 1.28.0 | conda pin missing |
  | `rhdf5filters` | not pinned at all | 1.18.1 | conda pin missing |


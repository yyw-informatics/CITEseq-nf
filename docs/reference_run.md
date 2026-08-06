# Reference run (2026-04-28)

This page records the full pipeline run that the baseline tables come from.
It exists so that a new owner has something to compare against after running
the pipeline for the first time.

Read the warning in the next section before you use any number on this page.

---

## A partial check of the generic code, 2026-08-04

The 2026-04-28 run used the older study-specific code. On **2026-08-04** the
generic code was given one real-data check. It was **not** a second full run,
and it did not compare against the 2026-04-28 numbers for most of what it
checked. Read this section carefully before treating it as a reproduction.

| Item | Value |
|---|---|
| Date | 2026-08-04 |
| Wall clock | 2 h 21 min 33 s, exit 0 |
| Tasks | 7 ran fresh, 12 came from cache |
| Libraries computed by the generic code | 2 of 16 |
| Libraries reused from the 2026-04-28 run | 14 of 16 |
| GPU steps | both ran, totalVI and scVI on separate queues |

The run was launched with `--n_libraries 2`, but a bug meant only the
parallel fan-out honoured it. Every step that reads the library list from the
metadata still used all 16 entries and picked up 14 per-library files left
over from 2026-04-28. The bug was fixed after this run.

The wall clock is not comparable with the 2 h 43 min above. That one ran 61
tasks with nothing cached. This one resumed after an earlier failed attempt,
and most of its time was scVI training.

### What the six checks were

`tests/validate_phase_a.sge` reported 6 PASS, 0 FAIL, 0 SKIP. Their baseline
is a Phase A output directory outside this repository, dated **2026-04-13**,
which is earlier than the 2026-04-28 run and is not `tests/baseline/`.

- Four h5ad checks compare structure only: shape, obs column names, layer
  names, obsm keys. No values.
- One CSV check compares row and column counts only.
- One CSV check compares values: `thresholdr_thresholds.csv` was an **exact
  match** after sorting rows.

### The embedding numbers

These were compared by eye against the committed 2026-04-28 table, and were
not one of the six checks.

| Metric | 2026-08-04 | 2026-04-28 | Difference |
|---|---|---|---|
| batch_silhouette, totalVI | -0.0570 | -0.0519 | 9.9 % |
| batch_silhouette, scVI | -0.0346 | -0.0367 | 5.6 % |
| bio_silhouette MARKER_C, totalVI | 0.1163 | 0.1219 | 4.6 % |
| bio_silhouette MARKER_C, scVI | 0.0882 | 0.0766 | 15.2 % |

Both runs chose `X_totalVI`. The automatic check in
`tests/validate_citeseq_nf_e2e.sge` uses a 2 percent tolerance, so these
differences would not pass it. Read them as the same sign, the same order of
magnitude, and the same decision, not as a match.

### The panel file was reconstructed

The panel YAML system was added after 2026-04-28, so no panel file existed
for this data. One was reconstructed on 2026-08-04 from the study's saved
panel annotations, the ADT feature names read from the CellRanger `.h5`
files, and the study's own RNA concordance map. The exact threshold match
above is the best evidence that the reconstruction is correct.

---

## Warning: these numbers are indicative, not a bit-exact target

The exact code revision that produced these numbers was Nextflow revision
`4d1b7954d4`. **That revision was never committed. It no longer exists.**
Nobody can check it out, and nobody can reproduce it exactly.

The pipeline in this repository is the same pipeline in structure and in
method, but the source is not byte-identical to what ran on 2026-04-28.

On top of that, two steps are stochastic (they use random numbers):

- **totalVI** (step 06) and **scVI** (step 07) are deep generative models.
  They start from random weights, shuffle mini-batches at random, and split
  train/validation at random. Two runs on the same input give different
  numbers.
- Anything computed *after* those models (silhouette scores, neighbour graph,
  UMAP, clustering) inherits that randomness.

So:

| Expect to match | Do not expect to match |
|---|---|
| The **shape** of every output table (same columns, same number of rows) | The exact decimal values of totalVI or scVI losses |
| Per-library cell counts and QC pass rates (these steps are deterministic) | Silhouette scores beyond the first two decimals |
| The **order of magnitude** and the **sign** of the embedding metrics | The exact epoch at which early stopping fires |
| The overall pass/fail verdicts of ADT gating fitness | Byte-identical CSV files |

If a new number is in the same range as the number here, that is a pass. If a
new number is 10x different, or has the opposite sign, or a whole column has
disappeared, that is a real problem worth investigating.

---

## What was run

| Item | Value |
|---|---|
| Date | 2026-04-28 |
| Start | 13:42:57 local time |
| End | 16:25:43 local time |
| Wall clock | 2 h 42 min 46 s (about 2 h 43 min) |
| Nextflow | 25.10.4 (DSL 2) |
| Revision reported by Nextflow | `4d1b7954d4` (never committed, see warning above) |
| Libraries processed | 16 |
| Donors | 6 (anonymised as `donor_01` … `donor_06`) |
| Executor | SGE (Sun Grid Engine), `-profile sge` |
| Entry point | default (full pipeline: Phase A QC, then Phase B embedding) |
| Cached tasks | none. All 61 tasks actually ran, so the wall clock is a real end-to-end time |

### Shape of the 16 libraries

Six donors, two groups, two timepoints, with technical replicates for some
donors. Labels in the baseline tables look like `donor_01_aged_d0_rep1`.

### Cluster setup

The driver job asked for 4 slots, 16 GB memory, and a 12 hour limit. Nextflow
then submitted each step as its own SGE job.

Two steps need a GPU. They were sent to two **different** GPU queues so they
could train at the same time instead of queueing behind each other:

| Step | Label | CPUs | Memory | Time limit | Queue supplied by |
|---|---|---|---|---|---|
| 06 totalVI | `gpu_totalvi` | 4 | 64 GB | 4 h | `--sge_gpu_totalvi_queue` |
| 07 scVI | `gpu_scvi` | 8 | 12 GB | 4 h | `--sge_gpu_scvi_queue` |

Both steps request `-l gpu=1`. The queue strings are site-specific, so they
are passed on the command line and are not stored in this repository.

Rough timing inside the run:

| Segment | Approximate time |
|---|---|
| Phase A (DecontX, QC, scDblFinder, DSB, CLR, RNA norm, ThresholdR, KDE, fitness) | 13:43 to about 14:11 |
| Step 06 totalVI training started | about 14:11 |
| Step 06 totalVI finished | 14:50 |
| Step 07 scVI finished | 16:19 |
| Step 08 embedding comparison finished | 16:25 |

totalVI and scVI ran in parallel. scVI was the long pole because its epoch
budget is much larger (400 versus 100).

---

## Verified software versions

These versions were read directly from the analysis environment as it was on
2026-04-28. This section is a record of that one run. It is not the build
recipe. Where a version here differs from [`../env/environment.yml`](../env/environment.yml),
the yml file is what the current code needs, so build your environment from the
yml. Two differences to know about: this page lists h5py 3.16.0 while the yml
pins h5py 3.15.1, and this page lists celda 1.22.1 while the yml pins 1.22.0.
`env/r_extra.R` raises celda to 1.22.1 during the build, so a built environment
does match.

### Core

| Tool | Version |
|---|---|
| Nextflow | 25.10.4 |
| Python | 3.10.19 |
| R | 4.4.3 |

### Python packages

| Package | Version |
|---|---|
| scanpy | 1.11.5 |
| anndata | 0.11.4 |
| scvi-tools | 1.3.3 |
| torch | 2.5.1 |
| lightning | 2.6.1 |
| jax | 0.6.2 |
| numpy | 2.2.6 |
| scipy | 1.15.2 |
| pandas | 2.3.3 |
| scikit-learn | 1.7.2 |
| numba | 0.63.1 |
| leidenalg | 0.11.0 |
| igraph | 1.0.0 |
| matplotlib | 3.10.8 |
| h5py | 3.16.0 |
| rpy2 | 3.6.6 |

Scrublet is not a separate package here. It is used through
`scanpy.pp.scrublet`, so its version follows scanpy.

### R packages

| Package | Version |
|---|---|
| celda (provides DecontX) | 1.22.1 |
| scDblFinder | 1.23.4 |
| dsb | 2.0.1 |
| ThresholdR | 0.1.0 |
| Seurat | 5.4.0 |
| SingleCellExperiment | 1.28.0 |
| Matrix | 1.7.4 |
| mclust | 6.1.2 |

The R package SoupX is not in this table on purpose. The R version of SoupX was
installed in that environment, but the pipeline never calls it. Step 00b uses
the Python port instead: it imports the `soupx` module, which comes from the
`soupx-python` package pinned in
[`../env/environment.yml`](../env/environment.yml).

Note: `numpy` 2.x and `scanpy` 1.11.x are a fairly new pairing. If you rebuild
the environment from scratch and hit import errors, this table is the
combination that is known to work.

---

## Key metrics from the run

### Cell QC (step 01)

| Metric | Value |
|---|---|
| Input cells across 16 libraries | 202,159 |
| Cells passing QC | 176,664 (87.4 %) |
| Per-library pass rate | 63.7 % to 96.4 % |
| Scrublet doublets removed | 139 total (0 to 64 per library) |
| Cells after doublet removal | 176,525 |

QC settings that produced this: at least 100 genes per cell, at most 20 %
mitochondrial reads, genes kept if seen in at least 3 cells, Scrublet expected
doublet rate 0.06.

One library stands out: the one at 63.7 % pass rate. That is normal for this
data set, not a bug. Expect one weak library.

### Doublet method comparison (step 01b)

Scrublet runs first and its doublets are already removed. So in
`doublet_comparison.csv` the Scrublet column is 0 everywhere by construction.
The interesting column is scDblFinder, which is run on the surviving cells:

| Metric | Value |
|---|---|
| scDblFinder calls | 9,447 cells across the cohort |
| Per-library rate | 2.35 % to 9.26 % |

The two methods disagree a lot. That disagreement is the point of the
comparison table. scDblFinder is far more aggressive than Scrublet on this
data.

### Ambient RNA removal, DecontX (step 00c)

| Metric | Value |
|---|---|
| Median contamination per library | 0.013 to 0.120 |
| Clusters found per library | 12 to 24 |
| Runtime per library | 66 s to 389 s (about 40 min of CPU time in total, run in parallel) |

Median contamination around 5 % is typical here. A library above 0.15 would be
worth a second look.

### ADT marker gating (steps 05a, 05b, 05c)

About 161,000 cells reached the ADT gating step. Four markers were gated.

**Important caveat about these two tables.** The gating steps in this run read
a stale `adt_counts_for_thresholdr.csv` left over from an earlier run, not the
file this run produced. That is why gating saw about 161,000 cells while the
QC step passed 176,664. So `thresholdr_fitness_metrics.csv` and
`kde_vs_thresholdr_comparison.csv` are dated 2026-04-02, not 2026-04-28, and
they describe a slightly different input. Treat their numbers as indicative of
shape and order of magnitude only. Every other table here comes from the
2026-04-28 run itself. A fresh run will regenerate all of them consistently.

Fitness verdicts (`thresholdr_fitness_metrics.csv`):

| Marker | Components chosen | Stain quality score | Verdict |
|---|---|---|---|
| MARKER_A | 2 | 3.73 | PASS |
| MARKER_C | 2 | 2.86 | FAIL |
| MARKER_B | 2 | 1.71 | FAIL |
| MARKER_D | 2 | 1.03 | FAIL |

Only MARKER_A passed. That is a real property of this antibody panel, not a
pipeline failure. The fitness step is doing its job: it is telling you that
three of four stains are weak or overlapping.

The KDE cross-check (`kde_vs_thresholdr_comparison.csv`) shows how much the
two threshold methods disagree:

| Marker | KDE positive | ThresholdR positive | Cells reclassified |
|---|---|---|---|
| MARKER_A | 33.5 % | 34.2 % | 1,072 |
| MARKER_B | 5.4 % | 9.1 % | 5,916 |
| MARKER_D | 7.0 % | 9.6 % | 4,161 |
| MARKER_C | 0.01 % | 2.5 % | 4,044 |

MARKER_C is the extreme case: the KDE valley method finds almost nothing (21
cells) while ThresholdR finds 4,065. When you see this pattern, trust neither
threshold blindly. This is exactly why the panel config supports an
`fgprob_override` list that switches a marker to the totalVI foreground
probability instead.

### totalVI training (step 06)

| Metric | Value |
|---|---|
| Latent dimensions | 15 |
| Max epochs allowed | 100 |
| Epochs actually run | 56 (early stopping fired) |
| Training ELBO, epoch 0 | 1365.8 |
| Training ELBO, final epoch | 844.7 |
| Learning rate schedule | 2e-3 down to 1e-4 |
| KL warm-up | 30 epochs |

The loss curve drops steeply for the first 5 epochs, then flattens. By epoch
30 it is essentially flat. A rerun should show the same curve shape. The exact
final ELBO will differ, and the early-stopping epoch will differ too.

### Embedding comparison (step 08)

The biological separation marker was MARKER_C.

| Metric | totalVI | scVI |
|---|---|---|
| Batch silhouette (lower is better mixing) | -0.0519 | -0.0367 |
| Biological silhouette on MARKER_C (higher is better) | 0.1219 | 0.0766 |
| Decision | `X_totalVI` | |

totalVI won on both counts: better batch mixing and clearly better separation
on the protein marker. The gap on the biological silhouette (0.12 versus 0.08)
is the meaningful one. Both batch silhouettes are small negative numbers,
which means batch mixing is good in both embeddings.

**What to check on a rerun:** the decision should still be `X_totalVI`, the
batch silhouettes should still be small and negative, and the totalVI
biological silhouette should still be clearly larger than the scVI one. The
third decimal place will move. That is fine.

---

## Where the baseline files are

All the small result tables from this run are saved in this repository under:

```
tests/baseline/
```

See `tests/baseline/README.md` for a one-line description of each file.

Donor labels in those files are pseudonyms, not anonymous. Real donor numbers
were replaced with `donor_01` through `donor_06`. Age group, timepoint,
replicate number and every numeric value are unchanged. The original sample
numbers were removed, because they would let anyone holding the study sample
sheet reverse the pseudonyms. The mapping is not stored in this repository,
but someone with access to the study records could rebuild it.

Marker names (MARKER_A, MARKER_D, MARKER_B, MARKER_C) are real study markers. They are not
personal data and were kept as-is.

The large outputs from this run (the `.h5ad` files, the trained model
checkpoints, the figures) are **not** in this repository. They were far too
large to commit and they no longer exist.

---

## What a new owner should do

1. **Run the pipeline once, end to end, on your own data.** Use the default
   entry point. See the main `README.md` for the command and for the required
   `--conda_env`, `--project_root`, `--outdir`, `--panel`, and GPU queue arguments.

2. **Compare shape first, values second.** Open your `qc_summary.csv` next to
   `tests/baseline/qc_summary.csv`. Do they have the same columns? Does yours
   have one row per library? If yes, the pipeline is wired correctly.

3. **Expect similar magnitudes, not identical values.** A QC pass rate in the
   80s or 90s, one weaker library, DecontX median contamination in the low
   percent range, scDblFinder in the low single-digit to ~10 % range. These
   are the ranges to sanity-check against.

4. **Do not chase exact totalVI or scVI numbers.** They are stochastic. If you
   rerun the same pipeline twice on the same input you will get two different
   ELBO curves and two slightly different silhouette scores. Only the
   qualitative conclusion should be stable.

5. **Treat a sign flip or a 10x change as a red flag.** For example, if the
   batch silhouette becomes clearly positive, or the totalVI biological
   silhouette drops below scVI, something changed in the data or in the
   method. Investigate before trusting downstream results.

6. **Re-record your own reference run.** Once you have a run you trust, copy
   its small CSVs into `tests/baseline/` (anonymising donor IDs the same way)
   and update this page with your date, your versions, and your numbers. Then
   commit the code revision that produced them, so the next person does not
   inherit the problem described in the warning at the top of this page.

### The one thing not to repeat

Nextflow printed a revision hash for this run, but that hash does not exist in
this repository's history. You can check for yourself:

```
git cat-file -t 4d1b7954d4
# fatal: Not a valid object name 4d1b7954d4
```

The code that produced these numbers was never committed anywhere that
survived. When the original machine went away, so did the exact source.

Before any run whose numbers you want to keep: commit first, run second, and
write the commit hash into your run log.

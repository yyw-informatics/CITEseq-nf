# Baseline results

These CSV files are the small result tables from the reference pipeline run of
**2026-04-28** (16 libraries, 6 donors, Nextflow 25.10.4). They are kept here
so a new owner can compare a fresh run against something known to work.

Full context, software versions, and the key numbers are in
[`../../docs/reference_run.md`](../../docs/reference_run.md).

## Read this first

**Donor labels are pseudonyms, not anonymous.** The real donor numbers were
replaced with `donor_01` through `donor_06`. Age group (`young` / `aged`),
timepoint (`d0` / `d7`), replicate (`rep1` / `rep2`) and **every numeric
value** are unchanged. The original sample numbers were removed, because they
would let anyone holding the study sample sheet reverse the pseudonyms.

The mapping back to the real donor numbers is not stored in this repository.
Someone with access to the study records could still rebuild it. Treat these
labels as pseudonyms and not as anonymous data.

Marker names (CD7, CD16, CD117, CRTH2) are real study markers. They are not
personal data, so they were left as-is.

**These are not a bit-exact target.** The code revision that produced them was
never committed and no longer exists, and totalVI and scVI are stochastic
(they use random numbers). Compare the *shape* of the tables and the
*magnitude* of the values, not the exact digits. See the warning at the top of
`docs/reference_run.md`.

## What each file is

| File | Produced by | What it holds |
|---|---|---|
| `qc_summary.csv` | step 01 (QC + Scrublet) | One row per library. Cells in, cells passing QC, pass percentage, the minimum-genes cutoff used, Scrublet doublets found, and cells left after doublet removal. |
| `doublet_comparison.csv` | step 01b (scDblFinder) | One row per library. Compares two doublet callers on the surviving cells. The Scrublet columns are 0 everywhere because Scrublet ran first in step 01 and its doublets were already removed. The useful columns are `scDblFinder_doublets` and `scDblFinder_pct`. |
| `doublet_comparison_decontx.csv` | step 01b | Same table for the DecontX decontamination branch. In this run the DecontX branch was the selected branch, so the two files are identical. Both are kept so a future run that uses the SoupX branch has both sides to compare. |
| `decontx_summary.csv` | step 00c (DecontX, R) | One row per library. Cells, clusters found, number of empty droplets used to estimate the ambient profile, and the median / mean / min / max ambient contamination fraction per cell, plus runtime in seconds. |
| `thresholdr_fitness_metrics.csv` | step 05c (fitness) | One row per gated ADT marker. Number of mixture components chosen, BIC margin, overlap coefficient, stain quality score, smallest component fraction, threshold-to-valley distance, RNA concordance ratio, and a PASS / FAIL verdict. This is the table that tells you whether a stain is usable. |
| `kde_vs_thresholdr_comparison.csv` | step 05b (KDE validation) | One row per gated ADT marker. Cross-checks the ThresholdR cutoff against a simple KDE-valley cutoff and against the totalVI foreground probability. Includes positive-cell counts, percentages, pairwise Jaccard overlaps, and how many cells the two methods classify differently. |
| `embedding_comparison.csv` | step 08 (embedding comparison) | Three rows. Batch silhouette and biological silhouette for totalVI versus scVI, plus the recorded decision of which embedding became the primary one. In this run the decision was `X_totalVI`. |
| `totalvi_training_history.csv` | step 06 (totalVI, GPU) | The training ELBO for each epoch. 56 epochs, because early stopping fired before the 100-epoch limit. Useful for checking the shape of the loss curve. |
| `totalvi_full_history.csv` | step 06 (totalVI, GPU) | The complete Lightning training log for the same 56 epochs: learning rate, KL weight, and the train and validation split of the loss into reconstruction and KL terms. Use this when the short history is not enough to explain a training problem. |

No file in this directory contains a cell barcode, a real donor ID, a file
path, a hostname, or a user name.

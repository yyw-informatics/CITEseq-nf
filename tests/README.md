# tests/: how to check the pipeline on your cluster

This directory holds four SGE job scripts. SGE (Sun Grid Engine) is the job
scheduler these scripts submit to. You start each one with `qsub`.

They exist so a new owner can answer one question: **does citeseq-nf still
work here, on my cluster, on my data?**

None of them hardcode a path, a user name, a host name, or a queue name. You
pass everything in as environment variables. The scripts find the repository
themselves, because they live inside it.

---

## 1. The four scripts

| Script | What it does | Needs a GPU? | Rough time |
|---|---|---|---|
| `run_full_pipeline.sge` | Runs the whole pipeline: Phase A (QC) then Phase B (embedding). Then checks that the five main output files exist. | Yes (the pipeline sends two steps to GPU queues) | About 3 hours for 16 libraries. Add queue waiting time. |
| `validate_citeseq_nf_e2e.sge` | End-to-end check of step 08. It builds a fresh reference result by running step 08 directly, then runs the same work through Nextflow, then compares the two. | No, in the default mode | 15 to 40 minutes. Hours if you set `RUN_FULL_DEFAULT=1`. |
| `validate_phase_a.sge` | Compares the Phase A output files against a known-good copy you put on disk yourself. | Only if you let it run the pipeline | About 10 minutes in compare-only mode (`SKIP_NEXTFLOW=1`). Hours otherwise. |
| `validate_phase_b_gpu.sge` | Records the numbers of the current totalVI and scVI outputs, retrains both models, then checks the new outputs against those numbers with a tolerance. | Yes | About 2 to 3 hours of training, plus queue waiting time. |

Some terms used above:

- **Phase A** is the QC part: ambient RNA removal, cell filtering, doublet
  calling, ADT normalisation, RNA normalisation, marker gating.
- **Phase B** is the embedding part: totalVI (step 06) and scVI (step 07)
  train on a GPU, then step 08 compares the two embeddings.
- **totalVI** and **scVI** are deep generative models. They learn a small set
  of numbers per cell (a "latent" representation) that summarises the cell.
- An **embedding** is that small set of numbers. It is what all the later
  plots and clustering use.
- **ADT** means antibody-derived tag: the protein measurement in CITE-seq.

### More detail on each script

**`run_full_pipeline.sge`** is not a test. It is the normal way to run the
pipeline as a batch job. It calls the pipeline with `-resume`, so a second run
reuses cached results for every step whose inputs and code did not change.
After the run it prints `OK` or `MISS` for five expected output files and ends
with `RESULT: PASS` or `RESULT: FAIL`.

**`validate_citeseq_nf_e2e.sge`** is the cheapest real test. By default it
uses the entry point `-entry from_totalvi_scvi`, which starts at step 08 and
skips all GPU training. It checks: cell and gene counts match exactly, the
embedding keys `X_totalVI`, `X_scVI` and `X_umap` are present, the ADT
`_denoised` and `_fg_prob` columns are present, the clustering column
`leiden_totalvi` is present, the chosen primary embedding matches, and the
silhouette scores agree within 2 percent. A silhouette score measures how well
separated groups of cells are.

**`validate_phase_a.sge`** compares files, not numbers you have to trust. For
each h5ad file it checks shape, `obs` column names, layer names and `obsm`
keys. For `thresholdr_thresholds.csv` it checks an exact match after sorting
rows. For `dsb_normalized_adt.csv` it checks the row and column count. It does
**not** create a baseline. You must put a known-good copy of the Phase A files
in `<PROJECT>/data/processed/phase_a_baseline/` first. If you have no such
copy, use `validate_citeseq_nf_e2e.sge` instead, because that one builds its
own reference.

Note one limit: citeseq-nf has no "Phase A only" entry point. If you let this
script run the pipeline, it runs Phase A **and** Phase B, so it needs the GPU
queue names. Only the Phase A files are compared at the end. Use
`SKIP_NEXTFLOW=1` to compare files that already exist and skip the run.

**`validate_phase_b_gpu.sge`** cannot compare files byte for byte, because
totalVI and scVI use random numbers. Instead it checks invariants: the h5ad
files and model directories exist and hold a `model.pt`, the cell count is
within 0.5 percent, the gene count is identical, the latent size matches
`scripts/config.py`, `protein_denoised` and `protein_fg_prob` have one column
per functional ADT marker in your panel, the `obsm` key set is unchanged, and
the last totalVI training ELBO is within 2 percent. ELBO is the training loss
of the model. Only totalVI has this ELBO check, because the script reads
`totalvi_training_history.csv` and there is no scVI equivalent.
This script runs the pipeline from the HVG checkpoint
(`-entry from_hvg`), so it retrains both models.

---

## 2. The variable contract

All four scripts read the same variable names. Export them in your shell, then
submit with `qsub -V` so the job inherits them.

### Always required

| Variable | Meaning |
|---|---|
| `PROJECT` | Your study root. This is the directory that holds `data/`, `results/`, `models/` and `figures/`. Outputs are written under it. |
| `CONDA_ENV` | Path to the conda environment built from `env/environment.yml`. It must contain `bin/python`, and `bin/Rscript` for the R steps. |
| `PANEL` | Path to your ADT panel YAML file. See `panels/example_panel.yaml`. The pipeline requires this. Use an **absolute** path. Each pipeline step runs in its own Nextflow work directory, so a relative path is not found there. Only `run_full_pipeline.sge` turns a relative path into an absolute one for you. |

### Required in some scripts

| Variable | Meaning | Required by |
|---|---|---|
| `DATA_ROOT` | CellRanger input root. One subdirectory per library. | Always in `run_full_pipeline.sge`, `validate_citeseq_nf_e2e.sge` and `validate_phase_b_gpu.sge`. The job stops at once if you leave it unset. The two validation scripts need it even when they do not run the pipeline, because `scripts/config.py` will not import without it. In `validate_phase_a.sge` only when the pipeline runs. |
| `GPU_QUEUE_TOTALVI` | SGE queue for step 06 totalVI, for example `gpu.q@your-gpu-host`. | Always in `run_full_pipeline.sge` and `validate_phase_b_gpu.sge`. In `validate_phase_a.sge` unless `SKIP_NEXTFLOW=1`. In `validate_citeseq_nf_e2e.sge` only with `RUN_FULL_DEFAULT=1`. |
| `GPU_QUEUE_SCVI` | SGE queue for step 07 scVI. | Same as above. |

Give the two GPU steps **different** queues if your cluster has more than one
GPU host. Then they train at the same time instead of waiting for each other.

### Optional, with safe defaults

| Variable | Default | Meaning |
|---|---|---|
| `SAMPLE_METADATA_JSON` | not set | JSON file mapping library id to sample metadata. Passed as `--sample_metadata_json`. |
| `FORCED_HVGS_FILE` | not set | Text file, one gene per line. Those genes are always kept as highly variable genes. Passed as `--forced_hvgs_file`. |
| `BIO_SEP_MARKER` | not set | The marker used to score biological separation. Pick a protein from your own panel that you expect to split your cells into two groups. Passed as `--biological_separation_marker`. There is no default marker, and none is assumed. If you leave it unset, the biological silhouette comparison in `validate_citeseq_nf_e2e.sge` is skipped. |
| `N_LIBRARIES` | not set | Run on only the first N libraries. Passed as `--n_libraries`. Every step honours it, so a pilot really does process N libraries and not the whole cohort. Set `SAMPLE_METADATA_JSON` as well: `N_LIBRARIES` selects the first N entries of that file. |
| `NEXTFLOW_MODULE` | `nextflow/25.10.4` | Module name to load if `nextflow` is not already on your `PATH`. Change it to whatever your cluster calls it. |
| `SCRATCH` | `${TMPDIR:-/tmp}` | Temporary directory for intermediate files. |
| `CITESEQ_NF` | derived | Path to this repository. The scripts work it out from their own location. Set it only if a script says it cannot find `main.nf`. |

`validate_phase_a.sge` has one extra rule: when it runs the pipeline, you must
set either `SAMPLE_METADATA_JSON` or `N_LIBRARIES`, otherwise the pipeline
cannot tell how many libraries to process.

### Switches

| Switch | Scripts | Effect |
|---|---|---|
| `SKIP_NEXTFLOW=1` | `validate_phase_a.sge` | Do not run the pipeline. Compare the Phase A files already on disk against your baseline copy. Nothing else runs. |
| `SKIP_NEXTFLOW=1` | `validate_phase_b_gpu.sge` | Do not run the pipeline. Read the current Phase B outputs and check them. See the warning below: this mode compares those files with themselves. |
| `SKIP_NEXTFLOW=1` | `validate_citeseq_nf_e2e.sge` | Do not run Nextflow. The script **still builds a fresh baseline** by running step 08 directly, then compares it against the pipeline outputs already on disk. Add `SKIP_BASELINE=1` if you want to skip that rebuild too. |
| `SKIP_BASELINE=1` | `validate_citeseq_nf_e2e.sge` | Reuse the saved baseline instead of building a fresh one. |
| `RUN_FULL_DEFAULT=1` | `validate_citeseq_nf_e2e.sge` | Run the full default workflow instead of step 08 only. This retrains on GPU and takes hours. |

Two warnings about `SKIP_NEXTFLOW=1`:

- In `validate_citeseq_nf_e2e.sge` this switch is **not** free. Building the
  baseline runs step 08 for real, which costs time and memory. Only
  `SKIP_BASELINE=1` skips that part.
- In `validate_phase_b_gpu.sge` the script always records its baseline numbers
  first, from the same files it later checks. So with `SKIP_NEXTFLOW=1` the
  cell count, gene count, `obsm` keys and ELBO are compared against themselves
  and always pass. Only three things are still meaningful in that mode: the
  files and `model.pt` exist, the latent size matches `scripts/config.py`, and
  the ADT column counts match your panel.

`validate_phase_b_gpu.sge` also asks for `GPU_QUEUE_TOTALVI` and
`GPU_QUEUE_SCVI` even with `SKIP_NEXTFLOW=1`. It checks them before it reads the
switch. Set both, or the job stops right away.

### About `CITESEQ_NF`

Under SGE, the scheduler copies your job script to a spool directory before it
runs. So `$0` inside the job is not the file you submitted. Each script first
tries to find the repository one level up from its own path. If that path has
no `main.nf`, it falls back to a `CITESEQ_NF` you set. If both fail, it exits
with a clear message. `run_full_pipeline.sge` tries two more places before it
gives up: the directory you submitted from, and the parent of that directory.
In practice `qsub -V` from the repository usually works without setting
anything.

---

## 3. How to submit

**Create `logs/` before you submit.** The `#$ -o` line in every script is a
relative path, `logs/<name>_$JOB_ID.log`. SGE will reject or lose the job if
that directory does not exist in the directory you submit from.

The examples below use placeholder paths. Replace them with your own.

### Common setup

```bash
mkdir -p logs

export PROJECT=/path/to/your/study
export CONDA_ENV=/path/to/your/envs/cite_seq
export PANEL=/path/to/your/panel.yaml   # absolute path, see the note above
export DATA_ROOT=/path/to/your/cellranger_outputs
export SAMPLE_METADATA_JSON=/path/to/your/study/sample_metadata.json
export BIO_SEP_MARKER=YOUR_MARKER  # a marker name from YOUR panel YAML
export GPU_QUEUE_TOTALVI='gpu.q@your-gpu-host-1'
export GPU_QUEUE_SCVI='gpu.q@your-gpu-host-2'
```

If your site blocks `qsub -V`, pass the variables one by one with `-v` instead:

```bash
qsub -v PROJECT=...,CONDA_ENV=...,PANEL=...,DATA_ROOT=... tests/<script>.sge
```

### Run the full pipeline

```bash
mkdir -p logs
qsub -V tests/run_full_pipeline.sge
```

### End-to-end check of step 08 (cheap, no GPU)

```bash
mkdir -p logs
qsub -V tests/validate_citeseq_nf_e2e.sge
```

Reuse the saved baseline instead of rebuilding it:

```bash
qsub -V -v SKIP_BASELINE=1 tests/validate_citeseq_nf_e2e.sge
```

Full workflow instead of step 08 only (slow, retrains on GPU):

```bash
qsub -V -v RUN_FULL_DEFAULT=1 tests/validate_citeseq_nf_e2e.sge
```

### Phase A comparison

Put a known-good copy of the Phase A files here first:

```
<PROJECT>/data/processed/phase_a_baseline/
```

Then compare without running the pipeline:

```bash
mkdir -p logs
qsub -V -v SKIP_NEXTFLOW=1 tests/validate_phase_a.sge
```

To run the pipeline as well (this also runs Phase B on GPU):

```bash
qsub -V tests/validate_phase_a.sge
```

### Phase B GPU check

```bash
mkdir -p logs
qsub -V tests/validate_phase_b_gpu.sge
```

Check the existing outputs without retraining:

```bash
qsub -V -v SKIP_NEXTFLOW=1 tests/validate_phase_b_gpu.sge
```

In this mode the script compares the files with themselves, so most checks pass
by construction. It is a quick way to confirm the files exist and still agree
with `scripts/config.py` and your panel. It is not a regression test. You still
need to set both GPU queue names, even though no GPU work is submitted.

### Where the output goes

Each job writes one log file to `logs/`, named after the job and the job id.
Every script prints a banner of the settings it actually resolved. If a job
fails, read that banner first. It tells you what the job really used, which is
often not what you thought you exported.

---

## 4. What to run first

A new owner starts with nothing on disk. Two scripts can create that first
output, because they run the pipeline themselves: `run_full_pipeline.sge`
always, and `validate_phase_a.sge` unless you set `SKIP_NEXTFLOW=1`. The other
two only read files that already exist. So the order is:

1. **`validate_phase_a.sge`, with `N_LIBRARIES` set to 2 or 3.** A small pilot
   on real data. This script runs the pipeline itself, so it creates the
   Phase A files it then checks. Setup problems appear here first, and on two
   or three libraries they appear in minutes.

   > **Not yet tested on real data.** `N_LIBRARIES` used to limit only the
   > parallel fan-out, so a pilot quietly processed the whole cohort. That was
   > fixed on 2026-08-05, after the last real run. The fix works in tests, but
   > nobody has run it against real data yet. If your pilot processes more
   > libraries than you asked for, that is the first thing to check. On this first pass there is no
   `phase_a_baseline/` on disk yet, so the compare stage has nothing to compare
   against. Read this pass as a check that Phase A runs at all. Note that the
   script also runs Phase B, so you need both GPU queue names, and the job
   continues for hours after Phase A is done. Cancel it with `qdel` once you
   have the Phase A result.

2. **`run_full_pipeline.sge`,** once the pilot is clean. It runs the whole
   pipeline on the full cohort and gives the later checks something to work
   with. If it fails early, you have a setup problem, not a science problem.

3. **`validate_citeseq_nf_e2e.sge`** in its default mode. This is the cheapest
   test and it needs no GPU. It checks that the Nextflow wiring gives the same
   answer as running the step directly.

4. **`validate_phase_a.sge` again, with `SKIP_NEXTFLOW=1`.** Copy a Phase A
   result you trust into `phase_a_baseline/` first. From then on this is a real
   regression check. It only compares files, so it is fast. Use it whenever you
   change a QC step.

5. **`validate_phase_b_gpu.sge`.** Run this last. It retrains both models on a
   GPU and takes hours. Run it when you change something in step 06 or 07, or
   when you move to a new GPU host and want to confirm training still works.

Short version: run the cheap checks before the expensive GPU one. The main
`README.md` gives the same order.

---

## 5. How to read the results

Every script ends with a summary line. Look for:

- `ALL CHECKS PASSED` or `ALL INVARIANT CHECKS PASSED`, or
- `RESULT: PASS` in `run_full_pipeline.sge`, or
- `VALIDATION FAILED` with a `FAIL:` line for every check that failed.

The exit code follows the same verdict, so you can chain these in a script.

### Numbers will not match exactly, and that is normal

totalVI and scVI are stochastic. They start from random weights, shuffle the
data at random, and split train and validation at random. Two runs on the same
input give different numbers. Everything computed after them (silhouette
scores, the neighbour graph, UMAP, clustering) inherits that randomness.

That is why the validation scripts compare with a tolerance instead of
comparing files:

| Check | Tolerance |
|---|---|
| Cell count (Phase B) | within 0.5 percent |
| Gene count | exact |
| Latent size | exact, from `scripts/config.py` |
| `obsm` key set | exact |
| Last totalVI training ELBO | within 2 percent |
| Silhouette scores (step 08) | within 2 percent |

Deterministic steps have no such excuse. QC cell counts, DSB output shape and
the ThresholdR cutoffs should reproduce exactly, and `validate_phase_a.sge`
checks them that way.

### Two different meanings of "baseline"

This is the easiest thing to confuse.

| Baseline | Where | Who makes it | Used for |
|---|---|---|---|
| Committed reference tables | `./baseline/` | The 2026-04-28 reference run. Already in the repository. | Reading by eye. Compare the shape and the size of your numbers. |
| Phase A baseline | `<PROJECT>/data/processed/phase_a_baseline/` | You, by hand | `validate_phase_a.sge` |
| Step 08 baseline | `<PROJECT>/data/processed/citeseq_nf_baseline/` | `validate_citeseq_nf_e2e.sge`, automatically | The same script |
| Phase B baseline | a JSON file in `SCRATCH`, deleted at the end | `validate_phase_b_gpu.sge`, automatically | The same script |

Only the first one is in git. The other three live in your study directory or
in scratch space.

### Comparing against the reference run

Read [`../docs/reference_run.md`](../docs/reference_run.md) first. It records
the one full run known to finish cleanly, with software versions and key
numbers. Read the warning at the top of that page: the exact code revision that
produced those numbers was never committed and no longer exists, so those
numbers are indicative, not a target.

[`./baseline/README.md`](./baseline/README.md) describes each committed CSV
file, one line per file.

Practical advice when you compare:

- **Check shape first, values second.** Same columns? One row per library?
  If yes, the pipeline is wired correctly.
- **Expect similar magnitudes, not identical values.**
- **Treat a sign flip or a 10x change as a red flag.** For example, a batch
  silhouette that becomes clearly positive means batch mixing got worse.
- **Do not chase exact totalVI or scVI numbers.** Only the qualitative
  conclusion should be stable.

Once you have a run you trust, replace the files in `./baseline/` with your
own and update `../docs/reference_run.md`. Commit the code first, run second,
and write the commit hash into your run log.

---

## 6. What these tests do NOT cover

Be honest with yourself about what a green result means here.

- **They are not unit tests.** No individual Python or R function is tested on
  its own. There are no small fixtures and no synthetic data. Every script
  needs a real study directory with real outputs.
- **They do not check that the science is right.** They check that the
  pipeline runs, that files exist, and that outputs stay the same shape. A
  wrong threshold or a wrong normalisation would still pass, as long as it
  produced the same shape.
- **They do not test steps 00 to 05 in isolation.** Only
  `validate_phase_a.sge` looks at Phase A outputs, and only for the files
  listed in it: DSB, CLR, RNA normalisation, ThresholdR thresholds and the QC
  h5ad. DecontX, SoupX, scDblFinder, the KDE validation and the fitness
  metrics are not compared to anything.
- **They do not validate figures.** No plot is checked.
- **They do not test the R steps directly.** If Rscript is missing,
  `validate_phase_a.sge` prints a warning and continues. The failure appears
  later, inside the pipeline.
- **They do not check the panel YAML contents.** They check that the file
  exists. A wrong marker name in the panel is found by the pipeline, not here.
- **They do not test error handling.** No script feeds bad input on purpose to
  confirm the pipeline fails cleanly.
- **They do not test any executor other than SGE.** Every script uses
  `-profile sge` and `qsub`. Moving to Slurm or to the cloud needs new
  profiles and new submission scripts.
- **They do not pin randomness.** No random seed is fixed, so a rerun always
  differs a little. That is why the tolerances exist.
- **The committed baseline in `./baseline/` is from a code revision that no
  longer exists.** It is a sanity reference, not a contract.

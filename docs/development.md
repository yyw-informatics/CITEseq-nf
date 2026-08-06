# Development

How citeseq-nf is put together, and how to extend it.

<details><summary>Design rules</summary>

1. **Params are override-only.** Defaults for scientific values live in
   `config.py` / `config.R`, not in `nextflow.config`. (Infrastructure paths
   like `python` / `rscript` are derived from `--conda_env` if not
   explicitly set.)
2. **Modules emit outputs through channels**, never by reading from
   `publishDir`.
3. **Distinct artifact per step.** No script overwrites another's output.
4. **Fixed staged file names.** Entry-point inputs use
   `path(file, name: ...)`.
5. **Software by label.** Every process gets `label 'env_cite_seq'`.
6. **Resources by label.** `cpus` / `memory` / `time` from
   `conf/resources.config`.
7. **Validation by invariant for stochastic steps** (totalVI, scVI,
   embedding comparison). Check cell counts to ±0.5%, ELBO and silhouette
   to ±2%, and the same set of embedding keys (the `obsm` slots). Do not
   compare the files exactly.

</details>

<details><summary>How to add a new generic step</summary>

1. Create the script under `scripts/python/` or `scripts/R/`.
2. Create the module in `modules/stepNN_name.nf`:
   - `input: path(file, name: 'staged_name.h5ad')`
   - `output: path('output_name.h5ad', emit: name)`
   - `publishDir` for `data/processed/`, `results/...`, `figures/`
   - `label 'cpu_small'` (or larger) + `label 'env_cite_seq'`
   - Script block: set all four `CITESEQ_*_DIR` env vars, `mkdir -p` output
     dirs, then call the script
3. Wire into `subworkflows/phase_a_qc.nf` or `subworkflows/phase_b_embed.nf`,
   or directly in `main.nf`.
4. Add a contract header (template below).
5. Add a resource label in `conf/resources.config` if a new tier is needed.

**Contract header template** (top-of-module comment):

```groovy
/*
 * Step NN: <name>
 *
 * Contract:
 *   Input:   <staged-input-file-name>
 *   Output:  <output-file-name>
 *   Adds:    <what gets added to adata>
 *
 * Stochastic: Yes|No (determines validation strategy)
 */
```

</details>

<details><summary>Repository layout</summary>

```
citeseq-nf/
├── run.sh                   # Launch wrapper, redirects .nextflow.log* to logs/nextflow/
├── main.nf                  # Top-level workflows + entry points
├── nextflow.config          # Params (override-only) + profiles
├── conf/
│   ├── resources.config     # Resource + environment labels
│   └── sge.config           # SGE executor + GPU queue settings (cluster-specific)
├── modules/                 # One .nf per pipeline step (14 modules)
│   ├── step00b_soupx.nf
│   ├── step00c_decontx.nf
│   ├── step01_qc.nf
│   ├── step01b_scdblfinder.nf
│   ├── step02_dsb.nf
│   ├── step03_clr.nf
│   ├── step04_rna_norm.nf
│   ├── step05a_thresholdr.nf
│   ├── step05b_kde.nf
│   ├── step05c_fitness.nf
│   ├── step06_totalvi.nf
│   ├── step07_scvi.nf
│   ├── step08_embedding.nf  # handoff point, pipeline ends here
│   └── select_branch.nf
├── subworkflows/
│   ├── phase_a_qc.nf        # decontam → QC → normalization → gating
│   └── phase_b_embed.nf     # totalVI + scVI + embedding comparison
├── scripts/
│   ├── config.py            # Generic config (Python) + panel-YAML loader
│   ├── config.R             # R mirror
│   ├── python/              # 9 analytical step scripts (Python)
│   ├── R/                   # 4 analytical step scripts (R)
│   └── admin/
│       └── select_branch.py # Used by modules/select_branch.nf
├── panels/
│   └── example_panel.yaml   # Schema template, copy and edit for your panel
├── env/                     # How to build the software environment
│   ├── README.md            # Build steps, timings, known problems
│   ├── environment.yml      # conda spec (Python + R + pip packages)
│   ├── r_extra.R            # R packages conda cannot supply
│   ├── r_session_info.txt   # Full package list of the working installation
│   └── vendor/              # Vendored ThresholdR source + its provenance
│       └── README.md        # How the tarball was rebuilt, and how it differs
├── tests/                   # 4 SGE scripts + committed reference tables
│   ├── README.md            # Variable contract per script + what is not covered
│   └── baseline/            # 9 anonymized result CSVs, mostly from 2026-04-28
│       └── README.md        # One line per baseline file
├── docs/
│   └── reference_run.md     # The one known-good full run, in detail
└── logs/                    # SGE job logs land here (kept by .gitkeep)
```

`conf/` holds cluster and resource configs. `modules/` and `subworkflows/`
together define the Nextflow DAG; module scripts invoke code from
`scripts/`. `env/`, `tests/` and `docs/` cover how to build the environment,
how to check the pipeline, and what the last known good run produced.

</details>


# Vendored R packages

This folder holds source copies of R packages that the pipeline needs but
that cannot be fetched from CRAN, Bioconductor, or conda. They are kept in
the repository so the pipeline stays installable even if the original
source disappears.

| file | package | version | size | sha256 |
|---|---|---|---|---|
| `ThresholdR_0.1.0.tar.gz` | ThresholdR | 0.1.0 | 7,636 bytes | `4d37abe4e17f96e29a2b3b108c073f6781872de371f21dab855950435cd7a204` |

Install with:

```bash
Rscript -e 'install.packages("env/vendor/ThresholdR_0.1.0.tar.gz", repos = NULL, type = "source")'
```

or, together with the rest of the R stack:

```bash
THRESHOLDR_SOURCE=vendor Rscript env/r_extra.R
```

---

## ThresholdR 0.1.0

ThresholdR fits Gaussian mixture models to ADT (antibody-derived tag)
marker distributions and turns them into positive/negative thresholds. The
pipeline uses it in `scripts/R/step05a_thresholdr_gating.R`.

* License: MIT
* Author and maintainer: Mohammad Oliaeimotlagh
* Upstream: <https://github.com/MDMotlagh/ThresholdR>
* Method paper: Oliaeimotlagh et al., *Automated denoising of CITE-seq data
  with ThresholdR*, Cell Reports Methods (2025), doi:10.1016/j.crmeth.2025.101088

### Why it is vendored

ThresholdR is not on CRAN and not on Bioconductor. The copy installed in
the analysis environment was built on 2026-03-29 under R 4.4.3, and its
`DESCRIPTION` carried **no** `Repository`, `RemoteUrl`, `RemoteType` or
`RemoteSha` field. R therefore kept no record of where the package came
from. Nothing in the project recorded the install command either.

If this copy were lost, step 5a would not be rebuildable with confidence.

### This copy is not the same as upstream

This is the important part. Upstream GitHub was checked and the repository
exists, but **no upstream commit was found that matches the installed
build.** The installed build is treated as a local patch of upstream.

Comparison against upstream commit
`9b639507d76e5348afda9ff633ca6c652bbd01b3` (the newest commit on `main`
that predates the 2026-03-29 build date):

| | this vendored copy | upstream `main` |
|---|---|---|
| `get_thresholds()` cut point | private helper `.intersect2mix()` | `AdaptGauss::Intersect2Mixtures(...)$CutX` |
| `DESCRIPTION` `Imports` | mixtools, mclust, foreach, ggplot2 | mixtools, mclust, foreach, ggplot2, **AdaptGauss** |
| `calculateBIC()`, `get_k_list()`, `fit_k()`, `publish_fittings_thr()` | same behaviour as upstream | same |

The two cut-point routines are different algorithms for the same quantity:

* `.intersect2mix()` solves `w1 * dnorm(x, m1, s1) == w2 * dnorm(x, m2, s2)`
  with `uniroot()` between the two component means. If that fails it falls
  back to a 10,000-point grid search over
  `[min(mu) - 3*max(sd), max(mu) + 3*max(sd)]`.
* `AdaptGauss::Intersect2Mixtures()` scans a 20,000-point grid over
  `[min(mu) - 2*sd, max(mu) + 2*sd]` and returns the point where the two
  weighted densities differ least.

**Consequence:** the `bi_cut`, `tri_cut1` and `tri_cut2` columns produced by
`get_thresholds()` can differ between this copy and a fresh upstream
install. The `bi_thr`, `tri_thr1` and `tri_thr2` columns are the noise
component's `mean + 3*sigma` and are computed the same way in both, so they
are unaffected. Step 5a gates on the `mean + 3*sigma` rule, so the gating
itself is robust, but the exported threshold table is not byte-identical
across the two.

Also note the trailing comma: the installed `DESCRIPTION` read
`Imports: mixtools, mclust, foreach, ggplot2,`. No upstream commit had
that form, which is further evidence of local editing.

### How this tarball was produced

No source tarball was archived, so the package was rebuilt from the
installed binary:

1. The R code was recovered by deparsing the package's lazy-load database
   (`R/ThresholdR.rdb`). All six objects were recovered, including the
   unexported `.intersect2mix()`.
2. Each recovered function was re-parsed and compared back against the
   object in the installed namespace. All six round-trip to identical
   formals and bodies.
3. `man/` pages were regenerated from the installed help database, so the
   documentation text is the original.
4. `NAMESPACE`, `LICENSE` and `data/readme.md` were copied unchanged.
5. `R CMD build` produced the tarball, which was then reinstalled into a
   clean library and compared function by function against the original
   binary. All six functions are identical.

Edits made during the rebuild, all recorded in the package's own
`inst/PROVENANCE.md`:

* Dropped the install-time `Built:`, `Packaged:`, `Author:` and
  `Maintainer:` fields from `DESCRIPTION`; R regenerates them on install.
* Removed the trailing comma in `Imports:`. The dependency list is
  unchanged.
* Replaced an absolute path in the `publish_fittings_thr` example with a
  relative one. Documentation text only.

What is **not** recoverable by this method: source comments, roxygen
blocks, and the original code formatting. The executed code is identical;
the human-readable source layout is not the original.

### Runtime dependencies

`ThresholdR` calls `require()` for its dependencies from inside its
functions rather than importing them at load time. That means a missing
dependency only fails once step 5a is already running. These must be
present:

`mixtools`, `mclust`, `foreach`, `ggplot2`.

`env/r_extra.R` installs all four at pinned versions.

If you install from upstream instead of from this tarball, you also need
`AdaptGauss`.

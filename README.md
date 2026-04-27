# nlmixr2boot

`nlmixr2boot` provides nonparametric bootstrap tools for `nlmixr2`
population PK/PD models, expanded from an initial implementation by
Matt Fidler and Vipul Mann in `nlmixr2extra`.

The package focuses on subject-level bootstrap workflows for `nlmixr2` fits,
including:

* `runBootstrap()` for cached, resumable bootstrap refits
* `sampling()` for subject-level resampling with replacement
* `bootstrapNormalizedData()` for population-level covariate normalization helpers
* `bootstrapFoldGen()` for stratified fold generation utilities
* `bootstrapOptimUnisampling()` for bounded uniform sampling with an approximate target median
* `bootstrapPlot()` for bootstrap delta-OFV diagnostics

`nlmixr2boot` works alongside `nlmixr2utils`, which provides the shared
`nlmixr2()`/`ini()`/`model()` helpers and worker-plan infrastructure used by
the split `nlmixr2` extension packages.

## Installation

The package is not on CRAN. Install it together with `nlmixr2utils`.

Using `pak`:

```r
pak::pkg_install(c(
  "kestrel99/nlmixr2utils",
  "kestrel99/nlmixr2boot"
))
```

Using `remotes`:

```r
remotes::install_github("kestrel99/nlmixr2utils")
remotes::install_github("kestrel99/nlmixr2boot")
```

## Basic Use

```r
library(nlmixr2data)
library(nlmixr2utils)
library(nlmixr2boot)

one_cmt <- function() {
  ini({
    tka <- 0.45
    tcl <- 1.00
    tv <- 3.45
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}

fit <- nlmixr2(
  one_cmt,
  data = nlmixr2data::theo_sd,
  est = "focei",
  control = list(print = 0L)
)

boot <- withr::with_tempdir({
  runBootstrap(fit, nboot = 5, restart = TRUE, workers = 1L)
})

boot$summary
```

The implementation follows standard nonparametric bootstrap workflows for
population PK/PD models, adapted for `nlmixr2` fits.

## How this differs from `nlmixr2extra`

The underlying statistical workflow is the same subject-level nonparametric
bootstrap that started in `nlmixr2extra`: resample subjects with replacement,
refit the model to each bootstrap dataset, and summarize the resulting fixed
effects and OMEGA estimates. The main differences in `nlmixr2boot` are in how
the bootstrap is executed and how results are stored.

Compared with `bootstrapFit()` in `nlmixr2extra`, `runBootstrap()`:

* writes each run to a dedicated bootstrap directory such as
  `<fitName>_boot_1` (or a user-supplied `outputDir`) using the shared
  numbered run-directory helpers from `nlmixr2utils`
* records extra run artifacts, including `boot_seed.rds`, `raw_results.csv`,
  and `bootstrap_summary.txt`, so the bootstrap is easier to inspect outside
  the fitted object
* resumes from shared on-disk task caches for sampled datasets, abbreviated
  model summaries, and fitted objects while preserving failed replicates in the
  canonical raw-results table
* supports explicit `workers` control through `nlmixr2utils` worker-plan
  helpers for sequential or parallel replicate refits
* returns a standalone result object with `seed`, `results`, `summary`,
  `outputDir`, and `timestamp`, instead of only returning the modified fit or
  the legacy fit/model lists
* stores bootstrap results in a single `fit$bootstrap` list when
  `updateFit = TRUE`, rather than rewriting several fit environment objects
  such as `parFixedDf`, `parFixed`, `bootSummary`, `bootCovMatrix`, and
  `covMethod`

One important practical consequence is that `nlmixr2boot` does not
automatically replace the active covariance method on the original fit.
Instead, the bootstrap covariance and correlation matrices are stored in
`fit$bootstrap`, so you can inspect or activate them later without the
bootstrap run mutating the rest of the fit summary in place.

For a worked example, see:
`vignette("runBootstrap", package = "nlmixr2boot")`.

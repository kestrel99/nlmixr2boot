# nlmixr2boot

`nlmixr2boot` provides nonparametric bootstrap tools for `nlmixr2`
population PK/PD models.

The package focuses on subject-level bootstrap workflows for `nlmixr2` fits,
including:

* `runBootstrap()` for cached, resumable bootstrap refits
* `sampling()` for subject-level resampling with replacement
* `normalizedData()` for population-level covariate normalization helpers
* `bootplot()` for bootstrap delta-OFV diagnostics

`nlmixr2boot` works alongside `nlmixr2utils`, which provides the shared
`nlmixr2()`/`ini()`/`model()` helpers and worker-plan infrastructure used by
the split `nlmixr2` extension packages.

## Installation

The package is not on CRAN. Install it together with `nlmixr2utils`.

Using `pak`:

```r
pak::pkg_install(c(
  "nlmixr2/nlmixr2utils",
  "nlmixr2/nlmixr2boot"
))
```

Using `remotes`:

```r
remotes::install_github("nlmixr2/nlmixr2utils")
remotes::install_github("nlmixr2/nlmixr2boot")
```

If you are working locally with the split repositories:

```r
devtools::install_local("../nlmixr2utils")
devtools::install_local("../nlmixr2boot")
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

## References

The implementation follows standard nonparametric bootstrap workflows for
population PK/PD models, adapted for `nlmixr2` fits.

For a worked example, see:
`vignette("runBootstrap", package = "nlmixr2boot")`.

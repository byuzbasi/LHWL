# Professional LHWR Model Report

Print a structured report for an LHWR fitted model. The report is
designed for applied analysis: it combines the model call, calibration
settings, a global Gaussian or binomial baseline, LHWR diagnostics,
local coefficient summaries, optional local inference, optional
diagnostic F-tests, and optional residual or collinearity diagnostics.

## Usage

``` r
lhwr_report(
  object,
  include_global = TRUE,
  include_inference = TRUE,
  include_f_tests = TRUE,
  include_moran = !is.null(moran),
  include_collinearity = !is.null(collinearity),
  local_inference = NULL,
  moran = NULL,
  collinearity = NULL,
  alpha = 0.05,
  adjust = "BH",
  digits = 4,
  ...
)
```

## Arguments

- object:

  An LHWR model object returned by
  [`lhwr()`](https://byuzbasi.github.io/LHWL/reference/lhwr.md).

- include_global:

  Logical. If `TRUE`, include the global linear or logistic baseline
  stored in `object$lm` when available.

- include_inference:

  Logical. If `TRUE`, include exploratory local coefficient inference
  when available or computable from the fit.

- include_f_tests:

  Logical. If `TRUE`, include diagnostic F-tests when stored in the fit
  or computable from the fit.

- include_moran:

  Logical. If `TRUE`, include the object supplied through `moran`. The
  function does not compute Moran's I automatically.

- include_collinearity:

  Logical. If `TRUE`, include the object supplied through
  `collinearity`. The function does not compute local collinearity
  diagnostics automatically.

- local_inference:

  Optional object returned by
  [`local_inference()`](https://byuzbasi.github.io/LHWL/reference/local_inference.md).
  If `NULL` and `include_inference = TRUE`, the function attempts to
  compute it from `object`.

- moran:

  Optional object returned by
  [`moran.lhwr()`](https://byuzbasi.github.io/LHWL/reference/moran.lhwr.md)
  or
  [`moran.hgwr()`](https://byuzbasi.github.io/LHWL/reference/moran.lhwr.md).

- collinearity:

  Optional object returned by
  [`local_collinearity()`](https://byuzbasi.github.io/LHWL/reference/local_collinearity.md).

- alpha:

  Significance level used when local inference must be computed.

- adjust:

  Multiplicity adjustment used when local inference must be computed.
  Passed to
  [`local_inference()`](https://byuzbasi.github.io/LHWL/reference/local_inference.md).

- digits:

  Number of significant digits used in printed numeric values.

- ...:

  Unused.

## Value

Invisibly returns an object of class `"lhwr_report"` containing the
tables printed in the report.

## Details

The function does not refit the model. Local inference is computed only
from quantities already stored in the fitted object, and therefore
requires a fit created with `diagnostics = "full"`. Moran and local
collinearity diagnostics can be passed through the `moran` and
`collinearity` arguments when they have already been computed
separately.

## Examples

``` r
set.seed(7)
H <- simulate_H2(80, r_max = 3)
Z <- H$Z
x1 <- rnorm(80)
x2 <- rnorm(80)
y <- 1 + 2 * x1 - 0.5 * x2 + rnorm(80, sd = 0.3)
dat <- data.frame(y, x1, x2)
fit <- lhwr(y ~ x1 + x2, data = dat, Z = Z, bw = 35,
            adaptive = TRUE, diagnostics = "full", f_tests = TRUE)
lhwr_report(fit, include_moran = FALSE, include_collinearity = FALSE)
```

# Local Inference for LHWR Coefficient Surfaces

Compute exploratory local t-tests and p-values for the coefficient
surfaces estimated by an LHWR model. The function uses the local
coefficient estimates, local standard errors, and local t-values stored
in the fitted object; it does not refit the model.

## Usage

``` r
local_inference(
  object,
  alpha = 0.05,
  df = c("edf", "residual", "normal"),
  adjust = c("BH", "none", "holm", "bonferroni", "BY", "fdr"),
  adjust_scope = c("coefficient", "global"),
  ...
)

local_inference.hgwr(
  object,
  alpha = 0.05,
  df = c("edf", "residual", "normal"),
  adjust = c("BH", "none", "holm", "bonferroni", "BY", "fdr"),
  adjust_scope = c("coefficient", "global"),
  ...
)
```

## Arguments

- object:

  An LHWR model object returned by
  [`lhwr()`](https://byuzbasi.github.io/LHWL/reference/lhwr.md).

- alpha:

  Significance level used to create the `significant` flag.

- df:

  Degrees-of-freedom rule. `"edf"` uses the model effective degrees of
  freedom stored in `object$LHWR.diagnostic$edf`; `"residual"` uses
  `n - tr(S)`; `"normal"` uses a standard normal reference distribution.

- adjust:

  Multiplicity adjustment passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). Use
  `"none"` for raw p-values.

- adjust_scope:

  Adjustment scope. `"coefficient"` adjusts p-values separately within
  each coefficient surface; `"global"` adjusts over all local
  coefficient tests at once.

- ...:

  Unused.

## Value

An object of class `"lhwr_local_inference"` with components: `local`, a
long data frame of local estimates, standard errors, t-values, raw
p-values, adjusted p-values, and significance flags; `summary`, a
coefficient-level summary table; and `settings`.

## Details

These p-values should be interpreted as local exploratory diagnostics.
Because LHWR produces many local tests over dependent neighborhoods,
multiplicity and local dependence matter. Use `adjust = "BH"`, `"holm"`,
or `"bonferroni"` when the goal is to identify robust local patterns.

Inference is conditional on the Lorentz-model coordinates and
neighborhood weights stored in `object`. If the representation was
estimated, the reported standard errors and p-values do not include
first-stage representation uncertainty or uncertainty from bandwidth
selection.

## Examples

``` r
set.seed(5)
H <- simulate_H2(80, r_max = 3)
Z <- H$Z
x1 <- rnorm(80)
y <- 1 + 2 * x1 + rnorm(80, sd = 0.3)
dat <- data.frame(y, x1)
fit <- lhwr(y ~ x1, dat, Z = Z, bw = 35, adaptive = TRUE, diagnostics = "full")
inf <- local_inference(fit, adjust = "BH")
inf
head(inf$local)
```

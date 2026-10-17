# Local Multicollinearity Diagnostics for LHWR

Diagnose local multicollinearity in a fitted LHWR model. The diagnostics
follow the classical GWR weighted linear algebra discussed for local
collinearity diagnostics, including Wheeler and Tiefelsdorf (2005), but
use the LHWR kernel weights induced by the fitted distance geometry.

## Usage

``` r
local_collinearity(object, ...)

# S3 method for class 'lhwr'
local_collinearity(
  object,
  threshold_cn = 30,
  threshold_vif = 10,
  center = TRUE,
  scale = TRUE,
  return_pairwise = FALSE,
  engine = c("auto", "cpp", "R"),
  digits = 3,
  ...
)
```

## Arguments

- object:

  A fitted model object.

- ...:

  Additional arguments passed to methods.

- threshold_cn:

  Threshold used to flag large local condition numbers.

- threshold_vif:

  Threshold used to flag large local VIF values.

- center:

  Logical. If `TRUE`, locally center predictor columns before computing
  condition numbers and local correlations.

- scale:

  Logical. If `TRUE`, locally scale predictor columns before computing
  condition numbers and local correlations.

- return_pairwise:

  Logical. If `TRUE`, return all pairwise local weighted correlations in
  long format.

- engine:

  Computational engine. `"cpp"` uses a compiled kNN implementation when
  cached kNN neighborhoods are stored in the fitted object. `"R"` uses
  the reference R implementation. `"auto"` selects `"cpp"` when possible
  and otherwise falls back to `"R"`.

- digits:

  Number of digits used by the print method.

## Value

A local collinearity diagnostic object.

## References

Wheeler, D. C. and Tiefelsdorf, M. (2005). Multicollinearity and
correlation among local regression coefficients in LHWR and related
local weighted regression models. *Journal of Geographical Systems*, 7,
161–187.

## Examples

``` r
set.seed(31)
H <- simulate_H2(60, r_max = 3)
x1 <- rnorm(60)
x2 <- 0.7 * x1 + rnorm(60, sd = 0.4)
y <- 1 + 2 * x1 - 0.5 * x2 + rnorm(60, sd = 0.3)
fit <- lhwr(y ~ x1 + x2, data.frame(y, x1, x2), Z = H$Z, bw = 25)
lc <- local_collinearity(fit)
lc
```

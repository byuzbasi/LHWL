# Fit the Low-Level LHWR Model

Fit local weighted least-squares regressions from a prepared response
vector, model matrix, distance matrix, and bandwidth vector. This
function is the computational engine behind
[`lhwr()`](https://byuzbasi.github.io/LHWL/reference/lhwr.md) and
assumes that all model inputs have already been constructed.

## Usage

``` r
hgwr.fit(y, X, D, bw, kernel = hgwr_kernel_names(), return_hat = TRUE)
```

## Arguments

- y:

  Numeric response vector of length `n`.

- X:

  Numeric `n x p` model matrix. The first column should usually be an
  intercept if an intercept is desired.

- D:

  Numeric `n x n` pairwise distance matrix used to construct kernel
  weights for each focal location.

- bw:

  Positive bandwidth value or numeric vector of length `n`. A scalar is
  recycled to all focal locations; a vector supplies one radius per
  location.

- kernel:

  Kernel used to convert distances into local weights. Available choices
  are `"bisquare"`, `"gaussian"`, `"exponential"`, `"tricube"`, and
  `"boxcar"`.

- return_hat:

  Logical. If `TRUE`, compute hat-matrix diagnostics and unscaled local
  standard-error components.

## Value

A list with local coefficient matrix `beta_hat`, local standard errors
`se_beta`, local t-values `t_beta`, fitted values `y_hat`, residuals,
residual sum of squares `RSS`, hat-matrix traces `trS` and `trStS`,
`AICc`, diagnostic summaries, effective bandwidth vector, hat-matrix
diagonal, local condition numbers, local effective sample sizes, and
local counts of positive kernel weights.

## Examples

``` r
set.seed(2)
H <- simulate_H2(40, r_max = 2)
D <- hyperbolic_distance_matrix(H$Z)
X <- cbind(Intercept = 1, x = rnorm(40))
y <- 1 + 2 * X[, "x"] + rnorm(40, sd = 0.2)
bw <- rep(1.5, 40)

fit <- hgwr.fit(y, X, D, bw, kernel = "gaussian")
fit$AICc
head(fit$beta_hat)
```

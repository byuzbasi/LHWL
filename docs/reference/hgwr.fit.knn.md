# Fit the Low-Level LHWR Model from Lorentz kNN Neighborhoods

Fit adaptive k-nearest-neighbor LHWR directly from Lorentz hyperboloid
coordinates without constructing a full pairwise distance matrix. This
is a low-level performance-oriented routine. For compact kernels such as
`"bisquare"`, `"tricube"`, and `"boxcar"`, it is algebraically
equivalent to
[`hgwr.fit()`](https://byuzbasi.github.io/LHWL/reference/hgwr.fit.md)
with bandwidths equal to the Lorentz radius of the `k`th nearest
neighbor at each focal point.

## Usage

``` r
hgwr.fit.knn(
  y,
  X,
  Z,
  k,
  kernel = c("bisquare", "tricube", "boxcar"),
  return_hat = TRUE,
  tol = 1e-06,
  check = TRUE,
  n_threads = 1L
)
```

## Arguments

- y:

  Numeric response vector of length `n`.

- X:

  Numeric `n x p` model matrix.

- Z:

  Numeric `n x 3` Lorentz hyperboloid coordinate matrix.

- k:

  Adaptive nearest-neighbor count. The local fit uses the focal
  observation and its `k` nearest neighbors under Lorentz distance.

- kernel:

  Compactly supported kernel. Available exact choices are `"bisquare"`,
  `"tricube"`, and `"boxcar"`.

- return_hat:

  Logical. If `TRUE`, compute hat-matrix diagnostics and unscaled local
  standard-error components.

- tol:

  Numerical tolerance used in Lorentz distance clipping.

- check:

  Logical. If `TRUE`, validate that `Z` lies on the upper Lorentz
  hyperboloid before fitting.

- n_threads:

  Number of OpenMP threads used by the kNN local fitting loop when
  OpenMP is available. The default `1` gives deterministic serial
  execution and is used automatically on systems without OpenMP support.

## Value

A list in the same format as
[`hgwr.fit()`](https://byuzbasi.github.io/LHWL/reference/hgwr.fit.md),
with an additional `bandwidth` vector containing the focal Lorentz kNN
radii.

## Examples

``` r
set.seed(7)
H <- simulate_H2(50, r_max = 3)
X <- cbind(Intercept = 1, x = rnorm(50))
y <- 1 + 2 * X[, "x"] + rnorm(50, sd = 0.2)
fit <- hgwr.fit.knn(y, X, H$Z, k = 20)
head(fit$beta_hat)
```

# Lorentz-Hyperbolic Moran's I Diagnostic for LHWR Residuals

Compute the Lorentz-hyperbolic Moran's I diagnostic (LH-Moran's I) for
hyperbolic residual autocorrelation from a fitted LHWR model. By
default, the neighborhood weights define hyperbolic neighborhoods from
Lorentz-hyperbolic distances stored in the fitted object.

## Usage

``` r
moran.lhwr(
  object,
  W = NULL,
  distance = c("lorentz", "euclidean"),
  kernel = c("bisquare", "gaussian", "exponential"),
  bandwidth = NULL,
  adaptive = TRUE,
  k = NULL,
  permutations = 999,
  row_standardize = TRUE,
  seed = NULL
)

moran.hgwr(
  object,
  W = NULL,
  distance = c("lorentz", "euclidean"),
  kernel = c("bisquare", "gaussian", "exponential"),
  bandwidth = NULL,
  adaptive = TRUE,
  k = NULL,
  permutations = 999,
  row_standardize = TRUE,
  seed = NULL
)
```

## Arguments

- object:

  An object returned by
  [`lhwr()`](https://byuzbasi.github.io/LHWL/reference/lhwr.md).

- W:

  Optional numeric `n x n` neighborhood weight matrix. If supplied, it
  is used directly after diagonal entries are set to zero.

- distance:

  Distance used when `W = NULL`. `"lorentz"` uses intrinsic
  Lorentz-hyperbolic distance on the Lorentz hyperboloid and is the
  default. `"euclidean"` is an explicit fallback using ordinary
  Euclidean distance on the projected `(z1, z2)` coordinates.

- kernel:

  Kernel used to convert distances into weights.

- bandwidth:

  Optional fixed bandwidth. If omitted and `adaptive = TRUE`, the fitted
  model bandwidth or `k` is used as the adaptive neighbor count.

- adaptive:

  Logical. If `TRUE`, use adaptive nearest-neighbor radii.

- k:

  Optional adaptive neighbor count. Overrides `bandwidth` for adaptive
  weights.

- permutations:

  Number of random permutations for the p-value. Use zero to skip
  permutation inference.

- row_standardize:

  Logical. If `TRUE`, row-standardize weights.

- seed:

  Optional random seed used for permutation inference.

## Value

An object of class `"hgwr_moran"` containing Moran's I, expected I,
permutation p-value, and settings.

## Details

These are descriptive fitted-residual diagnostics, not exact post-fit
tests: smoothing makes residuals nonexchangeable. The global permutation
expectation is -1/(n-1) for a fixed nonconstant vector and fixed
zero-diagonal weights, not the sampling expectation of estimated LHWR
residuals. Two-sided p-values use absolute deviations from that
expectation and a plus-one correction. A count of B permutations gives
p-value increments of 1/(B+1).

Supply nondefault diagnostic geometry, kernel and bandwidth settings
explicitly. They are not all inherited from the fitted model. The direct
kernel choices are bisquare, Gaussian and exponential; use an explicit W
for other kernels. A supplied W determines the diagnostic, so default
labels in the returned object do not describe how that matrix was
constructed. Validate its dimensions, row order, finite nonnegative
entries and isolated rows. Its diagonal is zeroed, then optional row
standardization is applied.

Stored kNN lists are reused only for compact bisquare diagnostic
weights. Gaussian and exponential diagnostics use the full distance
matrix, retaining their nonzero tails beyond the adaptive k-neighbour
radius. This requires quadratic memory in the number of observations. An
explicit W remains supported for user-defined diagnostic weights. See
[`vignette("lhwr-diagnostics-and-tuning", package = "LHWL")`](https://byuzbasi.github.io/LHWL/articles/lhwr-diagnostics-and-tuning.md)
for formulas, reproducible examples and explicit weight construction.

## Examples

``` r
set.seed(12)
H <- simulate_H2(60, r_max = 3)
x <- rnorm(60)
y <- 1 + 2 * x + rnorm(60, sd = 0.3)
fit <- lhwr(y ~ x, data.frame(y, x), Z = H$Z, bw = 25)
moran.lhwr(fit, permutations = 99)
```

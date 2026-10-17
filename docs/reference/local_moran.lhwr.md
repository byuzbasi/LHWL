# Lorentz-Hyperbolic LISA Diagnostic (LH-LISA) for LHWR Residuals

Compute the Lorentz-hyperbolic LISA diagnostic (LH-LISA), that is local
Moran's I under hyperbolic neighborhoods, for residuals from a fitted
LHWR model. The default weights use Lorentz-hyperbolic distance to
define each hyperbolic neighborhood.

## Usage

``` r
local_moran.lhwr(
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

local_moran.hgwr(
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

A data frame of class `"hgwr_local_moran"` with local Moran values,
permutation p-values, z-scores, and Lorentz-hyperbolic LISA cluster
labels.

## Details

The local statistic is the centred residual times its weighted residual
lag, divided by the mean squared centred residual. High/low
classifications refer to these two centred quantities, not response
levels or coefficient signs. The full residual vector is permuted,
rather than keeping the focal residual fixed. Two-sided absolute
deviations use each location's permutation mean; p-values use a plus-one
correction. These are descriptive post-fit checks, not exact tests for
fitted residuals.

`cluster_type` uses raw p-values \<= 0.05 and does not apply
multiplicity adjustment. Changing the plot's alpha only masks stored
classifications when requested; it does not reclassify or adjust
p-values. Keep raw classifications and any explicitly adjusted display
in separate columns. Adjustment alone does not validate the
fitted-residual permutation law. Record permutations and seed
separately, as not all settings are stored in the local result.

The diagnostic weight and kernel-specific cache rules described in
[`moran.lhwr()`](https://byuzbasi.github.io/LHWL/reference/moran.lhwr.md)
apply here too. See
[`vignette("lhwr-diagnostics-and-tuning", package = "LHWL")`](https://byuzbasi.github.io/LHWL/articles/lhwr-diagnostics-and-tuning.md)
and
[`vignette("lhwr-graphics", package = "LHWL")`](https://byuzbasi.github.io/LHWL/articles/lhwr-graphics.md)
for tables, cluster maps, Moran scatterplots and explicit adjusted
display rules.

## Examples

``` r
set.seed(13)
H <- simulate_H2(60, r_max = 3)
x <- rnorm(60)
y <- 1 + 2 * x + rnorm(60, sd = 0.3)
fit <- lhwr(y ~ x, data.frame(y, x), Z = H$Z, bw = 25)
lm <- local_moran.lhwr(fit, permutations = 99)
head(lm)
plot(lm)
```

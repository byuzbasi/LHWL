# Select LHWR Bandwidth

Select the kernel bandwidth used by LHWR local regressions. Candidate
Gaussian bandwidths are scored by corrected Akaike information criterion
(`AICc`) or leave-one-out squared-error cross-validation (`CV`).
Binomial bandwidths are scored by leave-one-out Bernoulli log loss. The
selected bandwidth is the candidate with the smallest score.

## Usage

``` r
bw.lhwr(
  formula,
  data,
  Z,
  approach = c("AICc", "CV"),
  kernel = hgwr_kernel_names(),
  adaptive = TRUE,
  bandwidth_type = c("knn", "volume"),
  search = NULL,
  geometry = c("lorentz", "tangent"),
  bw_grid = NULL,
  lower = NULL,
  upper = NULL,
  tol = 1e-04,
  max_iter = 200L,
  pilot_k = NULL,
  tangent_ref = NULL,
  dMat = NULL,
  engine = c("dense", "knn"),
  n_threads = 1L,
  family = c("gaussian", "binomial"),
  control = stats::glm.control(maxit = 100L)
)
```

## Arguments

- formula:

  Model formula of the form `response ~ predictors`.

- data:

  Data frame containing the response and predictor variables used in
  `formula`.

- Z:

  Numeric `n x 3` matrix of Lorentz hyperboloid coordinates. Row `i`
  gives the location of observation `i` and must satisfy
  `<z_i,z_i>_L = -1` with positive time coordinate.

- approach:

  Bandwidth selection criterion. `"AICc"` minimizes the corrected Akaike
  information criterion computed from the local fitted values and
  hat-matrix trace. `"CV"` minimizes the leave-one-out prediction error.

- kernel:

  Kernel used to convert distances into local weights. Available choices
  are `"bisquare"`, `"gaussian"`, `"exponential"`, `"tricube"`, and
  `"boxcar"`.

- adaptive:

  Logical. If `TRUE`, scalar bandwidths are interpreted as neighbor
  counts and converted to location-specific radii. If `FALSE`, scalar
  bandwidths are interpreted as fixed radii, or fixed hyperbolic volumes
  when `bandwidth_type = "volume"`.

- bandwidth_type:

  Bandwidth scale. `"knn"` uses the radius to the requested nearest
  neighbor. `"volume"` uses a hyperbolic ball-volume interpretation and
  is available only for `geometry = "lorentz"`.

- search:

  Bandwidth search method. `"golden"` uses golden-section optimization
  over an automatically chosen interval. `"grid"` evaluates the supplied
  candidate values in `bw_grid`. If omitted, `"golden"` is used when
  `bw_grid` is `NULL`, and `"grid"` is used when `bw_grid` is supplied.

- geometry:

  Distance geometry used to build kernel weights. `"lorentz"` uses
  intrinsic hyperbolic distances on the Lorentz hyperboloid. `"tangent"`
  maps locations to one tangent plane by the Lorentz logarithmic map and
  uses the positive-definite Lorentz metric induced on that tangent
  plane. Equivalently, these distances can be represented as Euclidean
  distances in a Lorentz-orthonormal tangent basis.

- bw_grid:

  Optional numeric vector of candidate bandwidths for `search = "grid"`.
  In adaptive fits these are neighbor counts; in fixed fits these are
  radii unless `bandwidth_type = "volume"`.

- lower:

  Optional lower search bound for `search = "golden"`.

- upper:

  Optional upper search bound for `search = "golden"`.

- tol:

  Golden-section stopping tolerance.

- max_iter:

  Maximum number of golden-section iterations.

- pilot_k:

  Pilot neighbor count used to estimate local point density for adaptive
  volume bandwidths. Ignored unless `adaptive = TRUE` and
  `bandwidth_type = "volume"`.

- tangent_ref:

  Optional reference point for `geometry = "tangent"`. If omitted, the
  normalized extrinsic Lorentz mean is used.

- dMat:

  Optional precomputed `n x n` distance matrix. If supplied, it is used
  directly and `geometry` is used only for labeling and validation of
  related arguments.

- engine:

  Computational engine. `"dense"` builds or uses a full distance matrix.
  `"knn"` uses exact Lorentz nearest-neighbor lists and avoids storing
  the full distance matrix; it is available for adaptive Lorentz kNN
  fits with compact kernels.

- n_threads:

  Number of OpenMP threads used by `engine = "knn"` when OpenMP is
  available. The default `1` preserves serial execution.

- family:

  Response family. Gaussian LHWR supports `"AICc"` and `"CV"`; binomial
  LHWR supports `"CV"` and uses leave-one-out Bernoulli log loss. When
  `family = "binomial"` and `approach` is omitted, `"CV"` is selected.

- control:

  Iteration controls for binomial LHWR, created by
  [`stats::glm.control()`](https://rdrr.io/r/stats/glm.control.html).
  Ignored for Gaussian bandwidth selection.

## Value

An object of class `"bw.lhwr"`, a list containing the selected bandwidth
(`bw`), selection criterion (`approach`), kernel, geometry, search
method, selected score, AICc/CV fields when applicable, hat-matrix trace
(`trS`), residual sum of squares (`RSS`), the evaluated
`selection.table`, and the original function call.

## Details

AICc and CV are selection criteria; post-fit AIC and BIC are not
additional choices for `approach`. Gaussian CV is a sum of squared
leave-one-out errors, using the supplied representation and candidate
radii. It does not relearn the representation or reselect bandwidth
within an outer validation loop. AICc is infinite when its denominator
n - tr(S) - 2 is nonpositive.

The adaptive neighbour rank is not a guarantee of exactly that many
positive kernel weights. Compact kernels have zero boundary weight; ties
and coincident locations matter. Inspect local positive-weight counts
and effective support. Grid search identifies the minimum among
evaluated candidates; review boundary selections before making a
scientifically justified change to the search. See
[`vignette("lhwr-diagnostics-and-tuning", package = "LHWL")`](https://byuzbasi.github.io/LHWL/articles/lhwr-diagnostics-and-tuning.md)
for the criterion formulas, kernel curves, volume rules and engine
compatibility.

## Examples

``` r
set.seed(20260615)
n <- 120
H <- simulate_H2(n, r_max = 5, uniform_area = TRUE)
Z <- H$Z
D <- hyperbolic_distance_matrix(Z, check = FALSE)

true_k <- 25
h <- apply(D, 1, function(d) sort(d[d > 0])[true_k])
centers <- round(seq(10, n - 10, length.out = 5))
Phi <- sapply(centers, function(j) {
  hgwr_kernel(D[, j], h[j], kernel = "bisquare")
})

beta0 <- 0.4 + 2 * (1.45 * Phi[, 1] - 1.25 * Phi[, 4])
beta1 <- 1.0 + 2 * (1.65 * Phi[, 2] - 1.20 * Phi[, 5])
beta2 <- -0.8 + 2 * (1.55 * Phi[, 3] - 1.10 * Phi[, 4])

x1 <- as.numeric(scale(rnorm(n)))
x2 <- as.numeric(scale(rnorm(n)))
eta <- beta0 + beta1 * x1 + beta2 * x2
y <- eta + rnorm(n, sd = sqrt(var(eta) / 8))
dat <- data.frame(y = y, x1 = x1, x2 = x2)

bw <- bw.lhwr(
  y ~ x1 + x2,
  data = dat,
  Z = Z,
  approach = "CV",
  kernel = "bisquare",
  adaptive = TRUE,
  geometry = "lorentz",
  bw_grid = seq(10, 80, by = 5)
)
bw
head(bw$selection.table)
```

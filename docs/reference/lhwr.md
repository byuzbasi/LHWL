# Fit a Lorentz Hyperbolic Weighted Regression Model

Estimate location-varying Gaussian or binomial regression coefficients
for observations located on the Lorentz hyperboloid. At each focal
point, the function fits an unpenalized weighted linear or logistic
regression using kernel weights derived from either intrinsic Lorentz
distances or tangent-space distances.

## Usage

``` r
lhwr(
  formula,
  data,
  Z,
  bw,
  kernel = hgwr_kernel_names(),
  adaptive = TRUE,
  bandwidth_type = c("knn", "volume"),
  geometry = c("lorentz", "tangent"),
  pilot_k = NULL,
  tangent_ref = NULL,
  dMat = NULL,
  engine = c("dense", "knn"),
  n_threads = 1L,
  diagnostics = c("full", "aicc", "none"),
  f_tests = FALSE,
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

- bw:

  Bandwidth used for local weighting. If `adaptive = TRUE`, a scalar is
  interpreted as a neighbor count and converted to a separate radius at
  each focal location. If `adaptive = FALSE`, a scalar is interpreted as
  a fixed radius, or as a fixed hyperbolic volume when
  `bandwidth_type = "volume"`. A numeric vector of length `n` may be
  used to supply one radius per location directly.

- kernel:

  Kernel used to convert distances into local weights. Available choices
  are `"bisquare"`, `"gaussian"`, `"exponential"`, `"tricube"`, and
  `"boxcar"`.

- adaptive:

  Logical. If `TRUE`, scalar `bw` is treated as a neighbor count. If
  `FALSE`, scalar `bw` is treated as a fixed distance or volume.

- bandwidth_type:

  Bandwidth scale. `"knn"` uses nearest-neighbor radii. `"volume"` uses
  hyperbolic ball volumes and is available only for
  `geometry = "lorentz"`.

- geometry:

  Distance geometry used to build kernel weights. `"lorentz"` uses
  intrinsic hyperbolic distances on the Lorentz hyperboloid. `"tangent"`
  maps locations to one tangent plane by the Lorentz logarithmic map and
  uses the positive-definite Lorentz metric induced on that tangent
  plane. Equivalently, these distances can be represented as Euclidean
  distances in a Lorentz-orthonormal tangent basis.

- pilot_k:

  Pilot neighbor count used to estimate local point density for adaptive
  volume bandwidths. Ignored unless `adaptive = TRUE` and
  `bandwidth_type = "volume"`.

- tangent_ref:

  Optional reference point for `geometry = "tangent"`. If omitted, the
  normalized extrinsic Lorentz mean is used.

- dMat:

  Optional precomputed `n x n` distance matrix. If supplied, it is used
  directly instead of computing a distance matrix from `Z`.

- engine:

  Computational engine. `"dense"` builds or uses a full distance matrix.
  `"knn"` uses exact Lorentz nearest-neighbor lists and avoids storing
  the full distance matrix; it is available for adaptive Lorentz kNN
  fits with compact kernels.

- n_threads:

  Number of OpenMP threads used by `engine = "knn"` when OpenMP is
  available. The default `1` preserves serial execution.

- diagnostics:

  Diagnostic detail level. `"full"` returns the traditional diagnostic
  table, local standard errors, t-values, and local R-squared values
  when available. `"aicc"` computes hat-matrix diagnostics for AICc but
  skips local standard-error and local R-squared columns in the output
  table. `"none"` skips hat-matrix diagnostics and returns only fitted
  values, residuals, local coefficients, RSS, and global R-squared.

- f_tests:

  Logical. If `TRUE`, compute LHWR diagnostic F-tests after model
  fitting. The default `FALSE` avoids the additional smoother
  reconstruction cost, which can be substantial for large data sets.
  Diagnostic F-tests are available only for Gaussian LHWR.

- family:

  Response family. `"gaussian"` fits local weighted least-squares
  regressions and `"binomial"` fits local weighted logistic regressions.
  The default preserves the original Gaussian LHWR behavior.

- control:

  Iteration controls for binomial LHWR, created by
  [`stats::glm.control()`](https://rdrr.io/r/stats/glm.control.html).
  `epsilon`, `maxit`, and `trace` are used by the native local IRLS
  solver. Ignored for Gaussian LHWR.

## Value

An object of class `"lhwr"`, a list containing model arguments,
diagnostic summaries, the global linear model fit, local coefficient
table `SDF`, coefficient names, local coefficient matrix `beta_hat`,
fitted values, residuals, distance matrix, selected geometry label,
effective bandwidth vector, hat-matrix diagonal, local condition
numbers, local effective sample sizes, local counts of positive kernel
weights, timing information, family-specific convergence information,
and the original function call. For binomial models, fitted values are
response probabilities and residuals are response residuals `y - p`.

## Details

Standard errors, local test statistics, and fitted-model diagnostics are
conditional on the supplied `Z` matrix and the fitted neighborhood
system. If `Z` is estimated from data, these quantities do not propagate
first-stage representation uncertainty. They also do not include
uncertainty from data-driven bandwidth selection.

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
  adaptive = TRUE,
  geometry = "lorentz",
  bw_grid = seq(10, 80, by = 5)
)

fit <- lhwr(
  y ~ x1 + x2,
  data = dat,
  Z = Z,
  bw = bw$bw,
  adaptive = TRUE,
  geometry = "lorentz"
)
fit
summary(fit)
head(fit$SDF)

beta_true <- cbind(Intercept = beta0, x1 = beta1, x2 = beta2)
sqrt(mean(rowSums((fit$beta_hat - beta_true)^2)))

fit_tangent <- lhwr(
  y ~ x1 + x2,
  data = dat,
  Z = Z,
  bw = bw$bw,
  adaptive = TRUE,
  geometry = "tangent"
)
fit_tangent$geometry

probability <- plogis(-0.4 + 1.1 * x1 - 0.7 * x2)
dat$event <- rbinom(n, 1, probability)
fit_binomial <- lhwr(
  event ~ x1 + x2,
  data = dat,
  Z = Z,
  bw = 40,
  family = "binomial",
  diagnostics = "full"
)
head(fitted(fit_binomial))
```

# LHWR Diagnostic F-Tests

Compute Leung-type, moment-matched diagnostic F-tests adapted to a
fitted LHWR model. The diagnostics use the LHWR smoother induced by the
model's kernel weights and distance geometry. Their reference
distributions use generally non-integer effective degrees of freedom.
The reference laws are approximate and condition on the fitted geometry
and selected bandwidth; they should not be read as post-selection
confirmatory tests. They are optional because coefficient-specific F3
diagnostics require reconstructing local coefficient operators and can
be costly for large data sets.

## Usage

``` r
diagnostic_tests(object, ...)

# S3 method for class 'lhwr'
diagnostic_tests(
  object,
  ols_model = NULL,
  digits = 3,
  engine = c("auto", "cpp", "R"),
  ...
)
```

## Arguments

- object:

  A fitted model object.

- ...:

  Additional arguments passed to methods.

- ols_model:

  Optional global linear model. If omitted, the `lm` object stored by
  [`lhwr()`](https://byuzbasi.github.io/LHWL/reference/lhwr.md) is used.

- digits:

  Number of digits used by the print method.

- engine:

  Computational engine for reconstructing F-test smoother components.
  `"cpp"` uses the compiled implementation, including the memory-saving
  kNN path when the model was fitted with `engine = "knn"`. If the
  fitted object stores kNN neighborhoods, the C++ engine reuses them
  instead of recomputing nearest neighbors; `"R"` uses the reference R
  implementation, and `"auto"` selects `"cpp"`.

## Value

A diagnostic test object.

## Details

F1 uses a lower-tail moment-matched F reference; F2 uses an upper-tail
reference for improvement over the stored global OLS model. F3 uses an
upper-tail reference for each coefficient's constancy across locations,
not for that coefficient being zero at one location. The quadratic forms
are generally dependent, so these are not exact nested-model ANOVA
tests. Independent homoskedastic Gaussian errors and fixed
stationary-model weights motivate the approximations; bandwidth
selection, correlated errors and varying nuisance coefficients limit
their interpretation.

Required moments and reference degrees of freedom must be positive; do
not interpret undefined p-values as zero. When supplying `ols_model`,
verify the same response, predictor columns and observation order, not
only row count. The default stored OLS model avoids a separate
sample/formula choice. See
[`vignette("lhwr-diagnostics-and-tuning", package = "LHWL")`](https://byuzbasi.github.io/LHWL/articles/lhwr-diagnostics-and-tuning.md)
for the statistics, reference degrees of freedom and executable
reporting tables.

## Examples

``` r
# Example 1: globally stationary coefficients.
# The diagnostics should usually not reject global adequacy or local
# coefficient constancy because the response does not depend on location.
set.seed(1)
H <- simulate_H2(80, r_max = 3)
Z <- H$Z
x1 <- rnorm(80)
x2 <- rnorm(80)
y <- 1 + 2 * x1 - 0.5 * x2 + rnorm(80, sd = 0.3)
dat <- data.frame(y, x1, x2)

fit_global <- lhwr(
  y ~ x1 + x2,
  data = dat,
  Z = Z,
  bw = 30,
  adaptive = TRUE,
  diagnostics = "aicc",
  f_tests = TRUE
)
fit_global$Ftests

# Example 2: hyperbolic coefficient heterogeneity.
# The coefficients are generated from localized bumps defined by
# Lorentz-hyperbolic distances, so the diagnostics may detect departures
# from the global model and coefficient-specific local variability.
set.seed(2)
H <- simulate_H2(80, r_max = 3, uniform_area = TRUE)
Z <- H$Z
D <- hyperbolic_distance_matrix(Z, check = FALSE)
x1 <- rnorm(80)
x2 <- rnorm(80)

phi1 <- exp(-D[, 15]^2 / 0.7^2)
phi2 <- exp(-D[, 45]^2 / 0.7^2)
phi3 <- exp(-D[, 70]^2 / 0.7^2)
beta0 <- 1 + 1.2 * phi1
beta1 <- 2 + 1.5 * phi2
beta2 <- -0.5 - 1.3 * phi3

eta <- beta0 + beta1 * x1 + beta2 * x2
y <- eta + rnorm(80, sd = 0.2)
dat <- data.frame(y, x1, x2)

fit_local <- lhwr(
  y ~ x1 + x2,
  data = dat,
  Z = Z,
  bw = 25,
  adaptive = TRUE,
  diagnostics = "aicc",
  f_tests = TRUE
)
diagnostic_tests(fit_local)
```

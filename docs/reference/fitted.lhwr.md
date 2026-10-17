# Extract and Predict from LHWR Models

Convenience S3 methods for fitted LHWR models. These methods make LHWR
objects behave like ordinary R model objects for common downstream
tasks. `predict.lhwr()` returns fitted values for the calibration data
when `newdata` is omitted. For new observations, supply both the
Euclidean covariates in `newdata` and their Lorentz-hyperboloid
coordinates in `Znew`. Prediction is then obtained by fitting a local
weighted Gaussian or binomial model at each new Lorentz location using
only the original calibration data.

## Usage

``` r
# S3 method for class 'lhwr'
fitted(object, ...)

# S3 method for class 'lhwr'
residuals(object, ...)

# S3 method for class 'lhwr'
coef(object, ...)

# S3 method for class 'lhwr'
predict(
  object,
  newdata = NULL,
  Znew = NULL,
  type = c("response", "link", "coefficients"),
  dMat = NULL,
  bw = NULL,
  kernel = NULL,
  adaptive = NULL,
  bandwidth_type = NULL,
  geometry = NULL,
  ...
)

# S3 method for class 'hgwr'
fitted(object, ...)

# S3 method for class 'hgwr'
residuals(object, ...)

# S3 method for class 'hgwr'
coef(object, ...)

# S3 method for class 'hgwr'
predict(
  object,
  newdata = NULL,
  Znew = NULL,
  type = c("response", "link", "coefficients"),
  dMat = NULL,
  bw = NULL,
  kernel = NULL,
  adaptive = NULL,
  bandwidth_type = NULL,
  geometry = NULL,
  ...
)
```

## Arguments

- object:

  An LHWR model object returned by
  [`lhwr()`](https://byuzbasi.github.io/LHWL/reference/lhwr.md).

- ...:

  Unused.

- newdata:

  Optional data frame containing predictor variables for new
  observations. If omitted,
  [`predict()`](https://rdrr.io/r/stats/predict.html) returns fitted
  values from the calibration data.

- Znew:

  Optional `m x 3` Lorentz-hyperboloid coordinate matrix for `newdata`.
  Required when `newdata` is supplied.

- type:

  Prediction type. `"response"` returns fitted responses (probabilities
  for binomial LHWR), `"link"` returns linear predictors, and
  `"coefficients"` returns the local coefficient matrix used for
  prediction.

- dMat:

  Optional `m x n` distance matrix from new observations to the original
  calibration observations. If supplied, it is used directly.

- bw:

  Optional prediction bandwidth. By default the bandwidth stored in the
  fitted model is used.

- kernel:

  Optional kernel name. By default the fitted model kernel is used.

- adaptive:

  Optional logical bandwidth flag. By default the fitted model setting
  is used.

- bandwidth_type:

  Optional bandwidth type. By default the fitted model setting is used.

- geometry:

  Optional prediction geometry. By default the fitted model geometry is
  used.

## Value

[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) returns the
fitted response vector;
[`residuals()`](https://rdrr.io/r/stats/residuals.html) returns
residuals; [`coef()`](https://rdrr.io/r/stats/coef.html) returns the
`n x p` calibration local coefficient matrix.
[`predict()`](https://rdrr.io/r/stats/predict.html) returns either a
numeric response vector or an `m x p` matrix of local coefficients for
`newdata`.

## Details

Tangent-E prediction reuses the stored training reference point for
every prediction batch. If an older object has no stored reference, it
is reconstructed from `object$Z`, never from `Znew`. Identical mapped
points have exactly zero distance; nearby distinct tangent coordinates
are not rounded to zero by a pairwise-distance tolerance. With unchanged
model settings, predicting the training rows reproduces the calibration
fit up to floating-point error.

## Examples

``` r
set.seed(11)
H <- simulate_H2(50, r_max = 3)
x <- rnorm(50)
y <- 1 + 2 * x + rnorm(50, sd = 0.3)
fit <- lhwr(y ~ x, data.frame(y, x), Z = H$Z, bw = 20)
head(fitted(fit))
head(residuals(fit))
head(coef(fit))

Znew <- simulate_H2(3, r_max = 3)$Z
newdata <- data.frame(x = c(-1, 0, 1))
predict(fit, newdata = newdata, Znew = Znew)
```

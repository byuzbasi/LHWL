# Lorentz-Hyperbolic Weighted Penalized Regression

Fit local Gaussian or binomial penalized regressions whose observation
weights are defined by Lorentz-hyperbolic neighborhoods.

## Usage

``` r
lhwpr(
  formula,
  data,
  Z,
  bw,
  family = c("gaussian", "binomial"),
  penalty = c("lasso", "mcp", "scad"),
  alpha = 1,
  lambda = NULL,
  nlambda = 100L,
  lambda_min_ratio = NULL,
  gamma = NULL,
  penalty_factor = NULL,
  case_weights = NULL,
  standardize = TRUE,
  kernel = hgwr_kernel_names(),
  adaptive = TRUE,
  bandwidth_type = c("knn", "volume"),
  geometry = c("lorentz", "tangent"),
  tangent_ref = NULL,
  dMat = NULL,
  engine = c("dense", "knn"),
  keep_path = c("auto", "all", "selected"),
  control = lhwpr_control()
)

# S3 method for class 'lhwpr'
print(x, ...)

# S3 method for class 'lhwpr'
summary(object, ...)

# S3 method for class 'summary.lhwpr'
print(x, ...)

# S3 method for class 'lhwpr'
coef(object, s = NULL, ...)

# S3 method for class 'lhwpr'
fitted(object, ...)

# S3 method for class 'lhwpr'
residuals(object, ...)

# S3 method for class 'lhwpr'
predict(
  object,
  newdata = NULL,
  Znew = NULL,
  type = c("response", "link", "coefficients"),
  s = NULL,
  ...
)

# S3 method for class 'lhwpr'
plot(
  x,
  focus = 1L,
  xvar = c("lambda", "lambda_ratio"),
  ...
)
```

## Arguments

- formula:

  A model formula with one response and one or more predictors.

- data:

  A data frame containing variables in `formula`.

- Z:

  An `n` by 3 matrix of Lorentz-hyperboloid coordinates.

- bw:

  Kernel bandwidth, or neighbor count for adaptive kNN fitting.

- family:

  Gaussian or binomial response family.

- penalty:

  Structural penalty: lasso, MCP, or SCAD. It is ignored when
  `alpha = 0`.

- alpha:

  Mixing value in `[0, 1]`. Zero is ridge; a value strictly between zero
  and one is elastic net and requires the lasso penalty; one gives
  lasso, MCP, or SCAD.

- lambda:

  Optional strictly decreasing common absolute lambda path.

- nlambda:

  Number of automatically generated lambda values.

- lambda_min_ratio:

  Ratio of the final to initial automatic lambda.

- gamma:

  MCP or SCAD concavity parameter.

- penalty_factor:

  Optional non-negative multiplier for each non-intercept coefficient.

- case_weights:

  Optional non-negative case weights.

- standardize:

  Whether predictors are centered and scaled once on the full analysis
  sample.

- kernel:

  Lorentz-distance kernel.

- adaptive:

  Whether `bw` is adaptive.

- bandwidth_type:

  Neighbor-count or hyperbolic-volume bandwidth scale.

- geometry:

  Intrinsic Lorentz or tangent-space distance geometry.

- tangent_ref:

  Optional Lorentz point defining the tangent plane.

- dMat:

  Optional precomputed distance matrix.

- engine:

  Dense or exact compact-kernel kNN computation.

- keep_path:

  Whether to retain all local coefficient paths.

- control:

  A control object from
  [`lhwpr_control()`](https://byuzbasi.github.io/LHWL/reference/lhwpr_control.md).

- x, object:

  A fitted `lhwpr` object, or its summary.

- s:

  A lambda value exactly represented in the fitted path.

- newdata:

  Optional predictor data for new focal locations.

- Znew:

  Lorentz coordinates corresponding to `newdata`.

- type:

  Prediction scale or local coefficients.

- focus:

  One-based focal-location index for a path plot.

- xvar:

  Horizontal path scale.

- ...:

  Additional arguments passed to the plotting method.

## Details

At every focal location, kernel and case weights are multiplied and
normalized before fitting. Predictor standardization is global rather
than focal. For a given bandwidth, every focal model uses the same
absolute lambda values.

Lasso and elastic-net paths use warm starts, sequential strong
screening, and a full KKT repair pass. MCP and SCAD use curvature-aware
scalar minimization and full stationarity scans. Binomial models use
proximal Newton/IRLS.

The intercept is never penalized. Group penalties are not part of LHWPR.
Classical unpenalized LHWR inference is not attached to penalized fits.

With `keep_path = "auto"`, a full dense coefficient cube is retained
only when it contains at most ten million values. Selected coefficients
and all fitted-value paths are retained regardless. Use
`keep_path = "all"` to override this memory guard.

## Value

`lhwpr()` returns an object of class `"lhwpr"` with selected local
coefficients and predictions, lambda metadata, optional paths,
convergence and KKT diagnostics, preprocessing metadata, and geometry
information. Extractor methods return the corresponding component.

## Examples

``` r
set.seed(1)
H <- simulate_H2(40, r_max = 2)
dat <- data.frame(
  y = rnorm(40),
  x1 = rnorm(40),
  x2 = rnorm(40)
)
fit <- lhwpr(
  y ~ x1 + x2,
  data = dat,
  Z = H$Z,
  bw = 15,
  alpha = 0.5,
  nlambda = 10
)
coef(fit)
```

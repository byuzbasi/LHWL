# LHWR

`LHWR` is an R package for **Lorentz Hyperbolic Weighted
Regression**, a local weighted regression framework on the two-dimensional
Lorentz hyperboloid model. The response and predictors are real-valued; the
local neighborhood structure is defined by Lorentz-hyperboloid coordinates
`Z` and intrinsic hyperbolic distances.

## Quick Start

```r
library(LHWR)

# Build Lorentz hyperboloid coordinates from any meaningful
# two-component representation of the observations.
Z <- make_lhwr_Z_features(dat[, c("feature_score_1", "feature_score_2")])

# Select an adaptive bandwidth by AICc or CV.
bw <- bw.lhwr(y ~ x1 + x2, data = dat, Z = Z, approach = "AICc")

# Fit the local weighted regression on the Lorentz hyperboloid.
fit <- lhwr(y ~ x1 + x2, data = dat, Z = Z, bw = bw$bw)
summary(fit)
plot(fit, coefficient = "x1")
```

The package currently provides:

- Lorentz inner products, norms, and hyperbolic distances
- Hyperboloid simulation, projection, and validation helpers
- Gaussian, exponential, bisquare, tricube, and boxcar LHWR kernels
- AICc and CV bandwidth selection
- Adaptive kNN and volume-style bandwidth handling
- Local weighted LHWR fitting with diagnostics
- Rcpp/RcppArmadillo compute kernels for Lorentz distances, local WLS, CV,
  and bandwidth grid evaluation
- Basic plotting, tests, and an introductory vignette


## Does LHWR require geographic coordinates?

No. Classical geographically weighted regression defines locality from
longitude-latitude or other spatial coordinates. LHWR defines locality from
Lorentz-hyperboloid neighborhoods. The required `Z` matrix is an `n x 3` matrix
whose rows lie on

```text
H^2 = { z in R^{1,2}: -z0^2 + z1^2 + z2^2 = -1, z0 > 0 }.
```

Thus `Z` may come from geographic coordinates, latent feature maps, similarity
embeddings, trade or social networks, clinical representations, functional data
scores, or shape-analysis scores. For non-geographic data, LHWR should be
interpreted as locally weighted regression over a latent hyperbolic geometry;
the model is useful only when the chosen representation captures meaningful
neighborhood relations among observations.

## Constructing Z for LHWR

LHWR provides helpers for turning raw representations into Lorentz hyperboloid
coordinates. If each observation has a meaningful two-component representation,
`make_lhwr_Z_features()` maps those two components to
`(sqrt(1 + z1^2 + z2^2), z1, z2)`:

```r
X_embed <- dat[, c("feature_score_1", "feature_score_2")]
Z <- make_lhwr_Z_features(X_embed, scale = TRUE)
check_lhwr_Z(Z)
```

For relation data, `make_lhwr_Z_network()` first embeds an adjacency,
similarity, or network matrix into two components, then maps the result to the
Lorentz hyperboloid:

```r
Z <- make_lhwr_Z_network(A, method = "spectral", scale = TRUE)
```

These functions are intentionally representation-agnostic. The two components
may summarize geographic coordinates, economic indicators, financial indicators,
clinical variables, functional data scores, shape/SRVF scores, similarity
matrices, trade networks, social networks, transportation networks, or general
graph embeddings. The list is not restrictive; the key requirement is that the
representation defines a scientifically meaningful neighborhood geometry.

```r
fit <- lhwr(
  y ~ x1 + x2,
  data = dat,
  Z = Z,
  bw = 40,
  kernel = "bisquare",
  adaptive = TRUE
)
```


## Plotting without confusion

LHWR plotting functions are designed to work with simple defaults. For most
analyses, users can call a single plotting function and obtain a clean figure:

```r
plot(fit, type = "coefficient", coefficient = "x1")
plot_lhwr_embedding(Z, values = y, labels = rownames(dat), title = NULL)
```

For manuscripts or presentations, the same figures can be adjusted without
rewriting the plotting code. Use `lhwr_plot_control()` to change label, legend,
colorbar, title, and point-size scaling in one place:

```r
ctl <- lhwr_plot_control(
  label_cex = 0.65,
  legend_cex = 0.65,
  colorbar_cex = 0.65,
  title_cex = 1.10
)

plot_lhwr_embedding(
  Z,
  values = y,
  labels = dat$country_code,
  value_label = "GDP growth",
  title = NULL,
  show_size_legend = FALSE,
  show_info = FALSE,
  plot_control = ctl,
  output_file = "lhwr_map.png"
)
```

Direct arguments such as `label_cex`, `legend_cex`, `colorbar_cex`, and
`title_cex` can also be supplied directly. These only change the visual display;
they do not alter Lorentz coordinates, hyperbolic distances, bandwidths, or
model estimates.

## Hyperbolic Residual Autocorrelation Diagnostics

Fitted LHWR residuals can be checked for remaining residual autocorrelation
under hyperbolic neighborhoods using the Lorentz-hyperbolic Moran's I
diagnostic (LH-Moran's I) and the Lorentz-hyperbolic LISA diagnostic
(LH-LISA). By default, the diagnostics construct weights from intrinsic
Lorentz-hyperbolic distances stored in the model object.

```r
fit <- lhwr(
  y ~ x1 + x2,
  data = dat,
  Z = Z,
  bw = 40,
  adaptive = TRUE
)

moran.lhwr(fit, permutations = 99)

lisa <- local_moran.lhwr(fit, permutations = 99)
head(lisa)
plot(lisa)
```

The diagnostics also accept a user-supplied weight matrix `W` or a Euclidean
fallback, requested explicitly, based on the projected `(z1, z2)` coordinates:

```r
moran.lhwr(fit, distance = "euclidean", k = 20)
local_moran.lhwr(fit, W = W)
```

## Diagnostic F-tests

Diagnostic F-tests are optional because they reconstruct smoother components
for the fitted LHWR model. They are off by default and can be requested
at fit time:

```r
fit <- lhwr(
  y ~ x1 + x2,
  data = dat,
  Z = Z,
  bw = 40,
  adaptive = TRUE,
  f_tests = TRUE
)

fit$Ftests
diagnostic_tests(fit)
```

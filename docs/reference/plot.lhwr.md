# Plot LHWR Objects

Plot local model outputs from a fitted LHWR model or plot a bandwidth
selection curve. Model plots can show local coefficients, standard
errors, local t-values, local p-values, local significance flags, fitted
values, residuals, or local R-squared values. The same quantity can be
drawn either on the projected `(z1, z2)` coordinates or on a perspective
view of the Lorentz hyperboloid.

## Usage

``` r
# S3 method for class 'lhwr'
plot(
  x,
  type = c("coefficient", "se", "t_value", "p_value", "p_adjusted", "significance",
    "residual", "fitted", "local_r2"),
  coefficient = NULL,
  Z = NULL,
  geometry = c("projection", "hyperboloid"),
  palette = "Viridis",
  point_size = 1.1,
  main = NULL,
  inference_df = c("edf", "residual", "normal"),
  p_adjust = c("BH", "none", "holm", "bonferroni", "BY", "fdr"),
  alpha = 0.05,
  ...
)
```

## Arguments

- x:

  An object returned by
  [`lhwr()`](https://byuzbasi.github.io/LHWL/reference/lhwr.md) or
  [`bw.lhwr()`](https://byuzbasi.github.io/LHWL/reference/bw.lhwr.md).

- type:

  Quantity to plot for model objects. `"coefficient"` plots a local
  coefficient estimate. `"se"` plots its local standard error.
  `"t_value"` plots its local t-value. `"p_value"` and `"p_adjusted"`
  plot local inference p-values. `"significance"` plots a 0/1 local
  significance flag based on adjusted p-values. `"residual"`,
  `"fitted"`, and `"local_r2"` plot model diagnostics stored in the
  fitted object.

- coefficient:

  Coefficient name to plot for coefficient-level plot types. Defaults to
  the first local coefficient.

- Z:

  Optional `n x 3` Lorentz coordinate matrix. If omitted for model
  objects, the coordinates stored in the fitted object are used.

- geometry:

  Plotting geometry. `"projection"` draws a fast two-dimensional map
  using `(z1, z2)`. `"hyperboloid"` draws the upper Lorentz hyperboloid
  surface `x0 = sqrt(1 + x1^2 + x2^2)` and overlays the observations.

- palette:

  Color palette passed to
  [`grDevices::hcl.colors()`](https://rdrr.io/r/grDevices/palettes.html).

- point_size:

  Point size.

- main:

  Optional plot title.

- inference_df:

  Degrees-of-freedom rule used when plotting local p-values; passed to
  [`local_inference()`](https://byuzbasi.github.io/LHWL/reference/local_inference.md).

- p_adjust:

  Multiplicity adjustment used when plotting adjusted p-values or
  significance; passed to
  [`local_inference()`](https://byuzbasi.github.io/LHWL/reference/local_inference.md).

- alpha:

  Significance level used when `type = "significance"`.

- ...:

  Additional arguments passed to base plotting functions.

## Value

Invisibly returns `x`.

## Examples

``` r
set.seed(4)
H <- simulate_H2(60, r_max = 3)
Z <- H$Z
x1 <- rnorm(60)
y <- 1 + 2 * x1 + rnorm(60, sd = 0.3)
dat <- data.frame(y = y, x1 = x1)
fit <- lhwr(y ~ x1, dat, Z = Z, bw = 30, adaptive = TRUE, diagnostics = "full")
plot(fit, coefficient = "x1")
plot(fit, type = "p_value", coefficient = "x1")
plot(fit, type = "significance", coefficient = "x1")
plot(fit, type = "residual")
plot(fit, coefficient = "x1", geometry = "hyperboloid")
```

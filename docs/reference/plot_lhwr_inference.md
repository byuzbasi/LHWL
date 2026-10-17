# Plot LHWR Local Inference

Plot local p-values, adjusted p-values, or local significance for a
fitted LHWR model. The Poincare geometry gives a clean disk map;
projection and hyperboloid modes reuse
[`plot.lhwr()`](https://byuzbasi.github.io/LHWL/reference/plot.lhwr.md).

## Usage

``` r
plot_lhwr_inference(
  object,
  coefficient = NULL,
  type = c("p_value", "p_adjusted", "significance"),
  geometry = c("poincare", "projection", "hyperboloid"),
  alpha = 0.05,
  p_adjust = c("BH", "none", "holm", "bonferroni", "BY", "fdr"),
  inference_df = c("edf", "residual", "normal"),
  show_info = FALSE,
  ...
)
```

## Arguments

- object:

  An LHWR model object fitted with `diagnostics = "full"`.

- coefficient:

  Coefficient name.

- type:

  Inference surface: `"p_value"`, `"p_adjusted"`, or `"significance"`.

- geometry:

  Plot geometry: `"poincare"`, `"projection"`, or `"hyperboloid"`.

- alpha:

  Significance level.

- p_adjust:

  Multiplicity adjustment.

- inference_df:

  Degrees-of-freedom rule.

- show_info:

  Logical; draw the explanatory information box in the Poincare view.

- ...:

  Additional arguments passed to the plotting function.

## Value

Invisibly returns the map object for Poincare plots or `object`.

## Examples

``` r
H <- simulate_H2(45, r_max = 2)
x <- rnorm(45)
y <- 1 + x + rnorm(45, sd = 0.2)
fit <- lhwr(y ~ x, data.frame(y, x), Z = H$Z, bw = 18, adaptive = TRUE,
            diagnostics = "full")
plot_lhwr_inference(fit, coefficient = "x", type = "significance")
```

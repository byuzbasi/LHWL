# Plot LHWR Model Surfaces

User-friendly wrapper around
[`plot.lhwr()`](https://byuzbasi.github.io/LHWL/reference/plot.lhwr.md)
for coefficient, fitted-value, residual, and local R-squared maps.

## Usage

``` r
plot_lhwr_model(
  object,
  type = c("coefficient", "se", "t_value", "fitted", "residual", "local_r2"),
  coefficient = NULL,
  ...
)
```

## Arguments

- object:

  An LHWR model object.

- type:

  Model surface type: `"coefficient"`, `"se"`, `"t_value"`, `"fitted"`,
  `"residual"`, or `"local_r2"`.

- coefficient:

  Coefficient name for coefficient-level surfaces.

- ...:

  Additional arguments passed to
  [`plot.lhwr()`](https://byuzbasi.github.io/LHWL/reference/plot.lhwr.md).

## Value

Invisibly returns `object`.

## Examples

``` r
H <- simulate_H2(35, r_max = 2)
x <- rnorm(35)
y <- 1 + x + rnorm(35, sd = 0.2)
fit <- lhwr(y ~ x, data.frame(y, x), Z = H$Z, bw = 15, adaptive = TRUE)
plot_lhwr_model(fit, type = "coefficient", coefficient = "x")
```

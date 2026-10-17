# Plot LHWR Diagnostics

Plot common diagnostics from bandwidth and fitted LHWR objects.

## Usage

``` r
plot_lhwr_diagnostics(
  object,
  type = c("bandwidth", "residual", "local_r2", "condition_number",
    "effective_sample_size"),
  ...
)
```

## Arguments

- object:

  An LHWR model object or a bandwidth object returned by
  [`bw.lhwr()`](https://byuzbasi.github.io/LHWL/reference/bw.lhwr.md).

- type:

  Diagnostic type. For bandwidth objects use `"bandwidth"`. For model
  objects use `"residual"`, `"local_r2"`, `"condition_number"`, or
  `"effective_sample_size"`.

- ...:

  Additional plotting arguments.

## Value

Invisibly returns `object`.

## Examples

``` r
H <- simulate_H2(35, r_max = 2)
x <- rnorm(35)
y <- 1 + x + rnorm(35, sd = 0.2)
fit <- lhwr(y ~ x, data.frame(y, x), Z = H$Z, bw = 15, adaptive = TRUE)
plot_lhwr_diagnostics(fit, type = "residual")
```

# Plot Lorentz-Hyperbolic LISA Diagnostics (LH-LISA)

Plot Lorentz-Hyperbolic LISA Diagnostics (LH-LISA)

## Usage

``` r
# S3 method for class 'hgwr_local_moran'
plot(x, significant_only = FALSE, alpha = 0.05, ...)
```

## Arguments

- x:

  An object returned by
  [`local_moran.lhwr()`](https://byuzbasi.github.io/LHWL/reference/local_moran.lhwr.md).

- significant_only:

  Logical. If `TRUE`, nonsignificant observations are plotted in grey.

- alpha:

  Significance level used when `significant_only = TRUE`.

- ...:

  Additional arguments passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisibly returns `x`.

# HGWR Kernel Interface

Dispatch to one of the supported kernel functions and return the
associated weights.

## Usage

``` r
hgwr_kernel(d, bw, kernel = hgwr_kernel_names())
```

## Arguments

- d:

  Numeric vector or matrix of nonnegative distances.

- bw:

  Positive bandwidth value, or a vector of bandwidths compatible with
  the columns of `d`.

- kernel:

  Kernel name. Supported choices are `"bisquare"`, `"gaussian"`,
  `"exponential"`, `"tricube"`, and `"boxcar"`.

## Value

Numeric weights with the same shape as `d`.

## Details

Kernel shape and bandwidth have different roles. Gaussian and
exponential kernels have noncompact support; bisquare, tricube and
boxcar are zero at and beyond the bandwidth. The exact kNN fitting
engine supports only those three compact kernels. The Moran/LISA direct
kernel choices form a smaller set; an explicit diagnostic weight matrix
can be supplied separately. See
[`vignette("lhwr-diagnostics-and-tuning", package = "LHWL")`](https://byuzbasi.github.io/LHWL/articles/lhwr-diagnostics-and-tuning.md)
for all five formulas, overlaid weight curves and supported option
combinations.

## Examples

``` r
d <- c(0, 0.5, 1, 2)
hgwr_kernel(d, bw = 1, kernel = "bisquare")
hgwr_kernel(d, bw = 1, kernel = "gaussian")
```

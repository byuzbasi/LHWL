# Tricube HGWR Kernel

Compute compact tricube kernel weights. Weights are `(1 - (d / bw)^3)^3`
for `d < bw` and zero otherwise.

## Usage

``` r
hgwr_tricube(d, bw)
```

## Arguments

- d:

  Numeric vector or matrix of nonnegative distances.

- bw:

  Positive bandwidth value, or a vector of bandwidths compatible with
  the columns of `d`.

## Value

Numeric weights with the same shape as `d`.

## Examples

``` r
d <- c(0, 0.5, 1, 2)
hgwr_tricube(d, bw = 1)
```

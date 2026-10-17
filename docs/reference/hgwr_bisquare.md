# Bisquare HGWR Kernel

Compute compact bisquare kernel weights. Weights are
`(1 - (d / bw)^2)^2` for `d < bw` and zero otherwise.

## Usage

``` r
hgwr_bisquare(d, bw)
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
hgwr_bisquare(d, bw = 1)
```

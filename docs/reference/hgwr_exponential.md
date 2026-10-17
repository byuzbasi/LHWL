# Exponential HGWR Kernel

Compute exponential kernel weights `exp(-d / bw)`.

## Usage

``` r
hgwr_exponential(d, bw)
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
hgwr_exponential(d, bw = 1)
```

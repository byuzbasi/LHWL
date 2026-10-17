# Lorentz Norm

Compute the signed Lorentz norm `<z,z>_L`.

## Usage

``` r
lorentz_norm(z)
```

## Arguments

- z:

  A numeric vector of length 3 or an `n x 3` matrix.

## Value

A numeric vector.

## Examples

``` r
z <- c(cosh(1), sinh(1), 0)
lorentz_norm(z)
```

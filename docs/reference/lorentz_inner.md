# Lorentz Inner Product

Compute `<z,w>_L = -z0 w0 + z1 w1 + z2 w2` row-wise for points in the
Lorentz model. Vectors and matrices with three columns are accepted.

## Usage

``` r
lorentz_inner(z, w)
```

## Arguments

- z:

  A numeric vector of length 3 or an `n x 3` matrix.

- w:

  A numeric vector of length 3 or an `n x 3` matrix.

## Value

A numeric vector of Lorentz inner products.

## Examples

``` r
z <- c(cosh(1), sinh(1), 0)
lorentz_inner(z, z)
```

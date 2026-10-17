# Hyperbolic Distance on the Lorentz Hyperboloid

Compute `dH(z,w) = acosh(-<z,w>_L)`. The compiled implementation uses an
equivalent stable hyperbolic identity when the direct Lorentz
contraction is susceptible to cancellation at large radial coordinates.

## Usage

``` r
hyperbolic_distance(z, w, check = TRUE, tol = 1e-06)
```

## Arguments

- z:

  A numeric vector of length 3 or an `n x 3` matrix.

- w:

  A numeric vector of length 3 or an `n x 3` matrix.

- check:

  Logical; if `TRUE`, validate that inputs are on the upper sheet of the
  hyperboloid.

- tol:

  Numerical tolerance used in hyperboloid checks.

## Value

A numeric vector of distances.

## Examples

``` r
z <- c(1, 0, 0)
w <- c(cosh(1), sinh(1), 0)
hyperbolic_distance(z, w)
```

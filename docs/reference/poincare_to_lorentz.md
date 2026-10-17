# Convert Poincare Disk Coordinates to Lorentz Coordinates

Convert two-dimensional Poincare disk coordinates to the upper sheet of
the Lorentz hyperboloid. Input rows must satisfy `rowSums(P^2) < 1`.

## Usage

``` r
poincare_to_lorentz(P, tol = 1e-10)
```

## Arguments

- P:

  An `n x 2` numeric matrix of Poincare disk coordinates.

- tol:

  Numerical tolerance used when checking the unit disk boundary.

## Value

An `n x 3` Lorentz coordinate matrix with columns `z0`, `z1`, `z2`.

## Examples

``` r
P <- matrix(c(0, 0, 0.2, -0.1), ncol = 2, byrow = TRUE)
Z <- poincare_to_lorentz(P)
check_hyperboloid(Z)
```

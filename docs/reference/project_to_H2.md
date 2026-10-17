# Project Spatial Coordinates to the Lorentz Hyperboloid

Map Euclidean coordinates `(x, y)` to `(sqrt(1 + x^2 + y^2), x, y)`.

## Usage

``` r
project_to_H2(S)
```

## Arguments

- S:

  An `n x 2` numeric matrix.

## Value

An `n x 3` matrix on the upper hyperboloid.

## Examples

``` r
S <- matrix(c(0, 0, 1, 0), ncol = 2, byrow = TRUE)
Z <- project_to_H2(S)
check_hyperboloid(Z)
```

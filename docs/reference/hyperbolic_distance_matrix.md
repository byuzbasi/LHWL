# Hyperbolic Distance Matrix

Compute pairwise hyperbolic distances between rows of `Z` and optionally
rows of `W`.

## Usage

``` r
hyperbolic_distance_matrix(Z, W = NULL, check = TRUE, tol = 1e-06)
```

## Arguments

- Z:

  An `n x 3` matrix of hyperboloid points.

- W:

  Optional `m x 3` matrix of hyperboloid points. If omitted, a symmetric
  `n x n` matrix is returned.

- check:

  Logical; validate hyperboloid constraints.

- tol:

  Numerical tolerance.

## Value

A numeric distance matrix.

## Examples

``` r
H <- simulate_H2(5, r_max = 2)
hyperbolic_distance_matrix(H$Z)
```

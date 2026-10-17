# Lorentz Hyperbolic k-Nearest Neighbors

Compute exact k-nearest neighbors under the intrinsic Lorentz hyperbolic
distance without returning the full pairwise distance matrix.

## Usage

``` r
lorentz_knn(Z, k, include_self = FALSE, check = TRUE, tol = 1e-06)
```

## Arguments

- Z:

  An `n x 3` matrix of Lorentz hyperboloid points.

- k:

  Number of nearest neighbors to return.

- include_self:

  If `TRUE`, each row includes itself as the first neighbor with
  distance zero. If `FALSE`, self-neighbors are excluded.

- check:

  Logical; validate hyperboloid constraints.

- tol:

  Numerical tolerance.

## Value

A list with `index`, an `n x k` integer matrix of one-based neighbor
indices, `distance`, an `n x k` matrix of corresponding hyperbolic
distances, `k`, and `include_self`.

## Examples

``` r
H <- simulate_H2(20, r_max = 3)
knn <- lorentz_knn(H$Z, k = 5)
str(knn)
```

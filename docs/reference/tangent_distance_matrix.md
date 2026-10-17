# Tangent-Space Distance Matrix

Compute a distance matrix after mapping Lorentz hyperboloid points to
`T_mu H^2` with the logarithmic map. Distances are computed using the
Lorentz inner product restricted to the tangent space, which is positive
definite. Equivalently, this is Euclidean distance in a
Lorentz-orthonormal tangent basis.

## Usage

``` r
tangent_distance_matrix(Z, mu = NULL, tol = 1e-10)
```

## Arguments

- Z:

  `n x 3` Lorentz hyperboloid coordinates.

- mu:

  Reference point. If `NULL`, the normalized extrinsic mean is used.

- tol:

  Intrinsic-distance threshold used by
  [`lorentz_log_map()`](https://byuzbasi.github.io/LHWL/reference/lorentz_log_map.md);
  points within this distance of `mu` are mapped to zero.

## Value

A symmetric `n x n` tangent-space distance matrix.

## Details

With \\v_i=\log\_\mu(z_i)\\ and \\v_j=\log\_\mu(z_j)\\, the distance is
\$\$d\_{T,\mu}(z_i,z_j)=\sqrt{\langle v_i-v_j,v_i-v_j\rangle_L}.\$\$
Here \\i,j\\ index input rows and
[`lorentz_inner()`](https://byuzbasi.github.io/LHWL/reference/lorentz_inner.md)
defines the inner product. Its restriction to the common tangent plane
is positive definite. These distances define the Tangent-E weighting
approximation: radial distances from `mu` are preserved in exact
arithmetic, but general pairwise hyperbolic distances are not. The
implementation uses coordinate differences and gives identical mapped
points exactly zero distance. The `tol` argument applies to each point's
distance from `mu` in the log map, not to arbitrary pairwise distances
between distinct tangent coordinates.

## Examples

``` r
H <- simulate_H2(5, r_max = 2)
tangent_distance_matrix(H$Z)
```

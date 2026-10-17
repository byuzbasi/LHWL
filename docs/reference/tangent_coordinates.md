# Tangent-Space Coordinates

Map Lorentz hyperboloid points to two-dimensional coordinates in a
Lorentz-orthonormal basis of `T_mu H^2`.

## Usage

``` r
tangent_coordinates(Z, mu = NULL)
```

## Arguments

- Z:

  `n x 3` Lorentz hyperboloid coordinates.

- mu:

  Reference point. If `NULL`, the normalized extrinsic mean is used.

## Value

An `n x 2` matrix of tangent coordinates.

## Details

The columns express
[`lorentz_log_map()`](https://byuzbasi.github.io/LHWL/reference/lorentz_log_map.md)
in an orthonormal tangent basis. Their ordinary Euclidean norm equals
the Lorentz norm of the tangent vector. This routine uses the log map's
default zero-distance threshold `1e-10`. To compare training and new
points, pass the same `mu` in both calls; separately estimated centers
define different coordinate systems.

## Examples

``` r
H <- simulate_H2(5, r_max = 2)
tangent_coordinates(H$Z)
```

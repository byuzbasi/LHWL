# Extrinsic Mean Projected to the Lorentz Hyperboloid

Compute the Euclidean coordinate mean of Lorentz points and project it
back to the upper sheet of the hyperboloid. This is used as the default
reference point for tangent-space approximations.

## Usage

``` r
extrinsic_lorentz_mean(Z)
```

## Arguments

- Z:

  `n x 3` Lorentz hyperboloid coordinates.

## Value

A length-3 Lorentz hyperboloid point.

## Details

For valid upper-sheet points, let \\m = n^{-1}\sum_i z_i\\ be the
ambient coordinate mean of the \\n\\ input rows and let
\\\langle\cdot,\cdot\rangle_L\\ denote
[`lorentz_inner()`](https://byuzbasi.github.io/LHWL/reference/lorentz_inner.md).
The reference is \$\$\mu = m / \sqrt{-\langle m,m\rangle_L}.\$\$ This
normalization puts the mean on the unit hyperboloid; it is not the
intrinsic Frechet mean. The implementation evaluates the normalization
stably and reconstructs the positive time coordinate. A zero spatial
mean gives \\\mu=(1,0,0)\\. Estimate the reference from training
locations and reuse it when mapping new locations.

## Examples

``` r
H <- simulate_H2(10, r_max = 2)
mu <- extrinsic_lorentz_mean(H$Z)
check_hyperboloid(mu)
```

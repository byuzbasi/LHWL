# Lorentz Logarithmic Map

Map points on the Lorentz hyperboloid to the tangent space at `mu`.

## Usage

``` r
lorentz_log_map(Z, mu = NULL, tol = 1e-10)
```

## Arguments

- Z:

  `n x 3` Lorentz hyperboloid coordinates.

- mu:

  Reference point on the Lorentz hyperboloid. If `NULL`, the normalized
  extrinsic mean is used.

- tol:

  Intrinsic-distance threshold for setting a log-map row to zero; the
  default is `1e-10` in unit-curvature distance units.

## Value

An `n x 3` matrix of tangent vectors in `T_mu H^2`.

## Details

For unit-hyperboloid points \\z\\ and \\\mu\\, define
\\a=-\langle\mu,z\rangle_L\\ and the intrinsic distance
\\d=\operatorname{acosh}(a)\\ using
[`lorentz_inner()`](https://byuzbasi.github.io/LHWL/reference/lorentz_inner.md)
and
[`hyperbolic_distance()`](https://byuzbasi.github.io/LHWL/reference/hyperbolic_distance.md).
The logarithmic map is
\$\$\log\_\mu(z)=\frac{d}{\sinh(d)}(z-a\mu),\qquad d\>0.\$\$ The result
is a tangent vector orthogonal to \\\mu\\ under the Lorentz inner
product. At \\z=\mu\\ the continuous extension is the zero vector.
Numerically, `tol` is compared with the distance \\d\\, not with
\\a-1\\: rows with \\d\leq\mathrm{tol}\\ are set to zero. Other rows use
stable, algebraically equivalent coordinate formulas.

## Examples

``` r
H <- simulate_H2(5, r_max = 2)
V <- lorentz_log_map(H$Z)
dim(V)
```

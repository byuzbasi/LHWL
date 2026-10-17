# Prepare Coordinates for LHWR

Convert common two-dimensional coordinate representations to the `n x 3`
Lorentz hyperboloid matrix used by
[`bw.lhwr()`](https://byuzbasi.github.io/LHWL/reference/bw.lhwr.md) and
[`lhwr()`](https://byuzbasi.github.io/LHWL/reference/lhwr.md).

## Usage

``` r
prepare_lorentz_coordinates(
  embedding,
  model = c("lorentz", "poincare", "spatial"),
  tol = 1e-06,
  stop_on_fail = TRUE
)
```

## Arguments

- embedding:

  Numeric coordinate matrix.

- model:

  Coordinate model of `embedding`.

- tol:

  Numerical tolerance used in validation.

- stop_on_fail:

  If `TRUE`, throw an error when the resulting coordinates are not valid
  Lorentz hyperboloid points. If `FALSE`, return the matrix and attach
  the row-wise validation result as attribute `valid_hyperboloid`.

## Value

An `n x 3` Lorentz coordinate matrix with attributes `input_model` and
`valid_hyperboloid`.

## Details

This is a data-preparation helper. It does not estimate an embedding.
Use `model = "lorentz"` when the input is already on the Lorentz
hyperboloid, `model = "poincare"` for coordinates inside the Poincare
disk, and `model = "spatial"` for ordinary two-dimensional coordinates
that should be embedded as `(sqrt(1 + x^2 + y^2), x, y)`.

## Examples

``` r
S <- matrix(c(0, 0, 1, 0), ncol = 2, byrow = TRUE)
Z_spatial <- prepare_lorentz_coordinates(S, model = "spatial")

P <- matrix(c(0, 0, 0.2, -0.1), ncol = 2, byrow = TRUE)
Z_poincare <- prepare_lorentz_coordinates(P, model = "poincare")

Z_checked <- prepare_lorentz_coordinates(Z_poincare, model = "lorentz")
all(check_hyperboloid(Z_checked))
```

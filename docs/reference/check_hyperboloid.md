# Check Hyperboloid Constraints

Verify `<z,z>_L = -1` and `z0 > 0` for all rows.

## Usage

``` r
check_hyperboloid(Z, tol = 1e-06, stop_on_fail = FALSE)
```

## Arguments

- Z:

  An `n x 3` matrix of Lorentz coordinates.

- tol:

  Numerical tolerance.

- stop_on_fail:

  If `TRUE`, throw an error when checks fail.

## Value

A logical vector indicating valid rows.

## Examples

``` r
H <- simulate_H2(5, r_max = 2)
check_hyperboloid(H$Z)
```

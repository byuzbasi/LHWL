# Safe Condition Number

Compute a condition number from a symmetric cross-product matrix.

## Usage

``` r
safe_condition_number(G, tol = 1e-10)
```

## Arguments

- G:

  Symmetric positive semidefinite matrix.

- tol:

  Eigenvalue tolerance.

## Value

Numeric condition number, or `Inf` if the matrix is singular.

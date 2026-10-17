# Weighted R-squared

Compute weighted R-squared for a weighted least-squares auxiliary
regression.

## Usage

``` r
weighted_r2(y, X, w)
```

## Arguments

- y:

  Numeric response vector.

- X:

  Numeric design matrix, usually including an intercept.

- w:

  Nonnegative weight vector.

## Value

Weighted R-squared, or `NA_real_` if it cannot be computed.

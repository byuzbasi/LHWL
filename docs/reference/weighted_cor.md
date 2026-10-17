# Weighted Correlation Matrix

Compute a weighted correlation matrix for numeric columns.

## Usage

``` r
weighted_cor(X, w, center = TRUE, scale = TRUE)
```

## Arguments

- X:

  Numeric matrix.

- w:

  Nonnegative weight vector.

- center:

  Logical. If `TRUE`, center columns by weighted means.

- scale:

  Logical. If `TRUE`, scale columns by weighted standard deviations.

## Value

Weighted correlation matrix.

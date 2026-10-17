# LHWR Cross-Validation Score

LHWR Cross-Validation Score

## Usage

``` r
CV.lhwr(object)

CV.hgwr(object)
```

## Arguments

- object:

  An LHWR bandwidth or model object.

## Value

A single numeric leave-one-out cross-validation score if it was computed
for `object`; otherwise `NA_real_`.

## Details

Gaussian CV is a sum of squared leave-one-out errors; binomial tuning
uses Bernoulli log loss. The representation and candidate radii are held
fixed for internal bandwidth scoring. This extractor does not run
cross-validation and a stored score is not an independent outer
evaluation after tuning.

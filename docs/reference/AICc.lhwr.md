# LHWR AICc

LHWR AICc

## Usage

``` r
AICc.lhwr(object)

AICc.hgwr(object)
```

## Arguments

- object:

  An LHWR model object returned by
  [`lhwr()`](https://byuzbasi.github.io/LHWL/reference/lhwr.md) or a
  low-level fit returned by
  [`hgwr.fit()`](https://byuzbasi.github.io/LHWL/reference/hgwr.fit.md).

## Value

A single numeric value containing the corrected Akaike information
criterion. Smaller values indicate a better fitted model under this
criterion.

## Details

This extractor does not refit or select a model. AICc requires
hat-matrix diagnostics; a nonpositive n - tr(S) - 2 denominator produces
infinity. Gaussian AIC and BIC are also stored in
`object$LHWR.diagnostic`, but
[`bw.lhwr()`](https://byuzbasi.github.io/LHWL/reference/bw.lhwr.md)
directly selects only by AICc or CV. See
[`vignette("lhwr-diagnostics-and-tuning", package = "LHWL")`](https://byuzbasi.github.io/LHWL/articles/lhwr-diagnostics-and-tuning.md)
for formulas.

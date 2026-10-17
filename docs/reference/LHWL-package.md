# LHWL Package

Lorentz Hyperbolic Weighted Learning on the Lorentz hyperboloid.

## Details

The LHWR module provides unpenalized Gaussian and binomial local
regression. The LHWPR module provides Gaussian and binomial ridge,
lasso, elastic net, MCP, and SCAD paths with joint
bandwidth-regularization tuning. Geometry, bandwidth selection, local
diagnostics, and plotting utilities support both modules where
statistically appropriate.

LHWR is the original unpenalized regression foundation. LHWPR-Lasso is a
separate penalized extension, not a replacement name for LHWR. Software
support for other penalties or response families does not imply that the
original Gaussian LHWR theory supplies guarantees for those extensions.

Use
[`make_lhwr_Z_features()`](https://byuzbasi.github.io/LHWL/reference/make_lhwr_Z.md)
or
[`make_lhwr_Z_network()`](https://byuzbasi.github.io/LHWL/reference/make_lhwr_Z.md)
to construct a representation,
[`bw.lhwr()`](https://byuzbasi.github.io/LHWL/reference/bw.lhwr.md) and
[`lhwr()`](https://byuzbasi.github.io/LHWL/reference/lhwr.md) for
unpenalized regression, and
[`tune.lhwpr()`](https://byuzbasi.github.io/LHWL/reference/tune.lhwpr.md)
for joint penalized tuning. For a complete introduction use
[`vignette("lorentz-hgwr-introduction", package = "LHWL")`](https://byuzbasi.github.io/LHWL/articles/lorentz-hgwr-introduction.md);
for Lasso use
[`vignette("lhwpr-lasso", package = "LHWL")`](https://byuzbasi.github.io/LHWL/articles/lhwpr-lasso.md).

The supplied representation defines locality; it is not automatically a
predictor or a validated latent geometry. Validation of a learned
representation requires rebuilding it within the relevant training
splits. Local Gaussian diagnostics condition on the fitted
representation and smoothing specification. Penalized activation is not
a significance test.

## See also

Useful links:

- <https://github.com/byuzbasi/LHWL>

- Report bugs at <https://github.com/byuzbasi/LHWL/issues>

## Author

**Maintainer**: Bahadır Yüzbaşı <b.yzb@hotmail.com>

Authors:

- Zühal Küçükarslan Yüzbaşı

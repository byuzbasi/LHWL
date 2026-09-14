#' LHWL Package
#'
#' Lorentz Hyperbolic Weighted Learning on the Lorentz hyperboloid.
#'
#' The LHWR module provides unpenalized Gaussian and binomial local regression.
#' The LHWPR module provides Gaussian and binomial ridge, lasso, elastic net,
#' MCP, and SCAD paths with joint bandwidth-regularization tuning. Geometry,
#' bandwidth selection, local diagnostics, and plotting utilities support both
#' modules where statistically appropriate.
#'
#' LHWR is the original unpenalized regression foundation. LHWPR-Lasso is a
#' separate penalized extension, not a replacement name for LHWR. Software
#' support for other penalties or response families does not imply that the
#' original Gaussian LHWR theory supplies guarantees for those extensions.
#'
#' Use [make_lhwr_Z_features()] or [make_lhwr_Z_network()] to construct a
#' representation, [bw.lhwr()] and [lhwr()] for unpenalized regression, and
#' [tune.lhwpr()] for joint penalized tuning. For a complete introduction use
#' `vignette("lorentz-hgwr-introduction", package = "LHWL")`; for Lasso use
#' `vignette("lhwpr-lasso", package = "LHWL")`.
#'
#' The supplied representation defines locality; it is not automatically a
#' predictor or a validated latent geometry. Validation of a learned
#' representation requires rebuilding it within the relevant training splits.
#' Local Gaussian diagnostics condition on the fitted representation and
#' smoothing specification. Penalized activation is not a significance test.
#'
#' @keywords internal
#' @importFrom Rcpp evalCpp
#' @importFrom MASS ginv
#' @importFrom graphics matplot par plot
#' @importFrom grDevices hcl.colors
#' @importFrom stats model.frame model.matrix model.response terms
#' @importFrom stats printCoefmat
#' @useDynLib LHWL, .registration = TRUE
"_PACKAGE"

#' LHWR Package
#'
#' Lorentz Hyperbolic Weighted Regression on the Lorentz hyperboloid.
#'
#' @keywords internal
#' @importFrom Rcpp evalCpp
#' @importFrom MASS ginv
#' @importFrom graphics matplot par plot
#' @importFrom grDevices hcl.colors
#' @importFrom stats model.frame model.matrix model.response terms
#' @importFrom stats printCoefmat
#' @useDynLib LHWR, .registration = TRUE
"_PACKAGE"

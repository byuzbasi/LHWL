#' Gaussian HGWR Kernel
#'
#' Compute Gaussian kernel weights `exp(-0.5 * (d / bw)^2)`.
#'
#' @param d Numeric vector or matrix of nonnegative distances.
#' @param bw Positive bandwidth value, or a vector of bandwidths compatible
#'   with the columns of `d`.
#' @return Numeric weights with the same shape as `d`.
#' @examples
#' d <- c(0, 0.5, 1, 2)
#' hgwr_gaussian(d, bw = 1)
#' @export
hgwr_gaussian <- function(d, bw) {
  check_bandwidth(bw)
  r <- scaled_distance(d, bw)
  exp(-0.5 * r^2)
}

#' Exponential HGWR Kernel
#'
#' Compute exponential kernel weights `exp(-d / bw)`.
#'
#' @param d Numeric vector or matrix of nonnegative distances.
#' @param bw Positive bandwidth value, or a vector of bandwidths compatible
#'   with the columns of `d`.
#' @return Numeric weights with the same shape as `d`.
#' @examples
#' d <- c(0, 0.5, 1, 2)
#' hgwr_exponential(d, bw = 1)
#' @export
hgwr_exponential <- function(d, bw) {
  check_bandwidth(bw)
  r <- scaled_distance(d, bw)
  exp(-r)
}

#' Bisquare HGWR Kernel
#'
#' Compute compact bisquare kernel weights. Weights are
#' `(1 - (d / bw)^2)^2` for `d < bw` and zero otherwise.
#'
#' @param d Numeric vector or matrix of nonnegative distances.
#' @param bw Positive bandwidth value, or a vector of bandwidths compatible
#'   with the columns of `d`.
#' @return Numeric weights with the same shape as `d`.
#' @examples
#' d <- c(0, 0.5, 1, 2)
#' hgwr_bisquare(d, bw = 1)
#' @export
hgwr_bisquare <- function(d, bw) {
  check_bandwidth(bw)
  r <- scaled_distance(d, bw)
  w <- (1 - r^2)^2
  w[r >= 1] <- 0
  w
}

#' Tricube HGWR Kernel
#'
#' Compute compact tricube kernel weights. Weights are
#' `(1 - (d / bw)^3)^3` for `d < bw` and zero otherwise.
#'
#' @param d Numeric vector or matrix of nonnegative distances.
#' @param bw Positive bandwidth value, or a vector of bandwidths compatible
#'   with the columns of `d`.
#' @return Numeric weights with the same shape as `d`.
#' @examples
#' d <- c(0, 0.5, 1, 2)
#' hgwr_tricube(d, bw = 1)
#' @export
hgwr_tricube <- function(d, bw) {
  check_bandwidth(bw)
  r <- scaled_distance(d, bw)
  w <- (1 - r^3)^3
  w[r >= 1] <- 0
  w
}

#' Boxcar HGWR Kernel
#'
#' Compute compact boxcar kernel weights. Weights are one for `d < bw` and zero
#' otherwise.
#'
#' @param d Numeric vector or matrix of nonnegative distances.
#' @param bw Positive bandwidth value, or a vector of bandwidths compatible
#'   with the columns of `d`.
#' @return Numeric weights with the same shape as `d`.
#' @examples
#' d <- c(0, 0.5, 1, 2)
#' hgwr_boxcar(d, bw = 1)
#' @export
hgwr_boxcar <- function(d, bw) {
  check_bandwidth(bw)
  r <- scaled_distance(d, bw)
  ifelse(r < 1, 1, 0)
}

#' HGWR Kernel Interface
#'
#' Dispatch to one of the supported kernel functions and return the associated
#' weights.
#'
#' @param d Numeric vector or matrix of nonnegative distances.
#' @param bw Positive bandwidth value, or a vector of bandwidths compatible
#'   with the columns of `d`.
#' @param kernel Kernel name. Supported choices are `"bisquare"`,
#'   `"gaussian"`, `"exponential"`, `"tricube"`, and `"boxcar"`.
#' @return Numeric weights with the same shape as `d`.
#' @examples
#' d <- c(0, 0.5, 1, 2)
#' hgwr_kernel(d, bw = 1, kernel = "bisquare")
#' hgwr_kernel(d, bw = 1, kernel = "gaussian")
#' @export
hgwr_kernel <- function(d, bw, kernel = hgwr_kernel_names()) {
  kernel <- match.arg(kernel)
  switch(
    kernel,
    gaussian = hgwr_gaussian(d, bw),
    exponential = hgwr_exponential(d, bw),
    bisquare = hgwr_bisquare(d, bw),
    tricube = hgwr_tricube(d, bw),
    boxcar = hgwr_boxcar(d, bw)
  )
}

hgwr_kernel_names <- function() {
  c("bisquare", "gaussian", "exponential", "tricube", "boxcar")
}

check_bandwidth <- function(bw) {
  if (!is.numeric(bw) || any(!is.finite(bw)) || any(bw <= 0)) {
    stop("Bandwidth values must be positive finite numbers.", call. = FALSE)
  }
  invisible(TRUE)
}

scaled_distance <- function(d, bw) {
  if (is.matrix(d) && length(bw) == ncol(d)) {
    sweep(d, 2L, bw, "/")
  } else {
    d / bw
  }
}

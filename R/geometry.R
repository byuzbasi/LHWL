#' Lorentz Inner Product
#'
#' Compute `<z,w>_L = -z0 w0 + z1 w1 + z2 w2` row-wise for points in the
#' Lorentz model. Vectors and matrices with three columns are accepted.
#'
#' @param z A numeric vector of length 3 or an `n x 3` matrix.
#' @param w A numeric vector of length 3 or an `n x 3` matrix.
#' @return A numeric vector of Lorentz inner products.
#' @examples
#' z <- c(cosh(1), sinh(1), 0)
#' lorentz_inner(z, z)
#' @export
lorentz_inner <- function(z, w) {
  z <- as_lorentz_matrix(z, "z")
  w <- as_lorentz_matrix(w, "w")
  
  .lorentz_inner_cpp(z, w)
}

#' Lorentz Norm
#'
#' Compute the signed Lorentz norm `<z,z>_L`.
#'
#' @param z A numeric vector of length 3 or an `n x 3` matrix.
#' @return A numeric vector.
#' @examples
#' z <- c(cosh(1), sinh(1), 0)
#' lorentz_norm(z)
#' @export
lorentz_norm <- function(z) {
  lorentz_inner(z, z)
}

#' Hyperbolic Distance on the Lorentz Hyperboloid
#'
#' Compute `dH(z,w) = acosh(-<z,w>_L)`. The compiled implementation uses an
#' equivalent stable hyperbolic identity when the direct Lorentz contraction
#' is susceptible to cancellation at large radial coordinates.
#'
#' @param z A numeric vector of length 3 or an `n x 3` matrix.
#' @param w A numeric vector of length 3 or an `n x 3` matrix.
#' @param check Logical; if `TRUE`, validate that inputs are on the
#'   upper sheet of the hyperboloid.
#' @param tol Numerical tolerance used in hyperboloid checks.
#' @return A numeric vector of distances.
#' @examples
#' z <- c(1, 0, 0)
#' w <- c(cosh(1), sinh(1), 0)
#' hyperbolic_distance(z, w)
#' @export
hyperbolic_distance <- function(z, w, check = TRUE, tol = 1e-6) {
  z <- as_lorentz_matrix(z, "z")
  w <- as_lorentz_matrix(w, "w")
  
  if (check) {
    check_hyperboloid(z, tol = tol, stop_on_fail = TRUE)
    check_hyperboloid(w, tol = tol, stop_on_fail = TRUE)
  }
  
  .hyperbolic_distance_cpp(z, w, tol = tol)
}

#' Hyperbolic Distance Matrix
#'
#' Compute pairwise hyperbolic distances between rows of `Z` and optionally
#' rows of `W`.
#'
#' @param Z An `n x 3` matrix of hyperboloid points.
#' @param W Optional `m x 3` matrix of hyperboloid points. If omitted, a
#'   symmetric `n x n` matrix is returned.
#' @param check Logical; validate hyperboloid constraints.
#' @param tol Numerical tolerance.
#' @return A numeric distance matrix.
#' @examples
#' H <- simulate_H2(5, r_max = 2)
#' hyperbolic_distance_matrix(H$Z)
#' @export
hyperbolic_distance_matrix <- function(Z, W = NULL, check = TRUE, tol = 1e-6) {
  Z <- as_lorentz_matrix(Z, "Z")
  symmetric <- is.null(W)
  W <- if (symmetric) Z else as_lorentz_matrix(W, "W")
  
  if (check) {
    check_hyperboloid(Z, tol = tol, stop_on_fail = TRUE)
    if (!symmetric) check_hyperboloid(W, tol = tol, stop_on_fail = TRUE)
  }
  
  .hyperbolic_distance_matrix_cpp(Z, W, symmetric = symmetric, tol = tol)
}

#' Lorentz Hyperbolic k-Nearest Neighbors
#'
#' Compute exact k-nearest neighbors under the intrinsic Lorentz hyperbolic
#' distance without returning the full pairwise distance matrix.
#'
#' @param Z An `n x 3` matrix of Lorentz hyperboloid points.
#' @param k Number of nearest neighbors to return.
#' @param include_self If `TRUE`, each row includes itself as the first
#'   neighbor with distance zero. If `FALSE`, self-neighbors are excluded.
#' @param check Logical; validate hyperboloid constraints.
#' @param tol Numerical tolerance.
#' @return A list with `index`, an `n x k` integer matrix of one-based neighbor
#'   indices, `distance`, an `n x k` matrix of corresponding hyperbolic
#'   distances, `k`, and `include_self`.
#' @examples
#' H <- simulate_H2(20, r_max = 3)
#' knn <- lorentz_knn(H$Z, k = 5)
#' str(knn)
#' @export
lorentz_knn <- function(Z, k, include_self = FALSE, check = TRUE, tol = 1e-6) {
  Z <- as_lorentz_matrix(Z, "Z")
  if (check) {
    check_hyperboloid(Z, tol = tol, stop_on_fail = TRUE)
  }
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k)) {
    stop("`k` must be a positive integer scalar.", call. = FALSE)
  }
  k <- as.integer(round(k))
  n <- nrow(Z)
  max_k <- if (isTRUE(include_self)) n else n - 1L
  if (k < 1L || k > max_k) {
    stop("`k` must be between 1 and the number of available neighbors.", call. = FALSE)
  }
  .lorentz_knn_cpp(Z, k = k, include_self = isTRUE(include_self), tol = tol)
}

as_lorentz_matrix <- function(x, name) {
  if (!is.numeric(x)) {
    stop("`", name, "` must be numeric.", call. = FALSE)
  }
  
  if (is.null(dim(x))) {
    if (length(x) != 3L) {
      stop("`", name, "` must have length 3.", call. = FALSE)
    }
    x <- matrix(x, nrow = 1L)
  } else {
    x <- as.matrix(x)
    if (ncol(x) != 3L) {
      stop("`", name, "` must have exactly 3 columns.", call. = FALSE)
    }
  }
  
  storage.mode(x) <- "double"
  x
}

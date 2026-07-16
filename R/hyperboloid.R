#' Simulate Points on the Two-Dimensional Hyperboloid
#'
#' @param n Number of points.
#' @param r_max Maximum hyperbolic radius.
#' @param uniform_area If `TRUE`, sample uniformly with respect to hyperbolic
#'   area inside the radius `r_max`; otherwise sample radius uniformly.
#' @return A list with Lorentz coordinates `Z`, spatial projection `S`, radius
#'   `r`, and angle `theta`.
#' @examples
#' set.seed(1)
#' H <- simulate_H2(10, r_max = 2)
#' head(H$Z)
#' check_hyperboloid(H$Z)
#' @export
simulate_H2 <- function(n, r_max = 8, uniform_area = TRUE) {
  if (!is.numeric(n) || length(n) != 1L || n < 1) {
    stop("`n` must be a positive scalar.", call. = FALSE)
  }
  if (!is.numeric(r_max) || length(r_max) != 1L || r_max <= 0) {
    stop("`r_max` must be a positive scalar.", call. = FALSE)
  }
  
  if (uniform_area) {
    u <- stats::runif(n)
    r <- acosh(1 + u * (cosh(r_max) - 1))
  } else {
    r <- stats::runif(n, 0, r_max)
  }
  theta <- stats::runif(n, 0, 2 * pi)
  
  Z <- cbind(
    cosh(r),
    sinh(r) * cos(theta),
    sinh(r) * sin(theta)
  )
  colnames(Z) <- c("z0", "z1", "z2")
  
  list(Z = Z, S = Z[, 2:3, drop = FALSE], r = r, theta = theta)
}

#' Project Spatial Coordinates to the Lorentz Hyperboloid
#'
#' Map Euclidean coordinates `(x, y)` to `(sqrt(1 + x^2 + y^2), x, y)`.
#'
#' @param S An `n x 2` numeric matrix.
#' @return An `n x 3` matrix on the upper hyperboloid.
#' @examples
#' S <- matrix(c(0, 0, 1, 0), ncol = 2, byrow = TRUE)
#' Z <- project_to_H2(S)
#' check_hyperboloid(Z)
#' @export
project_to_H2 <- function(S) {
  if (!is.numeric(S)) {
    stop("`S` must be numeric.", call. = FALSE)
  }
  S <- as.matrix(S)
  if (ncol(S) != 2L) {
    stop("`S` must have exactly 2 columns.", call. = FALSE)
  }
  
  Z <- cbind(sqrt(1 + rowSums(S^2)), S)
  colnames(Z) <- c("z0", "z1", "z2")
  Z
}

#' Convert Poincare Disk Coordinates to Lorentz Coordinates
#'
#' Convert two-dimensional Poincare disk coordinates to the upper sheet of the
#' Lorentz hyperboloid. Input rows must satisfy `rowSums(P^2) < 1`.
#'
#' @param P An `n x 2` numeric matrix of Poincare disk coordinates.
#' @param tol Numerical tolerance used when checking the unit disk boundary.
#' @return An `n x 3` Lorentz coordinate matrix with columns `z0`, `z1`, `z2`.
#' @examples
#' P <- matrix(c(0, 0, 0.2, -0.1), ncol = 2, byrow = TRUE)
#' Z <- poincare_to_lorentz(P)
#' check_hyperboloid(Z)
#' @export
poincare_to_lorentz <- function(P, tol = 1e-10) {
  if (!is.numeric(P)) {
    stop("`P` must be numeric.", call. = FALSE)
  }
  P <- as.matrix(P)
  if (ncol(P) != 2L) {
    stop("`P` must have exactly 2 columns.", call. = FALSE)
  }
  storage.mode(P) <- "double"
  
  norm2 <- rowSums(P^2)
  ok <- is.finite(norm2) & norm2 < (1 - tol)
  if (!all(ok)) {
    bad <- which(!ok)
    stop(
      "Rows are not valid Poincare disk coordinates: ",
      paste(utils::head(bad, 5L), collapse = ", "),
      if (length(bad) > 5L) ", ...",
      call. = FALSE
    )
  }
  
  den <- 1 - norm2
  Z <- cbind(
    (1 + norm2) / den,
    2 * P[, 1] / den,
    2 * P[, 2] / den
  )
  colnames(Z) <- c("z0", "z1", "z2")
  Z
}

#' Prepare Coordinates for LHWR
#'
#' Convert common two-dimensional coordinate representations to the `n x 3`
#' Lorentz hyperboloid matrix used by `bw.lhwr()` and `lhwr()`.
#'
#' This is a data-preparation helper. It does not estimate an embedding. Use
#' `model = "lorentz"` when the input is already on the Lorentz hyperboloid,
#' `model = "poincare"` for coordinates inside the Poincare disk, and
#' `model = "spatial"` for ordinary two-dimensional coordinates that should be
#' embedded as `(sqrt(1 + x^2 + y^2), x, y)`.
#'
#' @param embedding Numeric coordinate matrix.
#' @param model Coordinate model of `embedding`.
#' @param tol Numerical tolerance used in validation.
#' @param stop_on_fail If `TRUE`, throw an error when the resulting coordinates
#'   are not valid Lorentz hyperboloid points. If `FALSE`, return the matrix and
#'   attach the row-wise validation result as attribute `valid_hyperboloid`.
#' @return An `n x 3` Lorentz coordinate matrix with attributes `input_model`
#'   and `valid_hyperboloid`.
#' @examples
#' S <- matrix(c(0, 0, 1, 0), ncol = 2, byrow = TRUE)
#' Z_spatial <- prepare_lorentz_coordinates(S, model = "spatial")
#'
#' P <- matrix(c(0, 0, 0.2, -0.1), ncol = 2, byrow = TRUE)
#' Z_poincare <- prepare_lorentz_coordinates(P, model = "poincare")
#'
#' Z_checked <- prepare_lorentz_coordinates(Z_poincare, model = "lorentz")
#' all(check_hyperboloid(Z_checked))
#' @export
prepare_lorentz_coordinates <- function(embedding,
                                        model = c("lorentz", "poincare", "spatial"),
                                        tol = 1e-6,
                                        stop_on_fail = TRUE) {
  model <- match.arg(model)
  
  Z <- switch(
    model,
    lorentz = as_lorentz_matrix(embedding, "embedding"),
    poincare = poincare_to_lorentz(embedding, tol = min(tol, 1e-10)),
    spatial = project_to_H2(embedding)
  )
  
  valid <- check_hyperboloid(Z, tol = tol, stop_on_fail = FALSE)
  if (stop_on_fail && !all(valid)) {
    bad <- which(!valid)
    stop(
      "Prepared coordinates are not on the upper Lorentz hyperboloid: ",
      paste(utils::head(bad, 5L), collapse = ", "),
      if (length(bad) > 5L) ", ...",
      call. = FALSE
    )
  }
  
  attr(Z, "input_model") <- model
  attr(Z, "valid_hyperboloid") <- valid
  Z
}

#' Check Hyperboloid Constraints
#'
#' Verify `<z,z>_L = -1` and `z0 > 0` for all rows.
#'
#' @param Z An `n x 3` matrix of Lorentz coordinates.
#' @param tol Numerical tolerance.
#' @param stop_on_fail If `TRUE`, throw an error when checks fail.
#' @return A logical vector indicating valid rows.
#' @examples
#' H <- simulate_H2(5, r_max = 2)
#' check_hyperboloid(H$Z)
#' @export
check_hyperboloid <- function(Z, tol = 1e-6, stop_on_fail = FALSE) {
  Z <- as_lorentz_matrix(Z, "Z")
  residual <- abs(lorentz_norm(Z) + 1)
  scale <- pmax(1, abs(Z[, 1])^2 + rowSums(Z[, -1, drop = FALSE]^2))
  ok <- is.finite(rowSums(Z)) &
    (residual <= tol | residual / scale <= tol) &
    Z[, 1] > 0
  
  if (stop_on_fail && !all(ok)) {
    bad <- which(!ok)
    stop(
      "Rows are not on the upper Lorentz hyperboloid: ",
      paste(utils::head(bad, 5L), collapse = ", "),
      if (length(bad) > 5L) ", ...",
      call. = FALSE
    )
  }
  
  ok
}

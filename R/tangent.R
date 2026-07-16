#' Normalize a Vector to the Lorentz Hyperboloid
#'
#' Internal helper used by tangent-space approximations.
#'
#' @param x Numeric vector of length 3.
#' @return A point on the upper sheet of the Lorentz hyperboloid.
normalize_to_H2 <- function(x) {
  x <- as.numeric(x)
  if (length(x) != 3L) {
    stop("`x` must have length 3.", call. = FALSE)
  }
  
  norm_l <- as.numeric(lorentz_norm(matrix(x, nrow = 1)))
  if (!is.finite(norm_l) || norm_l >= 0) {
    x <- c(sqrt(1 + sum(x[-1]^2)), x[-1])
    norm_l <- as.numeric(lorentz_norm(matrix(x, nrow = 1)))
  }
  
  x <- x / sqrt(-norm_l)
  if (x[1] < 0) x <- -x
  x
}

#' Extrinsic Mean Projected to the Lorentz Hyperboloid
#'
#' Compute the Euclidean coordinate mean of Lorentz points and project it back
#' to the upper sheet of the hyperboloid. This is used as the default reference
#' point for tangent-space approximations.
#'
#' @param Z `n x 3` Lorentz hyperboloid coordinates.
#' @return A length-3 Lorentz hyperboloid point.
#' @examples
#' H <- simulate_H2(10, r_max = 2)
#' mu <- extrinsic_lorentz_mean(H$Z)
#' check_hyperboloid(mu)
#' @export
extrinsic_lorentz_mean <- function(Z) {
  Z <- as_lorentz_matrix(Z, "Z")
  check_hyperboloid(Z, stop_on_fail = TRUE)
  normalize_to_H2(colMeans(Z))
}

#' Lorentz Logarithmic Map
#'
#' Map points on the Lorentz hyperboloid to the tangent space at `mu`.
#'
#' @param Z `n x 3` Lorentz hyperboloid coordinates.
#' @param mu Reference point on the Lorentz hyperboloid. If `NULL`, the
#'   normalized extrinsic mean is used.
#' @param tol Numerical tolerance for zero-distance points.
#' @return An `n x 3` matrix of tangent vectors in `T_mu H^2`.
#' @examples
#' H <- simulate_H2(5, r_max = 2)
#' V <- lorentz_log_map(H$Z)
#' dim(V)
#' @export
lorentz_log_map <- function(Z, mu = NULL, tol = 1e-10) {
  Z <- as_lorentz_matrix(Z, "Z")
  check_hyperboloid(Z, stop_on_fail = TRUE)
  if (is.null(mu)) {
    mu <- extrinsic_lorentz_mean(Z)
  } else {
    mu <- as_lorentz_matrix(mu, "mu")
    if (nrow(mu) != 1L) {
      stop("`mu` must be a single Lorentz point.", call. = FALSE)
    }
    check_hyperboloid(mu, stop_on_fail = TRUE)
    mu <- as.numeric(mu)
  }
  
  ip <- as.numeric(lorentz_inner_matrix(Z, matrix(mu, nrow = 1)))
  arg <- pmax(-ip, 1)
  d <- acosh(arg)
  coef <- rep(1, length(d))
  nz <- d > tol
  coef[nz] <- d[nz] / sinh(d[nz])
  
  V <- (Z - matrix(cosh(d), nrow(Z), length(mu)) *
          matrix(mu, nrow(Z), length(mu), byrow = TRUE)) *
    matrix(coef, nrow(Z), length(mu))
  V[!nz, ] <- 0
  colnames(V) <- c("v0", "v1", "v2")
  V
}

lorentz_inner_matrix <- function(A, B = A) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  storage.mode(A) <- "double"
  storage.mode(B) <- "double"
  if (ncol(A) != 3L || ncol(B) != 3L) {
    stop("Lorentz inner-product matrices must have three columns.", call. = FALSE)
  }
  (-A[, 1, drop = FALSE]) %*% t(B[, 1, drop = FALSE]) +
    A[, 2:3, drop = FALSE] %*% t(B[, 2:3, drop = FALSE])
}

tangent_basis_at <- function(mu) {
  mu <- as.numeric(mu)
  candidates <- rbind(
    c(0, 1, 0),
    c(0, 0, 1),
    c(mu[2], mu[1], 0),
    c(mu[3], 0, mu[1])
  )
  basis <- matrix(0, 0, 3)
  for (a in seq_len(nrow(candidates))) {
    v <- candidates[a, ]
    v <- v + as.numeric(lorentz_inner(matrix(v, nrow = 1), matrix(mu, nrow = 1))) * mu
    if (nrow(basis)) {
      for (b in seq_len(nrow(basis))) {
        v <- v - as.numeric(lorentz_inner(matrix(v, nrow = 1), matrix(basis[b, ], nrow = 1))) *
          basis[b, ]
      }
    }
    nv <- as.numeric(lorentz_inner(matrix(v, nrow = 1), matrix(v, nrow = 1)))
    if (is.finite(nv) && nv > 1e-10) {
      basis <- rbind(basis, v / sqrt(nv))
    }
    if (nrow(basis) == 2L) break
  }
  if (nrow(basis) != 2L) {
    stop("Could not construct a Lorentz-orthonormal tangent basis.", call. = FALSE)
  }
  basis
}

#' Tangent-Space Coordinates
#'
#' Map Lorentz hyperboloid points to two-dimensional coordinates in a
#' Lorentz-orthonormal basis of `T_mu H^2`.
#'
#' @param Z `n x 3` Lorentz hyperboloid coordinates.
#' @param mu Reference point. If `NULL`, the normalized extrinsic mean is used.
#' @return An `n x 2` matrix of tangent coordinates.
#' @examples
#' H <- simulate_H2(5, r_max = 2)
#' tangent_coordinates(H$Z)
#' @export
tangent_coordinates <- function(Z, mu = NULL) {
  Z <- as_lorentz_matrix(Z, "Z")
  if (is.null(mu)) {
    mu <- extrinsic_lorentz_mean(Z)
  } else {
    mu <- as_lorentz_matrix(mu, "mu")
    if (nrow(mu) != 1L) {
      stop("`mu` must be a single Lorentz point.", call. = FALSE)
    }
    check_hyperboloid(mu, stop_on_fail = TRUE)
    mu <- as.numeric(mu)
  }
  
  V <- lorentz_log_map(Z, mu)
  B <- tangent_basis_at(mu)
  coords <- cbind(
    as.numeric(lorentz_inner_matrix(V, matrix(B[1, ], nrow = 1))),
    as.numeric(lorentz_inner_matrix(V, matrix(B[2, ], nrow = 1)))
  )
  colnames(coords) <- c("v1", "v2")
  coords
}

#' Tangent-Space Distance Matrix
#'
#' Compute a distance matrix after mapping Lorentz hyperboloid points to
#' `T_mu H^2` with the logarithmic map. Distances are computed using the
#' Lorentz inner product restricted to the tangent space, which is positive
#' definite. Equivalently, this is Euclidean distance in a Lorentz-orthonormal
#' tangent basis.
#'
#' @param Z `n x 3` Lorentz hyperboloid coordinates.
#' @param mu Reference point. If `NULL`, the normalized extrinsic mean is used.
#' @param tol Numerical tolerance used in the logarithmic map.
#' @return A symmetric `n x n` tangent-space distance matrix.
#' @examples
#' H <- simulate_H2(5, r_max = 2)
#' tangent_distance_matrix(H$Z)
#' @export
tangent_distance_matrix <- function(Z, mu = NULL, tol = 1e-10) {
  Z <- as_lorentz_matrix(Z, "Z")
  check_hyperboloid(Z, stop_on_fail = TRUE)
  if (is.null(mu)) {
    mu <- extrinsic_lorentz_mean(Z)
  } else {
    mu <- as_lorentz_matrix(mu, "mu")
    if (nrow(mu) != 1L) {
      stop("`mu` must be a single Lorentz point.", call. = FALSE)
    }
    check_hyperboloid(mu, stop_on_fail = TRUE)
    mu <- as.numeric(mu)
  }
  .tangent_distance_matrix_cpp(Z, matrix(mu, nrow = 1), tol = tol)
}

resolve_geometry_distance <- function(Z,
                                      dMat = NULL,
                                      geometry = c("lorentz", "tangent"),
                                      tangent_ref = NULL) {
  geometry <- match.arg(geometry)
  Z <- as_lorentz_matrix(Z, "Z")
  check_hyperboloid(Z, stop_on_fail = TRUE)
  
  if (!is.null(dMat)) {
    return(list(
      D = as.matrix(dMat),
      label = if (geometry == "lorentz") {
        "Lorentz hyperboloid"
      } else {
        "Tangent-space approximation (user-supplied distance matrix)"
      },
      tangent_ref = tangent_ref
    ))
  }
  
  if (geometry == "lorentz") {
    return(list(
      D = hyperbolic_distance_matrix(Z, check = FALSE),
      label = "Lorentz hyperboloid",
      tangent_ref = NULL
    ))
  }
  
  mu <- if (is.null(tangent_ref)) extrinsic_lorentz_mean(Z) else {
    tangent_ref <- as_lorentz_matrix(tangent_ref, "tangent_ref")
    if (nrow(tangent_ref) != 1L) {
      stop("`tangent_ref` must be a single Lorentz point.", call. = FALSE)
    }
    check_hyperboloid(tangent_ref, stop_on_fail = TRUE)
    as.numeric(tangent_ref)
  }
  
  list(
    D = tangent_distance_matrix(Z, mu = mu),
    label = "Tangent-space approximation",
    tangent_ref = mu
  )
}

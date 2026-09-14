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
  if (any(!is.finite(x))) {
    stop("`x` must contain only finite values.", call. = FALSE)
  }

  if (x[1] < 0) x <- -x
  spatial_norm <- sqrt(sum(x[-1]^2))
  norm_sq <- (x[1] - spatial_norm) * (x[1] + spatial_norm)
  if (is.finite(norm_sq) && norm_sq > 0) {
    spatial <- x[-1] / sqrt(norm_sq)
  } else {
    spatial <- x[-1]
  }
  c(sqrt(1 + sum(spatial^2)), spatial)
}

#' Extrinsic Mean Projected to the Lorentz Hyperboloid
#'
#' Compute the Euclidean coordinate mean of Lorentz points and project it back
#' to the upper sheet of the hyperboloid. This is used as the default reference
#' point for tangent-space approximations.
#'
#' For valid upper-sheet points, let \eqn{m = n^{-1}\sum_i z_i} be the ambient
#' coordinate mean of the \eqn{n} input rows and let \eqn{\langle\cdot,\cdot\rangle_L}
#' denote [lorentz_inner()]. The reference is
#' \deqn{\mu = m / \sqrt{-\langle m,m\rangle_L}.}
#' This normalization puts the mean on the unit hyperboloid; it is not the
#' intrinsic Frechet mean. The implementation evaluates the normalization
#' stably and reconstructs the positive time coordinate. A zero spatial mean
#' gives \eqn{\mu=(1,0,0)}. Estimate the reference from training locations and
#' reuse it when mapping new locations.
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

  spatial_mean <- colMeans(Z[, 2:3, drop = FALSE])
  spatial_mean_norm <- sqrt(sum(spatial_mean^2))
  if (spatial_mean_norm == 0) {
    return(c(1, 0, 0))
  }

  mean_direction <- spatial_mean / spatial_mean_norm
  spatial_norm <- sqrt(rowSums(Z[, 2:3, drop = FALSE]^2))
  radial <- asinh(spatial_norm)
  unit_direction <- matrix(0, nrow(Z), 2L)
  nonzero <- spatial_norm > 0
  unit_direction[nonzero, ] <-
    Z[nonzero, 2:3, drop = FALSE] / spatial_norm[nonzero]

  direction_gap <- rowSums(
    (unit_direction - matrix(mean_direction, nrow(Z), 2L, byrow = TRUE))^2
  ) / 2
  spatial_projection <- as.numeric(
    Z[, 2:3, drop = FALSE] %*% mean_direction
  )

  # These are stable evaluations of m0 - ||m_s|| and m0 + ||m_s||,
  # where m is the ambient coordinate mean.
  minus_part <- mean(exp(-radial) + spatial_norm * direction_gap)
  plus_part <- mean(sqrt(1 + spatial_norm^2) + spatial_projection)
  mean_lorentz_scale <- sqrt(minus_part * plus_part)
  if (!is.finite(mean_lorentz_scale) || mean_lorentz_scale <= 0) {
    stop("Could not normalize the extrinsic Lorentz mean.", call. = FALSE)
  }

  spatial <- spatial_mean / mean_lorentz_scale
  c(sqrt(1 + sum(spatial^2)), spatial)
}

stable_tangent_log_coordinates <- function(Z, mu, tol = 1e-10) {
  mu <- as.numeric(mu)
  mu_spatial_norm <- sqrt(sum(mu[2:3]^2))
  mu_radial <- asinh(mu_spatial_norm)

  if (mu_spatial_norm > 0) {
    mu_direction <- mu[2:3] / mu_spatial_norm
    mu_time <- sqrt(1 + mu_spatial_norm^2)
    basis <- rbind(
      c(mu_spatial_norm, mu_time * mu_direction),
      c(0, -mu_direction[2], mu_direction[1])
    )
  } else {
    mu_direction <- c(1, 0)
    mu_time <- 1
    basis <- rbind(c(0, 1, 0), c(0, 0, 1))
  }

  spatial_norm <- sqrt(rowSums(Z[, 2:3, drop = FALSE]^2))
  radial <- asinh(spatial_norm)
  unit_direction <- matrix(0, nrow(Z), 2L)
  nonzero <- spatial_norm > 0
  unit_direction[nonzero, ] <-
    Z[nonzero, 2:3, drop = FALSE] / spatial_norm[nonzero]

  if (mu_spatial_norm > 0) {
    direction_gap <- rowSums(
      (unit_direction - matrix(mu_direction, nrow(Z), 2L, byrow = TRUE))^2
    ) / 2
    radial_component <-
      sinh(radial - mu_radial) - spatial_norm * mu_time * direction_gap
    angular_component <-
      -Z[, 2] * mu_direction[2] + Z[, 3] * mu_direction[1]
  } else {
    radial_component <- Z[, 2]
    angular_component <- Z[, 3]
  }

  distance <- hyperbolic_distance(
    Z,
    matrix(mu, nrow = 1L),
    check = FALSE,
    tol = tol
  )
  multiplier <- rep(1, length(distance))
  nonzero_distance <- distance > tol
  ordinary <- nonzero_distance & distance < 350
  multiplier[ordinary] <- distance[ordinary] / sinh(distance[ordinary])
  large <- nonzero_distance & !ordinary
  multiplier[large] <- 2 * distance[large] * exp(-distance[large]) /
    (1 - exp(-2 * distance[large]))

  coordinates <- cbind(radial_component, angular_component) * multiplier
  coordinates[!nonzero_distance, ] <- 0
  colnames(coordinates) <- c("v1", "v2")
  list(coordinates = coordinates, basis = basis, distance = distance)
}

#' Lorentz Logarithmic Map
#'
#' Map points on the Lorentz hyperboloid to the tangent space at `mu`.
#'
#' For unit-hyperboloid points \eqn{z} and \eqn{\mu}, define
#' \eqn{a=-\langle\mu,z\rangle_L} and the intrinsic distance
#' \eqn{d=\operatorname{acosh}(a)} using [lorentz_inner()] and
#' [hyperbolic_distance()]. The logarithmic map is
#' \deqn{\log_\mu(z)=\frac{d}{\sinh(d)}(z-a\mu),\qquad d>0.}
#' The result is a tangent vector orthogonal to \eqn{\mu} under the Lorentz
#' inner product. At \eqn{z=\mu} the continuous extension is the zero vector.
#' Numerically, `tol` is compared with the distance \eqn{d}, not with
#' \eqn{a-1}: rows with \eqn{d\leq\mathrm{tol}} are set to zero. Other rows
#' use stable, algebraically equivalent coordinate formulas.
#'
#' @param Z `n x 3` Lorentz hyperboloid coordinates.
#' @param mu Reference point on the Lorentz hyperboloid. If `NULL`, the
#'   normalized extrinsic mean is used.
#' @param tol Intrinsic-distance threshold for setting a log-map row to zero;
#'   the default is `1e-10` in unit-curvature distance units.
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
  
  mapped <- stable_tangent_log_coordinates(Z, mu, tol = tol)
  V <- mapped$coordinates %*% mapped$basis
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
  spatial_norm <- sqrt(sum(mu[2:3]^2))
  if (spatial_norm == 0) {
    return(rbind(c(0, 1, 0), c(0, 0, 1)))
  }
  direction <- mu[2:3] / spatial_norm
  time <- sqrt(1 + spatial_norm^2)
  rbind(
    c(spatial_norm, time * direction),
    c(0, -direction[2], direction[1])
  )
}

#' Tangent-Space Coordinates
#'
#' Map Lorentz hyperboloid points to two-dimensional coordinates in a
#' Lorentz-orthonormal basis of `T_mu H^2`.
#'
#' The columns express [lorentz_log_map()] in an orthonormal tangent basis.
#' Their ordinary Euclidean norm equals the Lorentz norm of the tangent
#' vector. This routine uses the log map's default zero-distance threshold
#' `1e-10`. To compare training and new points, pass the same `mu` in both
#' calls; separately estimated centers define different coordinate systems.
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
  
  stable_tangent_log_coordinates(Z, mu)$coordinates
}

#' Tangent-Space Distance Matrix
#'
#' Compute a distance matrix after mapping Lorentz hyperboloid points to
#' `T_mu H^2` with the logarithmic map. Distances are computed using the
#' Lorentz inner product restricted to the tangent space, which is positive
#' definite. Equivalently, this is Euclidean distance in a Lorentz-orthonormal
#' tangent basis.
#'
#' With \eqn{v_i=\log_\mu(z_i)} and \eqn{v_j=\log_\mu(z_j)}, the distance is
#' \deqn{d_{T,\mu}(z_i,z_j)=\sqrt{\langle v_i-v_j,v_i-v_j\rangle_L}.}
#' Here \eqn{i,j} index input rows and [lorentz_inner()] defines the inner
#' product. Its restriction to the common tangent plane is positive definite.
#' These distances define the Tangent-E weighting approximation: radial
#' distances from `mu` are preserved in exact arithmetic, but general pairwise
#' hyperbolic distances are not. The implementation uses coordinate differences
#' and gives identical mapped points exactly zero distance. The `tol` argument
#' applies to each point's distance from `mu` in the log map, not to arbitrary
#' pairwise distances between distinct tangent coordinates.
#'
#' @param Z `n x 3` Lorentz hyperboloid coordinates.
#' @param mu Reference point. If `NULL`, the normalized extrinsic mean is used.
#' @param tol Intrinsic-distance threshold used by [lorentz_log_map()]; points
#'   within this distance of `mu` are mapped to zero.
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

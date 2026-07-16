#' Lorentz-Hyperbolic Moran's I Diagnostic for LHWR Residuals
#'
#' Compute the Lorentz-hyperbolic Moran's I diagnostic (LH-Moran's I) for
#' hyperbolic residual autocorrelation from a fitted LHWR model. By
#' default, the neighborhood weights define hyperbolic neighborhoods from
#' Lorentz-hyperbolic distances stored in the fitted object.
#'
#' @param object An object returned by [lhwr()].
#' @param W Optional numeric `n x n` neighborhood weight matrix. If supplied, it is
#'   used directly after diagonal entries are set to zero.
#' @param distance Distance used when `W = NULL`. `"lorentz"` uses intrinsic
#'   Lorentz-hyperbolic distance on the Lorentz hyperboloid and is the default.
#'   `"euclidean"` is an explicit fallback using ordinary Euclidean distance on
#'   the projected `(z1, z2)` coordinates.
#' @param kernel Kernel used to convert distances into weights.
#' @param bandwidth Optional fixed bandwidth. If omitted and `adaptive = TRUE`,
#'   the fitted model bandwidth or `k` is used as the adaptive neighbor count.
#' @param adaptive Logical. If `TRUE`, use adaptive nearest-neighbor radii.
#' @param k Optional adaptive neighbor count. Overrides `bandwidth` for
#'   adaptive weights.
#' @param permutations Number of random permutations for the p-value. Use zero
#'   to skip permutation inference.
#' @param row_standardize Logical. If `TRUE`, row-standardize weights.
#' @param seed Optional random seed used for permutation inference.
#' @return An object of class `"hgwr_moran"` containing Moran's I, expected I,
#'   permutation p-value, and settings.
#' @examples
#' set.seed(12)
#' H <- simulate_H2(60, r_max = 3)
#' x <- rnorm(60)
#' y <- 1 + 2 * x + rnorm(60, sd = 0.3)
#' fit <- lhwr(y ~ x, data.frame(y, x), Z = H$Z, bw = 25)
#' moran.lhwr(fit, permutations = 99)
#' @export
moran.lhwr <- function(object,
                       W = NULL,
                       distance = c("lorentz", "euclidean"),
                       kernel = c("bisquare", "gaussian", "exponential"),
                       bandwidth = NULL,
                       adaptive = TRUE,
                       k = NULL,
                       permutations = 999,
                       row_standardize = TRUE,
                       seed = NULL) {
  distance <- match.arg(distance)
  kernel <- match.arg(kernel)
  e <- hgwr_residual_vector(object)
  n <- length(e)
  W <- hgwr_moran_weights(
    object = object,
    W = W,
    distance = distance,
    kernel = kernel,
    bandwidth = bandwidth,
    adaptive = adaptive,
    k = k,
    row_standardize = row_standardize
  )
  I <- moran_statistic(e, W)
  expected_I <- -1 / (n - 1)
  p_value <- NA_real_
  permuted_I <- numeric()
  if (permutations > 0 && is.finite(I)) {
    if (!is.null(seed)) {
      old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) .Random.seed else NULL
      on.exit({
        if (is.null(old_seed)) {
          rm(".Random.seed", envir = .GlobalEnv)
        } else {
          assign(".Random.seed", old_seed, envir = .GlobalEnv)
        }
      }, add = TRUE)
      set.seed(seed)
    }
    permuted_I <- replicate(permutations, moran_statistic(sample(e), W))
    p_value <- permutation_p_value(I, permuted_I, center = expected_I)
  }

  out <- list(
    I = I,
    expected_I = expected_I,
    variance_I = NA_real_,
    p_value = p_value,
    permutations = permutations,
    n = n,
    distance = distance,
    kernel = kernel,
    bandwidth = attr(W, "bandwidth"),
    adaptive = adaptive,
    row_standardize = row_standardize,
    permuted_I = permuted_I
  )
  class(out) <- "hgwr_moran"
  out
}

#' Lorentz-Hyperbolic LISA Diagnostic (LH-LISA) for LHWR Residuals
#'
#' Compute the Lorentz-hyperbolic LISA diagnostic (LH-LISA), that is local
#' Moran's I under hyperbolic neighborhoods, for residuals from a fitted
#' LHWR model. The default weights use Lorentz-hyperbolic distance to
#' define each hyperbolic neighborhood.
#'
#' @inheritParams moran.lhwr
#' @return A data frame of class `"hgwr_local_moran"` with local Moran values,
#'   permutation p-values, z-scores, and Lorentz-hyperbolic LISA cluster labels.
#' @examples
#' set.seed(13)
#' H <- simulate_H2(60, r_max = 3)
#' x <- rnorm(60)
#' y <- 1 + 2 * x + rnorm(60, sd = 0.3)
#' fit <- lhwr(y ~ x, data.frame(y, x), Z = H$Z, bw = 25)
#' lm <- local_moran.lhwr(fit, permutations = 99)
#' head(lm)
#' plot(lm)
#' @export
local_moran.lhwr <- function(object,
                             W = NULL,
                             distance = c("lorentz", "euclidean"),
                             kernel = c("bisquare", "gaussian", "exponential"),
                             bandwidth = NULL,
                             adaptive = TRUE,
                             k = NULL,
                             permutations = 999,
                             row_standardize = TRUE,
                             seed = NULL) {
  distance <- match.arg(distance)
  kernel <- match.arg(kernel)
  e <- hgwr_residual_vector(object)
  n <- length(e)
  W <- hgwr_moran_weights(
    object = object,
    W = W,
    distance = distance,
    kernel = kernel,
    bandwidth = bandwidth,
    adaptive = adaptive,
    k = k,
    row_standardize = row_standardize
  )
  z <- e - mean(e)
  m2 <- sum(z^2) / n
  lag_z <- as.numeric(W %*% z)
  local_I <- if (m2 > .Machine$double.eps) z * lag_z / m2 else rep(NA_real_, n)
  p_value <- rep(NA_real_, n)
  z_score <- rep(NA_real_, n)

  if (permutations > 0 && all(is.finite(local_I))) {
    if (!is.null(seed)) {
      old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) .Random.seed else NULL
      on.exit({
        if (is.null(old_seed)) {
          rm(".Random.seed", envir = .GlobalEnv)
        } else {
          assign(".Random.seed", old_seed, envir = .GlobalEnv)
        }
      }, add = TRUE)
      set.seed(seed)
    }
    permuted <- matrix(NA_real_, nrow = permutations, ncol = n)
    for (b in seq_len(permutations)) {
      zb <- sample(z)
      permuted[b, ] <- zb * as.numeric(W %*% zb) / m2
    }
    p_value <- vapply(seq_len(n), function(i) {
      permutation_p_value(local_I[i], permuted[, i], center = mean(permuted[, i], na.rm = TRUE))
    }, numeric(1))
    perm_mean <- colMeans(permuted, na.rm = TRUE)
    perm_sd <- apply(permuted, 2L, stats::sd, na.rm = TRUE)
    z_score <- (local_I - perm_mean) / perm_sd
    z_score[!is.finite(z_score)] <- NA_real_
  }

  cluster_type <- lisa_cluster_type(z, lag_z, p_value)
  out <- data.frame(
    observation_id = seq_len(n),
    local_I = local_I,
    z_score = z_score,
    p_value = p_value,
    cluster_type = cluster_type,
    residual = e,
    residual_centered = z,
    spatial_lag = lag_z
  )
  attr(out, "distance") <- distance
  attr(out, "kernel") <- kernel
  attr(out, "bandwidth") <- attr(W, "bandwidth")
  attr(out, "adaptive") <- adaptive
  attr(out, "row_standardize") <- row_standardize
  class(out) <- c("hgwr_local_moran", "data.frame")
  out
}

#' @export
print.hgwr_moran <- function(x, ...) {
  cat("Lorentz-hyperbolic Moran's I diagnostic (LH-Moran's I)\n")
  cat("  I:          ", format(signif(x$I, 6)), "\n", sep = "")
  cat("  Expected I: ", format(signif(x$expected_I, 6)), "\n", sep = "")
  cat("  p-value:    ", format(signif(x$p_value, 6)), "\n", sep = "")
  cat("  Distance:   ", x$distance, "\n", sep = "")
  cat("  Kernel:     ", x$kernel, "\n", sep = "")
  cat("  n:          ", x$n, "\n", sep = "")
  invisible(x)
}

#' @export
summary.hgwr_moran <- function(object, ...) {
  print(object)
  invisible(object)
}

#' Plot Lorentz-Hyperbolic LISA Diagnostics (LH-LISA)
#'
#' @param x An object returned by `local_moran.lhwr()`.
#' @param significant_only Logical. If `TRUE`, nonsignificant observations are
#'   plotted in grey.
#' @param alpha Significance level used when `significant_only = TRUE`.
#' @param ... Additional arguments passed to [graphics::plot()].
#' @return Invisibly returns `x`.
#' @export
plot.hgwr_local_moran <- function(x, significant_only = FALSE, alpha = 0.05, ...) {
  values <- x$local_I
  cols <- lisa_cluster_colors(x$cluster_type)
  if (isTRUE(significant_only)) {
    cols[is.na(x$p_value) | x$p_value > alpha] <- "grey75"
  }
  graphics::plot(
    x$observation_id,
    values,
    pch = 19,
    col = cols,
    xlab = "Observation",
    ylab = "Lorentz-hyperbolic local Moran's I",
    main = "Lorentz-hyperbolic LISA diagnostic (LH-LISA)",
    ...
  )
  graphics::abline(h = 0, col = "grey60", lty = 2)
  graphics::legend(
    "topright",
    legend = names(lisa_cluster_palette()),
    col = lisa_cluster_palette(),
    pch = 19,
    bty = "n",
    cex = 0.8
  )
  invisible(x)
}

hgwr_residual_vector <- function(object) {
  if (!is_lhwr_object(object)) {
    stop("`object` must be an LHWR model object.", call. = FALSE)
  }
  e <- object$residuals
  if (is.null(e) && !is.null(object$SDF$residual)) {
    e <- object$SDF$residual
  }
  if (is.null(e)) {
    stop("No residuals found in `object`.", call. = FALSE)
  }
  e <- as.numeric(e)
  if (any(!is.finite(e))) {
    stop("Residuals must be finite.", call. = FALSE)
  }
  e
}

hgwr_moran_weights <- function(object, W, distance, kernel, bandwidth,
                               adaptive, k, row_standardize) {
  if (!is.null(W)) {
    W <- as.matrix(W)
    if (!is.numeric(W) || nrow(W) != ncol(W)) {
      stop("`W` must be a square numeric matrix.", call. = FALSE)
    }
    diag(W) <- 0
    attr(W, "bandwidth") <- NA_real_
    return(standardize_weights(W, row_standardize))
  }

  Z <- object$Z
  if (is.null(Z)) {
    stop("`object` does not contain Lorentz coordinates and `W` is NULL.", call. = FALSE)
  }
  if (isTRUE(adaptive) && distance == "lorentz") {
    W_knn <- hgwr_moran_weights_from_stored_knn(
      object = object,
      kernel = kernel,
      bandwidth = bandwidth,
      k = k,
      row_standardize = row_standardize
    )
    if (!is.null(W_knn)) {
      return(W_knn)
    }
  }
  Z <- as_lorentz_matrix(Z, "Z")
  D <- switch(
    distance,
    lorentz = hgwr_moran_lorentz_distance(object, Z),
    euclidean = as.matrix(stats::dist(Z[, 2:3, drop = FALSE]))
  )
  diag(D) <- 0
  n <- nrow(D)
  bw <- resolve_moran_bandwidth(object, D, bandwidth = bandwidth, adaptive = adaptive, k = k)
  W <- matrix(0, n, n)
  if (isTRUE(adaptive)) {
    for (i in seq_len(n)) {
      W[i, ] <- hgwr_kernel(D[i, ], bw[i], kernel = kernel)
    }
  } else {
    W[] <- hgwr_kernel(D, bw, kernel = kernel)
  }
  diag(W) <- 0
  attr(W, "bandwidth") <- bw
  standardize_weights(W, row_standardize)
}

hgwr_moran_lorentz_distance <- function(object, Z) {
  if (
    identical(object$geometry, "lorentz") &&
      !is.null(object$distance) &&
      is.matrix(object$distance) &&
      nrow(object$distance) == nrow(Z) &&
      ncol(object$distance) == nrow(Z)
  ) {
    return(object$distance)
  }
  hyperbolic_distance_matrix(Z, check = FALSE)
}

hgwr_moran_weights_from_stored_knn <- function(object, kernel, bandwidth, k,
                                               row_standardize) {
  index <- object$knn_index
  distance <- object$knn_distance
  if (is.null(index) || is.null(distance)) {
    return(NULL)
  }
  index <- as.matrix(index)
  distance <- as.matrix(distance)
  n <- nrow(index)
  if (
    n < 1L ||
      nrow(distance) != n ||
      ncol(distance) != ncol(index) ||
      is.null(object$Z) ||
      nrow(object$Z) != n
  ) {
    return(NULL)
  }

  k <- if (!is.null(k)) k else bandwidth
  if (is.null(k)) {
    k <- object$LHWR.arguments$bw
  }
  if (is.null(k) || length(k) != 1L || !is.finite(k)) {
    return(NULL)
  }
  k <- as.integer(round(k))
  if (k < 1L || k > ncol(index)) {
    return(NULL)
  }

  W <- matrix(0, n, n)
  for (i in seq_len(n)) {
    bw_i <- distance[i, k]
    if (!is.finite(bw_i) || bw_i <= 0) {
      return(NULL)
    }
    idx <- index[i, seq_len(k)]
    dij <- distance[i, seq_len(k)]
    ok <- is.finite(idx) & idx >= 1L & idx <= n & is.finite(dij)
    if (any(ok)) {
      W[i, idx[ok]] <- hgwr_kernel(dij[ok], bw_i, kernel = kernel)
    }
  }
  diag(W) <- 0
  attr(W, "bandwidth") <- distance[, k]
  standardize_weights(W, row_standardize)
}

resolve_moran_bandwidth <- function(object, D, bandwidth, adaptive, k) {
  if (isTRUE(adaptive)) {
    k <- if (!is.null(k)) k else bandwidth
    if (is.null(k)) {
      k <- object$LHWR.arguments$bw
    }
    if (is.null(k) || length(k) != 1L || !is.finite(k)) {
      stop("Adaptive Moran weights require scalar `k` or `bandwidth`.", call. = FALSE)
    }
    k <- as.integer(round(k))
    n <- nrow(D)
    if (k < 1L || k >= n) {
      stop("Adaptive neighbor count must be between 1 and n - 1.", call. = FALSE)
    }
    apply(D, 1L, function(d) sort(d[d > 0], partial = k)[k])
  } else {
    if (is.null(bandwidth)) {
      bandwidth <- object$LHWR.arguments$bw
    }
    if (is.null(bandwidth) || length(bandwidth) != 1L || !is.finite(bandwidth) || bandwidth <= 0) {
      stop("Fixed Moran weights require a positive scalar `bandwidth`.", call. = FALSE)
    }
    bandwidth
  }
}

standardize_weights <- function(W, row_standardize) {
  if (isTRUE(row_standardize)) {
    rs <- rowSums(W)
    ok <- rs > .Machine$double.eps
    W[ok, ] <- W[ok, , drop = FALSE] / rs[ok]
  }
  W
}

moran_statistic <- function(e, W) {
  z <- e - mean(e)
  denom <- sum(z^2)
  S0 <- sum(W)
  if (denom <= .Machine$double.eps || S0 <= .Machine$double.eps) {
    return(NA_real_)
  }
  length(e) / S0 * as.numeric(t(z) %*% W %*% z) / denom
}

permutation_p_value <- function(observed, permuted, center = 0) {
  permuted <- permuted[is.finite(permuted)]
  if (!length(permuted) || !is.finite(observed)) return(NA_real_)
  (sum(abs(permuted - center) >= abs(observed - center)) + 1) / (length(permuted) + 1)
}

lisa_cluster_type <- function(z, lag_z, p_value, alpha = 0.05) {
  significant <- is.finite(p_value) & p_value <= alpha
  out <- rep("Not Significant", length(z))
  out[significant & z > 0 & lag_z > 0] <- "High-High"
  out[significant & z < 0 & lag_z < 0] <- "Low-Low"
  out[significant & z > 0 & lag_z < 0] <- "High-Low"
  out[significant & z < 0 & lag_z > 0] <- "Low-High"
  factor(out, levels = names(lisa_cluster_palette()))
}

lisa_cluster_palette <- function() {
  c(
    "High-High" = "firebrick",
    "Low-Low" = "steelblue",
    "High-Low" = "orange",
    "Low-High" = "purple",
    "Not Significant" = "grey65"
  )
}

lisa_cluster_colors <- function(cluster_type) {
  pal <- lisa_cluster_palette()
  pal[as.character(cluster_type)]
}

#' @rdname moran.lhwr
#' @export
moran.hgwr <- moran.lhwr
#' @rdname local_moran.lhwr
#' @export
local_moran.hgwr <- local_moran.lhwr

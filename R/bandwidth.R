#' Select LHWR Bandwidth
#'
#' Select the kernel bandwidth used by LHWR local regressions. Candidate
#' bandwidths are scored by either corrected Akaike information criterion
#' (`AICc`) or leave-one-out cross-validation (`CV`). The selected bandwidth is
#' the candidate with the smallest score.
#'
#' @param formula Model formula of the form `response ~ predictors`.
#' @param data Data frame containing the response and predictor variables used
#'   in `formula`.
#' @param Z Numeric `n x 3` matrix of Lorentz hyperboloid coordinates. Row `i`
#'   gives the location of observation `i` and must satisfy
#'   `<z_i,z_i>_L = -1` with positive time coordinate.
#' @param approach Bandwidth selection criterion. `"AICc"` minimizes the
#'   corrected Akaike information criterion computed from the local fitted
#'   values and hat-matrix trace. `"CV"` minimizes the leave-one-out prediction
#'   error.
#' @param kernel Kernel used to convert distances into local weights. Available
#'   choices are `"bisquare"`, `"gaussian"`, `"exponential"`, `"tricube"`, and
#'   `"boxcar"`.
#' @param adaptive Logical. If `TRUE`, scalar bandwidths are interpreted as
#'   neighbor counts and converted to location-specific radii. If `FALSE`,
#'   scalar bandwidths are interpreted as fixed radii, or fixed hyperbolic
#'   volumes when `bandwidth_type = "volume"`.
#' @param bandwidth_type Bandwidth scale. `"knn"` uses the radius to the
#'   requested nearest neighbor. `"volume"` uses a hyperbolic ball-volume
#'   interpretation and is available only for `geometry = "lorentz"`.
#' @param search Bandwidth search method. `"golden"` uses golden-section
#'   optimization over an automatically chosen interval. `"grid"` evaluates the
#'   supplied candidate values in `bw_grid`. If omitted, `"golden"` is used when
#'   `bw_grid` is `NULL`, and `"grid"` is used when `bw_grid` is supplied.
#' @param geometry Distance geometry used to build kernel weights. `"lorentz"`
#'   uses intrinsic hyperbolic distances on the Lorentz hyperboloid.
#'   `"tangent"` maps locations to one tangent plane by the Lorentz logarithmic
#'   map and uses the positive-definite Lorentz metric induced on that tangent
#'   plane. Equivalently, these distances can be represented as Euclidean
#'   distances in a Lorentz-orthonormal tangent basis.
#' @param bw_grid Optional numeric vector of candidate bandwidths for
#'   `search = "grid"`. In adaptive fits these are neighbor counts; in fixed
#'   fits these are radii unless `bandwidth_type = "volume"`.
#' @param lower Optional lower search bound for `search = "golden"`.
#' @param upper Optional upper search bound for `search = "golden"`.
#' @param tol Golden-section stopping tolerance.
#' @param max_iter Maximum number of golden-section iterations.
#' @param pilot_k Pilot neighbor count used to estimate local point density for
#'   adaptive volume bandwidths. Ignored unless `adaptive = TRUE` and
#'   `bandwidth_type = "volume"`.
#' @param tangent_ref Optional reference point for `geometry = "tangent"`.
#'   If omitted, the normalized extrinsic Lorentz mean is used.
#' @param dMat Optional precomputed `n x n` distance matrix. If supplied, it is
#'   used directly and `geometry` is used only for labeling and validation of
#'   related arguments.
#' @param engine Computational engine. `"dense"` builds or uses a full distance
#'   matrix. `"knn"` uses exact Lorentz nearest-neighbor lists and avoids storing
#'   the full distance matrix; it is available for adaptive Lorentz kNN fits
#'   with compact kernels.
#' @param n_threads Number of OpenMP threads used by `engine = "knn"` when
#'   OpenMP is available. The default `1` preserves serial execution.
#' @return An object of class `"bw.lhwr"`, a list containing the selected
#'   bandwidth (`bw`), selection criterion (`approach`), kernel, geometry,
#'   search method, selected score, AICc/CV fields when applicable, hat-matrix
#'   trace (`trS`), residual sum of squares (`RSS`), the evaluated
#'   `selection.table`, and the original function call.
#' @examples
#' set.seed(20260615)
#' n <- 120
#' H <- simulate_H2(n, r_max = 5, uniform_area = TRUE)
#' Z <- H$Z
#' D <- hyperbolic_distance_matrix(Z, check = FALSE)
#'
#' true_k <- 25
#' h <- apply(D, 1, function(d) sort(d[d > 0])[true_k])
#' centers <- round(seq(10, n - 10, length.out = 5))
#' Phi <- sapply(centers, function(j) {
#'   hgwr_kernel(D[, j], h[j], kernel = "bisquare")
#' })
#'
#' beta0 <- 0.4 + 2 * (1.45 * Phi[, 1] - 1.25 * Phi[, 4])
#' beta1 <- 1.0 + 2 * (1.65 * Phi[, 2] - 1.20 * Phi[, 5])
#' beta2 <- -0.8 + 2 * (1.55 * Phi[, 3] - 1.10 * Phi[, 4])
#'
#' x1 <- as.numeric(scale(rnorm(n)))
#' x2 <- as.numeric(scale(rnorm(n)))
#' eta <- beta0 + beta1 * x1 + beta2 * x2
#' y <- eta + rnorm(n, sd = sqrt(var(eta) / 8))
#' dat <- data.frame(y = y, x1 = x1, x2 = x2)
#'
#' bw <- bw.lhwr(
#'   y ~ x1 + x2,
#'   data = dat,
#'   Z = Z,
#'   approach = "CV",
#'   kernel = "bisquare",
#'   adaptive = TRUE,
#'   geometry = "lorentz",
#'   bw_grid = seq(10, 80, by = 5)
#' )
#' bw
#' head(bw$selection.table)
#' @export
bw.lhwr <- function(formula, data, Z,
                    approach = c("AICc", "CV"),
                    kernel = hgwr_kernel_names(),
                    adaptive = TRUE,
                    bandwidth_type = c("knn", "volume"),
                    search = NULL,
                    geometry = c("lorentz", "tangent"),
                    bw_grid = NULL,
                    lower = NULL,
                    upper = NULL,
                    tol = 1e-4,
                    max_iter = 200L,
                    pilot_k = NULL,
                    tangent_ref = NULL,
                    dMat = NULL,
                    engine = c("dense", "knn"),
                    n_threads = 1L) {
  this.call <- match.call()
  approach <- match.arg(approach)
  kernel <- match.arg(kernel)
  bandwidth_type <- match.arg(bandwidth_type)
  geometry <- match.arg(geometry)
  engine <- match.arg(engine)
  n_threads <- resolve_n_threads_arg(n_threads)
  search <- resolve_bw_search(search, bw_grid)
  if (geometry == "tangent" && bandwidth_type == "volume") {
    stop("`bandwidth_type = \"volume\"` is only defined for `geometry = \"lorentz\"`.", call. = FALSE)
  }
  if (engine == "knn") {
    validate_knn_engine(geometry, adaptive, bandwidth_type, kernel, dMat, search)
  }
  
  mf <- hgwr_model_frame(formula, data)
  y <- mf$y
  X <- mf$X
  colnames(X) <- clean_coef_names(colnames(X))
  Z <- as_lorentz_matrix(Z, "Z")
  check_hyperboloid(Z, stop_on_fail = TRUE)
  distance <- if (engine == "dense") {
    resolve_geometry_distance(
      Z,
      dMat = dMat,
      geometry = geometry,
      tangent_ref = tangent_ref
    )
  } else {
    list(D = NULL, label = "Lorentz hyperboloid", tangent_ref = NULL)
  }
  D <- distance$D
  
  if (search == "grid" && is.null(bw_grid)) {
    bw_grid <- if (engine == "knn") {
      default_knn_bw_grid(length(y))
    } else {
      default_bw_grid(D, adaptive = adaptive, bandwidth_type = bandwidth_type)
    }
  }
  
  pilot_k_resolved <- if (adaptive && bandwidth_type == "volume") {
    resolve_pilot_k(length(y), pilot_k)
  } else {
    0L
  }
  selection <- if (engine == "knn") {
    if (search == "golden") {
      bw_hgwr_knn_golden_search(
        y = y,
        X = X,
        Z = Z,
        approach = approach,
        kernel = kernel,
        lower = lower,
        upper = upper,
        tol = tol,
        max_iter = max_iter,
        n_threads = n_threads
      )
    } else {
      score_bw_candidates_knn(
        y = y,
        X = X,
        Z = Z,
        bw = as.numeric(bw_grid),
        approach = approach,
        kernel = kernel,
        n_threads = n_threads
      )
    }
  } else if (search == "golden") {
    bw_hgwr_golden_search(
      y = y,
      X = X,
      D = D,
      approach = approach,
      kernel = kernel,
      adaptive = adaptive,
      bandwidth_type = bandwidth_type,
      pilot_k = pilot_k_resolved,
      lower = lower,
      upper = upper,
      tol = tol,
      max_iter = max_iter
    )
  } else {
    score_bw_candidates(
      y = y,
      X = X,
      D = D,
      bw = as.numeric(bw_grid),
      approach = approach,
      kernel = kernel,
      adaptive = adaptive,
      bandwidth_type = bandwidth_type,
      pilot_k = pilot_k_resolved
    )
  }
  best_id <- which.min(selection$score)
  if (!length(best_id) || !is.finite(selection$score[best_id])) {
    stop("No finite bandwidth selection score was found.", call. = FALSE)
  }
  result <- list(
    bw = selection$bw[best_id],
    approach = approach,
    kernel = kernel,
    adaptive = adaptive,
    bandwidth_type = bandwidth_type,
    search = search,
    lower = if (search == "golden") attr(selection, "lower") else lower,
    upper = if (search == "golden") attr(selection, "upper") else upper,
    tol = if (search == "golden") tol else NA_real_,
    geometry = geometry,
    engine = engine,
    n_threads = n_threads,
    pilot_k = if (adaptive && bandwidth_type == "volume") pilot_k_resolved else pilot_k,
    tangent_ref = distance$tangent_ref,
    distance = distance$label,
    selection.table = selection,
    score = selection$score[best_id],
    AICc = if (approach == "AICc") selection$score[best_id] else NA_real_,
    CV = if (approach == "CV") selection$score[best_id] else NA_real_,
    trS = selection$trS[best_id],
    RSS = selection$RSS[best_id],
    this.call = this.call
  )
  class(result) <- "bw.lhwr"
  result
}

#' Deprecated LHWR bandwidth-selection alias
#'
#' `bw.hgwr()` is retained temporarily for older scripts. New code should use
#' [bw.lhwr()].
#'
#' @param ... Arguments passed to [bw.lhwr()].
#' @return The result of [bw.lhwr()].
#' @export
bw.hgwr <- function(...) {
  .Deprecated("bw.lhwr", msg = "bw.hgwr() is deprecated. Please use bw.lhwr() for Lorentz Hyperbolic Weighted Regression.")
  bw.lhwr(...)
}

resolve_bandwidth <- function(D, bw, adaptive, bandwidth_type, pilot_k = NULL) {
  n <- nrow(D)
  if (length(bw) == n) {
    check_bandwidth(bw)
    return(as.numeric(bw))
  }
  
  if (length(bw) != 1L) {
    stop("`bw` must be scalar or have one value per observation.", call. = FALSE)
  }
  
  if (adaptive) {
    target_mass <- as.integer(round(bw))
    if (target_mass < 2L || target_mass >= n) {
      stop("Adaptive bandwidth must be an integer between 2 and n - 1.", call. = FALSE)
    }
    if (bandwidth_type == "knn") {
      out <- kth_neighbor_radius(D, target_mass)
    } else {
      k0 <- resolve_pilot_k(n, pilot_k)
      pilot_radius <- kth_neighbor_radius(D, k0)
      local_density <- k0 / hyperbolic_ball_volume(pilot_radius)
      target_volume <- target_mass / local_density
      out <- hyperbolic_radius_from_volume(target_volume)
    }
    return(as.numeric(out))
  }
  
  check_bandwidth(bw)
  if (bandwidth_type == "volume") {
    rep(hyperbolic_radius_from_volume(bw), n)
  } else {
    rep(bw, n)
  }
}

resolve_bw_search <- function(search, bw_grid) {
  if (is.null(search)) {
    return(if (is.null(bw_grid)) "golden" else "grid")
  }
  match.arg(search, c("golden", "grid"))
}

score_bw_candidates <- function(y, X, D, bw, approach, kernel, adaptive,
                                bandwidth_type, pilot_k) {
  .bw_hgwr_grid_cpp(
    y = y,
    X = X,
    D = D,
    bw_grid = as.numeric(bw),
    approach = approach,
    kernel = kernel,
    adaptive = adaptive,
    bandwidth_type = bandwidth_type,
    pilot_k = pilot_k
  )
}

score_bw_candidates_knn <- function(y, X, Z, bw, approach, kernel,
                                    n_threads = 1L) {
  .bw_hgwr_knn_grid_cpp(
    y = y,
    X = X,
    Z = Z,
    bw_grid = as.numeric(bw),
    approach = approach,
    kernel = kernel,
    tol = 1e-6,
    n_threads = n_threads
  )
}

bw_hgwr_knn_golden_search <- function(y, X, Z, approach, kernel,
                                      lower = NULL, upper = NULL,
                                      tol = 1e-4, max_iter = 200L,
                                      n_threads = 1L) {
  bounds <- resolve_knn_bw_bounds(nrow(Z), lower, upper)
  if (!is.finite(bounds$lower) || !is.finite(bounds$upper) || bounds$lower >= bounds$upper) {
    stop("Invalid adaptive kNN bandwidth search interval.", call. = FALSE)
  }
  selection <- .bw_hgwr_knn_golden_cpp(
    y = y,
    X = X,
    Z = Z,
    approach = approach,
    kernel = kernel,
    lower = bounds$lower,
    upper = bounds$upper,
    tol = tol,
    max_iter = as.integer(max_iter),
    dist_tol = 1e-6,
    n_threads = n_threads
  )
  selection <- unique_selection_table(selection)
  attr(selection, "lower") <- bounds$lower
  attr(selection, "upper") <- bounds$upper
  selection
}

bw_hgwr_golden_search <- function(y, X, D, approach, kernel, adaptive,
                                  bandwidth_type, pilot_k, lower = NULL,
                                  upper = NULL, tol = 1e-4,
                                  max_iter = 200L) {
  bounds <- resolve_bw_bounds(D, adaptive, bandwidth_type, lower, upper)
  if (!is.finite(bounds$lower) || !is.finite(bounds$upper) || bounds$lower >= bounds$upper) {
    stop("Invalid bandwidth search interval.", call. = FALSE)
  }
  selection <- .bw_hgwr_golden_cpp(
    y = y,
    X = X,
    D = D,
    approach = approach,
    kernel = kernel,
    adaptive = adaptive,
    bandwidth_type = bandwidth_type,
    pilot_k = pilot_k,
    lower = bounds$lower,
    upper = bounds$upper,
    tol = tol,
    max_iter = as.integer(max_iter)
  )
  selection <- unique_selection_table(selection)
  attr(selection, "lower") <- bounds$lower
  attr(selection, "upper") <- bounds$upper
  selection
}

resolve_bw_bounds <- function(D, adaptive, bandwidth_type,
                              lower = NULL, upper = NULL) {
  n <- nrow(D)
  if (adaptive) {
    default_lower <- if (n > 21L) 20 else 2
    default_upper <- n - 1L
  } else {
    d <- D[D > 0 & is.finite(D)]
    if (!length(d)) {
      stop("Distance matrix does not contain positive finite distances.", call. = FALSE)
    }
    max_dist <- max(d)
    default_lower <- max_dist / 5000
    default_upper <- max_dist
    if (bandwidth_type == "volume") {
      default_lower <- hyperbolic_ball_volume(default_lower)
      default_upper <- hyperbolic_ball_volume(default_upper)
    }
  }
  
  lower <- if (is.null(lower)) default_lower else lower
  upper <- if (is.null(upper)) default_upper else upper
  if (adaptive) {
    lower <- as.integer(round(lower))
    upper <- as.integer(round(upper))
    lower <- max(2L, lower)
    upper <- min(n - 1L, upper)
  }
  list(lower = as.numeric(lower), upper = as.numeric(upper))
}

resolve_knn_bw_bounds <- function(n, lower = NULL, upper = NULL) {
  default_lower <- if (n > 21L) 20 else 2
  default_upper <- n - 1L
  lower <- if (is.null(lower)) default_lower else lower
  upper <- if (is.null(upper)) default_upper else upper
  lower <- max(2L, as.integer(round(lower)))
  upper <- min(n - 1L, as.integer(round(upper)))
  list(lower = as.numeric(lower), upper = as.numeric(upper))
}

unique_selection_table <- function(selection) {
  if (is.null(selection) || !nrow(selection)) return(selection)
  selection <- selection[order(selection$iteration), , drop = FALSE]
  selection <- selection[!duplicated(selection$bw), , drop = FALSE]
  selection[order(selection$bw), , drop = FALSE]
}

resolve_pilot_k <- function(n, pilot_k = NULL) {
  if (is.null(pilot_k)) {
    pilot_k <- max(5L, round(sqrt(n)))
  }
  pilot_k <- as.integer(round(pilot_k))
  if (pilot_k < 2L || pilot_k >= n) {
    stop("`pilot_k` must be an integer between 2 and n - 1.", call. = FALSE)
  }
  pilot_k
}

kth_neighbor_radius <- function(D, k) {
  apply(D, 2, function(x) sort(x[x > 0])[k])
}

default_bw_grid <- function(D, adaptive, bandwidth_type) {
  n <- nrow(D)
  if (adaptive) {
    unique(round(seq(max(2, n * 0.05), max(2, n * 0.95), length.out = 20)))
  } else {
    d <- D[D > 0]
    radii <- seq(stats::quantile(d, 0.05), stats::quantile(d, 0.8), length.out = 20)
    if (bandwidth_type == "volume") hyperbolic_ball_volume(radii) else radii
  }
}

default_knn_bw_grid <- function(n) {
  unique(round(seq(max(2, n * 0.05), max(2, n * 0.95), length.out = 20)))
}

validate_knn_engine <- function(geometry, adaptive, bandwidth_type, kernel,
                                dMat = NULL, search = "grid") {
  if (geometry != "lorentz") {
    stop("`engine = \"knn\"` currently requires `geometry = \"lorentz\"`.", call. = FALSE)
  }
  if (!adaptive || bandwidth_type != "knn") {
    stop("`engine = \"knn\"` currently requires `adaptive = TRUE` and `bandwidth_type = \"knn\"`.", call. = FALSE)
  }
  if (!kernel %in% c("bisquare", "tricube", "boxcar")) {
    stop("`engine = \"knn\"` currently supports compact kernels: \"bisquare\", \"tricube\", and \"boxcar\".", call. = FALSE)
  }
  if (!is.null(dMat)) {
    stop("`dMat` cannot be supplied with `engine = \"knn\"`; distances are computed from `Z`.", call. = FALSE)
  }
  if (!search %in% c("grid", "golden")) {
    stop("`engine = \"knn\"` supports `search = \"grid\"` or `search = \"golden\"`.", call. = FALSE)
  }
  invisible(TRUE)
}

hgwr_cv_score <- function(y, X, D, bw, kernel = hgwr_kernel_names()) {
  kernel <- match.arg(kernel)
  y <- as.numeric(y)
  X <- as.matrix(X)
  D <- as.matrix(D)
  .hgwr_cv_score_cpp(y, X, D, bw, kernel = kernel)
}

hyperbolic_ball_volume <- function(r) {
  2 * pi * (cosh(r) - 1)
}

hyperbolic_radius_from_volume <- function(volume) {
  acosh(1 + volume / (2 * pi))
}

#' @export
print.bw.lhwr <- function(x, ...) {
  if (!inherits(x, "bw.lhwr")) {
    stop("`x` must inherit from class 'bw.lhwr'.", call. = FALSE)
  }
  hgwr_rule()
  cat(" LHWR bandwidth selection\n")
  cat(" Automatic calibration of the kernel bandwidth\n")
  hgwr_rule()
  info <- data.frame(
    Setting = c(
      "Criterion",
      "Geometry",
      "Distance",
      "Kernel",
      "Search",
      "Engine",
      "Bandwidth type",
      "Bandwidth scale",
      "Selected bandwidth",
      "Selected score",
      "tr(S)",
      "RSS",
      "Evaluated candidates"
    ),
    Value = c(
      hgwr_criterion_label(x$approach),
      hgwr_geometry_label(x$geometry),
      x$distance,
      hgwr_kernel_label(x$kernel),
      hgwr_search_label(x$search),
      hgwr_engine_label(x$engine),
      hgwr_bandwidth_type_label(x$bandwidth_type),
      if (isTRUE(x$adaptive)) "Adaptive" else "Fixed",
      hgwr_bandwidth_label(list(bw = x$bw, adaptive = x$adaptive)),
      hgwr_format_number(x$score),
      hgwr_format_number(x$trS),
      hgwr_format_number(x$RSS),
      nrow(x$selection.table)
    ),
    check.names = FALSE
  )
  print(info, row.names = FALSE, right = FALSE)
  hgwr_rule()
  invisible(x)
}

#' @export
print.bw.hgwr <- print.bw.lhwr

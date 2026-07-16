#' Fit the Low-Level LHWR Model
#'
#' Fit local weighted least-squares regressions from a prepared response vector,
#' model matrix, distance matrix, and bandwidth vector. This function is the
#' computational engine behind [lhwr()] and assumes that all model inputs
#' have already been constructed.
#'
#' @param y Numeric response vector of length `n`.
#' @param X Numeric `n x p` model matrix. The first column should usually be an
#'   intercept if an intercept is desired.
#' @param D Numeric `n x n` pairwise distance matrix used to construct kernel
#'   weights for each focal location.
#' @param bw Positive bandwidth value or numeric vector of length `n`. A scalar
#'   is recycled to all focal locations; a vector supplies one radius per
#'   location.
#' @param kernel Kernel used to convert distances into local weights. Available
#'   choices are `"bisquare"`, `"gaussian"`, `"exponential"`, `"tricube"`, and
#'   `"boxcar"`.
#' @param return_hat Logical. If `TRUE`, compute hat-matrix diagnostics and
#'   unscaled local standard-error components.
#' @return A list with local coefficient matrix `beta_hat`, local standard
#'   errors `se_beta`, local t-values `t_beta`, fitted values `y_hat`,
#'   residuals, residual sum of squares `RSS`, hat-matrix traces `trS` and
#'   `trStS`, `AICc`, diagnostic summaries, effective bandwidth vector,
#'   hat-matrix diagonal, local condition numbers, local effective sample sizes,
#'   and local counts of positive kernel weights.
#' @examples
#' set.seed(2)
#' H <- simulate_H2(40, r_max = 2)
#' D <- hyperbolic_distance_matrix(H$Z)
#' X <- cbind(Intercept = 1, x = rnorm(40))
#' y <- 1 + 2 * X[, "x"] + rnorm(40, sd = 0.2)
#' bw <- rep(1.5, 40)
#'
#' fit <- hgwr.fit(y, X, D, bw, kernel = "gaussian")
#' fit$AICc
#' head(fit$beta_hat)
#' @export
hgwr.fit <- function(y, X, D, bw, kernel = hgwr_kernel_names(),
                     return_hat = TRUE) {
  kernel <- match.arg(kernel)
  y <- as.numeric(y)
  X <- as.matrix(X)
  D <- as.matrix(D)
  
  n <- length(y)
  p <- ncol(X)
  if (nrow(X) != n || nrow(D) != n || ncol(D) != n) {
    stop("Dimensions of `y`, `X`, and `D` are not compatible.", call. = FALSE)
  }
  if (length(bw) == 1L) bw <- rep(bw, n)
  if (length(bw) != n) {
    stop("`bw` must be scalar or have one value per observation.", call. = FALSE)
  }
  
  fit_cpp <- .hgwr_fit_cpp(y, X, D, bw, kernel = kernel, return_hat = return_hat)
  beta_hat <- fit_cpp$beta_hat
  se_beta <- fit_cpp$se_unscaled
  y_hat <- as.numeric(fit_cpp$y_hat)
  s_diag <- as.numeric(fit_cpp$s_diag)
  s2_diag <- as.numeric(fit_cpp$s2_diag)
  local_kappa <- as.numeric(fit_cpp$local_kappa)
  local_neff <- as.numeric(fit_cpp$local_neff)
  local_n_positive <- as.integer(fit_cpp$local_n_positive)
  
  colnames(beta_hat) <- colnames(X)
  colnames(se_beta) <- paste0(colnames(X), "_SE")
  residuals <- y - y_hat
  diagnostic <- if (isTRUE(return_hat)) {
    hgwr_diagnostics(
      y,
      X,
      beta_hat,
      y_hat,
      residuals,
      s_diag,
      s2_diag,
      local_kappa = local_kappa,
      local_neff = local_neff,
      local_n_positive = local_n_positive
    )
  } else {
    hgwr_minimal_diagnostics(y, y_hat, residuals)
  }
  se_beta <- if (isTRUE(return_hat)) sqrt(diagnostic$sigma2 * se_beta) else se_beta
  t_beta <- beta_hat / se_beta
  colnames(t_beta) <- paste0(colnames(X), "_TV")
  
  list(
    beta_hat = beta_hat,
    se_beta = se_beta,
    t_beta = t_beta,
    y_hat = y_hat,
    residuals = residuals,
    RSS = diagnostic$RSS.gw,
    trS = diagnostic$tr.S,
    trStS = diagnostic$tr.StS,
    AICc = if (!is.null(diagnostic$AICc)) diagnostic$AICc else NA_real_,
    diagnostic = diagnostic,
    bandwidth = bw,
    hat_diag = s_diag,
    local_kappa = local_kappa,
    local_neff = local_neff,
    local_n_positive = local_n_positive
  )
}

#' Fit the Low-Level LHWR Model from Lorentz kNN Neighborhoods
#'
#' Fit adaptive k-nearest-neighbor LHWR directly from Lorentz hyperboloid
#' coordinates without constructing a full pairwise distance matrix. This is a
#' low-level performance-oriented routine. For compact kernels such as
#' `"bisquare"`, `"tricube"`, and `"boxcar"`, it is algebraically equivalent to
#' [hgwr.fit()] with bandwidths equal to the Lorentz radius of the `k`th nearest
#' neighbor at each focal point.
#'
#' @param y Numeric response vector of length `n`.
#' @param X Numeric `n x p` model matrix.
#' @param Z Numeric `n x 3` Lorentz hyperboloid coordinate matrix.
#' @param k Adaptive nearest-neighbor count. The local fit uses the focal
#'   observation and its `k` nearest neighbors under Lorentz distance.
#' @param kernel Compactly supported kernel. Available exact choices are
#'   `"bisquare"`, `"tricube"`, and `"boxcar"`.
#' @param return_hat Logical. If `TRUE`, compute hat-matrix diagnostics and
#'   unscaled local standard-error components.
#' @param tol Numerical tolerance used in Lorentz distance clipping.
#' @param check Logical. If `TRUE`, validate that `Z` lies on the upper Lorentz
#'   hyperboloid before fitting.
#' @param n_threads Number of OpenMP threads used by the kNN local fitting loop
#'   when OpenMP is available. The default `1` gives deterministic serial
#'   execution and is used automatically on systems without OpenMP support.
#' @return A list in the same format as [hgwr.fit()], with an additional
#'   `bandwidth` vector containing the focal Lorentz kNN radii.
#' @examples
#' set.seed(7)
#' H <- simulate_H2(50, r_max = 3)
#' X <- cbind(Intercept = 1, x = rnorm(50))
#' y <- 1 + 2 * X[, "x"] + rnorm(50, sd = 0.2)
#' fit <- hgwr.fit.knn(y, X, H$Z, k = 20)
#' head(fit$beta_hat)
#' @export
hgwr.fit.knn <- function(y, X, Z, k, kernel = c("bisquare", "tricube", "boxcar"),
                         return_hat = TRUE, tol = 1e-6, check = TRUE,
                         n_threads = 1L) {
  kernel <- match.arg(kernel)
  y <- as.numeric(y)
  X <- as.matrix(X)
  Z <- as_lorentz_matrix(Z, "Z")
  if (check) check_hyperboloid(Z, tol = tol, stop_on_fail = TRUE)

  n <- length(y)
  if (nrow(X) != n || nrow(Z) != n) {
    stop("Dimensions of `y`, `X`, and `Z` are not compatible.", call. = FALSE)
  }
  k <- as.integer(round(k))
  if (length(k) != 1L || k < 2L || k >= n) {
    stop("`k` must be an integer between 2 and n - 1.", call. = FALSE)
  }
  n_threads <- resolve_n_threads_arg(n_threads)

  fit_cpp <- .hgwr_fit_knn_cpp(
    y = y,
    X = X,
    Z = Z,
    k = k,
    kernel = kernel,
    return_hat = return_hat,
    tol = tol,
    n_threads = n_threads
  )
  beta_hat <- fit_cpp$beta_hat
  se_beta <- fit_cpp$se_unscaled
  y_hat <- as.numeric(fit_cpp$y_hat)
  s_diag <- as.numeric(fit_cpp$s_diag)
  s2_diag <- as.numeric(fit_cpp$s2_diag)
  local_kappa <- as.numeric(fit_cpp$local_kappa)
  local_neff <- as.numeric(fit_cpp$local_neff)
  local_n_positive <- as.integer(fit_cpp$local_n_positive)
  bw <- as.numeric(fit_cpp$bandwidth)
  knn_index <- fit_cpp$knn_index
  knn_distance <- fit_cpp$knn_distance

  colnames(beta_hat) <- colnames(X)
  colnames(se_beta) <- paste0(colnames(X), "_SE")
  residuals <- y - y_hat
  diagnostic <- if (isTRUE(return_hat)) {
    hgwr_diagnostics(
      y,
      X,
      beta_hat,
      y_hat,
      residuals,
      s_diag,
      s2_diag,
      local_kappa = local_kappa,
      local_neff = local_neff,
      local_n_positive = local_n_positive
    )
  } else {
    hgwr_minimal_diagnostics(y, y_hat, residuals)
  }
  se_beta <- if (isTRUE(return_hat)) sqrt(diagnostic$sigma2 * se_beta) else se_beta
  t_beta <- beta_hat / se_beta
  colnames(t_beta) <- paste0(colnames(X), "_TV")

  list(
    beta_hat = beta_hat,
    se_beta = se_beta,
    t_beta = t_beta,
    y_hat = y_hat,
    residuals = residuals,
    RSS = diagnostic$RSS.gw,
    trS = diagnostic$tr.S,
    trStS = diagnostic$tr.StS,
    AICc = if (!is.null(diagnostic$AICc)) diagnostic$AICc else NA_real_,
    diagnostic = diagnostic,
    bandwidth = bw,
    hat_diag = s_diag,
    local_kappa = local_kappa,
    local_neff = local_neff,
    local_n_positive = local_n_positive,
    knn_index = knn_index,
    knn_distance = knn_distance
  )
}

#' Fit a Lorentz Hyperbolic Weighted Regression Model
#'
#' Estimate location-varying regression coefficients for observations located
#' on the Lorentz hyperboloid. At each focal point, the function fits a weighted
#' linear regression using kernel weights derived from either intrinsic Lorentz
#' distances or tangent-space distances.
#'
#' @param formula Model formula of the form `response ~ predictors`.
#' @param data Data frame containing the response and predictor variables used
#'   in `formula`.
#' @param Z Numeric `n x 3` matrix of Lorentz hyperboloid coordinates. Row `i`
#'   gives the location of observation `i` and must satisfy
#'   `<z_i,z_i>_L = -1` with positive time coordinate.
#' @param bw Bandwidth used for local weighting. If `adaptive = TRUE`, a scalar
#'   is interpreted as a neighbor count and converted to a separate radius at
#'   each focal location. If `adaptive = FALSE`, a scalar is interpreted as a
#'   fixed radius, or as a fixed hyperbolic volume when
#'   `bandwidth_type = "volume"`. A numeric vector of length `n` may be used to
#'   supply one radius per location directly.
#' @param kernel Kernel used to convert distances into local weights. Available
#'   choices are `"bisquare"`, `"gaussian"`, `"exponential"`, `"tricube"`, and
#'   `"boxcar"`.
#' @param adaptive Logical. If `TRUE`, scalar `bw` is treated as a neighbor
#'   count. If `FALSE`, scalar `bw` is treated as a fixed distance or volume.
#' @param bandwidth_type Bandwidth scale. `"knn"` uses nearest-neighbor radii.
#'   `"volume"` uses hyperbolic ball volumes and is available only for
#'   `geometry = "lorentz"`.
#' @param geometry Distance geometry used to build kernel weights. `"lorentz"`
#'   uses intrinsic hyperbolic distances on the Lorentz hyperboloid.
#'   `"tangent"` maps locations to one tangent plane by the Lorentz logarithmic
#'   map and uses the positive-definite Lorentz metric induced on that tangent
#'   plane. Equivalently, these distances can be represented as Euclidean
#'   distances in a Lorentz-orthonormal tangent basis.
#' @param pilot_k Pilot neighbor count used to estimate local point density for
#'   adaptive volume bandwidths. Ignored unless `adaptive = TRUE` and
#'   `bandwidth_type = "volume"`.
#' @param tangent_ref Optional reference point for `geometry = "tangent"`.
#'   If omitted, the normalized extrinsic Lorentz mean is used.
#' @param dMat Optional precomputed `n x n` distance matrix. If supplied, it is
#'   used directly instead of computing a distance matrix from `Z`.
#' @param engine Computational engine. `"dense"` builds or uses a full distance
#'   matrix. `"knn"` uses exact Lorentz nearest-neighbor lists and avoids storing
#'   the full distance matrix; it is available for adaptive Lorentz kNN fits
#'   with compact kernels.
#' @param n_threads Number of OpenMP threads used by `engine = "knn"` when
#'   OpenMP is available. The default `1` preserves serial execution.
#' @param diagnostics Diagnostic detail level. `"full"` returns the traditional
#'   diagnostic table, local standard errors, t-values, and local R-squared
#'   values when available. `"aicc"` computes hat-matrix diagnostics for AICc
#'   but skips local standard-error and local R-squared columns in the output
#'   table. `"none"` skips hat-matrix diagnostics and returns only fitted
#'   values, residuals, local coefficients, RSS, and global R-squared.
#' @param f_tests Logical. If `TRUE`, compute LHWR diagnostic F-tests
#'   after model fitting. The default `FALSE` avoids the additional smoother
#'   reconstruction cost, which can be substantial for large data sets.
#' @return An object of class `"lhwr"`, a list containing model arguments,
#'   diagnostic summaries, the global linear model fit, local coefficient table
#'   `SDF`, coefficient names, local coefficient matrix `beta_hat`, fitted
#'   values, residuals, distance matrix, selected geometry label, effective
#'   bandwidth vector, hat-matrix diagonal, local condition numbers, local
#'   effective sample sizes, local counts of positive kernel weights, timing
#'   information, and the original function call.
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
#'   adaptive = TRUE,
#'   geometry = "lorentz",
#'   bw_grid = seq(10, 80, by = 5)
#' )
#'
#' fit <- lhwr(
#'   y ~ x1 + x2,
#'   data = dat,
#'   Z = Z,
#'   bw = bw$bw,
#'   adaptive = TRUE,
#'   geometry = "lorentz"
#' )
#' fit
#' summary(fit)
#' head(fit$SDF)
#'
#' beta_true <- cbind(Intercept = beta0, x1 = beta1, x2 = beta2)
#' sqrt(mean(rowSums((fit$beta_hat - beta_true)^2)))
#'
#' fit_tangent <- lhwr(
#'   y ~ x1 + x2,
#'   data = dat,
#'   Z = Z,
#'   bw = bw$bw,
#'   adaptive = TRUE,
#'   geometry = "tangent"
#' )
#' fit_tangent$geometry
#' @export
lhwr <- function(formula, data, Z, bw,
                       kernel = hgwr_kernel_names(),
                       adaptive = TRUE,
                       bandwidth_type = c("knn", "volume"),
                       geometry = c("lorentz", "tangent"),
                       pilot_k = NULL,
                       tangent_ref = NULL,
                       dMat = NULL,
                       engine = c("dense", "knn"),
                       n_threads = 1L,
                       diagnostics = c("full", "aicc", "none"),
                       f_tests = FALSE) {
  this.call <- match.call()
  timings <- list(start = Sys.time())
  kernel <- match.arg(kernel)
  bandwidth_type <- match.arg(bandwidth_type)
  geometry <- match.arg(geometry)
  engine <- match.arg(engine)
  diagnostics <- match.arg(diagnostics)
  f_tests <- isTRUE(f_tests)
  n_threads <- resolve_n_threads_arg(n_threads)
  if (geometry == "tangent" && bandwidth_type == "volume") {
    stop("`bandwidth_type = \"volume\"` is only defined for `geometry = \"lorentz\"`.", call. = FALSE)
  }
  if (engine == "knn") {
    validate_knn_engine(
      geometry = geometry,
      adaptive = adaptive,
      bandwidth_type = bandwidth_type,
      kernel = kernel,
      dMat = dMat,
      search = "grid"
    )
  }
  
  mf <- hgwr_model_frame(formula, data)
  y <- mf$y
  X <- mf$X
  colnames(X) <- clean_coef_names(colnames(X))
  
  Z <- as_lorentz_matrix(Z, "Z")
  check_hyperboloid(Z, stop_on_fail = TRUE)
  if (nrow(Z) != length(y)) {
    stop("`Z` must have the same number of rows as `data`.", call. = FALSE)
  }
  
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
  bw_vec <- if (engine == "knn") {
    if (length(bw) != 1L) {
      stop("`engine = \"knn\"` requires scalar adaptive neighbor-count `bw`.", call. = FALSE)
    }
    as.integer(round(bw))
  } else {
    resolve_bandwidth(D, bw, adaptive, bandwidth_type, pilot_k = pilot_k)
  }
  
  timings$calibration <- Sys.time()
  return_hat <- diagnostics != "none" || f_tests
  fit <- if (engine == "knn") {
    hgwr.fit.knn(
      y = y,
      X = X,
      Z = Z,
      k = bw_vec,
      kernel = kernel,
      return_hat = return_hat,
      check = FALSE,
      n_threads = n_threads
    )
  } else {
    hgwr.fit(y = y, X = X, D = D, bw = bw_vec, kernel = kernel, return_hat = return_hat)
  }
  timings$diagnostic <- Sys.time()
  
  beta_df <- as.data.frame(fit$beta_hat)
  SDF <- data.frame(
    y = y,
    yhat = fit$y_hat,
    residual = fit$residuals,
    beta_df,
    check.names = FALSE
  )
  if (diagnostics == "full") {
    SDF <- data.frame(
      SDF,
      as.data.frame(fit$se_beta),
      as.data.frame(fit$t_beta),
      Local_R2 = if (engine == "knn") {
        local_r2_knn(y, fit$y_hat, fit$knn_index, fit$knn_distance, fit$bandwidth, kernel)
      } else {
        local_r2(y, fit$y_hat, D, bw_vec, kernel)
      },
      check.names = FALSE
    )
  }
  
  lm_fit <- stats::lm(formula, data = data)
  result <- list(
    LHWR.arguments = list(
      formula = formula,
      bw = bw,
      kernel = kernel,
      adaptive = adaptive,
      bandwidth_type = bandwidth_type,
      geometry = geometry,
      engine = engine,
      n_threads = n_threads,
      diagnostics = diagnostics,
      f_tests = f_tests,
      pilot_k = if (adaptive && bandwidth_type == "volume") resolve_pilot_k(length(y), pilot_k) else pilot_k,
      tangent_ref = distance$tangent_ref,
      distance = distance$label
    ),
    LHWR.diagnostic = fit$diagnostic,
    lm = lm_fit,
    SDF = SDF,
    coef_names = colnames(fit$beta_hat),
    beta_hat = fit$beta_hat,
    y_hat = fit$y_hat,
    residuals = fit$residuals,
    Z = Z,
    distance = D,
    distance_label = distance$label,
    geometry = geometry,
    tangent_ref = distance$tangent_ref,
    bandwidth = if (engine == "knn") fit$bandwidth else bw_vec,
    hat_diag = fit$hat_diag,
    local_kappa = fit$local_kappa,
    local_neff = fit$local_neff,
    local_n_positive = fit$local_n_positive,
    knn_index = fit$knn_index,
    knn_distance = fit$knn_distance,
    timings = c(timings, list(stop = Sys.time())),
    this.call = this.call
  )
  
  result$HGWR.arguments <- result$LHWR.arguments
  result$HGWR.diagnostic <- result$LHWR.diagnostic
  class(result) <- "lhwr"
  if (f_tests) {
    result$Ftests <- diagnostic_tests(result)
  }
  result
}

#' Deprecated HGWR model-fitting alias
#'
#' `hgwr.basic()` is retained temporarily for older scripts. New code should use
#' [lhwr()].
#'
#' @param ... Arguments passed to [lhwr()].
#' @return The result of [lhwr()].
#' @export
hgwr.basic <- function(...) {
  .Deprecated("lhwr", msg = "hgwr.basic() is deprecated. Please use lhwr() for Lorentz Hyperbolic Weighted Regression.")
  lhwr(...)
}

hgwr_wls_one <- function(X, y, w, focus, return_hat = TRUE) {
  Xw <- X * w
  XtWX <- crossprod(X, Xw)
  XtWy <- crossprod(X, w * y)
  XtWX_inv <- tryCatch(solve(XtWX), error = function(e) MASS::ginv(XtWX))
  beta <- as.numeric(XtWX_inv %*% XtWy)
  n_positive <- sum(w > 0)
  neff <- sum(w)^2 / sum(w^2)
  neff <- if (is.finite(neff)) neff else 0
  kappa <- tryCatch(kappa(XtWX, exact = TRUE), error = function(e) Inf)
  
  se_unscaled <- rep(NA_real_, ncol(X))
  s_diag <- NA_real_
  s2_diag <- NA_real_
  
  if (return_hat) {
    C <- XtWX_inv %*% t(Xw)
    s_row <- as.numeric(X[focus, ] %*% C)
    s_diag <- s_row[focus]
    s2_diag <- sum(s_row^2)
    se_unscaled <- rowSums(C^2)
  }
  
  list(
    beta = beta,
    se_unscaled = se_unscaled,
    s_diag = s_diag,
    s2_diag = s2_diag,
    kappa = kappa,
    neff = neff,
    n_positive = n_positive
  )
}

hgwr_model_frame <- function(formula, data) {
  mf <- stats::model.frame(formula, data)
  y <- stats::model.response(mf)
  X <- stats::model.matrix(stats::terms(mf), mf)
  list(y = y, X = X)
}

resolve_n_threads_arg <- function(n_threads) {
  if (length(n_threads) != 1L || !is.finite(n_threads)) {
    stop("`n_threads` must be a positive integer.", call. = FALSE)
  }
  n_threads <- as.integer(round(n_threads))
  if (n_threads < 1L) {
    stop("`n_threads` must be a positive integer.", call. = FALSE)
  }
  n_threads
}

clean_coef_names <- function(x) {
  x[x == "(Intercept)"] <- "Intercept"
  x
}

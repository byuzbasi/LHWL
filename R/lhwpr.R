#' LHWPR Solver Controls
#'
#' Configure the numerical solver used by [lhwpr()].
#'
#' @param solver Coordinate-descent implementation. `"auto"` selects between
#'   residual-updating (`"naive"`) and Gaussian covariance-updating solvers
#'   from the problem dimensions.
#' @param screening Sequential strong screening for lasso and elastic net, or
#'   no screening. Every screened fit is followed by a full KKT repair scan.
#' @param tol Coordinate-descent convergence tolerance.
#' @param kkt_tol Full KKT/stationarity tolerance.
#' @param max_iter Maximum coordinate-descent sweeps per weighted least-squares
#'   problem.
#' @param irls_tol Binomial outer-loop tolerance.
#' @param irls_max_iter Maximum binomial proximal-Newton/IRLS iterations.
#' @param max_active Optional maximum active-set size. `NULL` imposes no limit.
#' @param early_stop Reserved for a later globally synchronized path solver;
#'   it must remain `FALSE` in this release.
#' @param ridge_path_tol Positive tolerance defining the ridge path's upper
#'   scale.
#' @param n_threads Number of focal-location threads.
#' @param trace Logical; retain verbose numerical tracing when available.
#' @return A validated control list.
#' @export
lhwpr_control <- function(
    solver = c("auto", "naive", "covariance"),
    screening = c("strong", "none"),
    tol = 1e-7,
    kkt_tol = 1e-7,
    max_iter = 100000L,
    irls_tol = 1e-7,
    irls_max_iter = 100L,
    max_active = NULL,
    early_stop = FALSE,
    ridge_path_tol = 1e-6,
    n_threads = 1L,
    trace = FALSE) {
  solver <- match.arg(solver)
  screening <- match.arg(screening)
  tol <- lhwpr_positive_scalar(tol, "tol")
  kkt_tol <- lhwpr_positive_scalar(kkt_tol, "kkt_tol")
  max_iter <- lhwpr_positive_integer(max_iter, "max_iter")
  irls_tol <- lhwpr_positive_scalar(irls_tol, "irls_tol")
  irls_max_iter <- lhwpr_positive_integer(irls_max_iter, "irls_max_iter")
  ridge_path_tol <- lhwpr_positive_scalar(ridge_path_tol, "ridge_path_tol")
  if (ridge_path_tol >= 1) {
    stop("`ridge_path_tol` must be smaller than one.", call. = FALSE)
  }
  n_threads <- resolve_n_threads_arg(n_threads)
  if (is.null(max_active)) {
    max_active <- -1L
  } else {
    max_active <- lhwpr_positive_integer(max_active, "max_active")
  }
  early_stop <- lhwpr_logical_scalar(early_stop, "early_stop")
  if (early_stop) {
    stop(
      "`early_stop = TRUE` is reserved for a later globally synchronized path solver; use `FALSE` in this release.",
      call. = FALSE
    )
  }
  trace <- lhwpr_logical_scalar(trace, "trace")
  structure(
    list(
      solver = solver,
      screening = screening,
      tol = tol,
      kkt_tol = kkt_tol,
      max_iter = max_iter,
      irls_tol = irls_tol,
      irls_max_iter = irls_max_iter,
      max_active = max_active,
      early_stop = early_stop,
      ridge_path_tol = ridge_path_tol,
      n_threads = n_threads,
      trace = trace
    ),
    class = "lhwpr_control"
  )
}

#' Lorentz-Hyperbolic Weighted Penalized Regression
#'
#' Fit local Gaussian or binomial regressions whose observation weights are
#' defined by Lorentz-hyperbolic kernels. The penalized numerical path is
#' solved in RcppArmadillo with warm starts and exact post-screening checks.
#'
#' @param formula,data Model formula and data frame.
#' @param Z Numeric `n x 3` Lorentz-hyperboloid coordinates.
#' @param bw Kernel bandwidth. With adaptive kNN bandwidths this is a neighbor
#'   count.
#' @param family `"gaussian"` or `"binomial"`.
#' @param penalty `"lasso"`, `"mcp"`, or `"scad"`. At `alpha = 0` this
#'   argument is accepted but ignored.
#' @param alpha Mixing value in `[0, 1]`: zero is ridge, values strictly
#'   between zero and one are elastic net with `penalty = "lasso"`, and one is
#'   lasso, MCP, or SCAD.
#' @param lambda Optional common absolute lambda path in decreasing order.
#' @param nlambda Number of automatically generated lambda values.
#' @param lambda_min_ratio Ratio between the last and first automatic lambda.
#' @param gamma MCP or SCAD concavity parameter.
#' @param penalty_factor Optional non-negative multiplier for each non-intercept
#'   coefficient.
#' @param case_weights Optional non-negative analysis weights.
#' @param standardize Logical; globally center/scale predictors once.
#' @param kernel Kernel name.
#' @param adaptive Logical; use adaptive bandwidths.
#' @param bandwidth_type `"knn"` or hyperbolic `"volume"`.
#' @param geometry `"lorentz"` or `"tangent"`.
#' @param tangent_ref Optional tangent-space reference point.
#' @param dMat Optional precomputed distance matrix.
#' @param engine Dense distance matrix or memory-saving exact kNN engine.
#' @param keep_path Store all coefficient paths, only the selected fit, or
#'   choose automatically with a ten-million-value dense-path memory guard.
#' @param control A value returned by [lhwpr_control()].
#' @return An object of class `"lhwpr"`.
#' @export
lhwpr <- function(
    formula,
    data,
    Z,
    bw,
    family = c("gaussian", "binomial"),
    penalty = c("lasso", "mcp", "scad"),
    alpha = 1,
    lambda = NULL,
    nlambda = 100L,
    lambda_min_ratio = NULL,
    gamma = NULL,
    penalty_factor = NULL,
    case_weights = NULL,
    standardize = TRUE,
    kernel = hgwr_kernel_names(),
    adaptive = TRUE,
    bandwidth_type = c("knn", "volume"),
    geometry = c("lorentz", "tangent"),
    tangent_ref = NULL,
    dMat = NULL,
    engine = c("dense", "knn"),
    keep_path = c("auto", "all", "selected"),
    control = lhwpr_control()) {
  this_call <- match.call()
  family <- match.arg(family)
  penalty <- match.arg(penalty)
  kernel <- match.arg(kernel)
  bandwidth_type <- match.arg(bandwidth_type)
  geometry <- match.arg(geometry)
  engine <- match.arg(engine)
  keep_path <- match.arg(keep_path)
  alpha <- lhwpr_alpha(alpha)
  lhwpr_validate_penalty_alpha(penalty, alpha)
  gamma <- lhwpr_gamma(penalty, alpha, gamma)
  control <- lhwpr_validate_control(control)
  lhwpr_validate_solver_family(control, family)
  nlambda <- lhwpr_positive_integer(nlambda, "nlambda")
  lambda <- lhwpr_validate_lambda(lambda)

  prepared <- lhwpr_prepare(
    formula = formula,
    data = data,
    Z = Z,
    family = family,
    penalty_factor = penalty_factor,
    case_weights = case_weights,
    standardize = standardize,
    geometry = geometry,
    tangent_ref = tangent_ref,
    dMat = dMat,
    engine = engine
  )
  if (is.null(lambda_min_ratio)) {
    lambda_min_ratio <- lhwpr_default_lambda_ratio(
      alpha, nrow(prepared$X), ncol(prepared$X), control$ridge_path_tol
    )
  }
  lambda_min_ratio <- lhwpr_ratio(lambda_min_ratio)
  neighborhood <- lhwpr_neighborhood(
    prepared = prepared,
    bw = bw,
    kernel = kernel,
    adaptive = adaptive,
    bandwidth_type = bandwidth_type,
    engine = engine
  )
  path_length <- if (is.null(lambda)) nlambda else length(lambda)
  path_elements <- nrow(prepared$X) *
    (ncol(prepared$X) + as.integer(prepared$intercept)) *
    path_length
  store_path <- switch(
    keep_path,
    all = TRUE,
    selected = FALSE,
    auto = path_length > 1L && path_elements <= 1e7
  )
  if (keep_path == "auto" && path_length > 1L && !store_path) {
    warning(
      "The dense coefficient path exceeds the automatic retention limit; ",
      "selected coefficients and the fitted-value path are retained. ",
      "Use `keep_path = \"all\"` to request the full coefficient cube.",
      call. = FALSE
    )
  }
  native <- lhwpr_native_fit(
    prepared = prepared,
    neighborhood = neighborhood,
    lambda = lambda,
    family = family,
    penalty = penalty,
    alpha = alpha,
    gamma = gamma,
    nlambda = nlambda,
    lambda_min_ratio = lambda_min_ratio,
    control = control,
    store_path = store_path,
    exclude_self = FALSE
  )
  out <- lhwpr_assemble(
    native = native,
    prepared = prepared,
    neighborhood = neighborhood,
    formula = formula,
    family = family,
    penalty = penalty,
    alpha = alpha,
    gamma = gamma,
    kernel = kernel,
    adaptive = adaptive,
    bandwidth_type = bandwidth_type,
    geometry = geometry,
    engine = engine,
    standardize = standardize,
    keep_path = keep_path,
    control = control,
    call = this_call,
    bw = bw
  )
  class(out) <- "lhwpr"
  out
}

lhwpr_prepare <- function(formula, data, Z, family, penalty_factor,
                          case_weights, standardize, geometry, tangent_ref,
                          dMat, engine) {
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a model formula.", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }
  mf <- stats::model.frame(
    formula, data = data, na.action = stats::na.fail
  )
  terms_object <- stats::terms(mf)
  if (length(stats::model.offset(mf))) {
    stop("Formula offsets are not yet supported by `lhwpr()`.", call. = FALSE)
  }
  response <- stats::model.response(mf)
  response_info <- lhwpr_response(response, family)
  y <- response_info$y
  X_full <- stats::model.matrix(terms_object, mf)
  intercept <- isTRUE(attr(terms_object, "intercept") == 1L)
  intercept_column <- match("(Intercept)", colnames(X_full), nomatch = 0L)
  if (intercept && intercept_column > 0L) {
    X <- X_full[, -intercept_column, drop = FALSE]
  } else {
    X <- X_full
  }
  if (!ncol(X)) {
    stop("`lhwpr()` requires at least one penalized predictor.", call. = FALSE)
  }
  storage.mode(X) <- "double"
  if (any(!is.finite(X))) {
    stop("The model matrix must contain only finite values.", call. = FALSE)
  }
  n <- nrow(X)
  case_weights <- lhwpr_case_weights(case_weights, n)
  penalty_factor <- lhwpr_penalty_factor(
    penalty_factor, ncol(X), colnames(X)
  )
  standardized <- lhwpr_standardize(
    X, case_weights, intercept = intercept,
    standardize = lhwpr_logical_scalar(standardize, "standardize")
  )

  Z <- as_lorentz_matrix(Z, "Z")
  check_hyperboloid(Z, stop_on_fail = TRUE)
  if (nrow(Z) != n) {
    stop(
      "`Z` must have the same number of rows as the analysis data.",
      call. = FALSE
    )
  }
  if (geometry == "tangent" && engine == "knn") {
    stop(
      "`engine = \"knn\"` currently requires `geometry = \"lorentz\"`.",
      call. = FALSE
    )
  }
  distance <- if (engine == "dense") {
    resolve_geometry_distance(
      Z,
      dMat = dMat,
      geometry = geometry,
      tangent_ref = tangent_ref
    )
  } else {
    if (!is.null(dMat)) {
      stop(
        "`dMat` cannot be supplied with `engine = \"knn\"`.",
        call. = FALSE
      )
    }
    list(
      D = NULL,
      label = "Lorentz hyperboloid",
      tangent_ref = NULL
    )
  }
  if (!is.null(distance$D) &&
      !identical(dim(distance$D), c(n, n))) {
    stop("`dMat` must be an `n x n` distance matrix.", call. = FALSE)
  }

  list(
    y = y,
    response_levels = response_info$levels,
    X = X,
    X_standardized = standardized$X,
    center = standardized$center,
    scale = standardized$scale,
    predictor_names = colnames(X),
    intercept = intercept,
    penalty_factor = penalty_factor,
    case_weights = case_weights,
    terms = terms_object,
    contrasts = attr(X_full, "contrasts"),
    xlevels = stats::.getXlevels(terms_object, mf),
    Z = Z,
    D = distance$D,
    distance_label = distance$label,
    tangent_ref = distance$tangent_ref,
    geometry = geometry,
    engine = engine
  )
}

lhwpr_neighborhood <- function(prepared, bw, kernel, adaptive,
                               bandwidth_type, engine) {
  adaptive <- lhwpr_logical_scalar(adaptive, "adaptive")
  if (adaptive && length(bw) == 1L &&
      (!is.numeric(bw) || !is.finite(bw) ||
       abs(bw - round(bw)) > sqrt(.Machine$double.eps))) {
    stop(
      "A scalar adaptive `bw` must be an integer neighbor count.",
      call. = FALSE
    )
  }
  if (prepared$geometry == "tangent" && bandwidth_type == "volume") {
    stop(
      "`bandwidth_type = \"volume\"` is only defined for Lorentz geometry.",
      call. = FALSE
    )
  }
  if (engine == "knn") {
    validate_knn_engine(
      geometry = prepared$geometry,
      adaptive = adaptive,
      bandwidth_type = bandwidth_type,
      kernel = kernel,
      dMat = NULL,
      search = "grid"
    )
    if (length(bw) != 1L || !is.finite(bw)) {
      stop(
        "`engine = \"knn\"` requires one finite neighbor count.",
        call. = FALSE
      )
    }
    k <- as.integer(round(bw))
    if (k < 2L || k >= nrow(prepared$X)) {
      stop(
        "Adaptive bandwidth must be an integer between 2 and n - 1.",
        call. = FALSE
      )
    }
    knn <- lorentz_knn(
      prepared$Z, k = k, include_self = FALSE, check = FALSE
    )
    return(list(
      engine = "knn",
      k = k,
      bandwidth = knn$distance[, k],
      index = knn$index,
      distance = knn$distance,
      kernel = kernel
    ))
  }
  bw_vector <- resolve_bandwidth(
    prepared$D,
    bw,
    adaptive = adaptive,
    bandwidth_type = bandwidth_type,
    pilot_k = NULL
  )
  list(
    engine = "dense",
    bandwidth = as.numeric(bw_vector),
    D = prepared$D,
    kernel = kernel
  )
}

lhwpr_native_options <- function(prepared, family, penalty, alpha, gamma,
                                 nlambda, lambda_min_ratio, control,
                                 store_path, exclude_self,
                                 selected_index = -1L, kernel,
                                 return_path_diagnostics = FALSE,
                                 return_active_frequency = FALSE,
                                 compute_lambda_diagnostics = TRUE) {
  list(
    family = family,
    penalty = penalty,
    alpha = alpha,
    gamma = gamma,
    intercept = prepared$intercept,
    solver = control$solver,
    screening = control$screening,
    tol = control$tol,
    kkt_tol = control$kkt_tol,
    max_iter = control$max_iter,
    irls_tol = control$irls_tol,
    irls_max_iter = control$irls_max_iter,
    max_active = control$max_active,
    early_stop = control$early_stop,
    ridge_path_tol = control$ridge_path_tol,
    n_threads = control$n_threads,
    trace = control$trace,
    store_path = store_path,
    return_inference = FALSE,
    return_path_diagnostics = return_path_diagnostics,
    return_active_frequency = return_active_frequency,
    compute_lambda_diagnostics = compute_lambda_diagnostics,
    selected_index = as.integer(selected_index),
    exclude_self = exclude_self,
    nlambda = as.integer(nlambda),
    lambda_min_ratio = lambda_min_ratio,
    kernel = kernel
  )
}

lhwpr_native_fit <- function(prepared, neighborhood, lambda, family,
                             penalty, alpha, gamma, nlambda,
                             lambda_min_ratio, control, store_path,
                             exclude_self, case_weights = prepared$case_weights,
                             X_focus = prepared$X_standardized,
                             return_path_diagnostics = FALSE,
                             return_active_frequency = FALSE,
                             compute_lambda_diagnostics = TRUE) {
  options <- lhwpr_native_options(
    prepared = prepared,
    family = family,
    penalty = penalty,
    alpha = alpha,
    gamma = gamma,
    nlambda = nlambda,
    lambda_min_ratio = lambda_min_ratio,
    control = control,
    store_path = store_path,
    exclude_self = exclude_self,
    selected_index = -1L,
    kernel = if (is.null(neighborhood$kernel)) {
      attr(neighborhood, "kernel")
    } else {
      neighborhood$kernel
    },
    return_path_diagnostics = return_path_diagnostics,
    return_active_frequency = return_active_frequency,
    compute_lambda_diagnostics = compute_lambda_diagnostics
  )
  if (is.null(options$kernel)) {
    stop("Internal LHWPR neighborhood is missing its kernel.", call. = FALSE)
  }
  lambda_native <- if (is.null(lambda)) numeric() else as.numeric(lambda)
  if (control$trace) {
    message(
      "LHWPR native fit: family=", family,
      ", penalty=", if (alpha == 0) "ridge" else penalty,
      ", locations=", nrow(X_focus),
      ", engine=", neighborhood$engine
    )
  }
  fit <- if (neighborhood$engine == "dense") {
    .lhwpr_fit_dense_cpp(
      y = prepared$y,
      X = prepared$X_standardized,
      X_focus = X_focus,
      D = neighborhood$D,
      bw = neighborhood$bandwidth,
      lambda = lambda_native,
      penalty_factor = prepared$penalty_factor,
      case_weights = case_weights,
      options = options
    )
  } else {
    .lhwpr_fit_knn_cpp(
      y = prepared$y,
      X = prepared$X_standardized,
      X_focus = X_focus,
      knn_index = neighborhood$index,
      knn_distance = neighborhood$distance,
      lambda = lambda_native,
      penalty_factor = prepared$penalty_factor,
      case_weights = case_weights,
      options = options
    )
  }
  if (control$trace) {
    message(
      "LHWPR native fit complete: ",
      sum(fit$converged != 0),
      "/",
      length(fit$converged),
      " focal-path solutions converged"
    )
  }
  fit
}

lhwpr_native_lambda_max <- function(prepared, neighborhood, family, penalty,
                                    alpha, gamma, nlambda,
                                    lambda_min_ratio, control, exclude_self,
                                    case_weights = prepared$case_weights) {
  options <- lhwpr_native_options(
    prepared = prepared,
    family = family,
    penalty = penalty,
    alpha = alpha,
    gamma = gamma,
    nlambda = nlambda,
    lambda_min_ratio = lambda_min_ratio,
    control = control,
    store_path = FALSE,
    exclude_self = exclude_self,
    selected_index = -1L,
    kernel = neighborhood$kernel
  )
  if (neighborhood$engine == "dense") {
    return(.lhwpr_lambda_max_dense_cpp(
      prepared$y,
      prepared$X_standardized,
      neighborhood$D,
      neighborhood$bandwidth,
      prepared$penalty_factor,
      case_weights,
      options
    ))
  }
  .lhwpr_lambda_max_knn_cpp(
    prepared$y,
    prepared$X_standardized,
    neighborhood$index,
    neighborhood$distance,
    prepared$penalty_factor,
    case_weights,
    options
  )
}

lhwpr_assemble <- function(native, prepared, neighborhood, formula, family,
                           penalty, alpha, gamma, kernel, adaptive,
                           bandwidth_type, geometry, engine, standardize,
                           keep_path, control, call, bw) {
  beta_hat <- lhwpr_unstandardize_matrix(
    native$beta_hat_standardized,
    prepared$center,
    prepared$scale,
    prepared$intercept,
    prepared$predictor_names
  )
  beta_path <- if (is.null(native$beta_path_standardized)) {
    NULL
  } else {
    lhwpr_unstandardize_path(
      native$beta_path_standardized,
      prepared$center,
      prepared$scale,
      prepared$intercept,
      prepared$predictor_names
    )
  }
  selected <- native$selected_index
  fitted_values <- native$fitted_path[, selected]
  eta <- native$linear_predictor_path[, selected]
  residual_values <- prepared$y - fitted_values
  coefficient_frame <- as.data.frame(beta_hat, check.names = FALSE)
  SDF <- data.frame(
    y = prepared$y,
    yhat = fitted_values,
    residual = residual_values,
    coefficient_frame,
    check.names = FALSE
  )
  effective_penalty <- if (alpha == 0) {
    "ridge"
  } else if (alpha < 1) {
    "elastic net"
  } else {
    penalty
  }
  neighborhood$kernel <- kernel
  list(
    call = call,
    formula = formula,
    terms = prepared$terms,
    contrasts = prepared$contrasts,
    xlevels = prepared$xlevels,
    family = family,
    response_levels = prepared$response_levels,
    penalty = penalty,
    effective_penalty = effective_penalty,
    alpha = alpha,
    gamma = gamma,
    lambda = native$lambda,
    lambda_ratio = native$lambda_ratio,
    lambda_max = native$lambda_max,
    lambda_local_max = native$lambda_local_max,
    selected_index = selected,
    selected_lambda = native$lambda[selected],
    beta_hat = beta_hat,
    beta_hat_standardized = native$beta_hat_standardized,
    beta_path = beta_path,
    beta_path_standardized = native$beta_path_standardized,
    fitted_path = native$fitted_path,
    linear_predictor_path = native$linear_predictor_path,
    y_hat = fitted_values,
    linear_predictor = eta,
    residuals = residual_values,
    SDF = SDF,
    iterations = native$iterations,
    irls_iterations = native$irls_iterations,
    converged = native$converged,
    kkt = native$kkt,
    n_active = native$n_active,
    local_neff = native$local_neff,
    local_n_positive = native$local_n_positive,
    y = prepared$y,
    X = prepared$X,
    X_standardized = prepared$X_standardized,
    standardization = list(
      center = prepared$center,
      scale = prepared$scale,
      enabled = isTRUE(standardize),
      scope = "global"
    ),
    penalty_factor = prepared$penalty_factor,
    case_weights = prepared$case_weights,
    Z = prepared$Z,
    D = if (engine == "dense") prepared$D else NULL,
    knn_index = neighborhood$index,
    knn_distance = neighborhood$distance,
    bandwidth = neighborhood$bandwidth,
    bw = bw,
    kernel = kernel,
    adaptive = adaptive,
    bandwidth_type = bandwidth_type,
    geometry = geometry,
    distance_label = prepared$distance_label,
    tangent_ref = prepared$tangent_ref,
    engine = engine,
    intercept = prepared$intercept,
    predictor_names = prepared$predictor_names,
    keep_path = keep_path,
    control = control,
    solver_used = if (family == "binomial") {
      "naive"
    } else if (control$solver == "auto") {
      "auto (focal-dimension selection)"
    } else {
      control$solver
    }
  )
}

lhwpr_unstandardize_matrix <- function(beta, center, scale, intercept,
                                       predictor_names) {
  beta <- as.matrix(beta)
  if (intercept) {
    slopes <- sweep(beta[, -1L, drop = FALSE], 2L, scale, "/")
    intercept_value <- beta[, 1L] - as.numeric(slopes %*% center)
    out <- cbind(Intercept = intercept_value, slopes)
  } else {
    out <- sweep(beta, 2L, scale, "/")
  }
  colnames(out) <- c(if (intercept) "Intercept", predictor_names)
  out
}

lhwpr_unstandardize_path <- function(beta_path, center, scale, intercept,
                                     predictor_names) {
  dimensions <- dim(beta_path)
  out <- array(
    NA_real_,
    dim = dimensions,
    dimnames = list(
      NULL,
      c(if (intercept) "Intercept", predictor_names),
      NULL
    )
  )
  for (index in seq_len(dimensions[3L])) {
    out[, , index] <- lhwpr_unstandardize_matrix(
      beta_path[, , index, drop = FALSE][, , 1L],
      center,
      scale,
      intercept,
      predictor_names
    )
  }
  out
}

lhwpr_standardize <- function(X, case_weights, intercept, standardize) {
  p <- ncol(X)
  if (!standardize) {
    return(list(X = X, center = numeric(p), scale = rep(1, p)))
  }
  normalized <- case_weights / sum(case_weights)
  center <- if (intercept) {
    colSums(X * normalized)
  } else {
    numeric(p)
  }
  centered <- sweep(X, 2L, center, "-")
  scale <- sqrt(colSums(centered^2 * normalized))
  bad <- !is.finite(scale) | scale <= sqrt(.Machine$double.eps)
  if (any(bad)) {
    stop(
      "Cannot standardize constant predictor(s): ",
      paste(colnames(X)[bad], collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  list(
    X = sweep(centered, 2L, scale, "/"),
    center = center,
    scale = scale
  )
}

lhwpr_response <- function(response, family) {
  if (family == "gaussian") {
    if (!is.numeric(response) || is.matrix(response)) {
      stop("Gaussian LHWPR requires a numeric response.", call. = FALSE)
    }
    y <- as.numeric(response)
    if (any(!is.finite(y))) {
      stop("The response must contain only finite values.", call. = FALSE)
    }
    return(list(y = y, levels = NULL))
  }
  if (is.factor(response)) {
    if (nlevels(response) != 2L) {
      stop("Binomial LHWPR requires exactly two response levels.", call. = FALSE)
    }
    levels <- levels(response)
    return(list(
      y = as.numeric(response == levels[2L]),
      levels = levels
    ))
  }
  if (!is.numeric(response) || is.matrix(response) ||
      any(!is.finite(response)) ||
      any(!response %in% c(0, 1))) {
    stop(
      "Binomial LHWPR requires a two-level factor or numeric 0/1 response.",
      call. = FALSE
    )
  }
  list(y = as.numeric(response), levels = c("0", "1"))
}

lhwpr_case_weights <- function(case_weights, n) {
  if (is.null(case_weights)) return(rep(1, n))
  if (!is.numeric(case_weights) || length(case_weights) != n ||
      any(!is.finite(case_weights)) || any(case_weights < 0) ||
      !any(case_weights > 0)) {
    stop(
      "`case_weights` must be a finite non-negative vector with at least one positive value.",
      call. = FALSE
    )
  }
  as.numeric(case_weights)
}

lhwpr_penalty_factor <- function(penalty_factor, p, predictor_names) {
  if (is.null(penalty_factor)) return(rep(1, p))
  if (!is.numeric(penalty_factor) || length(penalty_factor) != p ||
      any(!is.finite(penalty_factor)) || any(penalty_factor < 0)) {
    stop(
      "`penalty_factor` must contain one finite non-negative value per predictor.",
      call. = FALSE
    )
  }
  if (!is.null(names(penalty_factor))) {
    if (!setequal(names(penalty_factor), predictor_names)) {
      stop(
        "Named `penalty_factor` values must match model-matrix predictors.",
        call. = FALSE
      )
    }
    penalty_factor <- penalty_factor[predictor_names]
  }
  as.numeric(penalty_factor)
}

lhwpr_validate_control <- function(control) {
  if (!inherits(control, "lhwpr_control")) {
    stop("`control` must be created by `lhwpr_control()`.", call. = FALSE)
  }
  control
}

lhwpr_validate_solver_family <- function(control, family) {
  if (family == "binomial" && control$solver == "covariance") {
    stop(
      "`solver = \"covariance\"` is available only for Gaussian LHWPR.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

lhwpr_alpha <- function(alpha) {
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      !is.finite(alpha) || alpha < 0 || alpha > 1) {
    stop("`alpha` must be one finite value in [0, 1].", call. = FALSE)
  }
  as.numeric(alpha)
}

lhwpr_validate_penalty_alpha <- function(penalty, alpha) {
  if (alpha > 0 && alpha < 1 && penalty != "lasso") {
    stop(
      "`0 < alpha < 1` is elastic net and requires `penalty = \"lasso\"`; MCP and SCAD are available at `alpha = 1`.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

lhwpr_gamma <- function(penalty, alpha, gamma) {
  if (alpha < 1 || penalty == "lasso") {
    if (is.null(gamma)) return(NA_real_)
    warning("`gamma` is ignored unless `alpha = 1` with MCP or SCAD.")
    return(NA_real_)
  }
  if (is.null(gamma)) gamma <- if (penalty == "mcp") 3 else 3.7
  if (!is.numeric(gamma) || length(gamma) != 1L || !is.finite(gamma)) {
    stop("`gamma` must be one finite numeric value.", call. = FALSE)
  }
  if (penalty == "mcp" && gamma <= 1) {
    stop("MCP requires `gamma > 1`.", call. = FALSE)
  }
  if (penalty == "scad" && gamma <= 2) {
    stop("SCAD requires `gamma > 2`.", call. = FALSE)
  }
  as.numeric(gamma)
}

lhwpr_validate_lambda <- function(lambda) {
  if (is.null(lambda)) return(NULL)
  if (!is.numeric(lambda) || !length(lambda) ||
      any(!is.finite(lambda)) || any(lambda < 0)) {
    stop("`lambda` must contain finite non-negative values.", call. = FALSE)
  }
  lambda <- as.numeric(lambda)
  if (length(lambda) > 1L && any(diff(lambda) >= 0)) {
    stop("`lambda` must be strictly decreasing.", call. = FALSE)
  }
  lambda
}

lhwpr_default_lambda_ratio <- function(alpha, n, p, ridge_path_tol) {
  if (alpha == 0) return(ridge_path_tol^2)
  if (n > p) 1e-4 else 1e-2
}

lhwpr_ratio <- function(value) {
  if (!is.numeric(value) || length(value) != 1L ||
      !is.finite(value) || value <= 0 || value >= 1) {
    stop("`lambda_min_ratio` must lie strictly between zero and one.", call. = FALSE)
  }
  as.numeric(value)
}

lhwpr_positive_scalar <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L ||
      !is.finite(value) || value <= 0) {
    stop("`", name, "` must be one positive finite number.", call. = FALSE)
  }
  as.numeric(value)
}

lhwpr_positive_integer <- function(value, name) {
  value <- lhwpr_positive_scalar(value, name)
  rounded <- as.integer(round(value))
  if (abs(value - rounded) > sqrt(.Machine$double.eps)) {
    stop("`", name, "` must be an integer.", call. = FALSE)
  }
  rounded
}

lhwpr_logical_scalar <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop("`", name, "` must be TRUE or FALSE.", call. = FALSE)
  }
  isTRUE(value)
}

#' @export
print.lhwpr <- function(x, ...) {
  if (!inherits(x, "lhwpr")) {
    stop("`x` must inherit from class 'lhwpr'.", call. = FALSE)
  }
  cat("LHWPR fit\n")
  cat("  Family:   ", x$family, "\n", sep = "")
  cat("  Penalty:  ", x$effective_penalty, "\n", sep = "")
  cat("  Alpha:    ", format(x$alpha), "\n", sep = "")
  cat("  Bandwidth:", format(x$bw), "\n", sep = " ")
  cat("  Lambda:   ", format(x$selected_lambda, digits = 6), "\n", sep = "")
  converged <- x$converged[, x$selected_index] != 0
  cat(
    "  Converged: ",
    sum(converged),
    "/",
    length(converged),
    " focal fits\n",
    sep = ""
  )
  invisible(x)
}

#' @export
summary.lhwpr <- function(object, ...) {
  selected <- object$selected_index
  nonzero <- object$n_active[, selected]
  out <- list(
    call = object$call,
    family = object$family,
    penalty = object$effective_penalty,
    alpha = object$alpha,
    gamma = object$gamma,
    bandwidth = object$bw,
    lambda = object$selected_lambda,
    lambda_ratio = object$lambda_ratio[selected],
    n = length(object$y),
    p = ncol(object$X),
    converged = sum(object$converged[, selected] != 0),
    kkt_max = max(object$kkt[, selected]),
    active_min = min(nonzero),
    active_median = stats::median(nonzero),
    active_max = max(nonzero),
    neff_min = min(object$local_neff),
    neff_median = stats::median(object$local_neff)
  )
  class(out) <- "summary.lhwpr"
  out
}

#' @export
print.summary.lhwpr <- function(x, ...) {
  cat("Summary of LHWPR fit\n")
  cat("  Family / penalty: ", x$family, " / ", x$penalty, "\n", sep = "")
  cat("  Observations:      ", x$n, "\n", sep = "")
  cat("  Predictors:        ", x$p, "\n", sep = "")
  cat("  Bandwidth:         ", format(x$bandwidth), "\n", sep = "")
  cat("  Selected lambda:   ", format(x$lambda, digits = 6), "\n", sep = "")
  cat(
    "  Converged fits:    ",
    x$converged,
    "/",
    x$n,
    "\n",
    sep = ""
  )
  cat("  Maximum KKT error: ", format(x$kkt_max, digits = 4), "\n", sep = "")
  cat(
    "  Active coefficients (min/median/max): ",
    paste(c(x$active_min, x$active_median, x$active_max), collapse = " / "),
    "\n",
    sep = ""
  )
  invisible(x)
}

#' @export
coef.lhwpr <- function(object, s = NULL, ...) {
  if (is.null(s)) return(object$beta_hat)
  index <- lhwpr_lambda_index(object, s)
  if (is.null(object$beta_path)) {
    if (index != object$selected_index) {
      stop(
        "The requested coefficient path was not retained; refit with `keep_path = \"all\"`.",
        call. = FALSE
      )
    }
    return(object$beta_hat)
  }
  object$beta_path[, , index, drop = FALSE][, , 1L]
}

#' @export
fitted.lhwpr <- function(object, ...) {
  object$y_hat
}

#' @export
residuals.lhwpr <- function(object, ...) {
  object$residuals
}

lhwpr_lambda_index <- function(object, s) {
  if (!is.numeric(s) || length(s) != 1L || !is.finite(s) || s < 0) {
    stop("`s` must be one non-negative lambda value.", call. = FALSE)
  }
  difference <- abs(object$lambda - s)
  index <- which.min(difference)
  tolerance <- 1e-10 * (1 + abs(s))
  if (difference[index] > tolerance) {
    stop("`s` must equal a lambda value stored in the fit.", call. = FALSE)
  }
  index
}

#' @export
predict.lhwpr <- function(object, newdata = NULL, Znew = NULL,
                          type = c("response", "link", "coefficients"),
                          s = NULL, ...) {
  type <- match.arg(type)
  index <- if (is.null(s)) object$selected_index else {
    lhwpr_lambda_index(object, s)
  }
  if (is.null(newdata)) {
    if (type == "response") return(object$fitted_path[, index])
    if (type == "link") return(object$linear_predictor_path[, index])
    return(stats::coef(object, s = object$lambda[index]))
  }
  if (is.null(Znew)) {
    stop("`Znew` is required when `newdata` is supplied.", call. = FALSE)
  }
  Terms <- stats::delete.response(object$terms)
  mf <- stats::model.frame(
    Terms,
    newdata,
    na.action = stats::na.fail,
    xlev = object$xlevels
  )
  X_full <- stats::model.matrix(
    Terms, mf, contrasts.arg = object$contrasts
  )
  if (object$intercept) {
    intercept_column <- match("(Intercept)", colnames(X_full), nomatch = 0L)
    if (!intercept_column) {
      stop("Prediction model matrix is missing its intercept.", call. = FALSE)
    }
    X_new <- X_full[, -intercept_column, drop = FALSE]
  } else {
    X_new <- X_full
  }
  if (!identical(colnames(X_new), object$predictor_names)) {
    stop(
      "Prediction model-matrix columns do not match the fitted model.",
      call. = FALSE
    )
  }
  X_new_standardized <- sweep(
    sweep(X_new, 2L, object$standardization$center, "-"),
    2L,
    object$standardization$scale,
    "/"
  )
  Znew <- as_lorentz_matrix(Znew, "Znew")
  check_hyperboloid(Znew, stop_on_fail = TRUE)
  if (nrow(Znew) != nrow(X_new)) {
    stop(
      "`Znew` must have the same number of rows as `newdata`.",
      call. = FALSE
    )
  }
  D_new <- lhwpr_cross_distance(object, Znew)
  bw_new <- lhwpr_prediction_bandwidth(object, D_new)
  prepared <- list(
    y = object$y,
    X_standardized = object$X_standardized,
    penalty_factor = object$penalty_factor,
    case_weights = object$case_weights,
    intercept = object$intercept
  )
  neighborhood <- list(
    engine = "dense",
    D = D_new,
    bandwidth = bw_new,
    kernel = object$kernel
  )
  native <- lhwpr_native_fit(
    prepared = prepared,
    neighborhood = neighborhood,
    lambda = object$lambda[index],
    family = object$family,
    penalty = object$penalty,
    alpha = object$alpha,
    gamma = object$gamma,
    nlambda = 1L,
    lambda_min_ratio = object$lambda_ratio[index],
    control = object$control,
    store_path = FALSE,
    exclude_self = FALSE,
    X_focus = X_new_standardized,
    compute_lambda_diagnostics = FALSE
  )
  if (type == "response") return(as.numeric(native$fitted_path[, 1L]))
  if (type == "link") {
    return(as.numeric(native$linear_predictor_path[, 1L]))
  }
  lhwpr_unstandardize_matrix(
    native$beta_hat_standardized,
    object$standardization$center,
    object$standardization$scale,
    object$intercept,
    object$predictor_names
  )
}

lhwpr_cross_distance <- function(object, Znew) {
  if (object$geometry == "lorentz") {
    return(t(hyperbolic_distance_matrix(
      Znew, W = object$Z, check = FALSE
    )))
  }
  mu <- object$tangent_ref
  train <- tangent_coordinates(object$Z, mu = mu)
  target <- tangent_coordinates(Znew, mu = mu)
  squared <- outer(rowSums(train^2), rowSums(target^2), "+") -
    2 * train %*% t(target)
  sqrt(pmax(squared, 0))
}

lhwpr_prediction_bandwidth <- function(object, D_new) {
  if (length(object$bw) != 1L) {
    stop(
      "New-location prediction requires a scalar fitted `bw`; a location-specific bandwidth vector cannot be extrapolated.",
      call. = FALSE
    )
  }
  if (!object$adaptive) {
    if (object$bandwidth_type == "volume") {
      return(rep(hyperbolic_radius_from_volume(object$bw), ncol(D_new)))
    }
    return(rep(as.numeric(object$bw), ncol(D_new)))
  }
  target_k <- as.integer(round(object$bw))
  if (object$bandwidth_type == "knn") {
    return(apply(D_new, 2L, function(distance) {
      positive <- sort(distance[distance > 1e-6])
      all_distance <- sort(distance)
      candidates <- if (length(positive) >= target_k) positive else all_distance
      candidates[target_k]
    }))
  }
  pilot_k <- resolve_pilot_k(nrow(D_new), NULL)
  pilot_radius <- apply(D_new, 2L, function(distance) {
    sort(distance[distance > 1e-6])[pilot_k]
  })
  density <- pilot_k / hyperbolic_ball_volume(pilot_radius)
  hyperbolic_radius_from_volume(target_k / density)
}

#' @export
plot.lhwpr <- function(x, focus = 1L, xvar = c("lambda", "lambda_ratio"),
                       ...) {
  if (is.null(x$beta_path)) {
    stop(
      "Coefficient paths were not retained; refit with `keep_path = \"all\"`.",
      call. = FALSE
    )
  }
  focus <- lhwpr_positive_integer(focus, "focus")
  if (focus > dim(x$beta_path)[1L]) {
    stop("`focus` is outside the fitted locations.", call. = FALSE)
  }
  xvar <- match.arg(xvar)
  horizontal <- if (xvar == "lambda") x$lambda else x$lambda_ratio
  coefficient_index <- if (x$intercept) {
    seq.int(2L, dim(x$beta_path)[2L])
  } else {
    seq_len(dim(x$beta_path)[2L])
  }
  paths <- t(matrix(
    x$beta_path[focus, coefficient_index, , drop = FALSE],
    nrow = length(coefficient_index),
    ncol = dim(x$beta_path)[3L]
  ))
  graphics::matplot(
    log10(pmax(horizontal, .Machine$double.xmin)),
    paths,
    type = "l",
    lty = 1,
    xlab = paste0("log10(", xvar, ")"),
    ylab = "Local coefficient",
    ...
  )
  invisible(x)
}

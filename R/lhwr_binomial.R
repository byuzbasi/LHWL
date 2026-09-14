lhwr_response <- function(response, family) {
  family <- match.arg(family, c("gaussian", "binomial"))
  if (family == "gaussian") {
    if (!is.numeric(response) || is.matrix(response)) {
      stop("Gaussian LHWR requires a numeric response.", call. = FALSE)
    }
    y <- as.numeric(response)
    if (any(!is.finite(y))) {
      stop("The response must contain only finite values.", call. = FALSE)
    }
    return(list(y = y, levels = NULL))
  }

  if (is.factor(response)) {
    if (nlevels(response) != 2L) {
      stop("Binomial LHWR requires exactly two response levels.", call. = FALSE)
    }
    response_levels <- levels(response)
    y <- as.numeric(response == response_levels[2L])
    if (length(unique(y)) != 2L) {
      stop("Binomial LHWR requires both response classes in the analysis sample.", call. = FALSE)
    }
    return(list(
      y = y,
      levels = response_levels
    ))
  }
  if (is.logical(response)) {
    if (anyNA(response)) {
      stop("Binomial LHWR does not allow missing responses.", call. = FALSE)
    }
    y <- as.numeric(response)
    if (length(unique(y)) != 2L) {
      stop("Binomial LHWR requires both response classes in the analysis sample.", call. = FALSE)
    }
    return(list(y = y, levels = c("FALSE", "TRUE")))
  }
  if (!is.numeric(response) || is.matrix(response) ||
      any(!is.finite(response)) || any(!response %in% c(0, 1))) {
    stop(
      "Binomial LHWR requires a two-level factor, logical, or numeric 0/1 response.",
      call. = FALSE
    )
  }
  y <- as.numeric(response)
  if (length(unique(y)) != 2L) {
    stop("Binomial LHWR requires both response classes in the analysis sample.", call. = FALSE)
  }
  list(y = y, levels = c("0", "1"))
}

lhwr_validate_glm_control <- function(control) {
  if (!is.list(control) ||
      is.null(control$epsilon) ||
      is.null(control$maxit) ||
      is.null(control$trace)) {
    stop(
      "`control` must be a value returned by `stats::glm.control()`.",
      call. = FALSE
    )
  }
  if (!is.numeric(control$epsilon) || length(control$epsilon) != 1L ||
      !is.finite(control$epsilon) || control$epsilon <= 0) {
    stop("`control$epsilon` must be one positive finite number.", call. = FALSE)
  }
  if (!is.numeric(control$maxit) || length(control$maxit) != 1L ||
      !is.finite(control$maxit) || control$maxit < 1) {
    stop("`control$maxit` must be a positive integer.", call. = FALSE)
  }
  if (!is.logical(control$trace) || length(control$trace) != 1L ||
      is.na(control$trace)) {
    stop("`control$trace` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  list(
    epsilon = as.numeric(control$epsilon),
    maxit = as.integer(round(control$maxit)),
    trace = isTRUE(control$trace)
  )
}

lhwr_binomial_design <- function(X, intercept) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  if (any(!is.finite(X))) {
    stop("The model matrix must contain only finite values.", call. = FALSE)
  }
  intercept <- isTRUE(intercept)
  if (intercept) {
    intercept_column <- match("Intercept", colnames(X), nomatch = 0L)
    if (!intercept_column) {
      stop("Internal LHWR model matrix is missing its intercept.", call. = FALSE)
    }
    predictors <- X[, -intercept_column, drop = FALSE]
  } else {
    predictors <- X
  }
  if (!intercept && !ncol(predictors)) {
    stop("LHWR requires an intercept or at least one predictor.", call. = FALSE)
  }
  list(
    X = predictors,
    intercept = intercept,
    coefficient_names = c(
      if (intercept) "Intercept",
      colnames(predictors)
    )
  )
}

lhwr_binomial_native_options <- function(intercept, kernel, n_threads,
                                         control, return_inference,
                                         exclude_self) {
  control <- lhwr_validate_glm_control(control)
  list(
    family = "binomial",
    penalty = "lasso",
    alpha = 0,
    gamma = NA_real_,
    intercept = isTRUE(intercept),
    solver = "naive",
    screening = "none",
    tol = control$epsilon,
    kkt_tol = control$epsilon,
    max_iter = 100000L,
    irls_tol = control$epsilon,
    irls_max_iter = control$maxit,
    max_active = 0L,
    early_stop = FALSE,
    ridge_path_tol = 1e-6,
    n_threads = as.integer(n_threads),
    trace = control$trace,
    store_path = FALSE,
    return_inference = isTRUE(return_inference),
    selected_index = 0L,
    exclude_self = isTRUE(exclude_self),
    nlambda = 1L,
    lambda_min_ratio = 1,
    kernel = kernel
  )
}

lhwr_binomial_fit_native <- function(y, X, X_focus, neighborhood,
                                     intercept, kernel, n_threads, control,
                                     return_inference = FALSE,
                                     exclude_self = FALSE) {
  design <- lhwr_binomial_design(X, intercept)
  focus_design <- lhwr_binomial_design(X_focus, intercept)
  if (!identical(
    colnames(design$X),
    colnames(focus_design$X)
  )) {
    stop("Calibration and focal model matrices are not compatible.", call. = FALSE)
  }
  options <- lhwr_binomial_native_options(
    intercept = design$intercept,
    kernel = kernel,
    n_threads = n_threads,
    control = control,
    return_inference = return_inference,
    exclude_self = exclude_self
  )
  penalty_factor <- rep(0, ncol(design$X))
  case_weights <- rep(1, length(y))
  native <- if (neighborhood$engine == "dense") {
    .lhwpr_fit_dense_cpp(
      y = y,
      X = design$X,
      X_focus = focus_design$X,
      D = neighborhood$D,
      bw = neighborhood$bandwidth,
      lambda = 0,
      penalty_factor = penalty_factor,
      case_weights = case_weights,
      options = options
    )
  } else {
    .lhwpr_fit_knn_cpp(
      y = y,
      X = design$X,
      X_focus = focus_design$X,
      knn_index = neighborhood$index,
      knn_distance = neighborhood$distance,
      lambda = 0,
      penalty_factor = penalty_factor,
      case_weights = case_weights,
      options = options
    )
  }

  beta_hat <- as.matrix(native$beta_hat_standardized)
  colnames(beta_hat) <- design$coefficient_names
  se_beta <- as.matrix(native$standard_error)
  colnames(se_beta) <- paste0(design$coefficient_names, "_SE")
  z_beta <- beta_hat / se_beta
  colnames(z_beta) <- paste0(design$coefficient_names, "_ZV")
  selected <- as.integer(native$selected_index)
  fitted_values <- as.numeric(native$fitted_path[, selected])
  list(
    beta_hat = beta_hat,
    se_beta = se_beta,
    z_beta = z_beta,
    y_hat = fitted_values,
    linear_predictor = as.numeric(
      native$linear_predictor_path[, selected]
    ),
    residuals = if (length(fitted_values) == length(y)) {
      y - fitted_values
    } else {
      NULL
    },
    converged = as.logical(native$converged[, selected]),
    iterations = as.integer(native$iterations[, selected]),
    irls_iterations = as.integer(native$irls_iterations[, selected]),
    kkt = as.numeric(native$kkt[, selected]),
    local_kappa = as.numeric(native$local_kappa),
    local_rank = as.integer(native$local_rank),
    local_neff = as.numeric(native$local_neff),
    local_n_positive = as.integer(native$local_n_positive),
    local_weight_sum = as.numeric(native$local_weight_sum),
    local_class_variation = as.logical(native$local_class_variation),
    local_events = as.integer(native$local_events),
    local_nonevents = as.integer(native$local_nonevents),
    bandwidth = neighborhood$bandwidth,
    knn_index = neighborhood$index,
    knn_distance = neighborhood$distance
  )
}

lhwr_binomial_neighborhood <- function(Z, D, bw, kernel, engine,
                                       n_threads = 1L) {
  if (engine == "dense") {
    return(list(
      engine = "dense",
      D = D,
      bandwidth = as.numeric(bw),
      index = NULL,
      distance = NULL
    ))
  }
  k <- as.integer(round(bw))
  knn <- lorentz_knn(
    Z,
    k = k,
    include_self = FALSE,
    check = FALSE
  )
  list(
    engine = "knn",
    D = NULL,
    bandwidth = as.numeric(knn$distance[, k]),
    index = knn$index,
    distance = knn$distance
  )
}

lhwr_clip_probability <- function(probability) {
  pmin(pmax(as.numeric(probability), 1e-8), 1 - 1e-8)
}

lhwr_binomial_deviance <- function(y, probability) {
  probability <- lhwr_clip_probability(probability)
  -2 * sum(
    y * log(probability) +
      (1 - y) * log1p(-probability)
  )
}

lhwr_binomial_log_loss <- function(y, probability) {
  probability <- lhwr_clip_probability(probability)
  -mean(
    y * log(probability) +
      (1 - y) * log1p(-probability)
  )
}

lhwr_binomial_bw_selection <- function(
    y, X, Z, D, bw_grid, search, adaptive, bandwidth_type,
    kernel, engine, lower, upper, tol, max_iter, pilot_k,
    n_threads, intercept, control) {
  bounds <- if (engine == "knn") {
    resolve_knn_bw_bounds(length(y), lower, upper)
  } else {
    resolve_bw_bounds(D, adaptive, bandwidth_type, lower, upper)
  }
  if (!is.finite(bounds$lower) || !is.finite(bounds$upper) ||
      bounds$lower >= bounds$upper) {
    stop("Invalid binomial bandwidth search interval.", call. = FALSE)
  }

  candidates_needed <- if (search == "grid") {
    as.numeric(bw_grid)
  } else {
    c(bounds$lower, bounds$upper)
  }
  if (adaptive) {
    candidates_needed <- as.integer(round(candidates_needed))
  }
  if (any(!is.finite(candidates_needed)) ||
      any(candidates_needed <= 0)) {
    stop(
      "Binomial bandwidth candidates must be positive and finite.",
      call. = FALSE
    )
  }
  if (adaptive &&
      any(candidates_needed < 2L | candidates_needed >= length(y))) {
    stop(
      "Adaptive binomial bandwidth candidates must be integers between 2 and n - 1.",
      call. = FALSE
    )
  }

  knn_full <- NULL
  if (engine == "knn") {
    maximum_k <- if (search == "grid") {
      max(candidates_needed)
    } else {
      as.integer(round(bounds$upper))
    }
    knn_full <- lorentz_knn(
      Z,
      k = maximum_k,
      include_self = FALSE,
      check = FALSE
    )
  }

  trace <- list()
  evaluated <- numeric()
  evaluate <- function(candidate) {
    if (search == "golden") {
      candidate <- min(max(candidate, bounds$lower), bounds$upper)
    }
    if (adaptive) candidate <- as.integer(round(candidate))
    key <- format(candidate, scientific = FALSE, trim = TRUE, digits = 17)
    existing <- match(key, names(evaluated), nomatch = 0L)
    if (existing) return(unname(evaluated[existing]))

    neighborhood <- if (engine == "knn") {
      k <- as.integer(candidate)
      list(
        engine = "knn",
        D = NULL,
        bandwidth = as.numeric(knn_full$distance[, k]),
        index = knn_full$index[, seq_len(k), drop = FALSE],
        distance = knn_full$distance[, seq_len(k), drop = FALSE]
      )
    } else {
      list(
        engine = "dense",
        D = D,
        bandwidth = resolve_bandwidth(
          D,
          candidate,
          adaptive = adaptive,
          bandwidth_type = bandwidth_type,
          pilot_k = pilot_k
        ),
        index = NULL,
        distance = NULL
      )
    }

    fit <- lhwr_binomial_fit_native(
      y = y,
      X = X,
      X_focus = X,
      neighborhood = neighborhood,
      intercept = intercept,
      kernel = kernel,
      n_threads = n_threads,
      control = control,
      return_inference = FALSE,
      exclude_self = TRUE
    )
    raw_score <- lhwr_binomial_log_loss(y, fit$y_hat)
    failed <- sum(!fit$converged)
    score <- if (failed) Inf else raw_score
    trace[[length(trace) + 1L]] <<- data.frame(
      bw = as.numeric(candidate),
      score = score,
      trS = NA_real_,
      RSS = NA_real_,
      CV = score,
      raw_log_loss = raw_score,
      converged = sum(fit$converged),
      not_converged = failed,
      kappa_median = NA_real_,
      kappa_max = NA_real_,
      neff_min = min(fit$local_neff),
      neff_median = stats::median(fit$local_neff),
      n_positive_min = min(fit$local_n_positive),
      n_positive_median = stats::median(fit$local_n_positive),
      iteration = length(trace) + 1L,
      check.names = FALSE
    )
    evaluated[key] <<- score
    score
  }

  if (search == "grid") {
    invisible(vapply(candidates_needed, evaluate, numeric(1L)))
  } else {
    a <- bounds$lower
    b <- bounds$upper
    ratio <- (sqrt(5) - 1) / 2
    c_point <- b - ratio * (b - a)
    d_point <- a + ratio * (b - a)
    c_score <- evaluate(c_point)
    d_score <- evaluate(d_point)
    iteration <- 0L
    repeat {
      iteration <- iteration + 1L
      interval_small <- if (adaptive) {
        b - a <= 1
      } else {
        b - a <= tol
      }
      if (interval_small || iteration >= as.integer(max_iter)) break
      if (c_score <= d_score) {
        b <- d_point
        d_point <- c_point
        d_score <- c_score
        c_point <- b - ratio * (b - a)
        c_score <- evaluate(c_point)
      } else {
        a <- c_point
        c_point <- d_point
        c_score <- d_score
        d_point <- a + ratio * (b - a)
        d_score <- evaluate(d_point)
      }
    }
    final_candidates <- if (adaptive) {
      seq(
        max(bounds$lower, floor(a) - 1),
        min(bounds$upper, ceiling(b) + 1)
      )
    } else {
      c(a, b, (a + b) / 2)
    }
    invisible(vapply(final_candidates, evaluate, numeric(1L)))
  }

  selection <- do.call(rbind, trace)
  selection <- unique_selection_table(selection)
  attr(selection, "lower") <- bounds$lower
  attr(selection, "upper") <- bounds$upper
  selection
}

lhwr_binomial_diagnostics <- function(y, fitted, converged,
                                      local_class_variation,
                                      local_kappa = NULL,
                                      local_rank = NULL,
                                      coefficient_count = NULL,
                                      local_neff = NULL,
                                      local_n_positive = NULL) {
  probability <- lhwr_clip_probability(fitted)
  null_probability <- rep(mean(y), length(y))
  deviance <- lhwr_binomial_deviance(y, probability)
  null_deviance <- lhwr_binomial_deviance(y, null_probability)
  out <- list(
    family = "binomial",
    deviance = deviance,
    null.deviance = null_deviance,
    pseudo.R2 = if (null_deviance > 0) {
      1 - deviance / null_deviance
    } else {
      NA_real_
    },
    logLik = -0.5 * deviance,
    brier = mean((y - probability)^2),
    converged = sum(converged),
    not.converged = sum(!converged),
    class.variation.failures = sum(!local_class_variation)
  )
  finite_kappa <- local_kappa[is.finite(local_kappa)]
  if (length(finite_kappa)) {
    out$kappa.median <- stats::median(finite_kappa)
    out$kappa.max <- max(finite_kappa)
  }
  if (!is.null(local_rank) && !is.null(coefficient_count)) {
    out$rank.deficient <- sum(
      is.na(local_rank) | local_rank < coefficient_count
    )
  }
  if (!is.null(local_neff)) {
    out$neff.min <- min(local_neff, na.rm = TRUE)
    out$neff.median <- stats::median(local_neff, na.rm = TRUE)
  }
  if (!is.null(local_n_positive)) {
    out$n.positive.min <- min(local_n_positive, na.rm = TRUE)
    out$n.positive.median <- stats::median(local_n_positive, na.rm = TRUE)
  }
  out
}

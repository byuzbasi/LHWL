#' Joint Bandwidth and Regularization Tuning for LHWPR
#'
#' Evaluate bandwidth, lambda, and optional alpha/gamma grids using common
#' validation folds and a common absolute lambda path within each
#' bandwidth/penalty configuration. Gaussian fits also report conditional
#' effective degrees of freedom and information criteria.
#'
#' @param formula,data,Z Model specification; see [lhwpr()].
#' @param bw_grid Candidate bandwidths. When `NULL`, an adaptive kNN grid is
#'   generated automatically from `bw_lower`, `bw_upper`, and `n_bw`.
#' @param bw_lower,bw_upper Optional adaptive-kNN bounds used only when
#'   `bw_grid = NULL`.
#' @param n_bw Number of automatically generated adaptive-kNN candidates.
#' @param family,penalty Model family and structural penalty.
#' @param alpha_grid Candidate alpha values.
#' @param gamma_grid Optional MCP/SCAD gamma values.
#' @param lambda Optional common decreasing absolute lambda path.
#' @param nlambda,lambda_min_ratio Automatic path controls.
#' @param penalty_factor,case_weights,standardize Penalized preprocessing
#'   controls.
#' @param kernel,adaptive,bandwidth_type,geometry,tangent_ref,dMat,engine
#'   Locality and geometry controls.
#' @param resampling Leave-one-out, K-fold, or user-blocked validation.
#' @param nfolds Number of folds for automatically constructed K-fold
#'   validation.
#' @param fold_id Optional fold identifier per observation. It is required for
#'   blocked validation.
#'   Without `fold_id`, K-fold allocation cycles through the input row order
#'   deterministically; it is not randomly shuffled. Supply a frozen random
#'   or scientifically defined blocked assignment when that is intended.
#' @param rule Select the minimum score, the most regularized model within one
#'   standard error, or the sparsest model within one standard error. The
#'   one-standard-error rules are available only with `criterion = "cv"`.
#' @param control A value returned by [lhwpr_control()].
#' @param refit Logical; fit the selected model on all observations.
#' @param criterion Selection score. Cross-validation is the default; Gaussian
#'   models additionally support AIC, AICc, BIC, GCV, path-based EBIC, and a
#'   refit-based EBIC sensitivity criterion.
#' @param ebic_gamma EBIC model-space penalty in `[0, 1]`.
#' @param refine Logical or `NULL`. `NULL` refines automatically generated
#'   bandwidth grids and leaves supplied grids unchanged.
#' @param n_refine Number of candidates in the local bandwidth refinement.
#' @param refine_alpha Logical; add midpoints adjacent to the provisionally
#'   selected alpha, following the efficient two-stage idea used by `ensr`.
#' @param expand_bounds Logical; allow an automatically chosen lower bandwidth
#'   boundary to expand toward the model-dimension safety limit.
#' @param compare_lhwr Logical; evaluate an unpenalized LHWR reference on the
#'   same folds and bandwidths.
#' @param eligibility_kkt_tol Maximum recorded KKT violation allowed when a
#'   tuning row is eligible for selection. Candidate rows must also have all
#'   relevant focal fits converged. Cross-validation rows must pass the gate
#'   for both validation and full-data path fits.
#' @details Gaussian path-information criteria can select only rows with
#'   `edf_valid = TRUE`; refit EBIC additionally requires
#'   `post_edf_valid = TRUE`. Every selected row must also pass the recorded
#'   convergence and KKT gates. The one-standard-error rules preserve those
#'   same eligibility requirements. The selection table records eligibility,
#'   its reason, the logarithmic EBIC model-space term, and nonempty
#'   full-surface models for which that term is zero.
#'
#'   Tuning is conditional on the supplied representation and prepared design.
#'   Predictor standardization, the representation, and its tangent reference
#'   are prepared once for the supplied sample; the function does not learn an
#'   embedding or recompute these transformations separately within each fold.
#'   Use an outer training/validation workflow when assessing an entire learned
#'   representation pipeline or reporting performance after tuning. In that
#'   workflow, fit all data-dependent transformations on the outer training
#'   sample and apply its frozen transformations to validation observations.
#' @return An object of class `"tune.lhwpr"`.
#' @export
tune.lhwpr <- function(
    formula,
    data,
    Z,
    bw_grid = NULL,
    family = c("gaussian", "binomial"),
    penalty = c("lasso", "mcp", "scad"),
    alpha_grid = 1,
    gamma_grid = NULL,
    lambda = NULL,
    nlambda = 100L,
    lambda_min_ratio = NULL,
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
    resampling = c("kfold", "loo", "blocked"),
    nfolds = 5L,
    fold_id = NULL,
    rule = c("min", "1se", "sparse_1se"),
    control = lhwpr_control(),
    refit = TRUE,
    bw_lower = NULL,
    bw_upper = NULL,
    n_bw = 15L,
    criterion = c(
      "cv", "aicc", "aic", "bic", "gcv", "ebic", "ebic_refit"
    ),
    ebic_gamma = 0.5,
    refine = NULL,
    n_refine = 7L,
    refine_alpha = FALSE,
    expand_bounds = TRUE,
    compare_lhwr = TRUE,
    eligibility_kkt_tol = 1e-5) {
  this_call <- match.call()
  automatic_bandwidth <- is.null(bw_grid)
  family <- match.arg(family)
  penalty <- match.arg(penalty)
  kernel <- match.arg(kernel)
  bandwidth_type <- match.arg(bandwidth_type)
  geometry <- match.arg(geometry)
  engine <- match.arg(engine)
  resampling <- match.arg(resampling)
  rule <- match.arg(rule)
  criterion <- match.arg(criterion)
  control <- lhwpr_validate_control(control)
  lhwpr_validate_solver_family(control, family)
  refit <- lhwpr_logical_scalar(refit, "refit")
  refine_alpha <- lhwpr_logical_scalar(refine_alpha, "refine_alpha")
  expand_bounds <- lhwpr_logical_scalar(expand_bounds, "expand_bounds")
  compare_lhwr <- lhwpr_logical_scalar(compare_lhwr, "compare_lhwr")
  eligibility_kkt_tol <- lhwpr_positive_scalar(
    eligibility_kkt_tol, "eligibility_kkt_tol"
  )
  if (is.null(refine)) {
    refine <- automatic_bandwidth
  } else {
    refine <- lhwpr_logical_scalar(refine, "refine")
  }
  n_refine <- lhwpr_positive_integer(n_refine, "n_refine")
  if (n_refine < 3L) {
    stop("`n_refine` must be at least three.", call. = FALSE)
  }
  if (!is.numeric(ebic_gamma) || length(ebic_gamma) != 1L ||
      !is.finite(ebic_gamma) || ebic_gamma < 0 || ebic_gamma > 1) {
    stop("`ebic_gamma` must be one finite value in [0, 1].", call. = FALSE)
  }
  if (family != "gaussian" && criterion != "cv") {
    stop(
      "AIC, AICc, BIC, GCV, EBIC, and refit-EBIC tuning are currently defined only for Gaussian LHWPR; use `criterion = \"cv\"` for binomial models.",
      call. = FALSE
    )
  }
  if (criterion != "cv" && rule != "min") {
    stop(
      "The one-standard-error rules require `criterion = \"cv\"`.",
      call. = FALSE
    )
  }
  nlambda <- lhwpr_positive_integer(nlambda, "nlambda")
  lambda <- lhwpr_validate_lambda(lambda)
  alpha_grid <- vapply(alpha_grid, lhwpr_alpha, numeric(1))
  for (alpha in alpha_grid) {
    lhwpr_validate_penalty_alpha(penalty, alpha)
  }
  if (refine_alpha && penalty != "lasso") {
    stop(
      "`refine_alpha = TRUE` is available only for the lasso/elastic-net alpha continuum.",
      call. = FALSE
    )
  }
  gamma_grid <- lhwpr_tune_gamma(penalty, alpha_grid, gamma_grid)

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
  bandwidth_search <- lhwpr_tune_bandwidths(
    bw_grid = bw_grid,
    n = nrow(prepared$X),
    n_coef = ncol(prepared$X) + as.integer(prepared$intercept),
    adaptive = adaptive,
    bandwidth_type = bandwidth_type,
    bw_lower = bw_lower,
    bw_upper = bw_upper,
    n_bw = n_bw
  )
  bw_grid <- bandwidth_search$grid
  folds <- lhwpr_tune_folds(
    n = length(prepared$y),
    resampling = resampling,
    nfolds = nfolds,
    fold_id = fold_id
  )
  cache <- new.env(parent = emptyenv())
  lhwpr_tune_evaluate_grid(
    cache = cache,
    bw_grid = bw_grid,
    alpha_grid = alpha_grid,
    gamma_grid = gamma_grid,
    stage = "coarse",
    prepared = prepared,
    family = family,
    penalty = penalty,
    lambda = lambda,
    nlambda = nlambda,
    lambda_min_ratio = lambda_min_ratio,
    kernel = kernel,
    adaptive = adaptive,
    bandwidth_type = bandwidth_type,
    engine = engine,
    control = control,
    resampling = resampling,
    folds = folds,
    ebic_gamma = ebic_gamma
  )
  coarse_table <- lhwpr_tune_cache_table(cache)
  coarse_row <- lhwpr_select_tuning_row(
    coarse_table,
    criterion = criterion,
    rule = "min",
    kkt_tol = eligibility_kkt_tol
  )
  coarse_selected <- coarse_table[coarse_row, , drop = FALSE]

  refined_bw <- numeric()
  if (refine) {
    refined_bw <- lhwpr_tune_refined_bandwidths(
      selected_bw = coarse_selected$bw,
      evaluated_bw = bw_grid,
      bandwidth_search = bandwidth_search,
      adaptive = adaptive,
      bandwidth_type = bandwidth_type,
      n_refine = n_refine,
      expand_bounds = expand_bounds
    )
  }
  refined_alpha <- if (refine_alpha) {
    lhwpr_tune_refined_alpha(
      selected_alpha = coarse_selected$alpha,
      alpha_grid = alpha_grid
    )
  } else {
    numeric()
  }
  if (length(refined_bw) || length(refined_alpha)) {
    lhwpr_tune_evaluate_grid(
      cache = cache,
      bw_grid = sort(unique(c(bw_grid, refined_bw))),
      alpha_grid = sort(unique(c(alpha_grid, refined_alpha))),
      gamma_grid = gamma_grid,
      stage = "refined",
      prepared = prepared,
      family = family,
      penalty = penalty,
      lambda = lambda,
      nlambda = nlambda,
      lambda_min_ratio = lambda_min_ratio,
      kernel = kernel,
      adaptive = adaptive,
      bandwidth_type = bandwidth_type,
      engine = engine,
      control = control,
      resampling = resampling,
      folds = folds,
      ebic_gamma = ebic_gamma
    )
  }

  selection_table <- lhwpr_tune_cache_table(cache)
  selected_row <- lhwpr_select_tuning_row(
    selection_table,
    criterion = criterion,
    rule = rule,
    kkt_tol = eligibility_kkt_tol
  )
  eligibility <- lhwpr_tuning_eligibility(
    selection_table, criterion = criterion, kkt_tol = eligibility_kkt_tol
  )
  selection_table$eligible <- eligibility$eligible
  selection_table$ineligibility_reason <- eligibility$reason
  selection_table$selected <- FALSE
  selection_table$selected[selected_row] <- TRUE
  evaluated_bw <- sort(unique(selection_table$bw))
  selection_table$bw_boundary <- ifelse(
    selection_table$bw == min(evaluated_bw),
    "lower",
    ifelse(selection_table$bw == max(evaluated_bw), "upper", "none")
  )
  selected <- selection_table[selected_row, , drop = FALSE]
  bandwidth_search$coarse_grid <- bw_grid
  bandwidth_search$refined_grid <- sort(unique(refined_bw))
  bandwidth_search$evaluated_grid <- evaluated_bw
  bandwidth_search$expanded_lower <-
    length(refined_bw) &&
    bandwidth_search$mode == "automatic" &&
    min(refined_bw) < bandwidth_search$lower
  bandwidth_search$expanded_upper <-
    length(refined_bw) &&
    bandwidth_search$mode == "automatic" &&
    max(refined_bw) > bandwidth_search$upper
  bandwidth_search$selected_at_boundary <-
    selected$bw_boundary != "none"
  if (bandwidth_search$selected_at_boundary) {
    bandwidth_search$boundary_warning <- paste0(
      "The selected bandwidth lies on the ",
      selected$bw_boundary,
      " search boundary; inspect the bandwidth profile before final inference."
    )
  } else {
    bandwidth_search$boundary_warning <- NULL
  }

  configuration_keys <- unique(selection_table$configuration_key)
  selection_table$configuration <- match(
    selection_table$configuration_key, configuration_keys
  )
  selected <- selection_table[selected_row, , drop = FALSE]
  fold_scores <- lapply(
    configuration_keys,
    function(key) get(key, envir = cache, inherits = FALSE)$fold_scores
  )
  names(fold_scores) <- configuration_keys

  lhwr_comparison <- NULL
  if (compare_lhwr) {
    lhwr_comparison <- lhwpr_tune_lhwr_comparison(
      bw_grid = evaluated_bw,
      prepared = prepared,
      family = family,
      kernel = kernel,
      adaptive = adaptive,
      bandwidth_type = bandwidth_type,
      engine = engine,
      control = control,
      resampling = resampling,
      folds = folds,
      ebic_gamma = ebic_gamma
    )
  }

  selected_gamma <- if (penalty == "lasso" || selected$alpha == 0) {
    NULL
  } else {
    selected$gamma
  }

  selected_lambda_path <- lhwpr_selected_lambda_path(
    selection_table, selected
  )
  fit <- NULL
  refit_diagnostics <- NULL
  if (refit) {
    fit <- lhwpr(
      formula = formula,
      data = data,
      Z = Z,
      bw = selected$bw,
      family = family,
      penalty = penalty,
      alpha = selected$alpha,
      lambda = selected_lambda_path,
      nlambda = length(selected_lambda_path),
      gamma = selected_gamma,
      penalty_factor = penalty_factor,
      case_weights = case_weights,
      standardize = standardize,
      kernel = kernel,
      adaptive = adaptive,
      bandwidth_type = bandwidth_type,
      geometry = geometry,
      tangent_ref = tangent_ref,
      dMat = dMat,
      engine = engine,
      keep_path = "selected",
      control = control
    )
    refit_converged <- fit$converged[, fit$selected_index] != 0L
    refit_kkt_max <- max(fit$kkt[, fit$selected_index])
    refit_kkt_limit <- 5 * control$kkt_tol
    refit_passed <- all(refit_converged) && is.finite(refit_kkt_max) &&
      refit_kkt_max <= refit_kkt_limit
    refit_diagnostics <- data.frame(
      path_length = length(selected_lambda_path),
      selected_index = fit$selected_index,
      selected_lambda = fit$selected_lambda,
      convergence_rate = mean(refit_converged),
      kkt_max = refit_kkt_max,
      kkt_limit = refit_kkt_limit,
      passed = refit_passed,
      stringsAsFactors = FALSE
    )
    if (!refit_passed) {
      stop(
        "The selected LHWPR refit failed convergence/KKT checks: ",
        sum(refit_converged), "/", length(refit_converged),
        " focal fits converged; maximum KKT violation = ",
        format(refit_kkt_max, digits = 6),
        "; allowed = ", format(refit_kkt_limit, digits = 6), ".",
        call. = FALSE
      )
    }
  }
  out <- list(
    call = this_call,
    selection_table = selection_table,
    selected = selected,
    selected_index = selected_row,
    criterion = criterion,
    rule = rule,
    resampling = resampling,
    fold_id = folds,
    fold_scores = fold_scores,
    bandwidth_search = bandwidth_search,
    lhwr_comparison = lhwr_comparison,
    search_history = list(
      coarse = sort(unique(bw_grid)),
      refined = sort(unique(refined_bw)),
      alpha_refined = sort(unique(refined_alpha))
    ),
    ebic_gamma = ebic_gamma,
    eligibility_kkt_tol = eligibility_kkt_tol,
    selected_lambda_path = selected_lambda_path,
    refit_diagnostics = refit_diagnostics,
    family = family,
    penalty = penalty,
    kernel = kernel,
    fit = fit
  )
  class(out) <- "tune.lhwpr"
  out
}

lhwpr_tune_configurations <- function(bw_grid, alpha_grid, gamma_grid,
                                      penalty) {
  configurations <- expand.grid(
    bw = bw_grid,
    alpha = alpha_grid,
    gamma = gamma_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  if (penalty == "lasso") {
    configurations$gamma <- NA_real_
    configurations <- configurations[
      !duplicated(configurations[c("bw", "alpha")]),
      ,
      drop = FALSE
    ]
  } else if (any(configurations$alpha == 0)) {
    configurations$gamma[configurations$alpha == 0] <- NA_real_
    configurations <- configurations[
      !duplicated(configurations[c("bw", "alpha", "gamma")]),
      ,
      drop = FALSE
    ]
  }
  configurations
}

lhwpr_tune_key <- function(bw, alpha, gamma) {
  paste(
    format(bw, digits = 17, scientific = FALSE, trim = TRUE),
    format(alpha, digits = 17, scientific = FALSE, trim = TRUE),
    if (is.na(gamma)) {
      "NA"
    } else {
      format(gamma, digits = 17, scientific = FALSE, trim = TRUE)
    },
    sep = "|"
  )
}

lhwpr_tune_evaluate_grid <- function(cache, bw_grid, alpha_grid, gamma_grid,
                                     stage, prepared, family, penalty, lambda,
                                     nlambda, lambda_min_ratio, kernel,
                                     adaptive, bandwidth_type, engine, control,
                                     resampling, folds, ebic_gamma) {
  configurations <- lhwpr_tune_configurations(
    bw_grid, alpha_grid, gamma_grid, penalty
  )
  for (row in seq_len(nrow(configurations))) {
    configuration <- configurations[row, , drop = FALSE]
    key <- lhwpr_tune_key(
      configuration$bw, configuration$alpha, configuration$gamma
    )
    if (exists(key, envir = cache, inherits = FALSE)) next
    result <- lhwpr_tune_configuration(
      configuration = configuration,
      configuration_key = key,
      stage = stage,
      prepared = prepared,
      family = family,
      penalty = penalty,
      lambda = lambda,
      nlambda = nlambda,
      lambda_min_ratio = lambda_min_ratio,
      kernel = kernel,
      adaptive = adaptive,
      bandwidth_type = bandwidth_type,
      engine = engine,
      control = control,
      resampling = resampling,
      folds = folds,
      ebic_gamma = ebic_gamma
    )
    assign(key, result, envir = cache)
  }
  invisible(cache)
}

lhwpr_tune_configuration <- function(configuration, configuration_key, stage,
                                     prepared, family, penalty, lambda,
                                     nlambda, lambda_min_ratio, kernel,
                                     adaptive, bandwidth_type, engine, control,
                                     resampling, folds, ebic_gamma) {
  alpha <- configuration$alpha
  gamma <- if (penalty == "lasso" || alpha == 0) {
    NA_real_
  } else {
    configuration$gamma
  }
  ratio <- if (is.null(lambda_min_ratio)) {
    lhwpr_default_lambda_ratio(
      alpha,
      nrow(prepared$X),
      ncol(prepared$X),
      control$ridge_path_tol
    )
  } else {
    lhwpr_ratio(lambda_min_ratio)
  }
  neighborhood <- lhwpr_neighborhood(
    prepared = prepared,
    bw = configuration$bw,
    kernel = kernel,
    adaptive = adaptive,
    bandwidth_type = bandwidth_type,
    engine = engine
  )
  path_info <- lhwpr_tune_path(
    prepared = prepared,
    neighborhood = neighborhood,
    family = family,
    penalty = penalty,
    alpha = alpha,
    gamma = gamma,
    lambda = lambda,
    nlambda = nlambda,
    lambda_min_ratio = ratio,
    control = control,
    resampling = resampling,
    folds = folds
  )
  validation <- lhwpr_tune_predictions(
    prepared = prepared,
    neighborhood = neighborhood,
    family = family,
    penalty = penalty,
    alpha = alpha,
    gamma = gamma,
    lambda = path_info$lambda,
    lambda_min_ratio = ratio,
    control = control,
    resampling = resampling,
    folds = folds
  )
  loss <- lhwpr_validation_loss(
    prepared$y, validation$prediction, family
  )
  scores <- lhwpr_fold_scores(
    loss, folds, prepared$case_weights
  )
  cv_mean <- colMeans(scores)
  cv_se <- if (nrow(scores) > 1L) {
    apply(scores, 2L, stats::sd) / sqrt(nrow(scores))
  } else {
    rep(0, ncol(scores))
  }

  full_fit <- lhwpr_native_fit(
    prepared = prepared,
    neighborhood = neighborhood,
    lambda = path_info$lambda,
    family = family,
    penalty = penalty,
    alpha = alpha,
    gamma = gamma,
    nlambda = length(path_info$lambda),
    lambda_min_ratio = ratio,
    control = control,
    store_path = FALSE,
    exclude_self = FALSE,
    return_path_diagnostics = family == "gaussian",
    return_active_frequency = TRUE,
    compute_lambda_diagnostics = FALSE
  )
  metrics <- if (family == "gaussian") {
    metrics <- lhwpr_information_criteria(
      prepared = prepared,
      native = full_fit,
      ebic_gamma = ebic_gamma
    )
  } else {
    lhwpr_binomial_path_metrics(prepared, full_fit)
  }

  table <- data.frame(
    configuration_key = configuration_key,
    stage = stage,
    bw = configuration$bw,
    alpha = alpha,
    gamma = gamma,
    lambda = path_info$lambda,
    lambda_ratio = path_info$lambda / path_info$lambda_max,
    cv_mean = cv_mean,
    cv_se = cv_se,
    mean = cv_mean,
    se = cv_se,
    cv_convergence_rate = colMeans(validation$converged != 0),
    cv_kkt_max = apply(validation$kkt, 2L, max),
    metrics,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  list(table = table, fold_scores = scores)
}

lhwpr_tune_cache_table <- function(cache) {
  keys <- ls(cache, all.names = TRUE)
  tables <- lapply(
    keys,
    function(key) get(key, envir = cache, inherits = FALSE)$table
  )
  out <- do.call(rbind, tables)
  out <- out[order(
    out$bw, out$alpha, out$gamma,
    -out$lambda, na.last = TRUE
  ), , drop = FALSE]
  rownames(out) <- NULL
  out
}

lhwpr_tune_path <- function(prepared, neighborhood, family, penalty, alpha,
                            gamma, lambda, nlambda, lambda_min_ratio, control,
                            resampling, folds) {
  if (!is.null(lambda)) {
    maxima <- if (all(prepared$penalty_factor == 0) &&
                  all(lambda == 0)) {
      1
    } else {
      lhwpr_tune_lambda_maxima(
        prepared,
        neighborhood,
        family,
        penalty,
        alpha,
        gamma,
        nlambda,
        lambda_min_ratio,
        control,
        resampling,
        folds
      )
    }
    return(list(lambda = lambda, lambda_max = maxima))
  }
  maximum <- lhwpr_tune_lambda_maxima(
    prepared,
    neighborhood,
    family,
    penalty,
    alpha,
    gamma,
    nlambda,
    lambda_min_ratio,
    control,
    resampling,
    folds
  )
  path <- if (nlambda == 1L) {
    maximum
  } else {
    exp(seq(
      log(maximum),
      log(maximum * lambda_min_ratio),
      length.out = nlambda
    ))
  }
  list(lambda = path, lambda_max = maximum)
}

lhwpr_tune_lambda_maxima <- function(prepared, neighborhood, family, penalty,
                                     alpha, gamma, nlambda, lambda_min_ratio,
                                     control, resampling, folds) {
  if (resampling == "loo") {
    result <- lhwpr_native_lambda_max(
      prepared = prepared,
      neighborhood = neighborhood,
      family = family,
      penalty = penalty,
      alpha = alpha,
      gamma = gamma,
      nlambda = nlambda,
      lambda_min_ratio = lambda_min_ratio,
      control = control,
      exclude_self = TRUE
    )
    return(result$lambda_max)
  }
  maximum <- 0
  for (fold in sort(unique(folds))) {
    training_weight <- prepared$case_weights
    training_weight[folds == fold] <- 0
    result <- lhwpr_native_lambda_max(
      prepared = prepared,
      neighborhood = neighborhood,
      family = family,
      penalty = penalty,
      alpha = alpha,
      gamma = gamma,
      nlambda = nlambda,
      lambda_min_ratio = lambda_min_ratio,
      control = control,
      exclude_self = FALSE,
      case_weights = training_weight
    )
    maximum <- max(maximum, result$lambda_max)
  }
  maximum
}

lhwpr_tune_predictions <- function(prepared, neighborhood, family, penalty,
                                   alpha, gamma, lambda, lambda_min_ratio,
                                   control, resampling, folds) {
  n <- length(prepared$y)
  prediction <- matrix(NA_real_, nrow = n, ncol = length(lambda))
  converged <- matrix(0L, nrow = n, ncol = length(lambda))
  kkt <- matrix(Inf, nrow = n, ncol = length(lambda))
  if (resampling == "loo") {
    fit <- lhwpr_native_fit(
      prepared = prepared,
      neighborhood = neighborhood,
      lambda = lambda,
      family = family,
      penalty = penalty,
      alpha = alpha,
      gamma = gamma,
      nlambda = length(lambda),
      lambda_min_ratio = lambda_min_ratio,
      control = control,
      store_path = FALSE,
      exclude_self = TRUE,
      compute_lambda_diagnostics = FALSE
    )
    prediction[,] <- fit$fitted_path
    converged[,] <- fit$converged
    kkt[,] <- fit$kkt
    return(list(
      prediction = prediction,
      converged = converged,
      kkt = kkt
    ))
  }
  for (fold in sort(unique(folds))) {
    validation <- which(folds == fold)
    training_weight <- prepared$case_weights
    training_weight[validation] <- 0
    focal_neighborhood <- lhwpr_tune_focal_neighborhood(
      neighborhood, validation
    )
    fit <- lhwpr_native_fit(
      prepared = prepared,
      neighborhood = focal_neighborhood,
      lambda = lambda,
      family = family,
      penalty = penalty,
      alpha = alpha,
      gamma = gamma,
      nlambda = length(lambda),
      lambda_min_ratio = lambda_min_ratio,
      control = control,
      store_path = FALSE,
      exclude_self = identical(neighborhood$engine, "knn"),
      case_weights = training_weight,
      X_focus = prepared$X_standardized[validation, , drop = FALSE],
      compute_lambda_diagnostics = FALSE
    )
    prediction[validation, ] <- fit$fitted_path
    converged[validation, ] <- fit$converged
    kkt[validation, ] <- fit$kkt
  }
  list(prediction = prediction, converged = converged, kkt = kkt)
}

lhwpr_tune_focal_neighborhood <- function(neighborhood, focus) {
  out <- neighborhood
  if (identical(neighborhood$engine, "dense")) {
    out$D <- neighborhood$D[, focus, drop = FALSE]
  } else {
    out$index <- neighborhood$index[focus, , drop = FALSE]
    out$distance <- neighborhood$distance[focus, , drop = FALSE]
  }
  if (length(neighborhood$bandwidth) > 1L) {
    out$bandwidth <- neighborhood$bandwidth[focus]
  }
  out
}

lhwpr_validation_loss <- function(y, prediction, family) {
  if (family == "gaussian") {
    return((y - prediction)^2)
  }
  probability <- pmin(pmax(prediction, 1e-12), 1 - 1e-12)
  -(
    y * log(probability) +
      (1 - y) * log1p(-probability)
  )
}

lhwpr_fold_scores <- function(loss, folds, case_weights) {
  fold_values <- sort(unique(folds))
  out <- matrix(
    NA_real_,
    nrow = length(fold_values),
    ncol = ncol(loss)
  )
  for (index in seq_along(fold_values)) {
    use <- folds == fold_values[index]
    weights <- case_weights[use]
    total <- sum(weights)
    if (!is.finite(total) || total <= 0) {
      stop(
        "Every validation fold must contain positive total case weight.",
        call. = FALSE
      )
    }
    weights <- weights / total
    out[index, ] <- colSums(loss[use, , drop = FALSE] * weights)
  }
  out
}

lhwpr_criterion_column <- function(criterion) {
  switch(
    criterion,
    cv = "cv_mean",
    aic = "aic",
    aicc = "aicc",
    bic = "bic",
    gcv = "gcv",
    ebic = "ebic",
    ebic_refit = "ebic_refit"
  )
}

lhwpr_tuning_eligibility <- function(selection_table, criterion,
                                     kkt_tol = 1e-5) {
  score_column <- lhwpr_criterion_column(criterion)
  score <- selection_table[[score_column]]
  eligible <- is.finite(score)
  reasons <- rep("eligible", nrow(selection_table))
  reasons[!eligible] <- "non-finite score"

  require_gate <- function(rate_name, kkt_name, label) {
    if (!all(c(rate_name, kkt_name) %in% names(selection_table))) {
      stop(
        "The tuning table is missing `", rate_name, "` or `", kkt_name,
        "` for numerical eligibility.",
        call. = FALSE
      )
    }
    rate <- selection_table[[rate_name]]
    kkt <- selection_table[[kkt_name]]
    gate <- is.finite(rate) & rate >= 1 - 1e-12 &
      is.finite(kkt) & kkt <= kkt_tol
    newly_ineligible <- eligible & !gate
    reasons[newly_ineligible] <<- paste0(label, " convergence/KKT gate")
    eligible <<- eligible & gate
  }

  require_gate("convergence_rate", "kkt_max", "full-fit")
  if (criterion == "cv") {
    require_gate("cv_convergence_rate", "cv_kkt_max", "validation")
  }
  if (criterion %in% c("aic", "aicc", "bic", "gcv", "ebic")) {
    if (!"edf_valid" %in% names(selection_table)) {
      stop(
        "The tuning table is missing `edf_valid` for ",
        toupper(criterion), " selection.",
        call. = FALSE
      )
    }
    edf_gate <-
      !is.na(selection_table$edf_valid) &
      selection_table$edf_valid
    reasons[eligible & !edf_gate] <- "invalid effective degrees of freedom"
    eligible <- eligible & edf_gate
  }
  if (criterion == "ebic_refit") {
    if (!"post_edf_valid" %in% names(selection_table)) {
      stop(
        "The tuning table is missing `post_edf_valid` for refit-EBIC selection.",
        call. = FALSE
      )
    }
    post_gate <-
      !is.na(selection_table$post_edf_valid) &
      selection_table$post_edf_valid
    reasons[eligible & !post_gate] <-
      "invalid post-selection effective degrees of freedom"
    eligible <- eligible & post_gate
  }
  list(eligible = eligible, reason = reasons)
}

lhwpr_select_tuning_row <- function(selection_table, criterion, rule,
                                    kkt_tol = 1e-5) {
  score_column <- lhwpr_criterion_column(criterion)
  score <- selection_table[[score_column]]
  eligibility <- lhwpr_tuning_eligibility(
    selection_table, criterion = criterion, kkt_tol = kkt_tol
  )
  eligible_mask <- eligibility$eligible
  finite <- which(eligible_mask)
  if (!length(finite)) {
    stop(
      "No numerically eligible ", toupper(criterion),
      " score was produced by the candidate models.",
      call. = FALSE
    )
  }
  minimum_value <- min(score[finite])
  tied <- finite[
    abs(score[finite] - minimum_value) <=
      sqrt(.Machine$double.eps) * (1 + abs(minimum_value))
  ]
  sparsity <- ifelse(
    is.finite(selection_table$surface_size),
    selection_table$surface_size,
    selection_table$mean_active
  )
  minimum <- tied[order(
    sparsity[tied],
    -selection_table$lambda_ratio[tied],
    -selection_table$bw[tied],
    -selection_table$alpha[tied]
  )[1L]]
  if (rule == "min") return(minimum)
  eligible <- which(
    eligible_mask &
      selection_table$cv_mean <=
      selection_table$cv_mean[minimum] +
      selection_table$cv_se[minimum]
  )
  ordering <- if (rule == "sparse_1se") {
    order(
      sparsity[eligible],
      -selection_table$lambda_ratio[eligible],
      -selection_table$bw[eligible],
      -selection_table$alpha[eligible],
      selection_table$cv_mean[eligible]
    )
  } else {
    order(
      -selection_table$lambda_ratio[eligible],
      -selection_table$bw[eligible],
      -selection_table$alpha[eligible],
      selection_table$cv_mean[eligible]
    )
  }
  eligible[ordering[1L]]
}

lhwpr_selected_lambda_path <- function(selection_table, selected) {
  required <- c("bw", "lambda")
  if (!is.data.frame(selection_table) || !is.data.frame(selected) ||
      nrow(selected) != 1L ||
      !all(required %in% names(selection_table)) ||
      !all(required %in% names(selected))) {
    stop(
      "A one-row selection and its tuning table are required for exact-path refit.",
      call. = FALSE
    )
  }
  same_configuration <- selection_table$bw == selected$bw
  if (all(c("configuration_key") %in% names(selection_table)) &&
      "configuration_key" %in% names(selected)) {
    same_configuration <- same_configuration &
      selection_table$configuration_key == selected$configuration_key
  } else {
    for (name in intersect(
      c("alpha", "gamma"),
      intersect(names(selection_table), names(selected))
    )) {
      value <- selected[[name]][[1L]]
      same_configuration <- same_configuration & if (is.na(value)) {
        is.na(selection_table[[name]])
      } else {
        !is.na(selection_table[[name]]) &
          selection_table[[name]] == value
      }
    }
  }
  candidate <- as.numeric(selection_table$lambda[same_configuration])
  candidate <- candidate[is.finite(candidate) & candidate >= 0]
  candidate <- sort(unique(candidate), decreasing = TRUE)
  selected_lambda <- as.numeric(selected$lambda[[1L]])
  selected_index <- match(
    sprintf("%.17g", selected_lambda), sprintf("%.17g", candidate)
  )
  if (is.na(selected_index)) {
    stop(
      "The selected lambda is absent from its exact tuning configuration.",
      call. = FALSE
    )
  }
  path <- candidate[seq_len(selected_index)]
  path[[length(path)]] <- selected_lambda
  if (any(!is.finite(path)) || any(path < 0) ||
      (length(path) > 1L && any(diff(path) >= 0))) {
    stop(
      "The selected tuning prefix is not a finite strictly decreasing lambda path.",
      call. = FALSE
    )
  }
  path
}

lhwpr_tune_folds <- function(n, resampling, nfolds, fold_id) {
  if (resampling == "loo") return(seq_len(n))
  if (is.null(fold_id)) {
    if (resampling == "blocked") {
      stop(
        "`fold_id` is required for blocked validation.",
        call. = FALSE
      )
    }
    nfolds <- lhwpr_positive_integer(nfolds, "nfolds")
    if (nfolds < 2L || nfolds > n) {
      stop("`nfolds` must be between 2 and n.", call. = FALSE)
    }
    return(rep(seq_len(nfolds), length.out = n))
  }
  if (length(fold_id) != n || anyNA(fold_id)) {
    stop("`fold_id` must contain one non-missing value per observation.", call. = FALSE)
  }
  folds <- match(fold_id, unique(fold_id))
  if (length(unique(folds)) < 2L) {
    stop("At least two validation folds are required.", call. = FALSE)
  }
  folds
}

lhwpr_tune_bandwidths <- function(bw_grid, n, n_coef, adaptive,
                                  bandwidth_type, bw_lower, bw_upper, n_bw) {
  n_bw <- lhwpr_positive_integer(n_bw, "n_bw")
  if (!is.null(bw_grid)) {
    if (!is.numeric(bw_grid) || !length(bw_grid) ||
        any(!is.finite(bw_grid)) || any(bw_grid <= 0)) {
      stop("`bw_grid` must contain positive finite values.", call. = FALSE)
    }
    if (!is.null(bw_lower) || !is.null(bw_upper)) {
      stop("`bw_lower` and `bw_upper` require `bw_grid = NULL`.", call. = FALSE)
    }
    return(list(
      mode = "supplied",
      grid = unique(as.numeric(bw_grid)),
      lower = NA_real_,
      upper = NA_real_,
      hard_lower = NA_real_,
      hard_upper = NA_real_,
      explicit_lower = TRUE,
      explicit_upper = TRUE,
      n_bw = length(unique(bw_grid))
    ))
  }
  if (!isTRUE(adaptive) || bandwidth_type != "knn") {
    stop(
      "Automatic LHWPR bandwidth grids currently require adaptive kNN bandwidths; supply `bw_grid` otherwise.",
      call. = FALSE
    )
  }
  if (n < 4L) {
    stop("Automatic bandwidth selection requires at least four observations.", call. = FALSE)
  }
  hard_lower <- min(n - 2L, max(2L, n_coef + 2L))
  hard_upper <- n - 1L
  default_lower <- min(n - 2L, max(20L, 2L * n_coef + 2L))
  default_upper <- hard_upper
  lower <- if (is.null(bw_lower)) default_lower else bw_lower
  upper <- if (is.null(bw_upper)) default_upper else bw_upper
  if (!is.numeric(lower) || length(lower) != 1L || !is.finite(lower) ||
      !is.numeric(upper) || length(upper) != 1L || !is.finite(upper)) {
    stop("`bw_lower` and `bw_upper` must be finite scalar values.", call. = FALSE)
  }
  lower <- as.integer(round(lower))
  upper <- as.integer(round(upper))
  if (lower < 2L || upper > n - 1L || lower >= upper) {
    stop(
      "Automatic adaptive-kNN bounds must satisfy `2 <= bw_lower < bw_upper <= n - 1`.",
      call. = FALSE
    )
  }
  grid <- unique(as.numeric(as.integer(round(seq(lower, upper, length.out = n_bw)))))
  if (length(grid) < 2L) {
    stop("Automatic bandwidth bounds must generate at least two candidates.", call. = FALSE)
  }
  list(
    mode = "automatic",
    grid = grid,
    lower = lower,
    upper = upper,
    hard_lower = hard_lower,
    hard_upper = hard_upper,
    explicit_lower = !is.null(bw_lower),
    explicit_upper = !is.null(bw_upper),
    n_bw = length(grid)
  )
}

lhwpr_tune_refined_bandwidths <- function(
    selected_bw, evaluated_bw, bandwidth_search, adaptive, bandwidth_type,
    n_refine, expand_bounds) {
  values <- sort(unique(as.numeric(evaluated_bw)))
  if (length(values) < 2L) return(numeric())
  selected_index <- which.min(abs(values - selected_bw))
  left <- values[max(1L, selected_index - 1L)]
  right <- values[min(length(values), selected_index + 1L)]

  if (selected_index == 1L &&
      bandwidth_search$mode == "automatic" &&
      expand_bounds &&
      !bandwidth_search$explicit_lower) {
    step <- values[2L] - values[1L]
    left <- max(bandwidth_search$hard_lower, values[1L] - step)
  }
  if (selected_index == length(values) &&
      bandwidth_search$mode == "automatic" &&
      expand_bounds &&
      !bandwidth_search$explicit_upper) {
    step <- values[length(values)] - values[length(values) - 1L]
    right <- min(bandwidth_search$hard_upper, values[length(values)] + step)
  }
  candidates <- seq(left, right, length.out = n_refine)
  if (isTRUE(adaptive) && bandwidth_type == "knn") {
    candidates <- as.numeric(as.integer(round(candidates)))
  }
  candidates <- sort(unique(candidates))
  keep <- vapply(
    candidates,
    function(value) {
      !any(abs(values - value) <=
        sqrt(.Machine$double.eps) * (1 + abs(value)))
    },
    logical(1)
  )
  candidates[keep]
}

lhwpr_tune_refined_alpha <- function(selected_alpha, alpha_grid) {
  values <- sort(unique(as.numeric(alpha_grid)))
  if (length(values) < 2L) return(numeric())
  selected_index <- which.min(abs(values - selected_alpha))
  neighbor_index <- unique(c(selected_index - 1L, selected_index + 1L))
  neighbor_index <- neighbor_index[
    neighbor_index >= 1L & neighbor_index <= length(values)
  ]
  sort(unique((values[selected_index] + values[neighbor_index]) / 2))
}

lhwpr_empty_path_metrics <- function(nlambda) {
  data.frame(
    rss = rep(NA_real_, nlambda),
    deviance = rep(NA_real_, nlambda),
    logLik = rep(NA_real_, nlambda),
    trS = rep(NA_real_, nlambda),
    edf = rep(NA_real_, nlambda),
    aic = rep(NA_real_, nlambda),
    aicc = rep(NA_real_, nlambda),
    bic = rep(NA_real_, nlambda),
    gcv = rep(NA_real_, nlambda),
    ebic = rep(NA_real_, nlambda),
    ebic_refit = rep(NA_real_, nlambda),
    ebic_model_space_log = rep(NA_real_, nlambda),
    ebic_model_space_zero = rep(FALSE, nlambda),
    post_rss = rep(NA_real_, nlambda),
    post_deviance = rep(NA_real_, nlambda),
    post_edf = rep(NA_real_, nlambda),
    post_residual_df = rep(NA_real_, nlambda),
    post_edf_valid = rep(FALSE, nlambda),
    mean_active = rep(NA_real_, nlambda),
    median_active = rep(NA_real_, nlambda),
    surface_size = rep(NA_integer_, nlambda),
    convergence_rate = rep(NA_real_, nlambda),
    kkt_max = rep(NA_real_, nlambda),
    edf_valid = rep(FALSE, nlambda),
    stringsAsFactors = FALSE
  )
}

lhwpr_information_criteria <- function(prepared, native, ebic_gamma) {
  if (is.null(native$hat_diag_path) ||
      is.null(native$post_fitted_path) ||
      is.null(native$post_hat_diag_path) ||
      is.null(native$surface_active_count)) {
    stop("Internal Gaussian path diagnostics are incomplete.", call. = FALSE)
  }
  positive <- prepared$case_weights > 0
  n <- sum(positive)
  analysis_weights <- prepared$case_weights
  analysis_weights <- analysis_weights * n / sum(analysis_weights)
  weighted_rss <- function(fitted) {
    colSums(
      (prepared$y - fitted)^2 * analysis_weights
    )
  }
  gaussian_deviance <- function(rss) {
    n * (
      log(pmax(rss / n, .Machine$double.xmin)) +
        log(2 * pi) + 1
    )
  }

  rss <- weighted_rss(native$fitted_path)
  trS <- colSums(native$hat_diag_path)
  edf <- trS + 1
  deviance <- gaussian_deviance(rss)
  log_likelihood <- -0.5 * deviance
  aic <- deviance + 2 * edf
  denominator <- n - edf - 1
  aicc <- ifelse(
    denominator > 0,
    aic + 2 * edf * (edf + 1) / denominator,
    Inf
  )
  bic <- deviance + log(n) * edf
  gcv_denominator <- 1 - trS / n
  gcv <- ifelse(
    gcv_denominator > 0,
    (rss / n) / gcv_denominator^2,
    Inf
  )

  post_rss <- weighted_rss(native$post_fitted_path)
  post_trS <- colSums(native$post_hat_diag_path)
  post_edf <- post_trS + 1
  post_deviance <- gaussian_deviance(post_rss)
  active_count <- as.matrix(native$surface_active_count)
  surface_size <- colSums(active_count > 0)
  p_search <- sum(prepared$penalty_factor > 0)
  model_space <- vapply(
    surface_size,
    function(size) {
      if (size < 0 || size > p_search) return(Inf)
      if (p_search == 0L) return(0)
      lchoose(p_search, size)
    },
    numeric(1)
  )
  ebic <- bic + 2 * ebic_gamma * model_space
  ebic_model_space_zero <- is.finite(model_space) &
    model_space == 0 &
    surface_size > 0
  post_residual_df <- n - post_edf
  post_edf_valid <- is.finite(post_edf) &
    post_edf >= 1 &
    post_residual_df > 1
  ebic_refit <- post_deviance + log(n) * post_edf +
    2 * ebic_gamma * model_space
  ebic_refit[!post_edf_valid] <- Inf

  data.frame(
    rss = rss,
    deviance = deviance,
    logLik = log_likelihood,
    trS = trS,
    edf = edf,
    aic = aic,
    aicc = aicc,
    bic = bic,
    gcv = gcv,
    ebic = ebic,
    ebic_refit = ebic_refit,
    ebic_model_space_log = model_space,
    ebic_model_space_zero = ebic_model_space_zero,
    post_rss = post_rss,
    post_deviance = post_deviance,
    post_edf = post_edf,
    post_residual_df = post_residual_df,
    post_edf_valid = post_edf_valid,
    mean_active = colMeans(native$n_active),
    median_active = apply(native$n_active, 2L, stats::median),
    surface_size = as.integer(surface_size),
    convergence_rate = colMeans(native$converged != 0),
    kkt_max = apply(native$kkt, 2L, max),
    edf_valid = is.finite(edf) & edf >= 1 & edf < n - 1,
    stringsAsFactors = FALSE
  )
}

lhwpr_binomial_path_metrics <- function(prepared, native) {
  metrics <- lhwpr_empty_path_metrics(ncol(native$fitted_path))
  n <- sum(prepared$case_weights > 0)
  analysis_weights <- prepared$case_weights
  analysis_weights <- analysis_weights * n / sum(analysis_weights)
  probability <- pmin(
    pmax(native$fitted_path, 1e-12),
    1 - 1e-12
  )
  y <- matrix(
    prepared$y,
    nrow = length(prepared$y),
    ncol = ncol(probability)
  )
  weights <- matrix(
    analysis_weights,
    nrow = length(analysis_weights),
    ncol = ncol(probability)
  )
  negative_log_likelihood <- -colSums(
    weights * (
      y * log(probability) +
        (1 - y) * log1p(-probability)
    )
  )
  active_count <- as.matrix(native$surface_active_count)
  metrics$rss <- colSums(weights * (y - probability)^2)
  metrics$deviance <- 2 * negative_log_likelihood
  metrics$logLik <- -negative_log_likelihood
  metrics$mean_active <- colMeans(native$n_active)
  metrics$median_active <- apply(
    native$n_active, 2L, stats::median
  )
  metrics$surface_size <- as.integer(colSums(active_count > 0))
  metrics$convergence_rate <- colMeans(native$converged != 0)
  metrics$kkt_max <- apply(native$kkt, 2L, max)
  metrics
}

lhwpr_tune_lhwr_comparison <- function(
    bw_grid, prepared, family, kernel, adaptive, bandwidth_type, engine,
    control, resampling, folds, ebic_gamma) {
  reference <- prepared
  reference$penalty_factor[] <- 0
  rows <- lapply(
    bw_grid,
    function(bw) {
      result <- lhwpr_tune_configuration(
        configuration = data.frame(
          bw = bw, alpha = 1, gamma = NA_real_
        ),
        configuration_key = paste0("LHWR|", format(bw, digits = 17)),
        stage = "reference",
        prepared = reference,
        family = family,
        penalty = "lasso",
        lambda = 0,
        nlambda = 1L,
        lambda_min_ratio = 1e-4,
        kernel = kernel,
        adaptive = adaptive,
        bandwidth_type = bandwidth_type,
        engine = engine,
        control = control,
        resampling = resampling,
        folds = folds,
        ebic_gamma = ebic_gamma
      )
      result$table
    }
  )
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out$model <- "LHWR"
  out
}

lhwpr_tune_gamma <- function(penalty, alpha_grid, gamma_grid) {
  if (penalty == "lasso") return(NA_real_)
  if (any(!alpha_grid %in% c(0, 1))) {
    stop(
      "With MCP or SCAD, `alpha_grid` may contain only zero (ridge) and one.",
      call. = FALSE
    )
  }
  if (is.null(gamma_grid)) {
    gamma_grid <- if (penalty == "mcp") 3 else 3.7
  }
  vapply(
    gamma_grid,
    function(value) lhwpr_gamma(penalty, 1, value),
    numeric(1)
  )
}

#' @export
print.tune.lhwpr <- function(x, ...) {
  cat("LHWPR joint tuning\n")
  cat("  Resampling: ", x$resampling, "\n", sep = "")
  cat("  Criterion:  ", toupper(x$criterion), "\n", sep = "")
  cat("  Rule:       ", x$rule, "\n", sep = "")
  cat("  Bandwidth:  ", format(x$selected$bw), "\n", sep = "")
  cat("  Alpha:      ", format(x$selected$alpha), "\n", sep = "")
  if (is.finite(x$selected$gamma)) {
    cat("  Gamma:      ", format(x$selected$gamma), "\n", sep = "")
  }
  cat("  Lambda:     ", format(x$selected$lambda, digits = 6), "\n", sep = "")
  score_column <- lhwpr_criterion_column(x$criterion)
  cat(
    "  Score:      ",
    format(x$selected[[score_column]], digits = 6),
    "\n",
    sep = ""
  )
  if ("eligible" %in% names(x$selection_table)) {
    cat(
      "  Eligible:   ", sum(x$selection_table$eligible), "/",
      nrow(x$selection_table), " candidate rows\n", sep = ""
    )
  }
  if (!is.null(x$refit_diagnostics)) {
    cat(
      "  Refit:      exact path prefix (", x$refit_diagnostics$path_length,
      " levels), KKT certificate ",
      if (isTRUE(x$refit_diagnostics$passed)) "passed" else "failed",
      "\n", sep = ""
    )
  }
  if (isTRUE(x$bandwidth_search$selected_at_boundary)) {
    cat("  Boundary:   selected bandwidth is on the search boundary\n")
  }
  invisible(x)
}

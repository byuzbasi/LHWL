#' LHWR Diagnostic F-Tests
#'
#' Compute Leung-type, moment-matched diagnostic F-tests adapted to a fitted
#' LHWR model. The diagnostics use the LHWR smoother induced by the model's
#' kernel weights and distance geometry. Their reference distributions use
#' generally non-integer effective degrees of freedom. The reference laws are
#' approximate and condition on the fitted geometry and selected bandwidth;
#' they should not be read as post-selection confirmatory tests. They are
#' optional because coefficient-specific F3 diagnostics require reconstructing
#' local coefficient operators and can be costly for large data sets.
#'
#' @param object A fitted model object.
#' @param ... Additional arguments passed to methods.
#' @return A diagnostic test object.
#' @details
#' F1 uses a lower-tail moment-matched F reference; F2 uses an upper-tail
#' reference for improvement over the stored global OLS model. F3 uses an
#' upper-tail reference for each coefficient's constancy across locations,
#' not for that coefficient being zero at one location. The quadratic forms
#' are generally dependent, so these are not exact nested-model ANOVA tests.
#' Independent homoskedastic Gaussian errors and fixed stationary-model weights
#' motivate the approximations; bandwidth selection, correlated errors and
#' varying nuisance coefficients limit their interpretation.
#'
#' Required moments and reference degrees of freedom must be positive; do not
#' interpret undefined p-values as zero. When supplying `ols_model`, verify the
#' same response, predictor columns and observation order, not only row count.
#' The default stored OLS model avoids a separate sample/formula choice.
#' See `vignette("lhwr-diagnostics-and-tuning", package = "LHWL")` for the
#' statistics, reference degrees of freedom and executable reporting tables.
#' @examples
#' # Example 1: globally stationary coefficients.
#' # The diagnostics should usually not reject global adequacy or local
#' # coefficient constancy because the response does not depend on location.
#' set.seed(1)
#' H <- simulate_H2(80, r_max = 3)
#' Z <- H$Z
#' x1 <- rnorm(80)
#' x2 <- rnorm(80)
#' y <- 1 + 2 * x1 - 0.5 * x2 + rnorm(80, sd = 0.3)
#' dat <- data.frame(y, x1, x2)
#'
#' fit_global <- lhwr(
#'   y ~ x1 + x2,
#'   data = dat,
#'   Z = Z,
#'   bw = 30,
#'   adaptive = TRUE,
#'   diagnostics = "aicc",
#'   f_tests = TRUE
#' )
#' fit_global$Ftests
#'
#' # Example 2: hyperbolic coefficient heterogeneity.
#' # The coefficients are generated from localized bumps defined by
#' # Lorentz-hyperbolic distances, so the diagnostics may detect departures
#' # from the global model and coefficient-specific local variability.
#' set.seed(2)
#' H <- simulate_H2(80, r_max = 3, uniform_area = TRUE)
#' Z <- H$Z
#' D <- hyperbolic_distance_matrix(Z, check = FALSE)
#' x1 <- rnorm(80)
#' x2 <- rnorm(80)
#'
#' phi1 <- exp(-D[, 15]^2 / 0.7^2)
#' phi2 <- exp(-D[, 45]^2 / 0.7^2)
#' phi3 <- exp(-D[, 70]^2 / 0.7^2)
#' beta0 <- 1 + 1.2 * phi1
#' beta1 <- 2 + 1.5 * phi2
#' beta2 <- -0.5 - 1.3 * phi3
#'
#' eta <- beta0 + beta1 * x1 + beta2 * x2
#' y <- eta + rnorm(80, sd = 0.2)
#' dat <- data.frame(y, x1, x2)
#'
#' fit_local <- lhwr(
#'   y ~ x1 + x2,
#'   data = dat,
#'   Z = Z,
#'   bw = 25,
#'   adaptive = TRUE,
#'   diagnostics = "aicc",
#'   f_tests = TRUE
#' )
#' diagnostic_tests(fit_local)
#' @export
diagnostic_tests <- function(object, ...) {
  UseMethod("diagnostic_tests")
}

#' @rdname diagnostic_tests
#' @param ols_model Optional global linear model. If omitted, the `lm` object
#'   stored by [lhwr()] is used.
#' @param digits Number of digits used by the print method.
#' @param engine Computational engine for reconstructing F-test smoother
#'   components. `"cpp"` uses the compiled implementation, including the
#'   memory-saving kNN path when the model was fitted with `engine = "knn"`.
#'   If the fitted object stores kNN neighborhoods, the C++ engine reuses them
#'   instead of recomputing nearest neighbors;
#'   `"R"` uses the reference R implementation, and `"auto"` selects `"cpp"`.
#' @export
diagnostic_tests.lhwr <- function(object, ols_model = NULL, digits = 3,
                                  engine = c("auto", "cpp", "R"), ...) {
  engine <- match.arg(engine)
  if (engine == "auto") engine <- "cpp"
  if (!is_lhwr_object(object)) {
    stop("`object` must inherit from class 'lhwr'.", call. = FALSE)
  }
  if (identical(object$family %||% object$LHWR.arguments$family, "binomial")) {
    stop(
      "Diagnostic F-tests are defined only for Gaussian LHWR.",
      call. = FALSE
    )
  }
  if (is.null(ols_model)) {
    ols_model <- object$lm
  }
  if (is.null(ols_model) || !inherits(ols_model, "lm")) {
    stop("A global `lm` object is required for diagnostic F-tests.", call. = FALSE)
  }
  if (is.null(object$LHWR.diagnostic$tr.S) || is.null(object$LHWR.diagnostic$tr.StS)) {
    stop("Diagnostic F-tests require hat-matrix diagnostics. Refit with `diagnostics != \"none\"`.", call. = FALSE)
  }

  y <- as.numeric(object$SDF$y)
  X <- stats::model.matrix(ols_model)
  colnames(X) <- clean_coef_names(colnames(X))
  n <- length(y)
  p <- ncol(X)
  if (nrow(X) != n) {
    stop("The stored global model is not compatible with the LHWR fit.", call. = FALSE)
  }

  RSS_full <- RSS.lhwr(object)
  RSS_ols <- sum(stats::residuals(ols_model)^2)
  DF_ols <- stats::df.residual(ols_model)
  trS <- object$LHWR.diagnostic$tr.S
  trStS <- object$LHWR.diagnostic$tr.StS
  delta1 <- n - 2 * trS + trStS

  f_raw <- if (engine == "cpp") {
    hgwr_f_tests_cpp_components(
      object = object,
      y = y,
      X = X,
      RSS_full = RSS_full,
      RSS_ols = RSS_ols,
      DF_ols = DF_ols,
      trS = trS,
      trStS = trStS
    )
  } else {
    hgwr_f_tests_R_components(
      object = object,
      X = X,
      RSS_full = RSS_full,
      RSS_ols = RSS_ols,
      DF_ols = DF_ols,
      trS = trS,
      trStS = trStS
    )
  }

  global_tests <- data.frame(
    Test = c("F1", "F2"),
    Null_Hypothesis = c(
      "LHWR does not improve upon OLS",
      "Global OLS specification is adequate"
    ),
    F = c(f_raw$F1, f_raw$F2),
    df1 = c(f_raw$F1_df1, f_raw$F2_df1),
    df2 = c(f_raw$F1_df2, f_raw$F2_df2),
    p_value = c(f_raw$F1_p, f_raw$F2_p),
    check.names = FALSE
  )
  global_tests <- sanitize_f_tests(global_tests)
  coefficient_tests <- hgwr_f3_table(object$coef_names, f_raw)

  out <- list(
    global_tests = global_tests,
    coefficient_tests = coefficient_tests,
    diagnostics = list(
      n = n,
      p = p,
      RSS_ols = RSS_ols,
      RSS_full = RSS_full,
      trS = trS,
      trStS = trStS,
      delta1 = f_raw$L_delta1,
      delta2 = f_raw$L_delta2,
      delta1_hat = delta1,
      L_delta1 = f_raw$L_delta1,
      L_delta2 = f_raw$L_delta2,
      gamma1 = as.numeric(f_raw$gamma1),
      gamma2 = as.numeric(f_raw$gamma2),
      moment_method = f_raw$moment_method,
      bandwidth = object$bandwidth,
      kernel = object$LHWR.arguments$kernel,
      geometry = object$LHWR.arguments$geometry,
      distance = object$LHWR.arguments$distance,
      engine = engine
    ),
    digits = digits
  )
  class(out) <- "diagnostic_tests_lhwr"
  out
}

#' @export
print.diagnostic_tests_lhwr <- function(x, ...) {
  digits <- x$digits %||% 3L
  hgwr_rule()
  cat(" LHWR Diagnostic F-tests\n")
  hgwr_rule()
  cat("\nGlobal tests\n")
  hgwr_rule()
  global <- x$global_tests
  global_out <- data.frame(
    Test = global$Test,
    `Null Hypothesis` = global$Null_Hypothesis,
    F = format_f_number(global$F, digits),
    df1 = format_f_number(global$df1, digits),
    df2 = format_f_number(global$df2, digits),
    `p-value` = paste0(format_p_value(global$p_value, digits), signif_stars(global$p_value)),
    check.names = FALSE
  )
  print(global_out, row.names = FALSE, right = FALSE)

  cat("\nCoefficient-specific local variability tests\n")
  hgwr_rule()
  coef_tests <- x$coefficient_tests
  coef_out <- data.frame(
    Coefficient = coef_tests$Coefficient,
    F3 = format_f_number(coef_tests$F3, digits),
    df1 = format_f_number(coef_tests$df1, digits),
    df2 = format_f_number(coef_tests$df2, digits),
    `p-value` = paste0(format_p_value(coef_tests$p_value, digits), signif_stars(coef_tests$p_value)),
    check.names = FALSE
  )
  print(coef_out, row.names = FALSE, right = FALSE)
  hgwr_rule()
  cat(" Signif. codes:  *** p<0.001, ** p<0.01, * p<0.05\n")
  hgwr_rule()
  invisible(x)
}

hgwr_f_tests_cpp_components <- function(object, y, X, RSS_full, RSS_ols,
                                        DF_ols, trS, trStS) {
  if (identical(object$LHWR.arguments$engine, "knn")) {
    if (!is.null(object$knn_index) && !is.null(object$knn_distance)) {
      return(.hgwr_f_tests_knn_precomputed_cpp(
        y = y,
        X = X,
        knn_index = object$knn_index,
        knn_distance = object$knn_distance,
        kernel = object$LHWR.arguments$kernel,
        beta_hat = object$beta_hat,
        RSS_full = RSS_full,
        RSS_ols = RSS_ols,
        DF_ols = DF_ols,
        trS = trS,
        trStS = trStS
      ))
    }
    .hgwr_f_tests_knn_cpp(
      y = y,
      X = X,
      Z = object$Z,
      k = object$LHWR.arguments$bw,
      kernel = object$LHWR.arguments$kernel,
      beta_hat = object$beta_hat,
      RSS_full = RSS_full,
      RSS_ols = RSS_ols,
      DF_ols = DF_ols,
      trS = trS,
      trStS = trStS
    )
  } else {
    D <- hgwr_f_distance_matrix(object)
    .hgwr_f_tests_cpp(
      y = y,
      X = X,
      D = D,
      bw = object$bandwidth,
      kernel = object$LHWR.arguments$kernel,
      beta_hat = object$beta_hat,
      RSS_full = RSS_full,
      RSS_ols = RSS_ols,
      DF_ols = DF_ols,
      trS = trS,
      trStS = trStS
    )
  }
}

hgwr_f_tests_R_components <- function(object, X, RSS_full, RSS_ols,
                                      DF_ols, trS, trStS) {
  D <- hgwr_f_distance_matrix(object)
  bw <- object$bandwidth
  if (length(bw) == 1L) bw <- rep(bw, nrow(X))
  if (length(bw) != nrow(X)) {
    stop("Stored bandwidths are not compatible with the model size.", call. = FALSE)
  }
  kernel <- object$LHWR.arguments$kernel
  n <- nrow(X)
  p <- ncol(X)
  S <- matrix(0, n, n)
  B_ops <- lapply(seq_len(p), function(k) matrix(0, n, n))

  for (i in seq_len(n)) {
    w <- as.numeric(hgwr_kernel(D[, i], bw[i], kernel = kernel))
    Xw <- X * w
    XtWX <- crossprod(X, Xw)
    XtWX_inv <- tryCatch(solve(XtWX), error = function(e) MASS::ginv(XtWX))
    C <- XtWX_inv %*% t(Xw)
    s_row <- as.numeric(X[i, ] %*% C)
    S[i, ] <- s_row
    for (k in seq_len(p)) {
      B_ops[[k]][i, ] <- C[k, ]
    }
  }

  residual_operator <- diag(n) - S
  Q <- crossprod(residual_operator)
  q_diag <- diag(Q)
  L_delta1 <- sum(q_diag)
  L_delta2 <- sum(Q^2)
  delta1 <- L_delta1
  sigma2_delta1 <- RSS_full / L_delta1

  F1 <- (RSS_full / L_delta1) / (RSS_ols / DF_ols)
  F1_df1 <- L_delta1^2 / L_delta2
  F1_df2 <- DF_ols
  F1_p <- stats::pf(F1, F1_df1, F1_df2, lower.tail = TRUE)

  F2 <- ((RSS_ols - RSS_full) / (DF_ols - L_delta1)) / (RSS_ols / DF_ols)
  F2_df1 <- (DF_ols - L_delta1)^2 / (DF_ols - 2 * L_delta1 + L_delta2)
  F2_df2 <- DF_ols
  F2_p <- stats::pf(F2, F2_df1, F2_df2, lower.tail = FALSE)

  beta_hat <- as.matrix(object$beta_hat)
  F3 <- F3_df1 <- F3_df2 <- F3_p <- rep(NA_real_, p)
  gamma1 <- gamma2 <- rep(NA_real_, p)
  for (k in seq_len(p)) {
    beta_k <- beta_hat[, k]
    Vk2 <- (sum(beta_k^2) - sum(beta_k)^2 / n) / n
    centered_B <- sweep(B_ops[[k]], 2L, colMeans(B_ops[[k]]), FUN = "-")
    Ck <- crossprod(centered_B) / n
    gamma1[k] <- sum(diag(Ck))
    gamma2[k] <- sum(Ck^2)
    F3_df1[k] <- gamma1[k]^2 / gamma2[k]
    F3_df2[k] <- F1_df1
    F3[k] <- (Vk2 / gamma1[k]) / sigma2_delta1
    F3_p[k] <- stats::pf(F3[k], F3_df1[k], F3_df2[k], lower.tail = FALSE)
  }

  list(
    F1 = F1,
    F1_df1 = F1_df1,
    F1_df2 = F1_df2,
    F1_p = F1_p,
    F2 = F2,
    F2_df1 = F2_df1,
    F2_df2 = F2_df2,
    F2_p = F2_p,
    F3 = F3,
    F3_df1 = F3_df1,
    F3_df2 = F3_df2,
    F3_p = F3_p,
    q_diag = q_diag,
    L_delta1 = L_delta1,
    L_delta2 = L_delta2,
    delta1 = delta1,
    delta1_hat = n - 2 * trS + trStS,
    gamma1 = gamma1,
    gamma2 = gamma2,
    moment_method = "full_trace"
  )
}

hgwr_f_distance_matrix <- function(object) {
  if (!is.null(object$distance)) {
    return(object$distance)
  }
  if (identical(object$geometry, "tangent")) {
    return(tangent_distance_matrix(object$Z, mu = object$tangent_ref))
  }
  hyperbolic_distance_matrix(object$Z, check = FALSE)
}

hgwr_f3_table <- function(coefficient_names, f_raw) {
  out <- data.frame(
    Coefficient = coefficient_names,
    F3 = as.numeric(f_raw$F3),
    df1 = as.numeric(f_raw$F3_df1),
    df2 = as.numeric(f_raw$F3_df2),
    p_value = as.numeric(f_raw$F3_p),
    check.names = FALSE
  )
  sanitize_f_tests(out, statistic_col = "F3")
}

sanitize_f_tests <- function(x, statistic_col = "F") {
  bad <- !is.finite(x[[statistic_col]]) | x[[statistic_col]] < 0 |
    !is.finite(x$df1) | x$df1 <= 0 |
    !is.finite(x$df2) | x$df2 <= 0
  x[[statistic_col]][bad] <- NA_real_
  x$df1[bad] <- NA_real_
  x$df2[bad] <- NA_real_
  x$p_value[bad] <- NA_real_
  x$p_value[!is.na(x$p_value) & (x$p_value < 0 | x$p_value > 1)] <- NA_real_
  x
}

format_f_number <- function(x, digits) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}

format_p_value <- function(x, digits) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}

signif_stars <- function(p) {
  ifelse(
    is.na(p), "",
    ifelse(p < 0.001, " ***", ifelse(p < 0.01, " **", ifelse(p < 0.05, " *", "")))
  )
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' @export
diagnostic_tests.hgwr <- diagnostic_tests.lhwr
#' @export
print.diagnostic_tests_hgwr <- print.diagnostic_tests_lhwr

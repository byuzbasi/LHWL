test_that("LHWPR enforces the fixed alpha contract", {
  set.seed(100)
  H <- simulate_H2(24, r_max = 1.5)
  dat <- data.frame(y = rnorm(24), x = rnorm(24))

  expect_error(
    lhwpr(
      y ~ x, dat, H$Z, bw = 10,
      penalty = "mcp", alpha = 0.5, nlambda = 3
    ),
    "requires.*lasso"
  )
  expect_error(
    lhwpr(
      y ~ x, dat, H$Z, bw = 10,
      penalty = "scad", alpha = 0.5, nlambda = 3
    ),
    "requires.*lasso"
  )
  expect_error(
    lhwpr(
      y ~ x, dat, H$Z, bw = 10,
      penalty = "lasso", alpha = 1.1, nlambda = 3
    ),
    "\\[0, 1\\]"
  )
  expect_error(
    lhwpr(
      y ~ x, dat, H$Z, bw = 10,
      family = "binomial",
      control = lhwpr_control(solver = "covariance")
    ),
    "only for Gaussian"
  )
  expect_error(
    lhwpr_control(early_stop = TRUE),
    "reserved"
  )

  ridge <- lhwpr(
    y ~ x, dat, H$Z, bw = 10,
    penalty = "scad", alpha = 0, lambda = 0.2
  )
  expect_equal(ridge$effective_penalty, "ridge")
})

test_that("Gaussian ridge agrees with the local closed-form solution", {
  set.seed(101)
  n <- 32
  H <- simulate_H2(n, r_max = 2)
  dat <- data.frame(
    y = 0.7 + 1.3 * rnorm(n) + rnorm(n, sd = 0.2),
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  dat$y <- 0.7 + 1.3 * dat$x1 - 0.5 * dat$x2 + rnorm(n, sd = 0.2)
  lambda <- 0.35
  fit <- lhwpr(
    y ~ x1 + x2,
    dat,
    H$Z,
    bw = 15,
    alpha = 0,
    lambda = lambda,
    control = lhwpr_control(solver = "naive", tol = 1e-10, kkt_tol = 1e-9)
  )

  focus <- 7L
  D <- hyperbolic_distance_matrix(H$Z, check = FALSE)
  w <- hgwr_kernel(
    D[, focus], fit$bandwidth[focus], kernel = fit$kernel
  )
  w <- w / sum(w)
  A <- cbind(Intercept = 1, fit$X_standardized)
  penalty <- diag(c(0, rep(lambda, ncol(fit$X_standardized))))
  reference <- solve(
    crossprod(A, A * w) + penalty,
    crossprod(A, w * fit$y)
  )

  expect_equal(
    as.numeric(fit$beta_hat_standardized[focus, ]),
    as.numeric(reference),
    tolerance = 2e-7
  )
  expect_true(all(fit$converged[, fit$selected_index] != 0))
})

test_that("naive and covariance Gaussian solvers agree", {
  set.seed(102)
  n <- 38
  H <- simulate_H2(n, r_max = 2)
  dat <- data.frame(
    y = rnorm(n),
    x1 = rnorm(n),
    x2 = rnorm(n),
    x3 = rnorm(n)
  )
  lambda <- c(0.5, 0.15, 0.04)
  naive <- lhwpr(
    y ~ x1 + x2 + x3,
    dat,
    H$Z,
    bw = 17,
    lambda = lambda,
    alpha = 0.6,
    control = lhwpr_control(
      solver = "naive", tol = 1e-9, kkt_tol = 1e-8
    )
  )
  covariance <- lhwpr(
    y ~ x1 + x2 + x3,
    dat,
    H$Z,
    bw = 17,
    lambda = lambda,
    alpha = 0.6,
    control = lhwpr_control(
      solver = "covariance", tol = 1e-9, kkt_tol = 1e-8
    )
  )

  expect_equal(
    covariance$beta_path_standardized,
    naive$beta_path_standardized,
    tolerance = 2e-7
  )
  expect_equal(covariance$fitted_path, naive$fitted_path, tolerance = 2e-7)
})

test_that("lasso and elastic-net paths satisfy KKT conditions", {
  set.seed(103)
  n <- 42
  H <- simulate_H2(n, r_max = 2)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  dat <- data.frame(
    y = 1 + 1.5 * x1 + rnorm(n, sd = 0.3),
    x1 = x1,
    x2 = x2
  )
  lasso <- lhwpr(
    y ~ x1 + x2, dat, H$Z,
    bw = 19, alpha = 1, nlambda = 7
  )
  elastic <- lhwpr(
    y ~ x1 + x2, dat, H$Z,
    bw = 19, alpha = 0.4, nlambda = 7
  )

  expect_true(all(abs(lasso$beta_path_standardized[, -1, 1]) < 2e-7))
  expect_true(all(abs(elastic$beta_path_standardized[, -1, 1]) < 2e-7))
  expect_lt(max(lasso$kkt), 3e-7)
  expect_lt(max(elastic$kkt), 3e-7)
  expect_true(all(lasso$converged != 0))
  expect_true(all(elastic$converged != 0))
})

test_that("Gaussian MCP and SCAD paths are stationary", {
  set.seed(104)
  n <- 40
  H <- simulate_H2(n, r_max = 2)
  x1 <- rnorm(n)
  dat <- data.frame(
    y = 0.5 + 2 * x1 + rnorm(n, sd = 0.25),
    x1 = x1,
    x2 = rnorm(n)
  )

  for (penalty in c("mcp", "scad")) {
    fit <- lhwpr(
      y ~ x1 + x2,
      dat,
      H$Z,
      bw = 18,
      alpha = 1,
      penalty = penalty,
      nlambda = 6,
      control = lhwpr_control(tol = 1e-8, kkt_tol = 1e-7)
    )
    expect_true(all(is.finite(fit$beta_hat)))
    expect_lt(max(fit$kkt), 5e-7)
    expect_true(all(fit$converged != 0))
  }
})

test_that("binomial LHWPR returns stable probabilities", {
  set.seed(105)
  n <- 55
  H <- simulate_H2(n, r_max = 2)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  probability <- plogis(-0.3 + 1.2 * x1 - 0.8 * x2)
  dat <- data.frame(
    y = factor(rbinom(n, 1, probability), levels = c(0, 1)),
    x1 = x1,
    x2 = x2
  )
  fit <- lhwpr(
    y ~ x1 + x2,
    dat,
    H$Z,
    bw = 28,
    family = "binomial",
    alpha = 0.7,
    nlambda = 6
  )

  expect_true(all(fit$fitted_path > 0 & fit$fitted_path < 1))
  expect_true(all(is.finite(fit$linear_predictor_path)))
  expect_true(all(fit$converged != 0))
  expect_lt(max(fit$kkt), 5e-7)
  expect_equal(fitted(fit), fit$y_hat)
  expect_equal(residuals(fit), fit$y - fit$y_hat)

  for (specification in list(
    list(alpha = 0, penalty = "scad"),
    list(alpha = 1, penalty = "lasso"),
    list(alpha = 1, penalty = "mcp"),
    list(alpha = 1, penalty = "scad")
  )) {
    candidate <- lhwpr(
      y ~ x1 + x2,
      dat,
      H$Z,
      bw = 28,
      family = "binomial",
      alpha = specification$alpha,
      penalty = specification$penalty,
      nlambda = 4
    )
    expect_true(all(candidate$converged != 0))
    expect_lt(max(candidate$kkt), 5e-7)
  }
})

test_that("global standardization and zero penalty factors are preserved", {
  set.seed(106)
  n <- 36
  H <- simulate_H2(n, r_max = 2)
  dat <- data.frame(
    y = rnorm(n),
    x1 = rnorm(n, mean = 5, sd = 2),
    x2 = rnorm(n, mean = -3, sd = 4)
  )
  case_weights <- seq_len(n)
  fit <- lhwpr(
    y ~ x1 + x2,
    dat,
    H$Z,
    bw = 16,
    lambda = 1e3,
    penalty_factor = c(x1 = 0, x2 = 1),
    case_weights = case_weights
  )
  normalized <- case_weights / sum(case_weights)

  expect_equal(
    fit$standardization$center,
    colSums(fit$X * normalized)
  )
  expect_gt(max(abs(fit$beta_hat[, "x1"])), 1e-6)
  expect_lt(max(abs(fit$beta_hat[, "x2"])), 1e-6)

  no_intercept <- lhwpr(
    y ~ 0 + x1 + x2,
    dat,
    H$Z,
    bw = 16,
    lambda = 0.1
  )
  expect_equal(no_intercept$standardization$center, c(0, 0))
  expect_false("Intercept" %in% colnames(no_intercept$beta_hat))
})

test_that("dense and kNN LHWPR engines agree for compact kernels", {
  set.seed(107)
  n <- 44
  H <- simulate_H2(n, r_max = 2)
  dat <- data.frame(
    y = rnorm(n),
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  dense <- lhwpr(
    y ~ x1 + x2,
    dat,
    H$Z,
    bw = 18,
    lambda = c(0.3, 0.08),
    engine = "dense"
  )
  knn <- lhwpr(
    y ~ x1 + x2,
    dat,
    H$Z,
    bw = 18,
    lambda = c(0.3, 0.08),
    engine = "knn"
  )
  threaded <- lhwpr(
    y ~ x1 + x2,
    dat,
    H$Z,
    bw = 18,
    lambda = c(0.3, 0.08),
    engine = "dense",
    control = lhwpr_control(n_threads = 2)
  )

  expect_equal(knn$beta_path, dense$beta_path, tolerance = 1e-9)
  expect_equal(knn$fitted_path, dense$fitted_path, tolerance = 1e-9)
  expect_identical(threaded$beta_path, dense$beta_path)
  expect_identical(threaded$fitted_path, dense$fitted_path)
})

test_that("path retention and prediction methods use the selected lambda", {
  set.seed(108)
  n <- 34
  H <- simulate_H2(n, r_max = 2)
  dat <- data.frame(
    y = rnorm(n),
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  fit <- lhwpr(
    y ~ x1 + x2,
    dat,
    H$Z,
    bw = 15,
    lambda = c(0.3, 0.1, 0.03),
    keep_path = "selected"
  )

  expect_null(fit$beta_path)
  expect_equal(coef(fit), fit$beta_hat)
  expect_error(coef(fit, s = 0.3), "not retained")
  predicted <- predict(
    fit,
    newdata = dat[1:5, , drop = FALSE],
    Znew = H$Z[1:5, , drop = FALSE]
  )
  expect_equal(predicted, fit$y_hat[1:5], tolerance = 5e-8)
})

test_that("joint tuning uses reproducible folds and returns a refit", {
  set.seed(109)
  n <- 33
  H <- simulate_H2(n, r_max = 2)
  x <- rnorm(n)
  dat <- data.frame(y = 1 + x + rnorm(n, sd = 0.3), x = x)
  tuned <- tune.lhwpr(
    y ~ x,
    dat,
    H$Z,
    bw_grid = c(13, 17),
    alpha_grid = c(0, 0.5, 1),
    nlambda = 4,
    resampling = "kfold",
    nfolds = 3,
    rule = "1se"
  )

  expect_s3_class(tuned, "tune.lhwpr")
  expect_s3_class(tuned$fit, "lhwpr")
  expect_equal(tuned$fold_id, rep(1:3, length.out = n))
  expect_equal(sum(tuned$selection_table$selected), 1)
  expect_equal(tuned$fit$selected_lambda, tuned$selected$lambda)
  expect_identical(
    tuned$selected_lambda_path[[length(tuned$selected_lambda_path)]],
    as.numeric(tuned$selected$lambda)
  )
  expect_true(length(tuned$selected_lambda_path) >= 1L)
  expect_true(tuned$refit_diagnostics$passed)
  expect_equal(tuned$refit_diagnostics$convergence_rate, 1)
  expect_lte(
    tuned$refit_diagnostics$kkt_max,
    tuned$refit_diagnostics$kkt_limit
  )
  expect_true(all(c(
    "cv_mean", "cv_se", "rss", "deviance", "logLik", "trS", "edf",
    "aic", "aicc", "bic", "gcv", "ebic", "ebic_refit",
    "ebic_model_space_log", "ebic_model_space_zero", "edf_valid",
    "post_residual_df", "post_edf_valid", "surface_size",
    "convergence_rate", "eligible", "ineligibility_reason", "bw_boundary"
  ) %in% names(tuned$selection_table)))
  expect_true(tuned$selection_table$eligible[tuned$selected_index])
  expect_equal(nrow(tuned$lhwr_comparison), 2)

  expect_error(
    tune.lhwpr(
      y ~ x, dat, H$Z,
      bw_grid = 15,
      nlambda = 2,
      resampling = "blocked"
    ),
    "fold_id"
  )
})

test_that("joint tuning generates a safe automatic adaptive-kNN grid", {
  default_grid <- lhwpr_tune_bandwidths(
    bw_grid = NULL,
    n = 120,
    n_coef = 9,
    adaptive = TRUE,
    bandwidth_type = "knn",
    bw_lower = NULL,
    bw_upper = NULL,
    n_bw = 15
  )
  expect_identical(default_grid$mode, "automatic")
  expect_equal(default_grid$lower, 20)
  expect_equal(default_grid$upper, 119)
  expect_equal(length(default_grid$grid), 15)

  set.seed(111)
  n <- 42
  H <- simulate_H2(n, r_max = 2)
  dat <- data.frame(
    y = rnorm(n),
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  tuned <- tune.lhwpr(
    y ~ x1 + x2,
    dat,
    H$Z,
    bw_grid = NULL,
    bw_lower = 18,
    bw_upper = 30,
    n_bw = 4,
    nlambda = 3,
    resampling = "kfold",
    nfolds = 2,
    refine = FALSE
  )

  expect_identical(tuned$bandwidth_search$mode, "automatic")
  expect_equal(tuned$bandwidth_search$lower, 18)
  expect_equal(tuned$bandwidth_search$upper, 30)
  expect_equal(tuned$bandwidth_search$grid, c(18, 22, 26, 30))
  expect_equal(sort(unique(tuned$selection_table$bw)), c(18, 22, 26, 30))

  expect_error(
    tune.lhwpr(
      y ~ x1 + x2,
      dat,
      H$Z,
      bw_grid = NULL,
      adaptive = FALSE,
      nlambda = 2
    ),
    "adaptive kNN"
  )
})

test_that("information-criterion selection excludes invalid EDF rows", {
  selection_table <- data.frame(
    cv_mean = c(2, 1),
    cv_se = c(0.1, 0.1),
    aic = c(-100, 8),
    aicc = c(-100, 9),
    bic = c(-100, 10),
    gcv = c(-100, 11),
    ebic = c(-100, 12),
    ebic_refit = c(-200, 13),
    edf_valid = c(FALSE, TRUE),
    post_edf_valid = c(FALSE, TRUE),
    convergence_rate = c(1, 1),
    kkt_max = c(1e-8, 1e-8),
    cv_convergence_rate = c(1, 1),
    cv_kkt_max = c(1e-8, 1e-8),
    surface_size = c(3L, 1L),
    mean_active = c(3, 1),
    lambda_ratio = c(1e-4, 0.1),
    bw = c(12, 20),
    alpha = 1
  )

  for (criterion in c("aic", "aicc", "bic", "gcv", "ebic")) {
    expect_identical(
      lhwpr_select_tuning_row(selection_table, criterion, "min"),
      2L
    )
  }
  expect_identical(
    lhwpr_select_tuning_row(selection_table, "ebic_refit", "min"),
    2L
  )
  expect_identical(
    lhwpr_select_tuning_row(selection_table, "cv", "min"),
    2L
  )
})

test_that("all tuning rules preserve convergence and KKT eligibility", {
  selection_table <- data.frame(
    cv_mean = c(0.90, 1.00, 1.02, 1.03),
    cv_se = c(0.20, 0.20, 0.20, 0.20),
    convergence_rate = c(1, 1, 1, 0.99),
    kkt_max = c(2e-5, 1e-8, 1e-8, 1e-8),
    cv_convergence_rate = c(1, 1, 1, 1),
    cv_kkt_max = c(1e-8, 1e-8, 2e-5, 1e-8),
    surface_size = c(1L, 3L, 1L, 1L),
    mean_active = c(1, 3, 1, 1),
    lambda_ratio = c(0.9, 0.2, 0.8, 0.7),
    bw = c(14, 18, 22, 26),
    alpha = 1
  )

  eligibility <- lhwpr_tuning_eligibility(
    selection_table, "cv", kkt_tol = 1e-5
  )
  expect_identical(eligibility$eligible, c(FALSE, TRUE, FALSE, FALSE))
  for (rule in c("min", "1se", "sparse_1se")) {
    expect_identical(
      lhwpr_select_tuning_row(
        selection_table, "cv", rule, kkt_tol = 1e-5
      ),
      2L
    )
  }
})

test_that("selected refit uses the exact decreasing tuning prefix", {
  selection_table <- data.frame(
    configuration_key = c(rep("20|1|NA", 4), rep("30|1|NA", 2)),
    bw = c(rep(20, 4), rep(30, 2)),
    alpha = 1,
    gamma = NA_real_,
    lambda = c(0.8, 0.3, 0.1, 0.03, 0.7, 0.2),
    stringsAsFactors = FALSE
  )
  selected <- selection_table[3, , drop = FALSE]
  path <- lhwpr_selected_lambda_path(selection_table, selected)
  expect_identical(path, c(0.8, 0.3, 0.1))
  expect_true(all(diff(path) < 0))
  expect_identical(path[[length(path)]], selected$lambda)

  expect_error(
    lhwpr_selected_lambda_path(
      selection_table,
      transform(selected, lambda = 0.11)
    ),
    "absent"
  )
})

test_that("Gaussian tuning criteria and EBIC follow their documented formulas", {
  set.seed(112)
  n <- 40
  H <- simulate_H2(n, r_max = 1.8)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rnorm(n)
  dat <- data.frame(
    y = 0.5 + 1.4 * x1 - 0.7 * x2 + rnorm(n, sd = 0.35),
    x1 = x1,
    x2 = x2,
    x3 = x3
  )
  tuned <- tune.lhwpr(
    y ~ x1 + x2 + x3,
    dat,
    H$Z,
    bw_grid = c(22, 30),
    alpha_grid = c(0.5, 1),
    nlambda = 4,
    nfolds = 2,
    criterion = "ebic",
    ebic_gamma = 0.5,
    refine = FALSE,
    refit = FALSE
  )

  expect_true(is.finite(tuned$selected$ebic))
  expect_true(tuned$selected$edf_valid)
  finite_aicc <- which(is.finite(tuned$selection_table$aicc))[1L]
  row <- tuned$selection_table[finite_aicc, , drop = FALSE]
  expect_equal(
    row$aicc,
    row$aic + 2 * row$edf * (row$edf + 1) /
      (n - row$edf - 1),
    tolerance = 1e-10
  )

  selected <- tuned$selected
  expected_ebic <- selected$bic +
    2 * 0.5 * lchoose(3, selected$surface_size)
  expect_equal(selected$ebic, expected_ebic, tolerance = 1e-9)
  expect_equal(
    selected$ebic_model_space_log,
    lchoose(3, selected$surface_size),
    tolerance = 1e-12
  )
  expect_identical(
    selected$ebic_model_space_zero,
    selected$surface_size > 0 &&
      lchoose(3, selected$surface_size) == 0
  )
  expect_equal(
    tuned$selection_table$post_residual_df,
    n - tuned$selection_table$post_edf,
    tolerance = 1e-12
  )
  expected_post_valid <- is.finite(tuned$selection_table$post_edf) &
    tuned$selection_table$post_edf >= 1 &
    tuned$selection_table$post_residual_df > 1
  expect_identical(
    tuned$selection_table$post_edf_valid,
    expected_post_valid
  )
  valid <- which(expected_post_valid)[1L]
  expect_false(is.na(valid))
  expected_refit <- tuned$selection_table$post_deviance[[valid]] +
    log(n) * tuned$selection_table$post_edf[[valid]] +
    2 * 0.5 * lchoose(3, tuned$selection_table$surface_size[[valid]])
  expect_equal(
    tuned$selection_table$ebic_refit[[valid]],
    expected_refit,
    tolerance = 1e-9
  )
  expect_true(all(is.infinite(
    tuned$selection_table$ebic_refit[!expected_post_valid]
  )))
  expect_true(all(
    tuned$selection_table$surface_size >= 0 &
      tuned$selection_table$surface_size <= 3
  ))

  tuned_refit <- tune.lhwpr(
    y ~ x1 + x2 + x3,
    dat,
    H$Z,
    bw_grid = c(22, 30),
    alpha_grid = c(0.5, 1),
    nlambda = 4,
    nfolds = 2,
    criterion = "ebic_refit",
    ebic_gamma = 0.5,
    refine = FALSE,
    refit = FALSE
  )
  expect_true(is.finite(tuned_refit$selected$ebic_refit))
  expect_true(tuned_refit$selected$post_edf_valid)

  output <- tempfile(fileext = ".pdf")
  grDevices::pdf(output)
  expect_silent(plot(tuned, type = "all"))
  grDevices::dev.off()
  expect_true(file.exists(output))
})

test_that("ridge path EDF equals the exact conditional smoother trace", {
  set.seed(114)
  n <- 30
  H <- simulate_H2(n, r_max = 1.7)
  dat <- data.frame(
    y = rnorm(n),
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  lambda <- 0.25
  tuned <- tune.lhwpr(
    y ~ x1 + x2,
    dat,
    H$Z,
    bw_grid = 18,
    alpha_grid = 0,
    lambda = lambda,
    nfolds = 2,
    refine = FALSE,
    compare_lhwr = FALSE,
    refit = FALSE
  )
  fit <- lhwpr(
    y ~ x1 + x2,
    dat,
    H$Z,
    bw = 18,
    alpha = 0,
    lambda = lambda,
    keep_path = "selected"
  )
  distance <- hyperbolic_distance_matrix(H$Z, check = FALSE)
  design <- cbind(1, fit$X_standardized)
  penalty <- diag(c(0, rep(lambda, ncol(fit$X_standardized))))
  leverage <- numeric(n)
  for (focus in seq_len(n)) {
    weight <- hgwr_kernel(
      distance[, focus],
      fit$bandwidth[focus],
      kernel = fit$kernel
    )
    weight <- weight / sum(weight)
    hessian <- crossprod(design, design * weight) + penalty
    leverage[focus] <- weight[focus] * drop(
      design[focus, , drop = FALSE] %*%
        solve(hessian, design[focus, ])
    )
  }
  expect_equal(
    tuned$selection_table$trS,
    sum(leverage),
    tolerance = 2e-8
  )
  expect_equal(
    tuned$selection_table$edf,
    sum(leverage) + 1,
    tolerance = 2e-8
  )
})

test_that("tuning refinement and family guards are explicit", {
  search <- lhwpr_tune_bandwidths(
    bw_grid = NULL,
    n = 80,
    n_coef = 5,
    adaptive = TRUE,
    bandwidth_type = "knn",
    bw_lower = NULL,
    bw_upper = NULL,
    n_bw = 5
  )
  refined <- lhwpr_tune_refined_bandwidths(
    selected_bw = search$grid[1L],
    evaluated_bw = search$grid,
    bandwidth_search = search,
    adaptive = TRUE,
    bandwidth_type = "knn",
    n_refine = 7,
    expand_bounds = TRUE
  )
  expect_true(length(refined) > 0)
  expect_true(min(refined) >= search$hard_lower)
  expect_true(any(refined < search$lower))
  expect_equal(
    lhwpr_tune_refined_alpha(0.5, c(0, 0.5, 1)),
    c(0.25, 0.75)
  )

  set.seed(113)
  n <- 30
  H <- simulate_H2(n, r_max = 1.5)
  x <- rnorm(n)
  dat <- data.frame(
    y = rbinom(n, 1, stats::plogis(x)),
    x = x
  )
  expect_error(
    tune.lhwpr(
      y ~ x,
      dat,
      H$Z,
      bw_grid = c(16, 22),
      family = "binomial",
      criterion = "aicc",
      nlambda = 2
    ),
    "only for Gaussian"
  )

  binomial_tuned <- tune.lhwpr(
    y ~ x,
    dat,
    H$Z,
    bw_grid = 22,
    family = "binomial",
    alpha_grid = 1,
    nlambda = 2,
    nfolds = 2,
    refine = FALSE,
    compare_lhwr = FALSE,
    refit = FALSE
  )
  expect_true(all(is.finite(
    binomial_tuned$selection_table$surface_size
  )))
  expect_true(all(is.finite(
    binomial_tuned$selection_table$deviance
  )))
  expect_true(all(is.na(binomial_tuned$selection_table$aicc)))
})

test_that("adding LHWPR does not change existing LHWR fits", {
  set.seed(110)
  n <- 30
  H <- simulate_H2(n, r_max = 2)
  dat <- data.frame(y = rnorm(n), x = rnorm(n))
  before <- lhwr(
    y ~ x, dat, H$Z, bw = 14,
    diagnostics = "none"
  )
  invisible(lhwpr(
    y ~ x, dat, H$Z, bw = 14,
    lambda = 0.1
  ))
  after <- lhwr(
    y ~ x, dat, H$Z, bw = 14,
    diagnostics = "none"
  )

  expect_identical(after$beta_hat, before$beta_hat)
  expect_identical(after$y_hat, before$y_hat)
})

test_that("Gaussian LHWPR lasso is stable when p exceeds n", {
  set.seed(20260906)
  n <- 28L
  p <- 36L
  H <- simulate_H2(n, r_max = 2)
  X <- replicate(p, rnorm(n))
  colnames(X) <- sprintf("x%02d", seq_len(p))
  beta <- c(1.4, -1, 0.7, rep(0, p - 3L))
  dat <- data.frame(
    y = 0.5 + drop(X %*% beta) + rnorm(n, sd = 0.25),
    X,
    check.names = FALSE
  )

  fit <- lhwpr(
    reformulate(colnames(X), response = "y"),
    dat,
    H$Z,
    bw = 20,
    penalty = "lasso",
    alpha = 1,
    lambda = c(0.8, 0.25, 0.08),
    keep_path = "selected",
    control = lhwpr_control(
      n_threads = 1,
      tol = 1e-8,
      kkt_tol = 1e-7
    )
  )

  expect_gt(p, n)
  expect_equal(dim(fit$beta_hat), c(n, p + 1L))
  expect_true(all(is.finite(fit$beta_hat)))
  expect_true(all(is.finite(fit$y_hat)))
  expect_true(all(fit$converged[, fit$selected_index] != 0))
  expect_lt(max(fit$kkt[, fit$selected_index]), 5e-7)
})

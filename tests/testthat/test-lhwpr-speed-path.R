test_that("fold tuning solves only validation focal locations without changing results", {
  set.seed(20260820)
  n <- 35L
  geometry <- simulate_H2(n, r_max = 1.7)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rnorm(n)
  dat <- data.frame(
    y = 0.4 + 1.2 * x1 - 0.6 * x2 + rnorm(n, sd = 0.35),
    x1 = x1,
    x2 = x2,
    x3 = x3
  )
  folds <- rep(c(1L, 1L, 2L, 3L), length.out = n)
  lambda <- c(0.45, 0.16, 0.05)
  control <- lhwpr_control(
    solver = "naive", tol = 1e-9, kkt_tol = 1e-8,
    n_threads = 1L
  )

  for (engine in c("dense", "knn")) {
    prepared <- lhwpr_prepare(
      y ~ x1 + x2 + x3,
      dat,
      geometry$Z,
      family = "gaussian",
      penalty_factor = NULL,
      case_weights = NULL,
      standardize = TRUE,
      geometry = "lorentz",
      tangent_ref = NULL,
      dMat = NULL,
      engine = engine
    )
    neighborhood <- lhwpr_neighborhood(
      prepared,
      bw = 18L,
      kernel = "bisquare",
      adaptive = TRUE,
      bandwidth_type = "knn",
      engine = engine
    )

    reference_prediction <- matrix(NA_real_, n, length(lambda))
    reference_converged <- matrix(0L, n, length(lambda))
    reference_kkt <- matrix(Inf, n, length(lambda))
    for (fold in sort(unique(folds))) {
      validation <- which(folds == fold)
      training_weight <- prepared$case_weights
      training_weight[validation] <- 0
      reference <- lhwpr_native_fit(
        prepared = prepared,
        neighborhood = neighborhood,
        lambda = lambda,
        family = "gaussian",
        penalty = "lasso",
        alpha = 1,
        gamma = NA_real_,
        nlambda = length(lambda),
        lambda_min_ratio = min(lambda) / max(lambda),
        control = control,
        store_path = FALSE,
        exclude_self = FALSE,
        case_weights = training_weight,
        compute_lambda_diagnostics = FALSE
      )
      reference_prediction[validation, ] <- reference$fitted_path[
        validation, , drop = FALSE
      ]
      reference_converged[validation, ] <- reference$converged[
        validation, , drop = FALSE
      ]
      reference_kkt[validation, ] <- reference$kkt[
        validation, , drop = FALSE
      ]
    }

    optimized <- lhwpr_tune_predictions(
      prepared = prepared,
      neighborhood = neighborhood,
      family = "gaussian",
      penalty = "lasso",
      alpha = 1,
      gamma = NA_real_,
      lambda = lambda,
      lambda_min_ratio = min(lambda) / max(lambda),
      control = control,
      resampling = "kfold",
      folds = folds
    )

    expect_equal(optimized$prediction, reference_prediction, tolerance = 1e-12)
    expect_identical(optimized$converged, reference_converged)
    expect_equal(optimized$kkt, reference_kkt, tolerance = 1e-12)
  }
})

test_that("explicit-lambda internal fits can skip unused lambda diagnostics", {
  set.seed(20260821)
  n <- 30L
  geometry <- simulate_H2(n, r_max = 1.5)
  dat <- data.frame(y = rnorm(n), x1 = rnorm(n), x2 = rnorm(n))
  prepared <- lhwpr_prepare(
    y ~ x1 + x2,
    dat,
    geometry$Z,
    family = "gaussian",
    penalty_factor = NULL,
    case_weights = NULL,
    standardize = TRUE,
    geometry = "lorentz",
    tangent_ref = NULL,
    dMat = NULL,
    engine = "dense"
  )
  neighborhood <- lhwpr_neighborhood(
    prepared,
    bw = 16L,
    kernel = "bisquare",
    adaptive = TRUE,
    bandwidth_type = "knn",
    engine = "dense"
  )
  control <- lhwpr_control(solver = "naive", n_threads = 1L)
  arguments <- list(
    prepared = prepared,
    neighborhood = neighborhood,
    lambda = c(0.3, 0.08),
    family = "gaussian",
    penalty = "lasso",
    alpha = 1,
    gamma = NA_real_,
    nlambda = 2L,
    lambda_min_ratio = 0.08 / 0.3,
    control = control,
    store_path = FALSE,
    exclude_self = FALSE
  )
  with_diagnostics <- do.call(
    lhwpr_native_fit,
    c(arguments, list(compute_lambda_diagnostics = TRUE))
  )
  without_diagnostics <- do.call(
    lhwpr_native_fit,
    c(arguments, list(compute_lambda_diagnostics = FALSE))
  )

  expect_identical(without_diagnostics$fitted_path, with_diagnostics$fitted_path)
  expect_identical(without_diagnostics$converged, with_diagnostics$converged)
  expect_identical(without_diagnostics$kkt, with_diagnostics$kkt)
  expect_true(is.na(without_diagnostics$lambda_max))
  expect_true(all(is.na(without_diagnostics$lambda_local_max)))
  expect_true(all(is.na(without_diagnostics$lambda_ratio)))

  retained_path <- do.call(
    lhwpr_native_fit,
    c(
      arguments[c(
        "prepared", "neighborhood", "lambda", "family", "penalty",
        "alpha", "gamma", "nlambda", "lambda_min_ratio", "control",
        "exclude_self"
      )],
      list(store_path = TRUE, compute_lambda_diagnostics = FALSE)
    )
  )
  expect_identical(
    without_diagnostics$beta_hat_standardized,
    retained_path$beta_hat_standardized
  )
  expect_identical(without_diagnostics$fitted_path, retained_path$fitted_path)
  expect_identical(without_diagnostics$converged, retained_path$converged)
  expect_identical(without_diagnostics$kkt, retained_path$kkt)
})

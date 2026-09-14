test_that("diagnostic_tests returns global and coefficient tables", {
  set.seed(20260620)
  H <- simulate_H2(45, r_max = 2.5)
  x1 <- rnorm(45)
  x2 <- rnorm(45)
  y <- 1 + 2 * x1 - 0.5 * x2 + rnorm(45, sd = 0.25)
  dat <- data.frame(y = y, x1 = x1, x2 = x2)

  fit <- lhwr(
    y ~ x1 + x2,
    data = dat,
    Z = H$Z,
    bw = 20,
    adaptive = TRUE,
    diagnostics = "aicc",
    f_tests = FALSE
  )
  out <- diagnostic_tests(fit)

  expect_s3_class(out, "diagnostic_tests_lhwr")
  expect_true(is.data.frame(out$global_tests))
  expect_true(is.data.frame(out$coefficient_tests))
  expect_equal(out$global_tests$Test, c("F1", "F2"))
  expect_equal(
    out$global_tests$Null_Hypothesis,
    c("LHWR does not improve upon OLS", "Global OLS specification is adequate")
  )
  expect_equal(out$coefficient_tests$Coefficient, fit$coef_names)
})

test_that("lhwr can attach F-tests when requested", {
  set.seed(20260621)
  H <- simulate_H2(40, r_max = 2.5)
  x1 <- rnorm(40)
  y <- 1 + 2 * x1 + rnorm(40, sd = 0.25)
  dat <- data.frame(y = y, x1 = x1)

  fit <- lhwr(
    y ~ x1,
    data = dat,
    Z = H$Z,
    bw = 18,
    adaptive = TRUE,
    diagnostics = "aicc",
    f_tests = TRUE
  )

  expect_s3_class(fit$Ftests, "diagnostic_tests_lhwr")
  pvals <- c(fit$Ftests$global_tests$p_value, fit$Ftests$coefficient_tests$p_value)
  expect_true(all(is.na(pvals) | (pvals >= 0 & pvals <= 1)))
})

test_that("diagnostic_tests print output includes expected headers", {
  set.seed(20260622)
  H <- simulate_H2(35, r_max = 2)
  x1 <- rnorm(35)
  y <- 1 + 2 * x1 + rnorm(35, sd = 0.25)
  dat <- data.frame(y = y, x1 = x1)

  fit <- lhwr(
    y ~ x1,
    data = dat,
    Z = H$Z,
    bw = 16,
    adaptive = TRUE,
    diagnostics = "aicc"
  )
  out <- diagnostic_tests(fit)

  printed <- capture.output(print(out))
  expect_true(any(grepl("LHWR Diagnostic F-tests", printed, fixed = TRUE)))
  expect_true(any(grepl("Global tests", printed, fixed = TRUE)))
  expect_true(any(grepl("Coefficient-specific local variability tests", printed, fixed = TRUE)))
})

test_that("diagnostic_tests C++ engine matches R reference engine", {
  set.seed(20260623)
  H <- simulate_H2(35, r_max = 2)
  x1 <- rnorm(35)
  x2 <- rnorm(35)
  y <- 1 + 2 * x1 - 0.5 * x2 + rnorm(35, sd = 0.25)
  dat <- data.frame(y = y, x1 = x1, x2 = x2)

  fit <- lhwr(
    y ~ x1 + x2,
    data = dat,
    Z = H$Z,
    bw = 16,
    adaptive = TRUE,
    diagnostics = "aicc"
  )
  out_cpp <- diagnostic_tests(fit, engine = "cpp")
  out_R <- diagnostic_tests(fit, engine = "R")

  expect_equal(out_cpp$global_tests$F, out_R$global_tests$F, tolerance = 1e-8)
  expect_equal(out_cpp$global_tests$p_value, out_R$global_tests$p_value, tolerance = 1e-8)
  expect_equal(out_cpp$coefficient_tests$F3, out_R$coefficient_tests$F3, tolerance = 1e-8)
  expect_equal(out_cpp$coefficient_tests$p_value, out_R$coefficient_tests$p_value, tolerance = 1e-8)
  expect_equal(out_cpp$diagnostics$L_delta2, out_R$diagnostics$L_delta2, tolerance = 1e-8)
  expect_equal(out_cpp$diagnostics$gamma1, out_R$diagnostics$gamma1, tolerance = 1e-8)
  expect_equal(out_cpp$diagnostics$gamma2, out_R$diagnostics$gamma2, tolerance = 1e-8)
})

test_that("diagnostic_tests kNN C++ engine matches dense C++ engine", {
  set.seed(20260624)
  H <- simulate_H2(45, r_max = 2.5)
  x1 <- rnorm(45)
  x2 <- rnorm(45)
  y <- 1 + 2 * x1 - 0.5 * x2 + rnorm(45, sd = 0.25)
  dat <- data.frame(y = y, x1 = x1, x2 = x2)

  fit_dense <- lhwr(
    y ~ x1 + x2,
    data = dat,
    Z = H$Z,
    bw = 18,
    adaptive = TRUE,
    kernel = "bisquare",
    diagnostics = "aicc",
    engine = "dense"
  )
  fit_knn <- lhwr(
    y ~ x1 + x2,
    data = dat,
    Z = H$Z,
    bw = 18,
    adaptive = TRUE,
    kernel = "bisquare",
    diagnostics = "aicc",
    engine = "knn"
  )

  expect_equal(dim(fit_knn$knn_index), c(45L, 18L))
  expect_equal(dim(fit_knn$knn_distance), c(45L, 18L))

  out_dense <- diagnostic_tests(fit_dense, engine = "cpp")
  out_knn <- diagnostic_tests(fit_knn, engine = "cpp")

  expect_equal(out_knn$global_tests$F, out_dense$global_tests$F, tolerance = 1e-8)
  expect_equal(out_knn$global_tests$p_value, out_dense$global_tests$p_value, tolerance = 1e-8)
  expect_equal(out_knn$coefficient_tests$F3, out_dense$coefficient_tests$F3, tolerance = 1e-8)
  expect_equal(out_knn$coefficient_tests$p_value, out_dense$coefficient_tests$p_value, tolerance = 1e-8)
  expect_equal(out_knn$diagnostics$L_delta2, out_dense$diagnostics$L_delta2, tolerance = 1e-8)
  expect_equal(out_knn$diagnostics$gamma1, out_dense$diagnostics$gamma1, tolerance = 1e-8)
  expect_equal(out_knn$diagnostics$gamma2, out_dense$diagnostics$gamma2, tolerance = 1e-8)
})

test_that("diagnostic_tests kNN C++ engine reuses stored neighborhoods", {
  set.seed(20260625)
  H <- simulate_H2(45, r_max = 2.5)
  x1 <- rnorm(45)
  x2 <- rnorm(45)
  y <- 1 + 2 * x1 - 0.5 * x2 + rnorm(45, sd = 0.25)
  dat <- data.frame(y = y, x1 = x1, x2 = x2)

  fit_knn <- lhwr(
    y ~ x1 + x2,
    data = dat,
    Z = H$Z,
    bw = 18,
    adaptive = TRUE,
    kernel = "bisquare",
    diagnostics = "aicc",
    engine = "knn"
  )
  fit_recomputed <- fit_knn
  fit_recomputed$knn_index <- NULL
  fit_recomputed$knn_distance <- NULL

  out_cached <- diagnostic_tests(fit_knn, engine = "cpp")
  out_recomputed <- diagnostic_tests(fit_recomputed, engine = "cpp")

  expect_equal(out_cached$global_tests$F, out_recomputed$global_tests$F, tolerance = 1e-8)
  expect_equal(out_cached$global_tests$p_value, out_recomputed$global_tests$p_value, tolerance = 1e-8)
  expect_equal(out_cached$coefficient_tests$F3, out_recomputed$coefficient_tests$F3, tolerance = 1e-8)
  expect_equal(out_cached$coefficient_tests$p_value, out_recomputed$coefficient_tests$p_value, tolerance = 1e-8)
})

test_that("diagnostic tests use full trace moments and correct F tails", {
  set.seed(20260808)
  n <- 32
  H <- simulate_H2(n, r_max = 2.5)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1 + 1.5 * x1 - 0.75 * x2 + rnorm(n, sd = 0.3)
  dat <- data.frame(y = y, x1 = x1, x2 = x2)

  fit <- lhwr(
    y ~ x1 + x2,
    data = dat,
    Z = H$Z,
    bw = 15,
    adaptive = TRUE,
    kernel = "bisquare",
    diagnostics = "aicc",
    engine = "dense"
  )
  ols <- fit$lm
  X <- stats::model.matrix(ols)
  D <- hyperbolic_distance_matrix(H$Z, check = FALSE)
  bw <- fit$bandwidth
  if (length(bw) == 1L) bw <- rep(bw, n)

  S <- matrix(0, n, n)
  B_ops <- lapply(seq_len(ncol(X)), function(k) matrix(0, n, n))
  for (i in seq_len(n)) {
    w <- as.numeric(hgwr_kernel(D[, i], bw[i], kernel = "bisquare"))
    Xw <- X * w
    C <- solve(crossprod(X, Xw), t(Xw))
    S[i, ] <- as.numeric(X[i, ] %*% C)
    for (k in seq_len(ncol(X))) B_ops[[k]][i, ] <- C[k, ]
  }

  residual_operator <- diag(n) - S
  Q <- crossprod(residual_operator)
  expected_delta1 <- sum(diag(Q))
  expected_delta2 <- sum(Q^2)
  M <- diag(n) - matrix(1 / n, n, n)
  expected_gamma1 <- expected_gamma2 <- numeric(ncol(X))
  for (k in seq_len(ncol(X))) {
    Ck <- crossprod(B_ops[[k]], M %*% B_ops[[k]]) / n
    expected_gamma1[k] <- sum(diag(Ck))
    expected_gamma2[k] <- sum(Ck^2)
  }

  raw <- hgwr_f_tests_R_components(
    object = fit,
    X = X,
    RSS_full = RSS.lhwr(fit),
    RSS_ols = sum(stats::residuals(ols)^2),
    DF_ols = stats::df.residual(ols),
    trS = fit$LHWR.diagnostic$tr.S,
    trStS = fit$LHWR.diagnostic$tr.StS
  )

  expect_equal(raw$L_delta1, expected_delta1, tolerance = 1e-10)
  expect_equal(raw$L_delta2, expected_delta2, tolerance = 1e-10)
  expect_gt(raw$L_delta2, sum(diag(Q)^2))
  expect_equal(raw$gamma1, expected_gamma1, tolerance = 1e-10)
  expect_equal(raw$gamma2, expected_gamma2, tolerance = 1e-10)
  expect_equal(raw$F1_p, stats::pf(raw$F1, raw$F1_df1, raw$F1_df2, lower.tail = TRUE))
  expect_equal(raw$F2_p, stats::pf(raw$F2, raw$F2_df1, raw$F2_df2, lower.tail = FALSE))
  expect_equal(
    raw$F3_p,
    stats::pf(raw$F3, raw$F3_df1, raw$F3_df2, lower.tail = FALSE)
  )
  expect_identical(raw$moment_method, "full_trace")
})

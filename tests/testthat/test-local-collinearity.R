test_that("local_collinearity returns expected diagnostic tables", {
  set.seed(20260627)
  H <- simulate_H2(45, r_max = 2.5)
  x1 <- rnorm(45)
  x2 <- 0.6 * x1 + rnorm(45, sd = 0.5)
  y <- 1 + 2 * x1 - 0.5 * x2 + rnorm(45, sd = 0.25)
  dat <- data.frame(y = y, x1 = x1, x2 = x2)

  fit <- lhwr(
    y ~ x1 + x2,
    data = dat,
    Z = H$Z,
    bw = 18,
    adaptive = TRUE,
    kernel = "bisquare",
    diagnostics = "aicc"
  )
  out <- local_collinearity(fit)

  expect_s3_class(out, "local_collinearity_lhwr")
  expect_true(is.data.frame(out$condition_number))
  expect_true(is.data.frame(out$local_vif))
  expect_true(is.data.frame(out$local_correlation_summary))
  expect_equal(nrow(out$condition_number), 45)
  expect_equal(nrow(out$local_vif), 90)
})

test_that("local condition numbers and VIF flags are well formed", {
  set.seed(20260628)
  H <- simulate_H2(40, r_max = 2.5)
  x1 <- rnorm(40)
  x2 <- rnorm(40)
  y <- 1 + 2 * x1 - 0.5 * x2 + rnorm(40, sd = 0.25)
  dat <- data.frame(y = y, x1 = x1, x2 = x2)

  fit <- lhwr(y ~ x1 + x2, data = dat, Z = H$Z, bw = 16, diagnostics = "aicc")
  out <- local_collinearity(fit)

  expect_true(all(out$condition_number$local_CN >= 0 | is.infinite(out$condition_number$local_CN)))
  expect_type(out$condition_number$flag_CN, "logical")
  finite_vif <- out$local_vif$local_VIF[is.finite(out$local_vif$local_VIF)]
  expect_true(all(finite_vif >= 1 - 1e-10))
  expect_type(out$local_vif$flag_VIF, "logical")
})

test_that("local_collinearity print output includes expected headers", {
  set.seed(20260629)
  H <- simulate_H2(35, r_max = 2)
  x1 <- rnorm(35)
  x2 <- rnorm(35)
  y <- 1 + 2 * x1 - 0.5 * x2 + rnorm(35, sd = 0.25)
  dat <- data.frame(y = y, x1 = x1, x2 = x2)

  fit <- lhwr(y ~ x1 + x2, data = dat, Z = H$Z, bw = 15, diagnostics = "aicc")
  out <- local_collinearity(fit)
  printed <- capture.output(print(out))

  expect_true(any(grepl("LHWR Local Multicollinearity Diagnostics", printed, fixed = TRUE)))
  expect_true(any(grepl("Neighbor count", printed, fixed = TRUE)))
  expect_true(any(grepl("Median radius", printed, fixed = TRUE)))
  expect_true(any(grepl("Condition number summary", printed, fixed = TRUE)))
  expect_true(any(grepl("Local VIF summary", printed, fixed = TRUE)))
  expect_true(any(grepl("# Locations VIF", printed, fixed = TRUE)))
  expect_true(any(grepl("Local weighted correlation summary", printed, fixed = TRUE)))
})

test_that("local_collinearity kNN cached neighborhoods match recomputed weights", {
  set.seed(20260630)
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

  cached <- local_collinearity(fit_knn)
  recomputed <- local_collinearity(fit_recomputed)

  expect_equal(cached$condition_number$local_CN, recomputed$condition_number$local_CN, tolerance = 1e-10)
  expect_equal(cached$local_vif$local_VIF, recomputed$local_vif$local_VIF, tolerance = 1e-10)
  expect_equal(
    cached$local_correlation_summary$max_abs_local_correlation,
    recomputed$local_correlation_summary$max_abs_local_correlation,
    tolerance = 1e-10
  )
})

test_that("local_collinearity C++ engine matches R reference engine", {
  set.seed(20260702)
  H <- simulate_H2(45, r_max = 2.5)
  x1 <- rnorm(45)
  x2 <- 0.5 * x1 + rnorm(45, sd = 0.6)
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
  out_cpp <- local_collinearity(fit_knn, engine = "cpp")
  out_R <- local_collinearity(fit_knn, engine = "R")

  expect_equal(out_cpp$condition_number$local_CN, out_R$condition_number$local_CN, tolerance = 1e-10)
  expect_equal(out_cpp$local_vif$local_VIF, out_R$local_vif$local_VIF, tolerance = 1e-10)
  expect_equal(
    out_cpp$local_correlation_summary$max_abs_local_correlation,
    out_R$local_correlation_summary$max_abs_local_correlation,
    tolerance = 1e-10
  )
  expect_equal(out_cpp$diagnostics$collinearity_engine, "cpp")
  expect_equal(out_R$diagnostics$collinearity_engine, "R")
})

test_that("local_collinearity can return pairwise local correlations", {
  set.seed(20260701)
  H <- simulate_H2(30, r_max = 2)
  x1 <- rnorm(30)
  x2 <- rnorm(30)
  y <- 1 + x1 + x2 + rnorm(30, sd = 0.2)
  dat <- data.frame(y = y, x1 = x1, x2 = x2)

  fit <- lhwr(y ~ x1 + x2, data = dat, Z = H$Z, bw = 12, diagnostics = "aicc")
  out <- local_collinearity(fit, return_pairwise = TRUE)

  expect_true(is.data.frame(out$local_correlation_pairs))
  expect_equal(nrow(out$local_correlation_pairs), 30)
  expect_true(all(c("variable_1", "variable_2", "local_correlation") %in% names(out$local_correlation_pairs)))
})

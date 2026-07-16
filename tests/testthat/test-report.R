test_that("lhwr_report prints a structured report", {
  set.seed(20260721)
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
    diagnostics = "full",
    f_tests = TRUE
  )
  txt <- capture.output(out <- lhwr_report(fit, include_moran = FALSE, include_collinearity = FALSE))

  expect_s3_class(out, "lhwr_report")
  expect_true(any(grepl("LHWR Model Report", txt)))
  expect_true(any(grepl("Global OLS baseline", txt)))
  expect_true(any(grepl("LHWR diagnostics", txt)))
  expect_true(any(grepl("Local coefficient inference", txt)))
  expect_true(any(grepl("Diagnostic F-tests", txt)))
  expect_equal(out$call, fit$this.call)
})

test_that("lhwr_report can suppress optional sections", {
  set.seed(20260722)
  H <- simulate_H2(35, r_max = 2)
  x1 <- rnorm(35)
  y <- 1 + 2 * x1 + rnorm(35, sd = 0.25)
  dat <- data.frame(y = y, x1 = x1)

  fit <- lhwr(y ~ x1, data = dat, Z = H$Z, bw = 16, adaptive = TRUE, diagnostics = "aicc")
  txt <- capture.output(out <- lhwr_report(
    fit,
    include_global = FALSE,
    include_inference = FALSE,
    include_f_tests = FALSE,
    include_moran = FALSE,
    include_collinearity = FALSE
  ))

  expect_s3_class(out, "lhwr_report")
  expect_false(any(grepl("Global OLS baseline", txt)))
  expect_false(any(grepl("Local coefficient inference", txt)))
  expect_false(any(grepl("Diagnostic F-tests", txt)))
})

test_that("lhwr_report accepts precomputed optional diagnostics", {
  set.seed(20260723)
  H <- simulate_H2(40, r_max = 2)
  x1 <- rnorm(40)
  x2 <- rnorm(40)
  y <- 1 + 2 * x1 - x2 + rnorm(40, sd = 0.25)
  dat <- data.frame(y = y, x1 = x1, x2 = x2)

  fit <- lhwr(y ~ x1 + x2, data = dat, Z = H$Z, bw = 18, adaptive = TRUE, diagnostics = "full")
  inf <- local_inference(fit)
  mi <- moran.lhwr(fit, k = 18, permutations = 0)
  txt <- capture.output(out <- lhwr_report(
    fit,
    local_inference = inf,
    moran = mi,
    include_collinearity = FALSE,
    include_f_tests = FALSE
  ))

  expect_s3_class(out, "lhwr_report")
  expect_true(any(grepl("Residual hyperbolic autocorrelation", txt)))
  expect_s3_class(out$local_inference, "lhwr_local_inference")
})

test_that("local_inference returns finite local p-values", {
  set.seed(20260713)
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
    diagnostics = "full"
  )
  out <- local_inference(fit, adjust = "BH")

  expect_s3_class(out, "lhwr_local_inference")
  expect_equal(nrow(out$local), length(fit$coef_names) * nrow(fit$SDF))
  expect_true(all(is.na(out$local$p_value) | (out$local$p_value >= 0 & out$local$p_value <= 1)))
  expect_true(all(is.na(out$local$p_adjusted) | (out$local$p_adjusted >= 0 & out$local$p_adjusted <= 1)))
  expect_true(all(out$local$coefficient %in% fit$coef_names))
})

test_that("local_inference supports degrees-of-freedom choices", {
  set.seed(20260714)
  H <- simulate_H2(40, r_max = 2)
  x1 <- rnorm(40)
  y <- 1 + 2 * x1 + rnorm(40, sd = 0.25)
  dat <- data.frame(y = y, x1 = x1)

  fit <- lhwr(y ~ x1, data = dat, Z = H$Z, bw = 18, adaptive = TRUE, diagnostics = "full")
  edf <- local_inference(fit, df = "edf")
  residual <- local_inference(fit, df = "residual")
  normal <- local_inference(fit, df = "normal")

  expect_true(is.finite(edf$settings$reference_df))
  expect_true(is.finite(residual$settings$reference_df))
  expect_true(is.infinite(normal$settings$reference_df))
})

test_that("local_inference adjustment scope is valid", {
  set.seed(20260715)
  H <- simulate_H2(42, r_max = 2)
  x1 <- rnorm(42)
  x2 <- rnorm(42)
  y <- 1 + 2 * x1 - x2 + rnorm(42, sd = 0.25)
  dat <- data.frame(y = y, x1 = x1, x2 = x2)

  fit <- lhwr(y ~ x1 + x2, data = dat, Z = H$Z, bw = 18, adaptive = TRUE, diagnostics = "full")
  by_coef <- local_inference(fit, adjust_scope = "coefficient")
  global <- local_inference(fit, adjust_scope = "global")

  expect_equal(nrow(by_coef$local), nrow(global$local))
  expect_true(all(is.na(global$local$p_adjusted) | global$local$p_adjusted >= global$local$p_value - 1e-12))
})

test_that("plot.lhwr supports p-value and significance maps", {
  set.seed(20260716)
  H <- simulate_H2(35, r_max = 2)
  x1 <- rnorm(35)
  y <- 1 + 2 * x1 + rnorm(35, sd = 0.25)
  dat <- data.frame(y = y, x1 = x1)

  fit <- lhwr(y ~ x1, data = dat, Z = H$Z, bw = 16, adaptive = TRUE, diagnostics = "full")
  f1 <- tempfile(fileext = ".png")
  f2 <- tempfile(fileext = ".png")
  grDevices::png(f1)
  expect_silent(plot(fit, type = "p_value", coefficient = "x1"))
  grDevices::dev.off()
  grDevices::png(f2)
  expect_silent(plot(fit, type = "significance", coefficient = "x1"))
  grDevices::dev.off()
  expect_true(all(file.exists(c(f1, f2))))
})

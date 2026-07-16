test_that("plotting grammar wrappers run", {
  set.seed(20260717)
  H <- simulate_H2(40, r_max = 2.5)
  Z <- H$Z
  x <- rnorm(40)
  y <- 1 + x + rnorm(40, sd = 0.2)
  dat <- data.frame(y = y, x = x)
  fit <- lhwr(y ~ x, dat, Z = Z, bw = 16, adaptive = TRUE, diagnostics = "full")
  edges <- data.frame(from = 1:5, to = 6:10)

  f1 <- tempfile(fileext = ".png")
  f2 <- tempfile(fileext = ".png")
  f3 <- tempfile(fileext = ".png")
  f4 <- tempfile(fileext = ".png")
  f5 <- tempfile(fileext = ".png")
  f6 <- tempfile(fileext = ".png")

  expect_silent(plot_lhwr_geometry(Z, type = "poincare", output_file = f1))
  expect_silent(plot_lhwr_embedding(Z, values = H$r, output_file = f2))
  expect_silent(plot_lhwr_network_map(Z, edges = edges, output_file = f3))
  expect_silent(plot_lhwr_inference(fit, coefficient = "x", type = "p_value", output_file = f4))
  expect_silent(plot_lhwr_inference(fit, coefficient = "x", type = "significance", output_file = f5))

  grDevices::png(f6)
  expect_silent(plot_lhwr_model(fit, type = "coefficient", coefficient = "x"))
  grDevices::dev.off()

  expect_true(all(file.exists(c(f1, f2, f3, f4, f5, f6))))
})

test_that("geometry radial and diagnostics plots run", {
  set.seed(20260718)
  H <- simulate_H2(35, r_max = 2)
  Z <- H$Z
  x <- rnorm(35)
  y <- 1 + x + rnorm(35, sd = 0.25)
  dat <- data.frame(y = y, x = x)
  fit <- lhwr(y ~ x, dat, Z = Z, bw = 15, adaptive = TRUE, diagnostics = "full")

  f1 <- tempfile(fileext = ".png")
  f2 <- tempfile(fileext = ".png")
  f3 <- tempfile(fileext = ".png")

  grDevices::png(f1)
  expect_silent(plot_lhwr_geometry(Z, type = "radial"))
  grDevices::dev.off()

  grDevices::png(f2)
  expect_silent(plot_lhwr_diagnostics(fit, type = "residual"))
  grDevices::dev.off()

  grDevices::png(f3)
  expect_silent(plot_lhwr_diagnostics(fit, type = "condition_number"))
  grDevices::dev.off()

  expect_true(all(file.exists(c(f1, f2, f3))))
})

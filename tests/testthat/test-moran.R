test_that("global Moran's I is finite for HGWR residuals", {
  set.seed(101)
  H <- simulate_H2(45, r_max = 3)
  x <- rnorm(45)
  y <- 1 + 2 * x + rnorm(45, sd = 0.2)
  fit <- lhwr(y ~ x, data.frame(y, x), Z = H$Z, bw = 20)

  out <- moran.lhwr(fit, permutations = 19, seed = 1)

  expect_s3_class(out, "hgwr_moran")
  expect_true(is.finite(out$I))
})

test_that("global Moran's I expected value is -1/(n - 1)", {
  set.seed(102)
  H <- simulate_H2(35, r_max = 3)
  x <- rnorm(35)
  y <- 1 + x + rnorm(35, sd = 0.2)
  fit <- lhwr(y ~ x, data.frame(y, x), Z = H$Z, bw = 18)

  out <- moran.lhwr(fit, permutations = 0)

  expect_equal(out$expected_I, -1 / (out$n - 1))
})

test_that("global Moran permutation p-value is in [0, 1]", {
  set.seed(103)
  H <- simulate_H2(40, r_max = 3)
  x <- rnorm(40)
  y <- 1 + x + rnorm(40, sd = 0.2)
  fit <- lhwr(y ~ x, data.frame(y, x), Z = H$Z, bw = 20)

  out <- moran.lhwr(fit, permutations = 19, seed = 2)

  expect_true(is.na(out$p_value) || (out$p_value >= 0 && out$p_value <= 1))
})

test_that("supplied W can reproduce internally constructed weights", {
  set.seed(104)
  H <- simulate_H2(36, r_max = 3)
  D <- hyperbolic_distance_matrix(H$Z, check = FALSE)
  k <- 14
  bw <- apply(D, 1L, function(d) sort(d[d > 0], partial = k)[k])
  W <- matrix(0, nrow(D), ncol(D))
  for (i in seq_len(nrow(D))) W[i, ] <- hgwr_kernel(D[i, ], bw[i], kernel = "bisquare")
  diag(W) <- 0
  W <- W / rowSums(W)

  x <- rnorm(36)
  y <- 1 + x + rnorm(36, sd = 0.2)
  fit <- lhwr(y ~ x, data.frame(y, x), Z = H$Z, bw = k)

  internal <- moran.lhwr(fit, k = k, permutations = 0)
  supplied <- moran.lhwr(fit, W = W, permutations = 0)

  expect_equal(supplied$I, internal$I, tolerance = 1e-10)
})

test_that("zero residual variance does not crash Moran diagnostics", {
  set.seed(105)
  H <- simulate_H2(25, r_max = 3)
  x <- rnorm(25)
  y <- 1 + 2 * x
  fit <- lhwr(y ~ x, data.frame(y, x), Z = H$Z, bw = 15)

  out <- moran.lhwr(fit, permutations = 0)

  expect_true(is.na(out$I) || is.finite(out$I))
})

test_that("missing residuals produce a clear error", {
  set.seed(106)
  H <- simulate_H2(30, r_max = 3)
  x <- rnorm(30)
  y <- 1 + x + rnorm(30, sd = 0.2)
  fit <- lhwr(y ~ x, data.frame(y, x), Z = H$Z, bw = 15)
  fit$residuals <- NULL
  fit$SDF$residual <- NULL

  expect_error(moran.lhwr(fit), "No residuals")
})

test_that("local Moran output has one row per observation", {
  set.seed(107)
  H <- simulate_H2(38, r_max = 3)
  x <- rnorm(38)
  y <- 1 + x + rnorm(38, sd = 0.2)
  fit <- lhwr(y ~ x, data.frame(y, x), Z = H$Z, bw = 18)

  out <- local_moran.lhwr(fit, permutations = 19, seed = 3)

  expect_s3_class(out, "hgwr_local_moran")
  expect_equal(nrow(out), length(fit$residuals))
  expect_true(all(c("observation_id", "local_I", "z_score", "p_value", "cluster_type") %in% names(out)))
})

test_that("local Moran cluster labels are valid", {
  set.seed(108)
  H <- simulate_H2(38, r_max = 3)
  x <- rnorm(38)
  y <- 1 + x + rnorm(38, sd = 0.2)
  fit <- lhwr(y ~ x, data.frame(y, x), Z = H$Z, bw = 18)

  out <- local_moran.lhwr(fit, permutations = 19, seed = 4)
  valid <- c("High-High", "Low-Low", "High-Low", "Low-High", "Not Significant")

  expect_true(all(as.character(out$cluster_type) %in% valid))
})

test_that("Moran diagnostics reuse stored kNN neighborhoods without changing results", {
  set.seed(20260626)
  H <- simulate_H2(50, r_max = 2.5)
  x1 <- rnorm(50)
  x2 <- rnorm(50)
  y <- 1 + 2 * x1 - 0.5 * x2 + rnorm(50, sd = 0.25)
  dat <- data.frame(y = y, x1 = x1, x2 = x2)

  fit_knn <- lhwr(
    y ~ x1 + x2,
    data = dat,
    Z = H$Z,
    bw = 20,
    adaptive = TRUE,
    kernel = "bisquare",
    diagnostics = "aicc",
    engine = "knn"
  )
  fit_recomputed <- fit_knn
  fit_recomputed$knn_index <- NULL
  fit_recomputed$knn_distance <- NULL

  cached_global <- moran.lhwr(fit_knn, k = 20, permutations = 0)
  recomputed_global <- moran.lhwr(fit_recomputed, k = 20, permutations = 0)
  expect_equal(cached_global$I, recomputed_global$I, tolerance = 1e-12)

  cached_local <- local_moran.lhwr(fit_knn, k = 20, permutations = 0)
  recomputed_local <- local_moran.lhwr(fit_recomputed, k = 20, permutations = 0)
  expect_equal(cached_local$local_I, recomputed_local$local_I, tolerance = 1e-12)
  expect_equal(cached_local$spatial_lag, recomputed_local$spatial_lag, tolerance = 1e-12)
})

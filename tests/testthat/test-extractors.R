test_that("standard model extractors work for LHWR objects", {
  set.seed(17)
  H <- simulate_H2(60, r_max = 3, uniform_area = TRUE)
  x1 <- rnorm(60)
  x2 <- rnorm(60)
  y <- 1 + 0.7 * x1 - 0.2 * x2 + rnorm(60, sd = 0.3)
  dat <- data.frame(y = y, x1 = x1, x2 = x2)
  fit <- lhwr(y ~ x1 + x2, data = dat, Z = H$Z, bw = 25, adaptive = TRUE)

  expect_equal(fitted(fit), fit$y_hat)
  expect_equal(predict(fit), fit$y_hat)
  expect_equal(residuals(fit), fit$residuals)
  expect_equal(coef(fit), fit$beta_hat)
  expect_equal(dim(coef(fit)), c(60L, 3L))

  Hnew <- simulate_H2(4, r_max = 3, uniform_area = TRUE)
  newdata <- data.frame(x1 = rnorm(4), x2 = rnorm(4))
  pred <- predict(fit, newdata = newdata, Znew = Hnew$Z)
  pred_coef <- predict(fit, newdata = newdata, Znew = Hnew$Z, type = "coefficients")
  expect_length(pred, 4L)
  expect_equal(dim(pred_coef), c(4L, 3L))
  expect_true(all(is.finite(pred)))
  expect_true(all(is.finite(pred_coef)))
  expect_error(predict(fit, newdata = newdata), "Znew")
})

test_that("predict.lhwr supports tangent geometry for new observations", {
  set.seed(18)
  H <- simulate_H2(50, r_max = 3, uniform_area = TRUE)
  x <- rnorm(50)
  y <- 1 + 1.5 * x + rnorm(50, sd = 0.3)
  dat <- data.frame(y = y, x = x)
  fit <- lhwr(y ~ x, data = dat, Z = H$Z, bw = 20, adaptive = TRUE, geometry = "tangent")

  Hnew <- simulate_H2(3, r_max = 3, uniform_area = TRUE)
  pred <- predict(fit, newdata = data.frame(x = rnorm(3)), Znew = Hnew$Z)
  expect_length(pred, 3L)
  expect_true(all(is.finite(pred)))
})

test_that("newdata prediction reproduces adaptive calibration fits", {
  set.seed(181)
  H <- simulate_H2(70, r_max = 4, uniform_area = TRUE)
  x1 <- rnorm(70)
  x2 <- rnorm(70)
  y <- 0.5 + 1.2 * x1 - 0.8 * x2 + rnorm(70, sd = 0.2)
  dat <- data.frame(y, x1, x2)
  fit <- lhwr(
    y ~ x1 + x2,
    data = dat,
    Z = H$Z,
    bw = 24,
    adaptive = TRUE,
    kernel = "bisquare",
    engine = "dense"
  )

  pred <- predict(fit, newdata = dat, Znew = H$Z)
  pred_coef <- predict(fit, newdata = dat, Znew = H$Z, type = "coefficients")

  expect_equal(pred, fitted(fit), tolerance = 1e-8)
  expect_equal(unname(pred_coef), unname(coef(fit)), tolerance = 1e-8)
})

test_that("Tangent-E prediction preserves zero distances and calibration fits", {
  set.seed(20260909)
  Z <- simulate_H2(30, r_max = 3)$Z
  dat <- data.frame(x = rnorm(30), y = rnorm(30))
  for (adaptive in c(TRUE, FALSE)) {
    for (kernel in c("bisquare", "gaussian", "exponential", "tricube")) {
      fit <- lhwr(
        y ~ x, data = dat, Z = Z, bw = if (adaptive) 15 else 6,
        adaptive = adaptive, kernel = kernel, geometry = "tangent",
        diagnostics = "none", engine = "dense"
      )
      D <- LHWL:::lhwr_predict_distance(fit, Z, "tangent")
      expect_identical(diag(D), rep(0, nrow(Z)))
      expect_equal(D, tangent_distance_matrix(Z, fit$tangent_ref),
                   tolerance = 1e-12)
      expect_equal(predict(fit, newdata = dat, Znew = Z), fitted(fit),
                   tolerance = 1e-8)
      expect_equal(unname(predict(fit, newdata = dat, Znew = Z,
                                 type = "coefficients")),
                   unname(coef(fit)), tolerance = 1e-8)
      rows <- c(7L, 3L, 7L, 1L)
      expect_equal(predict(fit, newdata = dat[rows, , drop = FALSE],
                           Znew = Z[rows, , drop = FALSE]),
                   unname(fitted(fit)[rows]), tolerance = 1e-8)
    }
  }
})

test_that("Tangent-E calibration prediction supports duplicated locations", {
  set.seed(20260911)
  base_Z <- simulate_H2(30, r_max = 3)$Z
  Z <- rbind(base_Z, base_Z[c(1, 7), , drop = FALSE])
  dat <- data.frame(x = rnorm(nrow(Z)), y = rnorm(nrow(Z)))
  fit <- lhwr(y ~ x, data = dat, Z = Z, bw = 15, geometry = "tangent",
              diagnostics = "none", engine = "dense")
  D <- LHWL:::lhwr_predict_distance(fit, Z, "tangent")
  expect_identical(D[1, 31], 0)
  expect_identical(D[7, 32], 0)
  expect_equal(predict(fit, newdata = dat, Znew = Z), fitted(fit),
               tolerance = 1e-8)
})

test_that("Tangent-E prediction retains close distinct distances", {
  r <- 3
  eps <- 1e-9
  Z <- rbind(c(cosh(r), sinh(r), 0), c(cosh(r + eps), sinh(r + eps), 0))
  object <- list(Z = Z, tangent_ref = c(1, 0, 0))
  D <- LHWL:::lhwr_predict_distance(object, Z, "tangent")
  expect_identical(diag(D), c(0, 0))
  expect_gt(D[1, 2], 0)
  expect_equal(D, tangent_distance_matrix(Z, object$tangent_ref),
               tolerance = 1e-14)
  # cosh/sinh and the log map round coordinates at r = 3 on an absolute
  # floating-point scale, not relative to this 1e-9 separation. A relative
  # tolerance of 1e-14 would require an unattainable error around 1e-23.
  expect_lt(abs(D[1, 2] - eps), 1e-14)
})

test_that("Tangent-E prediction reuses its training reference in every batch", {
  set.seed(20260910)
  Z <- simulate_H2(35, r_max = 3)$Z
  dat <- data.frame(x = rnorm(35), y = rnorm(35))
  fit <- lhwr(y ~ x, data = dat, Z = Z, bw = 18, geometry = "tangent",
              diagnostics = "none", engine = "dense")
  Znew <- simulate_H2(3, r_max = 3)$Z
  newdata <- data.frame(x = rnorm(3))
  Q <- tangent_coordinates(Z, fit$tangent_ref)
  Qnew <- tangent_coordinates(Znew, fit$tangent_ref)
  Dref <- t(vapply(seq_len(nrow(Qnew)), function(i) {
    sqrt(rowSums(sweep(Q, 2L, Qnew[i, ], "-")^2))
  }, numeric(nrow(Q))))
  expect_equal(LHWL:::lhwr_predict_distance(fit, Znew, "tangent"), Dref,
               tolerance = 1e-12)
  batch <- predict(fit, newdata = newdata, Znew = Znew)
  single <- vapply(seq_len(nrow(Znew)), function(i) {
    predict(fit, newdata = newdata[i, , drop = FALSE],
            Znew = Znew[i, , drop = FALSE])
  }, numeric(1))
  expect_equal(batch, single, tolerance = 1e-12)
  expect_equal(batch, predict(fit, newdata = newdata, dMat = Dref),
               tolerance = 1e-12)
})

test_that("prediction preserves the selected kernel under rank deficiency", {
  X <- cbind(Intercept = 1, x = c(-1, 0, 1))
  y <- c(1, 2, 4)
  D_new <- matrix(c(0.1, 1, 2), nrow = 1)
  bw <- 0.2
  w <- hgwr_kernel(D_new[1, ], bw, kernel = "bisquare")
  XtWX <- crossprod(X, X * w)
  XtWy <- crossprod(X, w * y)
  expected <- MASS::ginv(XtWX) %*% XtWy

  observed <- LHWL:::lhwr_predict_local_coefficients(
    y, X, D_new, bw_new = bw, kernel = "bisquare"
  )
  expect_equal(as.numeric(observed), as.numeric(expected), tolerance = 1e-10)
})

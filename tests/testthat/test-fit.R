test_that("bandwidth selection and HGWR fitting run", {
  set.seed(1)
  H <- simulate_H2(35, r_max = 2)
  x <- rnorm(35)
  y <- 1 + 2 * x + rnorm(35, sd = 0.05)
  dat <- data.frame(y = y, x = x)
  
  bw <- bw.lhwr(
    y ~ x,
    data = dat,
    Z = H$Z,
    approach = "CV",
    bw_grid = c(8, 12),
    adaptive = TRUE
  )
  expect_s3_class(bw, "bw.lhwr")
  expect_true(bw$bw %in% c(8, 12))
  
  fit <- lhwr(y ~ x, dat, H$Z, bw = bw$bw)
  expect_s3_class(fit, "lhwr")
  expect_equal(nrow(fit$SDF), 35)
  expect_true(is.finite(AICc.lhwr(fit)))
  expect_true(RSS.lhwr(fit) >= 0)
})

test_that("tangent geometry bandwidth selection and fitting run", {
  set.seed(13)
  H <- simulate_H2(40, r_max = 3)
  x <- rnorm(40)
  y <- 1 + 2 * x + rnorm(40, sd = 0.1)
  dat <- data.frame(y = y, x = x)
  
  bw <- bw.lhwr(
    y ~ x,
    data = dat,
    Z = H$Z,
    approach = "AICc",
    bw_grid = c(10, 15),
    adaptive = TRUE,
    geometry = "tangent"
  )
  expect_s3_class(bw, "bw.lhwr")
  expect_equal(bw$geometry, "tangent")
  expect_match(bw$distance, "Tangent-space")
  
  fit <- lhwr(
    y ~ x,
    data = dat,
    Z = H$Z,
    bw = bw$bw,
    adaptive = TRUE,
    geometry = "tangent"
  )
  expect_s3_class(fit, "lhwr")
  expect_equal(fit$geometry, "tangent")
  expect_match(fit$distance_label, "Tangent-space")
  expect_equal(nrow(fit$SDF), 40)
  expect_true(is.finite(AICc.lhwr(fit)))
})

test_that("tangent geometry does not use hyperbolic volume bandwidths", {
  set.seed(14)
  H <- simulate_H2(25, r_max = 2)
  dat <- data.frame(y = rnorm(25), x = rnorm(25))
  
  expect_error(
    bw.lhwr(y ~ x, dat, H$Z, geometry = "tangent", bandwidth_type = "volume"),
    "only defined"
  )
  expect_error(
    lhwr(y ~ x, dat, H$Z, bw = 8, geometry = "tangent", bandwidth_type = "volume"),
    "only defined"
  )
})

test_that("volume adaptive bandwidths use local density", {
  set.seed(2)
  H <- simulate_H2(40, r_max = 3)
  D <- hyperbolic_distance_matrix(H$Z)
  
  bw_knn <- resolve_bandwidth(D, 10, adaptive = TRUE, bandwidth_type = "knn")
  bw_volume <- resolve_bandwidth(
    D,
    10,
    adaptive = TRUE,
    bandwidth_type = "volume",
    pilot_k = 5
  )
  
  expect_length(bw_volume, 40)
  expect_true(all(is.finite(bw_volume)))
  expect_false(isTRUE(all.equal(bw_knn, bw_volume)))
})

test_that("AICc diagnostics follow GWmodel-style formula", {
  set.seed(4)
  H <- simulate_H2(45, r_max = 3)
  x <- rnorm(45)
  y <- 1 + x + rnorm(45, sd = 0.1)
  X <- model.matrix(y ~ x)
  colnames(X) <- clean_coef_names(colnames(X))
  D <- hyperbolic_distance_matrix(H$Z)
  bw <- resolve_bandwidth(D, 16, adaptive = TRUE, bandwidth_type = "knn")
  fit <- hgwr.fit(y, X, D, bw)
  
  n <- length(y)
  rss <- sum((y - fit$y_hat)^2)
  trS <- sum(fit$hat_diag)
  dev <- n * log(rss / n) + n * log(2 * pi) + n
  aicc <- dev + 2 * n * (trS + 1) / (n - trS - 2)
  
  expect_equal(fit$RSS, rss)
  expect_equal(fit$trS, trS)
  expect_true(is.finite(fit$AICc))
  expect_equal(fit$AICc, aicc)
})

test_that("adaptive kNN bandwidth excludes self-distance", {
  set.seed(5)
  H <- simulate_H2(20, r_max = 3)
  D <- hyperbolic_distance_matrix(H$Z)
  bw <- resolve_bandwidth(D, 5, adaptive = TRUE, bandwidth_type = "knn")
  manual <- apply(D, 2, function(x) sort(x[x > 0])[5])
  
  expect_equal(bw, manual)
  expect_true(all(bw > 0))
})

test_that("hat matrix diagonal is the influence of y_i on fitted y_i", {
  set.seed(6)
  H <- simulate_H2(30, r_max = 3)
  x <- rnorm(30)
  y <- 1 + 2 * x + rnorm(30, sd = 0.1)
  X <- model.matrix(y ~ x)
  colnames(X) <- clean_coef_names(colnames(X))
  D <- hyperbolic_distance_matrix(H$Z)
  bw <- resolve_bandwidth(D, 12, adaptive = TRUE, bandwidth_type = "knn")
  i <- 7
  w <- hgwr_kernel(D[, i], bw[i], kernel = "bisquare")
  
  fit_i <- hgwr_wls_one(X, y, w, focus = i, return_hat = TRUE)
  Xw <- X * w
  C <- solve(crossprod(X, Xw)) %*% t(Xw)
  s_row <- as.numeric(X[i, , drop = FALSE] %*% C)
  
  expect_equal(fit_i$s_diag, s_row[i])
  expect_equal(fit_i$s2_diag, sum(s_row^2))
  expect_equal(sum(X[i, ] * fit_i$beta), sum(s_row * y))
})

test_that("LOOCV sets focal self-weight to zero", {
  set.seed(7)
  H <- simulate_H2(32, r_max = 3)
  x <- rnorm(32)
  y <- 1 + x + rnorm(32, sd = 0.1)
  X <- model.matrix(y ~ x)
  colnames(X) <- clean_coef_names(colnames(X))
  D <- hyperbolic_distance_matrix(H$Z)
  bw <- resolve_bandwidth(D, 14, adaptive = TRUE, bandwidth_type = "knn")
  
  cv_pkg <- hgwr_cv_score(y, X, D, bw, kernel = "bisquare")
  cv_manual <- 0
  for (i in seq_along(y)) {
    w <- hgwr_kernel(D[, i], bw[i], kernel = "bisquare")
    w[i] <- 0
    beta <- solve(crossprod(X, X * w), crossprod(X, w * y))
    pred <- sum(X[i, ] * as.numeric(beta))
    cv_manual <- cv_manual + (y[i] - pred)^2
  }
  
  expect_equal(cv_pkg, cv_manual)
})

test_that("C++ LOOCV score matches the R reference implementation", {
  ref_file <- test_path("../../../R_reference/bandwidth_R_reference.R")
  skip_if_not(file.exists(ref_file), "R reference files are local development artifacts")
  ref <- new.env(parent = asNamespace("LHWL"))
  sys.source(ref_file, envir = ref)
  
  set.seed(10)
  H <- simulate_H2(70, r_max = 4)
  x1 <- rnorm(70)
  x2 <- rnorm(70)
  y <- 1 + 2 * x1 - 0.5 * x2 + rnorm(70, sd = 0.2)
  X <- model.matrix(y ~ x1 + x2)
  colnames(X) <- clean_coef_names(colnames(X))
  D <- hyperbolic_distance_matrix(H$Z, check = FALSE)
  bw <- resolve_bandwidth(D, 24, adaptive = TRUE, bandwidth_type = "knn")
  
  expect_equal(
    hgwr_cv_score(y, X, D, bw, kernel = "bisquare"),
    ref$hgwr_cv_score(y, X, D, bw, kernel = "bisquare"),
    tolerance = 1e-10
  )
})

test_that("bandwidth selection table reports CV explicitly", {
  set.seed(8)
  H <- simulate_H2(35, r_max = 3)
  x <- rnorm(35)
  y <- 1 + 2 * x + rnorm(35, sd = 0.1)
  dat <- data.frame(y = y, x = x)
  
  bw <- bw.lhwr(
    y ~ x,
    data = dat,
    Z = H$Z,
    approach = "CV",
    bw_grid = c(10, 14),
    adaptive = TRUE
  )
  
  expect_true("CV" %in% names(bw$selection.table))
  expect_equal(bw$selection.table$CV, bw$selection.table$score)
  expect_equal(bw$bw, bw$selection.table$bw[which.min(bw$selection.table$score)])
})

test_that("golden bandwidth search runs without an explicit grid", {
  set.seed(81)
  H <- simulate_H2(45, r_max = 3)
  x <- rnorm(45)
  y <- 1 + 2 * x + rnorm(45, sd = 0.15)
  dat <- data.frame(y = y, x = x)
  
  bw <- bw.lhwr(
    y ~ x,
    data = dat,
    Z = H$Z,
    approach = "AICc",
    adaptive = TRUE,
    search = "golden"
  )
  
  expect_s3_class(bw, "bw.lhwr")
  expect_equal(bw$search, "golden")
  expect_true(bw$bw >= 2)
  expect_true(bw$bw < nrow(H$Z))
  expect_true(nrow(bw$selection.table) >= 2)
  expect_equal(bw$bw, bw$selection.table$bw[which.min(bw$selection.table$score)])
})

test_that("bandwidth grid keeps the explicit simulation search path", {
  set.seed(82)
  H <- simulate_H2(35, r_max = 2)
  dat <- data.frame(y = rnorm(35), x = rnorm(35))
  
  bw <- bw.lhwr(
    y ~ x,
    data = dat,
    Z = H$Z,
    approach = "CV",
    bw_grid = c(8, 12, 16),
    adaptive = TRUE
  )
  
  expect_equal(bw$search, "grid")
  expect_true(bw$bw %in% c(8, 12, 16))
  expect_equal(sort(bw$selection.table$bw), c(8, 12, 16))
})

test_that("C++ bandwidth grid selection matches the R reference implementation", {
  ref_file <- test_path("../../../R_reference/bandwidth_R_reference.R")
  skip_if_not(file.exists(ref_file), "R reference files are local development artifacts")
  ref <- new.env(parent = asNamespace("LHWL"))
  sys.source(ref_file, envir = ref)
  
  set.seed(11)
  H <- simulate_H2(60, r_max = 4)
  x1 <- rnorm(60)
  x2 <- rnorm(60)
  y <- 1 + 2 * x1 - 0.5 * x2 + rnorm(60, sd = 0.2)
  dat <- data.frame(y = y, x1 = x1, x2 = x2)
  D <- hyperbolic_distance_matrix(H$Z, check = FALSE)
  grid <- c(14, 22, 30)
  
  for (approach in c("AICc", "CV")) {
    for (bandwidth_type in c("knn", "volume")) {
      bw_r <- ref$bw.hgwr(
        y ~ x1 + x2,
        data = dat,
        Z = H$Z,
        approach = approach,
        bandwidth_type = bandwidth_type,
        bw_grid = grid,
        dMat = D,
        pilot_k = 8
      )
      bw_cpp <- bw.lhwr(
        y ~ x1 + x2,
        data = dat,
        Z = H$Z,
        approach = approach,
        bandwidth_type = bandwidth_type,
        bw_grid = grid,
        dMat = D,
        pilot_k = 8
      )
      
      expect_equal(bw_cpp$bw, bw_r$bw)
      expect_equal(bw_cpp$selection.table$score, bw_r$selection.table$score, tolerance = 1e-8)
      expect_equal(bw_cpp$selection.table$trS, bw_r$selection.table$trS, tolerance = 1e-8)
      expect_equal(bw_cpp$selection.table$RSS, bw_r$selection.table$RSS, tolerance = 1e-8)
    }
  }
})

test_that("C++ hgwr.fit matches the R reference implementation", {
  ref_file <- test_path("../../../R_reference/fit_R_reference.R")
  skip_if_not(file.exists(ref_file), "R reference files are local development artifacts")
  ref <- new.env(parent = asNamespace("LHWL"))
  sys.source(ref_file, envir = ref)
  
  set.seed(9)
  H <- simulate_H2(80, r_max = 4)
  x1 <- rnorm(80)
  x2 <- rnorm(80)
  y <- 1 + 2 * x1 - 0.5 * x2 + rnorm(80, sd = 0.2)
  X <- model.matrix(y ~ x1 + x2)
  colnames(X) <- clean_coef_names(colnames(X))
  D <- hyperbolic_distance_matrix(H$Z, check = FALSE)
  bw <- resolve_bandwidth(D, 24, adaptive = TRUE, bandwidth_type = "knn")
  
  fit_r <- ref$hgwr.fit(y, X, D, bw, kernel = "bisquare")
  fit_cpp <- hgwr.fit(y, X, D, bw, kernel = "bisquare")
  
  expect_equal(fit_cpp$beta_hat, fit_r$beta_hat, tolerance = 1e-10)
  expect_equal(fit_cpp$y_hat, fit_r$y_hat, tolerance = 1e-10)
  expect_equal(fit_cpp$hat_diag, fit_r$hat_diag, tolerance = 1e-10)
  expect_equal(fit_cpp$RSS, fit_r$RSS, tolerance = 1e-10)
  expect_equal(fit_cpp$AICc, fit_r$AICc, tolerance = 1e-8)
})

test_that("kNN HGWR fit matches dense adaptive fit for compact kernels", {
  set.seed(11)
  H <- simulate_H2(90, r_max = 4, uniform_area = TRUE)
  x1 <- rnorm(90)
  x2 <- rnorm(90)
  y <- 0.5 + 1.5 * x1 - 0.8 * x2 + rnorm(90, sd = 0.25)
  X <- model.matrix(y ~ x1 + x2)
  colnames(X) <- clean_coef_names(colnames(X))

  D <- hyperbolic_distance_matrix(H$Z, check = FALSE)
  k <- 28
  bw <- resolve_bandwidth(D, k, adaptive = TRUE, bandwidth_type = "knn")

  fit_dense <- hgwr.fit(y, X, D, bw, kernel = "bisquare")
  fit_knn <- hgwr.fit.knn(y, X, H$Z, k = k, kernel = "bisquare", check = FALSE)

  expect_equal(fit_knn$bandwidth, bw, tolerance = 1e-10)
  expect_equal(fit_knn$beta_hat, fit_dense$beta_hat, tolerance = 1e-10)
  expect_equal(fit_knn$y_hat, fit_dense$y_hat, tolerance = 1e-10)
  expect_equal(fit_knn$hat_diag, fit_dense$hat_diag, tolerance = 1e-10)
  expect_equal(fit_knn$RSS, fit_dense$RSS, tolerance = 1e-10)
  expect_equal(fit_knn$AICc, fit_dense$AICc, tolerance = 1e-8)
})

test_that("kNN engine matches dense engine for bandwidth selection and basic fit", {
  set.seed(12)
  H <- simulate_H2(85, r_max = 4, uniform_area = TRUE)
  x1 <- rnorm(85)
  x2 <- rnorm(85)
  y <- 0.2 + 1.2 * x1 - 0.4 * x2 + rnorm(85, sd = 0.3)
  dat <- data.frame(y = y, x1 = x1, x2 = x2)
  grid <- c(18, 25, 32)

  bw_dense <- bw.lhwr(
    y ~ x1 + x2,
    data = dat,
    Z = H$Z,
    approach = "AICc",
    kernel = "bisquare",
    adaptive = TRUE,
    bandwidth_type = "knn",
    search = "grid",
    bw_grid = grid,
    engine = "dense"
  )
  bw_knn <- bw.lhwr(
    y ~ x1 + x2,
    data = dat,
    Z = H$Z,
    approach = "AICc",
    kernel = "bisquare",
    adaptive = TRUE,
    bandwidth_type = "knn",
    search = "grid",
    bw_grid = grid,
    engine = "knn"
  )

  expect_equal(bw_knn$selection.table$score, bw_dense$selection.table$score, tolerance = 1e-8)
  expect_equal(bw_knn$selection.table$trS, bw_dense$selection.table$trS, tolerance = 1e-8)
  expect_equal(bw_knn$selection.table$RSS, bw_dense$selection.table$RSS, tolerance = 1e-10)
  expect_equal(bw_knn$bw, bw_dense$bw)

  bw_dense_cv <- bw.lhwr(
    y ~ x1 + x2,
    data = dat,
    Z = H$Z,
    approach = "CV",
    kernel = "bisquare",
    adaptive = TRUE,
    bandwidth_type = "knn",
    search = "grid",
    bw_grid = grid,
    engine = "dense"
  )
  bw_knn_cv <- bw.lhwr(
    y ~ x1 + x2,
    data = dat,
    Z = H$Z,
    approach = "CV",
    kernel = "bisquare",
    adaptive = TRUE,
    bandwidth_type = "knn",
    search = "grid",
    bw_grid = grid,
    engine = "knn"
  )
  expect_equal(bw_knn_cv$selection.table$score, bw_dense_cv$selection.table$score, tolerance = 1e-8)
  expect_equal(bw_knn_cv$selection.table$CV, bw_dense_cv$selection.table$CV, tolerance = 1e-8)
  expect_equal(bw_knn_cv$bw, bw_dense_cv$bw)

  fit_dense <- lhwr(
    y ~ x1 + x2,
    data = dat,
    Z = H$Z,
    bw = bw_dense$bw,
    kernel = "bisquare",
    adaptive = TRUE,
    bandwidth_type = "knn",
    engine = "dense"
  )
  fit_knn <- lhwr(
    y ~ x1 + x2,
    data = dat,
    Z = H$Z,
    bw = bw_knn$bw,
    kernel = "bisquare",
    adaptive = TRUE,
    bandwidth_type = "knn",
    engine = "knn"
  )

  expect_equal(fit_knn$beta_hat, fit_dense$beta_hat, tolerance = 1e-10)
  expect_equal(fit_knn$y_hat, fit_dense$y_hat, tolerance = 1e-10)
  expect_equal(fit_knn$LHWR.diagnostic$AICc, fit_dense$LHWR.diagnostic$AICc, tolerance = 1e-8)
})

test_that("kNN golden search matches dense golden search", {
  set.seed(13)
  H <- simulate_H2(90, r_max = 4, uniform_area = TRUE)
  x1 <- rnorm(90)
  x2 <- rnorm(90)
  y <- 0.8 + 1.1 * x1 - 0.7 * x2 + rnorm(90, sd = 0.25)
  dat <- data.frame(y = y, x1 = x1, x2 = x2)

  bw_dense <- bw.lhwr(
    y ~ x1 + x2,
    data = dat,
    Z = H$Z,
    approach = "AICc",
    kernel = "bisquare",
    adaptive = TRUE,
    bandwidth_type = "knn",
    search = "golden",
    lower = 20,
    upper = 60,
    engine = "dense"
  )
  bw_knn <- bw.lhwr(
    y ~ x1 + x2,
    data = dat,
    Z = H$Z,
    approach = "AICc",
    kernel = "bisquare",
    adaptive = TRUE,
    bandwidth_type = "knn",
    search = "golden",
    lower = 20,
    upper = 60,
    engine = "knn"
  )

  dense_tab <- bw_dense$selection.table[order(bw_dense$selection.table$bw), ]
  knn_tab <- bw_knn$selection.table[order(bw_knn$selection.table$bw), ]
  expect_equal(knn_tab$bw, dense_tab$bw)
  expect_equal(knn_tab$score, dense_tab$score, tolerance = 1e-8)
  expect_equal(knn_tab$RSS, dense_tab$RSS, tolerance = 1e-10)
  expect_equal(knn_tab$trS, dense_tab$trS, tolerance = 1e-8)
  expect_equal(bw_knn$bw, bw_dense$bw)
})

test_that("kNN thread argument preserves numerical results", {
  set.seed(14)
  H <- simulate_H2(95, r_max = 4, uniform_area = TRUE)
  x1 <- rnorm(95)
  x2 <- rnorm(95)
  y <- 0.5 + 1.3 * x1 - 0.2 * x2 + rnorm(95, sd = 0.3)
  X <- model.matrix(y ~ x1 + x2)
  colnames(X) <- clean_coef_names(colnames(X))

  fit_1 <- hgwr.fit.knn(y, X, H$Z, k = 30, kernel = "bisquare", check = FALSE, n_threads = 1)
  fit_2 <- hgwr.fit.knn(y, X, H$Z, k = 30, kernel = "bisquare", check = FALSE, n_threads = 2)

  expect_equal(fit_2$beta_hat, fit_1$beta_hat, tolerance = 1e-10)
  expect_equal(fit_2$y_hat, fit_1$y_hat, tolerance = 1e-10)
  expect_equal(fit_2$AICc, fit_1$AICc, tolerance = 1e-8)
})

test_that("lhwr diagnostic modes control output detail", {
  set.seed(15)
  H <- simulate_H2(70, r_max = 3, uniform_area = TRUE)
  x1 <- rnorm(70)
  x2 <- rnorm(70)
  y <- 0.4 + 1.2 * x1 - 0.6 * x2 + rnorm(70, sd = 0.3)
  dat <- data.frame(y = y, x1 = x1, x2 = x2)

  fit_full <- lhwr(
    y ~ x1 + x2,
    data = dat,
    Z = H$Z,
    bw = 30,
    kernel = "bisquare",
    adaptive = TRUE,
    diagnostics = "full"
  )
  fit_aicc <- lhwr(
    y ~ x1 + x2,
    data = dat,
    Z = H$Z,
    bw = 30,
    kernel = "bisquare",
    adaptive = TRUE,
    diagnostics = "aicc"
  )
  fit_none <- lhwr(
    y ~ x1 + x2,
    data = dat,
    Z = H$Z,
    bw = 30,
    kernel = "bisquare",
    adaptive = TRUE,
    diagnostics = "none"
  )

  expect_true("x1_SE" %in% names(fit_full$SDF))
  expect_true("Local_R2" %in% names(fit_full$SDF))
  expect_false("x1_SE" %in% names(fit_aicc$SDF))
  expect_false("Local_R2" %in% names(fit_aicc$SDF))
  expect_false("x1_SE" %in% names(fit_none$SDF))
  expect_false("Local_R2" %in% names(fit_none$SDF))

  expect_true(is.finite(AICc.lhwr(fit_full)))
  expect_equal(AICc.lhwr(fit_aicc), AICc.lhwr(fit_full), tolerance = 1e-8)
  expect_true(is.na(AICc.lhwr(fit_none)))
  expect_equal(fit_none$beta_hat, fit_full$beta_hat, tolerance = 1e-10)
  expect_equal(fit_none$y_hat, fit_full$y_hat, tolerance = 1e-10)
  expect_equal(RSS.lhwr(fit_none), RSS.lhwr(fit_full), tolerance = 1e-10)
})

test_that("kNN engine computes finite local R2 in full diagnostics", {
  set.seed(16)
  H <- simulate_H2(80, r_max = 3, uniform_area = TRUE)
  x1 <- rnorm(80)
  x2 <- rnorm(80)
  y <- 0.5 + 1.0 * x1 - 0.4 * x2 + rnorm(80, sd = 0.25)
  dat <- data.frame(y = y, x1 = x1, x2 = x2)

  fit <- lhwr(
    y ~ x1 + x2,
    data = dat,
    Z = H$Z,
    bw = 30,
    kernel = "bisquare",
    adaptive = TRUE,
    geometry = "lorentz",
    engine = "knn",
    diagnostics = "full"
  )

  expect_true("Local_R2" %in% names(fit$SDF))
  expect_true(any(is.finite(fit$SDF$Local_R2)))
  expect_false(all(is.na(fit$SDF$Local_R2)))
})

test_that("dense and kNN engines agree for nearby high-radius locations", {
  set.seed(161)
  n <- 40
  r <- 35
  theta <- seq(0, by = 1e-10, length.out = n)
  Z <- cbind(
    z0 = cosh(r),
    z1 = sinh(r) * cos(theta),
    z2 = sinh(r) * sin(theta)
  )
  x <- rnorm(n)
  y <- 1 + 2 * x + rnorm(n, sd = 0.1)
  dat <- data.frame(y, x)

  fit_dense <- lhwr(
    y ~ x, dat, Z = Z, bw = 12, adaptive = TRUE,
    kernel = "bisquare", engine = "dense"
  )
  fit_knn <- lhwr(
    y ~ x, dat, Z = Z, bw = 12, adaptive = TRUE,
    kernel = "bisquare", engine = "knn"
  )

  expect_equal(fit_knn$beta_hat, fit_dense$beta_hat, tolerance = 1e-10)
  expect_equal(fit_knn$y_hat, fit_dense$y_hat, tolerance = 1e-10)
})

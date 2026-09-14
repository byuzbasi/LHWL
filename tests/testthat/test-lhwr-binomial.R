binomial_lhwr_fixture <- function(n = 80L, seed = 2401L) {
  set.seed(seed)
  H <- simulate_H2(n, r_max = 2.5, uniform_area = TRUE)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  probability <- plogis(-0.35 + 0.9 * x1 - 0.55 * x2)
  y <- rbinom(n, 1, probability)
  if (length(unique(y)) != 2L) y[c(1L, 2L)] <- c(0, 1)
  list(
    data = data.frame(y = y, x1 = x1, x2 = x2),
    Z = H$Z,
    D = hyperbolic_distance_matrix(H$Z, check = FALSE)
  )
}

test_that("Gaussian LHWR remains the default and positional kernel calls work", {
  fixture <- binomial_lhwr_fixture(n = 45, seed = 1)
  dat <- transform(
    fixture$data,
    response = 1 + 1.2 * x1 - 0.4 * x2 + rnorm(45, sd = 0.2)
  )

  original <- lhwr(
    response ~ x1 + x2, dat, fixture$Z, bw = 20,
    diagnostics = "none"
  )
  explicit <- lhwr(
    response ~ x1 + x2, dat, fixture$Z, bw = 20,
    diagnostics = "none", family = "gaussian"
  )
  positional <- lhwr(
    response ~ x1 + x2, dat, fixture$Z, 20, "tricube",
    diagnostics = "none"
  )

  expect_equal(original$beta_hat, explicit$beta_hat, tolerance = 0)
  expect_equal(original$y_hat, explicit$y_hat, tolerance = 0)
  expect_equal(original$LHWR.diagnostic, explicit$LHWR.diagnostic, tolerance = 0)
  expect_equal(positional$LHWR.arguments$kernel, "tricube")
})

test_that("Binomial LHWR accepts numeric, logical, and two-level factor responses", {
  fixture <- binomial_lhwr_fixture(seed = 2)
  dat <- fixture$data
  dat$logical_y <- dat$y == 1
  dat$factor_y <- factor(
    ifelse(dat$y == 1, "event", "nonevent"),
    levels = c("nonevent", "event")
  )
  bandwidth <- max(fixture$D) + 1

  numeric_fit <- lhwr(
    y ~ x1 + x2, dat, fixture$Z, bandwidth,
    kernel = "boxcar", adaptive = FALSE,
    family = "binomial", diagnostics = "none"
  )
  logical_fit <- lhwr(
    logical_y ~ x1 + x2, dat, fixture$Z, bandwidth,
    kernel = "boxcar", adaptive = FALSE,
    family = "binomial", diagnostics = "none"
  )
  factor_fit <- lhwr(
    factor_y ~ x1 + x2, dat, fixture$Z, bandwidth,
    kernel = "boxcar", adaptive = FALSE,
    family = "binomial", diagnostics = "none"
  )

  expect_equal(numeric_fit$beta_hat, logical_fit$beta_hat, tolerance = 1e-10)
  expect_equal(numeric_fit$beta_hat, factor_fit$beta_hat, tolerance = 1e-10)
  expect_equal(factor_fit$response_levels, c("nonevent", "event"))

  invalid <- transform(dat, invalid_y = ifelse(y == 1, 2, 0))
  three_level <- transform(
    dat,
    three = factor(rep(c("a", "b", "c", "a"), length.out = nrow(dat)))
  )
  one_class <- transform(dat, constant = 0)
  grouped <- transform(dat, failure = 1 - y)

  expect_error(
    lhwr(invalid_y ~ x1, invalid, fixture$Z, bandwidth,
         adaptive = FALSE, family = "binomial"),
    "numeric 0/1"
  )
  expect_error(
    lhwr(three ~ x1, three_level, fixture$Z, bandwidth,
         adaptive = FALSE, family = "binomial"),
    "exactly two"
  )
  expect_error(
    lhwr(constant ~ x1, one_class, fixture$Z, bandwidth,
         adaptive = FALSE, family = "binomial"),
    "both response classes"
  )
  expect_error(
    lhwr(cbind(y, failure) ~ x1, grouped, fixture$Z, bandwidth,
         adaptive = FALSE, family = "binomial"),
    "requires a two-level factor"
  )
  expect_error(
    lhwr(y ~ x1 + offset(x2), dat, fixture$Z, bandwidth,
         adaptive = FALSE, family = "binomial"),
    "offsets"
  )
})

test_that("Full-support binomial LHWR matches global glm and Fisher information", {
  fixture <- binomial_lhwr_fixture(n = 90, seed = 3)
  dat <- fixture$data
  bandwidth <- max(fixture$D) + 1
  fit <- lhwr(
    y ~ x1 + x2, dat, fixture$Z, bandwidth,
    kernel = "boxcar", adaptive = FALSE,
    family = "binomial", diagnostics = "full"
  )
  global <- glm(y ~ x1 + x2, dat, family = binomial())
  expected_beta <- matrix(
    coef(global), nrow(dat), length(coef(global)), byrow = TRUE
  )
  design <- model.matrix(global)
  probability <- plogis(as.numeric(design %*% coef(global)))
  expected_se_vector <- sqrt(diag(solve(
    crossprod(design, design * as.numeric(probability * (1 - probability)))
  )))
  expected_se <- matrix(
    expected_se_vector,
    nrow(dat), length(coef(global)), byrow = TRUE
  )
  se_columns <- paste0(fit$coef_names, "_SE")
  z_columns <- paste0(fit$coef_names, "_ZV")

  expect_equal(unname(fit$beta_hat), expected_beta, tolerance = 2e-8)
  expect_equal(unname(fitted(fit)), unname(fitted(global)), tolerance = 2e-8)
  expect_equal(
    unname(as.matrix(fit$SDF[, se_columns])),
    expected_se,
    tolerance = 1e-7
  )
  expect_equal(
    unname(as.matrix(fit$SDF[, z_columns])),
    unname(fit$beta_hat / as.matrix(fit$SDF[, se_columns])),
    tolerance = 1e-10
  )
  expect_true(all(fit$converged))
  expect_true(all(fitted(fit) > 0 & fitted(fit) < 1))
  expect_equal(residuals(fit), dat$y - fitted(fit), tolerance = 0)
  expect_false(any(grepl("_TV$", names(fit$SDF))))
  expect_false("Local_R2" %in% names(fit$SDF))
})

test_that("Binomial LHWR supports no-intercept and intercept-only models", {
  fixture <- binomial_lhwr_fixture(n = 70, seed = 4)
  dat <- fixture$data
  bandwidth <- max(fixture$D) + 1

  no_intercept <- lhwr(
    y ~ 0 + x1 + x2, dat, fixture$Z, bandwidth,
    kernel = "boxcar", adaptive = FALSE,
    family = "binomial", diagnostics = "none"
  )
  global_no_intercept <- glm(
    y ~ 0 + x1 + x2, dat, family = binomial()
  )
  expect_equal(
    unname(no_intercept$beta_hat),
    matrix(coef(global_no_intercept), nrow(dat), 2, byrow = TRUE),
    tolerance = 1e-8
  )

  intercept_only <- lhwr(
    y ~ 1, dat, fixture$Z, bandwidth,
    kernel = "boxcar", adaptive = FALSE,
    family = "binomial", diagnostics = "full"
  )
  expected <- qlogis(mean(dat$y))
  expect_equal(
    as.numeric(intercept_only$beta_hat[, "Intercept"]),
    rep(expected, nrow(dat)),
    tolerance = 1e-8
  )
  expect_true(all(is.finite(intercept_only$SDF$Intercept_SE)))
})

test_that("Binomial dense and exact-kNN engines agree and are deterministic", {
  fixture <- binomial_lhwr_fixture(n = 100, seed = 5)
  dat <- fixture$data
  arguments <- list(
    formula = y ~ x1 + x2,
    data = dat,
    Z = fixture$Z,
    bw = 45,
    kernel = "bisquare",
    adaptive = TRUE,
    family = "binomial",
    diagnostics = "none"
  )
  dense <- suppressWarnings(do.call(lhwr, c(arguments, list(engine = "dense"))))
  knn_one <- suppressWarnings(do.call(
    lhwr, c(arguments, list(engine = "knn", n_threads = 1L))
  ))
  knn_two <- suppressWarnings(do.call(
    lhwr, c(arguments, list(engine = "knn", n_threads = 2L))
  ))

  expect_equal(dense$bandwidth, knn_one$bandwidth, tolerance = 1e-10)
  expect_equal(dense$beta_hat, knn_one$beta_hat, tolerance = 1e-9)
  expect_equal(knn_one$beta_hat, knn_two$beta_hat, tolerance = 0)
  expect_equal(knn_one$y_hat, knn_two$y_hat, tolerance = 0)
  expect_equal(knn_one$converged, knn_two$converged)
})

test_that("Binomial prediction supports response, link, and coefficients", {
  fixture <- binomial_lhwr_fixture(n = 75, seed = 6)
  dat <- fixture$data
  fit <- lhwr(
    y ~ x1 + x2, dat, fixture$Z, bw = 100,
    kernel = "boxcar", adaptive = FALSE,
    family = "binomial", diagnostics = "none"
  )
  set.seed(600)
  Znew <- simulate_H2(4, r_max = 2)$Z
  newdata <- data.frame(
    x1 = c(-1, 0, 0.5, 1),
    x2 = c(0.2, -0.3, 0.1, 0.8)
  )
  coefficients <- predict(
    fit, newdata = newdata, Znew = Znew, type = "coefficients"
  )
  link <- predict(fit, newdata = newdata, Znew = Znew, type = "link")
  response <- predict(
    fit, newdata = newdata, Znew = Znew, type = "response"
  )
  expected <- coef(glm(y ~ x1 + x2, dat, family = binomial()))

  expect_equal(
    unname(coefficients),
    matrix(expected, nrow(newdata), length(expected), byrow = TRUE),
    tolerance = 1e-8
  )
  expect_equal(
    link,
    as.numeric(model.matrix(~ x1 + x2, newdata) %*% expected),
    tolerance = 1e-8
  )
  expect_equal(response, plogis(link), tolerance = 0)
  expect_equal(predict(fit), fitted(fit), tolerance = 0)
  expect_equal(predict(fit, type = "link"), fit$linear_predictor, tolerance = 0)
})

test_that("Binomial diagnostics keep Gaussian-only quantities out of scope", {
  fixture <- binomial_lhwr_fixture(n = 75, seed = 7)
  dat <- fixture$data
  fit <- lhwr(
    y ~ x1 + x2, dat, fixture$Z, bw = max(fixture$D) + 1,
    kernel = "boxcar", adaptive = FALSE,
    family = "binomial", diagnostics = "full"
  )

  expect_named(
    fit$LHWR.diagnostic,
    c(
      "family", "deviance", "null.deviance", "pseudo.R2", "logLik",
      "brier", "converged", "not.converged",
      "class.variation.failures", "kappa.median", "kappa.max",
      "rank.deficient",
      "neff.min", "neff.median", "n.positive.min", "n.positive.median"
    )
  )
  expect_error(AICc.lhwr(fit), "not currently defined")
  expect_error(RSS.lhwr(fit), "Gaussian LHWR")
  expect_error(diagnostic_tests(fit), "only for Gaussian")
  expect_error(local_inference(fit), "only for Gaussian")
  expect_error(
    lhwr(
      y ~ x1, dat, fixture$Z, bw = 30,
      family = "binomial", diagnostics = "aicc"
    ),
    "not defined"
  )
  expect_error(
    lhwr(
      y ~ x1, dat, fixture$Z, bw = 30,
      family = "binomial", f_tests = TRUE
    ),
    "only for Gaussian"
  )
  report_output <- capture.output(
    report <- lhwr_report(
      fit,
      include_inference = FALSE,
      include_f_tests = FALSE,
      include_moran = FALSE,
      include_collinearity = FALSE
    )
  )
  expect_s3_class(report, "lhwr_report")
  expect_true(any(report$global$Metric == "Deviance"))
  expect_true(any(grepl("Global logistic baseline", report_output)))
})

test_that("Binomial bandwidth CV equals manual leave-one-out log loss", {
  fixture <- binomial_lhwr_fixture(n = 55, seed = 8)
  dat <- fixture$data
  bandwidth <- max(fixture$D) + 1
  selected <- bw.lhwr(
    y ~ x1 + x2,
    data = dat,
    Z = fixture$Z,
    kernel = "boxcar",
    adaptive = FALSE,
    bw_grid = c(bandwidth, bandwidth + 1),
    family = "binomial"
  )

  X <- model.matrix(y ~ x1 + x2, dat)
  manual_probability <- vapply(seq_len(nrow(dat)), function(i) {
    global_without_i <- glm.fit(
      x = X[-i, , drop = FALSE],
      y = dat$y[-i],
      family = binomial()
    )
    plogis(sum(X[i, ] * global_without_i$coefficients))
  }, numeric(1L))
  manual_probability <- pmin(
    pmax(manual_probability, 1e-8),
    1 - 1e-8
  )
  manual_score <- -mean(
    dat$y * log(manual_probability) +
      (1 - dat$y) * log1p(-manual_probability)
  )

  expect_equal(selected$family, "binomial")
  expect_equal(selected$approach, "CV")
  expect_equal(selected$selection.table$score, rep(manual_score, 2), tolerance = 1e-7)
  expect_true(all(selected$selection.table$not_converged == 0))
  expect_true(all(is.na(selected$selection.table$trS)))
  expect_true(all(is.na(selected$selection.table$RSS)))
  expect_error(
    bw.lhwr(
      y ~ x1, dat, fixture$Z,
      family = "binomial", approach = "AICc"
    ),
    "supports only"
  )
})

test_that("Single-class local neighborhoods are exposed without hidden ridge", {
  n <- 36
  coordinate <- seq(-3, 3, length.out = n)
  Z <- cbind(cosh(coordinate), sinh(coordinate), 0)
  set.seed(340)
  dat <- data.frame(
    y = c(rep(0, n / 2), rep(1, n / 2)),
    x = rnorm(n)
  )

  expect_warning(
    fit <- lhwr(
      y ~ x,
      data = dat,
      Z = Z,
      bw = 5,
      kernel = "boxcar",
      adaptive = TRUE,
      engine = "knn",
      family = "binomial",
      diagnostics = "full"
    ),
    "did not meet the convergence criterion"
  )
  expect_true(any(!fit$local_class_variation))
  expect_true(all(!fit$converged[!fit$local_class_variation]))
  expect_true(anyNA(fit$SDF$Intercept_SE[!fit$local_class_variation]))
  expect_true(all(is.finite(fitted(fit))))
})

test_that("tuning landscape exposes evaluated and eligible candidates", {
  set.seed(20260906)
  n <- 32L
  H <- simulate_H2(n, r_max = 1.6)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  dat <- data.frame(
    y = 0.4 + x1 - 0.5 * x2 + rnorm(n, sd = 0.35),
    x1 = x1,
    x2 = x2
  )
  tuned <- tune.lhwpr(
    y ~ x1 + x2,
    dat,
    H$Z,
    bw_grid = c(18, 24),
    nlambda = 4,
    nfolds = 2,
    refine = FALSE,
    compare_lhwr = FALSE,
    refit = FALSE
  )

  output <- tempfile(fileext = ".pdf")
  grDevices::pdf(output)
  plotted <- plot(tuned, type = "landscape")
  grDevices::dev.off()

  expect_true(file.exists(output))
  expect_equal(nrow(plotted), nrow(tuned$selection_table))
  expect_true(all(c("plot_eligible", "score_delta") %in% names(plotted)))
  expect_true(plotted$plot_eligible[plotted$selected])
})

test_that("LHWPR focal surfaces distinguish selection from significance", {
  set.seed(20260907)
  n <- 36L
  H <- simulate_H2(n, r_max = 1.8)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  dat <- data.frame(
    y = 0.5 + 1.1 * x1 + rnorm(n, sd = 0.4),
    x1 = x1,
    x2 = x2
  )
  fit <- lhwpr(
    y ~ x1 + x2,
    dat,
    H$Z,
    bw = 20,
    lambda = c(0.25, 0.08),
    keep_path = "selected"
  )

  coefficient_file <- tempfile(fileext = ".pdf")
  selection_file <- tempfile(fileext = ".pdf")
  coefficient <- plot_lhwpr_surface(
    fit, "x1", type = "coefficient", output_file = coefficient_file
  )
  selection <- plot_lhwpr_surface(
    fit, 1L, type = "selection", output_file = selection_file
  )

  expect_true(all(file.exists(c(coefficient_file, selection_file))))
  expect_equal(coefficient$coefficient, fit$beta_hat[, "x1"])
  expect_identical(coefficient$active, selection$active)
  expect_identical(unique(selection$type), "selection")
  expect_false(any(grepl("signific", names(selection), ignore.case = TRUE)))
  expect_error(
    plot_lhwpr_surface(fit, "Intercept", type = "selection"),
    "intercept"
  )
})

###############################################################################
# Reproducible LHWPR microbenchmark
#
# This development helper compares solver and neighborhood engines without
# adding a benchmarking package dependency. It does not write output files.
###############################################################################

benchmark_lhwpr <- function(n = 180L,
                            p = 12L,
                            nlambda = 30L,
                            repetitions = 3L,
                            n_threads = 2L,
                            seed = 20260726L) {
  stopifnot(
    n > 30L,
    p >= 2L,
    p < n,
    nlambda >= 2L,
    repetitions >= 1L,
    n_threads >= 1L
  )
  set.seed(seed)
  geometry <- simulate_H2(n, r_max = 2)
  predictors <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(predictors) <- paste0("x", seq_len(p))
  signal <- seq(1.2, 0.2, length.out = min(5L, p))
  response <- 0.5 + predictors[, seq_along(signal), drop = FALSE] %*% signal +
    rnorm(n, sd = 0.5)
  data <- data.frame(y = as.numeric(response), predictors, check.names = FALSE)
  formula <- stats::reformulate(colnames(predictors), response = "y")
  bandwidth <- max(20L, min(n - 1L, as.integer(round(0.3 * n))))

  configurations <- list(
    dense_naive_1 = list(
      engine = "dense",
      control = lhwpr_control(solver = "naive", n_threads = 1L)
    ),
    dense_covariance_1 = list(
      engine = "dense",
      control = lhwpr_control(solver = "covariance", n_threads = 1L)
    ),
    dense_auto_parallel = list(
      engine = "dense",
      control = lhwpr_control(solver = "auto", n_threads = n_threads)
    ),
    knn_auto_parallel = list(
      engine = "knn",
      control = lhwpr_control(solver = "auto", n_threads = n_threads)
    )
  )

  fits <- vector("list", length(configurations))
  names(fits) <- names(configurations)
  timings <- vector("list", length(configurations))
  for (configuration in seq_along(configurations)) {
    setting <- configurations[[configuration]]
    elapsed <- numeric(repetitions)
    for (iteration in seq_len(repetitions)) {
      gc(FALSE)
      timing <- system.time({
        fit <- lhwpr(
          formula = formula,
          data = data,
          Z = geometry$Z,
          bw = bandwidth,
          alpha = 0.7,
          nlambda = nlambda,
          engine = setting$engine,
          keep_path = "selected",
          control = setting$control
        )
      })
      elapsed[iteration] <- unname(timing["elapsed"])
    }
    fits[[configuration]] <- fit
    timings[[configuration]] <- data.frame(
      configuration = names(configurations)[configuration],
      minimum_seconds = min(elapsed),
      median_seconds = stats::median(elapsed),
      maximum_seconds = max(elapsed),
      repetitions = repetitions,
      stringsAsFactors = FALSE
    )
  }

  reference <- fits$dense_naive_1$beta_hat
  equality <- vapply(
    fits,
    function(fit) max(abs(fit$beta_hat - reference)),
    numeric(1)
  )
  list(
    timings = do.call(rbind, timings),
    maximum_coefficient_difference = equality,
    dimensions = c(n = n, p = p, nlambda = nlambda),
    bandwidth = bandwidth,
    seed = seed
  )
}

message("Loaded benchmark_lhwpr(). Example:")
message("benchmark_lhwpr(n = 180, p = 12, nlambda = 30)")

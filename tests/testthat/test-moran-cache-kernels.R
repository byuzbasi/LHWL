test_that("noncompact Moran and LISA weights never lose positive kernel tails", {
  set.seed(20260914)
  H <- simulate_H2(36, r_max = 2.5)
  x <- rnorm(36)
  y <- 1 + 1.3 * x + rnorm(36, sd = 0.4)
  fit <- lhwr(y ~ x, data.frame(y, x), Z = H$Z, bw = 12,
              kernel = "bisquare", adaptive = TRUE, engine = "knn",
              diagnostics = "aicc")
  full <- fit
  full$knn_index <- NULL
  full$knn_distance <- NULL
  weight_fun <- getFromNamespace("hgwr_moran_weights", "LHWL")
  cache_fun <- getFromNamespace("hgwr_moran_weights_from_stored_knn", "LHWL")
  D <- hyperbolic_distance_matrix(fit$Z, check = FALSE)
  diag(D) <- 0

  for (kernel in c("bisquare", "gaussian", "exponential")) {
    for (row_standardize in c(FALSE, TRUE)) {
      for (k in c(8, 12)) {
        weights <- function(object) weight_fun(object, W = NULL,
          distance = "lorentz", kernel = kernel, bandwidth = NULL,
          adaptive = TRUE, k = k, row_standardize = row_standardize)
        actual <- weights(fit)
        expected <- weights(full)
        radii <- apply(D, 1L, function(d) sort(d[d > 0], partial = k)[k])
        explicit <- t(vapply(seq_len(nrow(D)), function(i) {
          hgwr_kernel(D[i, ], radii[i], kernel = kernel)
        }, numeric(nrow(D))))
        diag(explicit) <- 0
        if (row_standardize) explicit <- explicit / rowSums(explicit)
        attr(explicit, "bandwidth") <- radii
        expect_equal(actual, expected, tolerance = 1e-12)
        expect_equal(actual, explicit, tolerance = 1e-12)
        expect_equal(diag(actual), rep(0, nrow(D)))
        if (kernel != "bisquare") {
          expect_true(all(actual[row(D) != col(D)] > 0))
          expect_null(cache_fun(fit, kernel, bandwidth = NULL, k = k,
                                row_standardize = row_standardize))
        } else {
          expect_equal(cache_fun(fit, kernel, bandwidth = NULL, k = k,
                                 row_standardize = row_standardize), actual)
        }
        global <- function(object) moran.lhwr(object, kernel = kernel, k = k,
          row_standardize = row_standardize, permutations = 9, seed = 413)
        local <- function(object) local_moran.lhwr(object, kernel = kernel, k = k,
          row_standardize = row_standardize, permutations = 9, seed = 413)
        expect_equal(global(fit), global(full), tolerance = 1e-12)
        expect_equal(local(fit), local(full), tolerance = 1e-12)
        expect_equal(global(fit)$I, moran.lhwr(fit, W = explicit,
          row_standardize = row_standardize, permutations = 0)$I, tolerance = 1e-12)
      }
    }
  }
})

test_that("fixed bandwidth and supplied diagnostic weights bypass the kNN cache", {
  set.seed(20260915)
  H <- simulate_H2(28, r_max = 2)
  x <- rnorm(28)
  y <- x + rnorm(28)
  fit <- lhwr(y ~ x, data.frame(y, x), Z = H$Z, bw = 10,
              kernel = "bisquare", adaptive = TRUE, engine = "knn",
              diagnostics = "aicc")
  full <- fit
  full$knn_index <- NULL
  full$knn_distance <- NULL
  W <- matrix(1, 28, 28)
  diag(W) <- 0
  for (kernel in c("bisquare", "gaussian", "exponential")) {
    expect_equal(moran.lhwr(fit, kernel = kernel, adaptive = FALSE,
                           bandwidth = 1.4, permutations = 0),
                 moran.lhwr(full, kernel = kernel, adaptive = FALSE,
                            bandwidth = 1.4, permutations = 0), tolerance = 1e-12)
    expect_equal(local_moran.lhwr(fit, kernel = kernel, adaptive = FALSE,
                                 bandwidth = 1.4, permutations = 0),
                 local_moran.lhwr(full, kernel = kernel, adaptive = FALSE,
                                 bandwidth = 1.4, permutations = 0), tolerance = 1e-12)
    expect_equal(moran.lhwr(fit, W = W, kernel = kernel, permutations = 0)$I,
                 moran.lhwr(full, W = W, permutations = 0)$I)
  }
})

test_that("kernels return expected weights", {
  d <- c(0, 0.5, 1, 2)
  expect_equal(hgwr_gaussian(d, 1)[1], 1)
  expect_equal(hgwr_exponential(d, 1), exp(-d))
  expect_equal(hgwr_bisquare(d, 1)[c(3, 4)], c(0, 0))
  expect_equal(hgwr_tricube(d, 1)[c(3, 4)], c(0, 0))
  expect_equal(hgwr_boxcar(d, 1), c(1, 1, 0, 0))
  expect_equal(hgwr_kernel(d, 1, kernel = "tricube"), hgwr_tricube(d, 1))
  expect_equal(hgwr_kernel(d, 1, kernel = "boxcar"), hgwr_boxcar(d, 1))
  expect_error(hgwr_kernel(d, 0), "positive")
})

test_that("all supported kernels run through the C++ fit path", {
  set.seed(123)
  H <- simulate_H2(20, r_max = 1.5)
  D <- hyperbolic_distance_matrix(H$Z)
  X <- cbind(Intercept = 1, x = rnorm(20))
  y <- rnorm(20)
  bw <- rep(1.5, 20)
  
  for (kernel in c("gaussian", "exponential", "bisquare", "tricube", "boxcar")) {
    fit <- hgwr.fit(y, X, D, bw, kernel = kernel)
    expect_equal(dim(fit$beta_hat), c(20L, 2L))
    expect_true(all(is.finite(fit$y_hat)))
  }
})

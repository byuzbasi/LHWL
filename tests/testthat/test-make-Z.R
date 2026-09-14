test_that("feature-map construction creates valid Lorentz coordinates", {
  X <- data.frame(
    score_a = c(10, 12, 14, 16),
    score_b = c(1, 4, 9, 16)
  )
  Z <- make_lhwr_Z_features(X, scale = TRUE)

  expect_equal(dim(Z), c(4L, 3L))
  expect_equal(colnames(Z), c("z0", "z1", "z2"))
  expect_true(check_lhwr_Z(Z))
  expect_true(all(check_hyperboloid(Z)))
  expect_equal(attr(Z, "lhwr_Z_source"), "features")
  expect_equal(attr(Z, "feature_names"), names(X))
})

test_that("feature-map construction can keep original component scale", {
  X <- matrix(c(0, 0, 1, 2, -3, 4), ncol = 2, byrow = TRUE)
  Z <- make_lhwr_Z_features(X, scale = FALSE)

  expect_equal(unname(Z[, 2:3]), X, tolerance = 1e-12)
  expect_true(check_lhwr_Z(Z))
})

test_that("feature-map construction rejects invalid raw inputs", {
  expect_error(make_lhwr_Z_features(matrix(1:9, ncol = 3)), "exactly 2 columns")
  expect_error(make_lhwr_Z_features(data.frame(a = 1:3, b = letters[1:3])), "numeric")
  expect_error(make_lhwr_Z_features(matrix(c(1, NA, 2, 3), ncol = 2)), "missing")
})

test_that("constant feature components are handled explicitly", {
  X <- cbind(constant = rep(2, 4), varying = 1:4)
  expect_warning(
    Z <- make_lhwr_Z_features(X, scale = TRUE),
    "Constant feature component"
  )
  expect_equal(unname(Z[, "z1"]), rep(0, 4))
  expect_true(check_lhwr_Z(Z))
})

test_that("network spectral construction creates valid Lorentz coordinates", {
  A <- matrix(c(
    0, 1, 1, 0, 0,
    1, 0, 1, 1, 0,
    1, 1, 0, 1, 1,
    0, 1, 1, 0, 1,
    0, 0, 1, 1, 0
  ), nrow = 5, byrow = TRUE)
  Z <- make_lhwr_Z_network(A, method = "spectral")

  expect_equal(dim(Z), c(5L, 3L))
  expect_true(check_lhwr_Z(Z))
  expect_equal(attr(Z, "lhwr_Z_source"), "network")
  expect_equal(attr(Z, "network_method"), "spectral")
})

test_that("network shortest-path MDS construction creates valid Lorentz coordinates", {
  A <- matrix(c(
    0, 2, 0, 0,
    2, 0, 1, 0,
    0, 1, 0, 3,
    0, 0, 3, 0
  ), nrow = 4, byrow = TRUE)
  Z <- make_lhwr_Z_network(A, method = "shortest_path_mds", scale = TRUE)

  expect_equal(dim(Z), c(4L, 3L))
  expect_true(check_lhwr_Z(Z))
  expect_equal(attr(Z, "network_method"), "shortest_path_mds")
})

test_that("network construction handles asymmetry and disconnection explicitly", {
  A_asym <- matrix(c(
    0, 1, 0,
    0, 0, 1,
    1, 0, 0
  ), nrow = 3, byrow = TRUE)
  expect_warning(
    Z <- make_lhwr_Z_network(A_asym, method = "spectral", symmetrize = TRUE),
    "not symmetric"
  )
  expect_true(check_lhwr_Z(Z))

  A_disc <- matrix(c(
    0, 1, 0, 0,
    1, 0, 0, 0,
    0, 0, 0, 1,
    0, 0, 1, 0
  ), nrow = 4, byrow = TRUE)
  expect_error(
    make_lhwr_Z_network(A_disc, method = "shortest_path_mds"),
    "disconnected"
  )
  expect_warning(
    Z_imp <- make_lhwr_Z_network(
      A_disc,
      method = "shortest_path_mds",
      disconnected = "impute"
    ),
    "disconnected"
  )
  expect_true(check_lhwr_Z(Z_imp))
})

test_that("generic constructor and compatibility aliases work", {
  X <- matrix(c(0, 0, 1, 1, -1, 2), ncol = 2, byrow = TRUE)
  Z1 <- make_lhwr_Z(X, method = "features", scale = FALSE)
  Z2 <- make_lorentz_Z_features(X, scale = FALSE)
  Z3 <- make_lorentz_Z(X, method = "features", scale = FALSE)

  expect_equal(unname(Z1), unname(Z2), tolerance = 1e-12)
  expect_equal(unname(Z1), unname(Z3), tolerance = 1e-12)
  expect_true(check_lorentz_Z(Z1))
})

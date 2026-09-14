test_that("project_to_H2 creates valid hyperboloid points", {
  S <- matrix(c(0, 0, 1, 0, 0, 2), ncol = 2, byrow = TRUE)
  Z <- project_to_H2(S)
  expect_true(all(check_hyperboloid(Z)))
  expect_equal(as.numeric(lorentz_norm(Z)), rep(-1, nrow(Z)), tolerance = 1e-10)
})

test_that("Poincare coordinates convert to valid Lorentz points", {
  P <- matrix(c(0, 0, 0.2, -0.1, -0.3, 0.4), ncol = 2, byrow = TRUE)
  Z <- poincare_to_lorentz(P)
  
  expect_equal(dim(Z), c(3L, 3L))
  expect_equal(unname(Z[1, ]), c(1, 0, 0), tolerance = 1e-12)
  expect_true(all(check_hyperboloid(Z)))
  expect_equal(as.numeric(lorentz_norm(Z)), rep(-1, nrow(Z)), tolerance = 1e-10)
})

test_that("prepare_lorentz_coordinates supports lorentz, poincare, and spatial input", {
  S <- matrix(c(0, 0, 1, 0, 0, 2), ncol = 2, byrow = TRUE)
  P <- matrix(c(0, 0, 0.2, -0.1, -0.3, 0.4), ncol = 2, byrow = TRUE)
  
  Z_spatial <- prepare_lorentz_coordinates(S, model = "spatial")
  Z_poincare <- prepare_lorentz_coordinates(P, model = "poincare")
  Z_lorentz <- prepare_lorentz_coordinates(Z_poincare, model = "lorentz")
  
  expect_equal(
    unname(as.vector(Z_spatial)),
    unname(as.vector(project_to_H2(S))),
    tolerance = 1e-12
  )
  expect_equal(
    unname(as.vector(Z_lorentz)),
    unname(as.vector(Z_poincare)),
    tolerance = 1e-12
  )
  expect_true(all(check_hyperboloid(Z_spatial)))
  expect_true(all(check_hyperboloid(Z_poincare)))
  expect_equal(attr(Z_spatial, "input_model"), "spatial")
  expect_equal(attr(Z_poincare, "input_model"), "poincare")
  expect_true(all(attr(Z_lorentz, "valid_hyperboloid")))
})

test_that("coordinate preparation rejects invalid inputs", {
  bad_poincare <- matrix(c(1, 0, 0.2, 0.1), ncol = 2, byrow = TRUE)
  bad_lorentz <- matrix(c(1, 1, 0), nrow = 1)
  
  expect_error(poincare_to_lorentz(bad_poincare), "Poincare disk")
  expect_error(
    prepare_lorentz_coordinates(bad_lorentz, model = "lorentz"),
    "upper Lorentz hyperboloid"
  )
  
  checked <- prepare_lorentz_coordinates(
    bad_lorentz,
    model = "lorentz",
    stop_on_fail = FALSE
  )
  expect_false(all(attr(checked, "valid_hyperboloid")))
})

test_that("hyperbolic distances are stable and symmetric", {
  H <- simulate_H2(8, r_max = 2)
  D <- hyperbolic_distance_matrix(H$Z)
  expect_equal(D, t(D), tolerance = 1e-10)
  expect_equal(diag(D), rep(0, 8), tolerance = 1e-10)
  expect_true(all(D >= 0))
})

test_that("Lorentz kNN matches dense distance ordering", {
  set.seed(19)
  H <- simulate_H2(35, r_max = 4)
  D <- hyperbolic_distance_matrix(H$Z, check = FALSE)
  k <- 7
  knn <- lorentz_knn(H$Z, k = k, include_self = FALSE)
  
  expect_equal(dim(knn$index), c(35L, k))
  expect_equal(dim(knn$distance), c(35L, k))
  for (i in seq_len(nrow(D))) {
    ord <- order(D[i, ], seq_along(D[i, ]))
    ord <- ord[ord != i][seq_len(k)]
    expect_equal(unname(knn$index[i, ]), ord)
    expect_equal(unname(knn$distance[i, ]), unname(D[i, ord]), tolerance = 1e-10)
  }
})

test_that("Lorentz kNN can include self as the first neighbor", {
  set.seed(20)
  H <- simulate_H2(15, r_max = 3)
  knn <- lorentz_knn(H$Z, k = 4, include_self = TRUE)
  
  expect_equal(knn$index[, 1], seq_len(15))
  expect_equal(knn$distance[, 1], rep(0, 15), tolerance = 1e-12)
  expect_true(all(knn$distance[, -1, drop = FALSE] >= 0))
})

test_that("hyperbolic distance uses acosh of negative Lorentz inner product", {
  H <- simulate_H2(10, r_max = 2)
  z <- H$Z[3, ]
  w <- H$Z[8, ]
  arg <- -lorentz_inner(z, w)
  
  expect_true(arg >= 1)
  expect_equal(
    hyperbolic_distance(z, w),
    acosh(arg),
    tolerance = 1e-12
  )
})

test_that("hyperboloid checks tolerate r_max 10 simulation roundoff", {
  set.seed(3)
  H <- simulate_H2(200, r_max = 10)
  expect_true(all(check_hyperboloid(H$Z)))
})

test_that("hyperboloid checks use scale-aware tolerance for large radii", {
  set.seed(4)
  H <- simulate_H2(300, r_max = 12)
  expect_true(all(check_hyperboloid(H$Z)))
})

test_that("distance matrix is stable for large hyperbolic radii", {
  set.seed(5)
  H <- simulate_H2(120, r_max = 15)
  D <- hyperbolic_distance_matrix(H$Z, check = TRUE)
  expect_equal(D, t(D), tolerance = 1e-10)
  expect_equal(diag(D), rep(0, 120), tolerance = 1e-10)
  expect_true(all(is.finite(D)))
  expect_true(all(D >= 0))
})

test_that("Lorentz distances remain accurate for nearby high-radius points", {
  for (r in c(20, 35)) {
    delta <- 1e-10
    Z <- rbind(
      c(cosh(r), sinh(r), 0),
      c(cosh(r), sinh(r) * cos(delta), sinh(r) * sin(delta))
    )
    expected <- 2 * asinh(sinh(r) * abs(sin(delta / 2)))

    expect_true(all(check_hyperboloid(Z)))
    expect_equal(hyperbolic_distance(Z[1, ], Z[2, ]), expected, tolerance = 1e-8)
    expect_equal(
      hyperbolic_distance_matrix(Z, check = FALSE)[1, 2],
      expected,
      tolerance = 1e-8
    )
    expect_equal(
      lorentz_knn(Z, k = 1, check = FALSE)$distance[1, 1],
      expected,
      tolerance = 1e-8
    )
  }
})

test_that("tangent-space coordinates and distances are well defined", {
  set.seed(12)
  H <- simulate_H2(30, r_max = 3)
  mu <- extrinsic_lorentz_mean(H$Z)
  V <- lorentz_log_map(H$Z, mu = mu)
  Q <- tangent_coordinates(H$Z, mu = mu)
  D <- tangent_distance_matrix(H$Z, mu = mu)
  D_ref <- as.matrix(stats::dist(Q))
  
  expect_equal(as.numeric(lorentz_norm(matrix(mu, nrow = 1))), -1, tolerance = 1e-8)
  tangent_ip <- -V[, 1] * mu[1] + V[, 2] * mu[2] + V[, 3] * mu[3]
  expect_equal(tangent_ip, rep(0, nrow(V)), tolerance = 1e-8)
  expect_equal(dim(Q), c(30L, 2L))
  expect_equal(dim(D), c(30L, 30L))
  expect_equal(unname(D), unname(D_ref), tolerance = 1e-8)
  expect_true(all(is.finite(D)))
  expect_equal(D, t(D), tolerance = 1e-12)
  expect_equal(unname(diag(D)), rep(0, 30), tolerance = 1e-12)
})

test_that("tangent geometry remains finite at large hyperbolic radii", {
  for (r in c(20, 35)) {
    delta <- 1e-10
    Z <- rbind(
      c(cosh(r), sinh(r), 0),
      c(cosh(r), sinh(r) * cos(delta), sinh(r) * sin(delta))
    )
    mu <- extrinsic_lorentz_mean(Z)
    V <- lorentz_log_map(Z, mu = mu)
    Q <- tangent_coordinates(Z, mu = mu)
    D <- tangent_distance_matrix(Z, mu = mu)

    expect_true(all(is.finite(mu)))
    expect_true(all(check_hyperboloid(mu)))
    expect_true(all(is.finite(V)))
    expect_true(all(is.finite(Q)))
    expect_true(all(is.finite(D)))
    expect_equal(unname(D), unname(as.matrix(stats::dist(Q))), tolerance = 1e-6)
    expect_equal(
      D[1, 2],
      hyperbolic_distance(Z[1, ], Z[2, ], check = FALSE),
      tolerance = 1e-8
    )
    expect_equal(D, t(D), tolerance = 1e-12)
  }
})

test_that("C++ distance matrix matches the R reference implementation", {
  ref_file <- test_path("../../../R_reference/geometry_R_reference.R")
  skip_if_not(file.exists(ref_file), "R reference files are local development artifacts")
  ref <- new.env(parent = asNamespace("LHWL"))
  sys.source(ref_file, envir = ref)
  
  set.seed(6)
  H <- simulate_H2(80, r_max = 10)
  D_r <- ref$hyperbolic_distance_matrix(H$Z, check = FALSE)
  D_cpp <- hyperbolic_distance_matrix(H$Z, check = FALSE)
  
  expect_equal(D_cpp, D_r, tolerance = 1e-8)
})

test_that("C++ Lorentz inner product matches the R reference implementation", {
  ref_file <- test_path("../../../R_reference/geometry_R_reference.R")
  skip_if_not(file.exists(ref_file), "R reference files are local development artifacts")
  ref <- new.env(parent = asNamespace("LHWL"))
  sys.source(ref_file, envir = ref)
  
  set.seed(7)
  H <- simulate_H2(40, r_max = 5)
  z <- H$Z[1:40, ]
  w <- H$Z[40:1, ]
  
  expect_equal(lorentz_inner(z, w), ref$lorentz_inner(z, w), tolerance = 1e-12)
})

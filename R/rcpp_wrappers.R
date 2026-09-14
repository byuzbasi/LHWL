.lorentz_inner_cpp <- function(z, w) {
  .Call(`_LHWL_lorentz_inner_cpp`, z, w)
}

.hyperbolic_distance_cpp <- function(z, w, tol = 1e-6) {
  .Call(`_LHWL_hyperbolic_distance_cpp`, z, w, tol)
}

.hyperbolic_distance_matrix_cpp <- function(Z, W = NULL, symmetric = is.null(W),
                                            tol = 1e-6) {
  if (is.null(W)) W <- Z
  .Call(`_LHWL_hyperbolic_distance_matrix_cpp`, Z, W, symmetric, tol)
}

.lorentz_knn_cpp <- function(Z, k, include_self = FALSE, tol = 1e-6) {
  .Call(`_LHWL_lorentz_knn_cpp`, Z, as.integer(k), include_self, tol)
}

.tangent_distance_matrix_cpp <- function(Z, mu, tol = 1e-10) {
  .Call(`_LHWL_tangent_distance_matrix_cpp`, Z, mu, tol)
}

.hgwr_fit_cpp <- function(y, X, D, bw, kernel = "bisquare", return_hat = TRUE) {
  .Call(`_LHWL_hgwr_fit_cpp`, y, X, D, bw, kernel, return_hat)
}

.hgwr_fit_knn_cpp <- function(y, X, Z, k, kernel = "bisquare",
                              return_hat = TRUE, tol = 1e-6,
                              n_threads = 1L) {
  .Call(
    `_LHWL_hgwr_fit_knn_cpp`,
    y,
    X,
    Z,
    as.integer(k),
    kernel,
    return_hat,
    tol,
    as.integer(n_threads)
  )
}

.hgwr_cv_score_cpp <- function(y, X, D, bw, kernel = "bisquare") {
  .Call(`_LHWL_hgwr_cv_score_cpp`, y, X, D, bw, kernel)
}

.hgwr_f_tests_cpp <- function(y, X, D, bw, kernel, beta_hat, RSS_full,
                              RSS_ols, DF_ols, trS, trStS) {
  .Call(
    `_LHWL_hgwr_f_tests_cpp`,
    y, X, D, bw, kernel, beta_hat, RSS_full, RSS_ols, DF_ols, trS, trStS
  )
}

.hgwr_f_tests_knn_cpp <- function(y, X, Z, k, kernel, beta_hat, RSS_full,
                                  RSS_ols, DF_ols, trS, trStS, tol = 1e-6) {
  .Call(
    `_LHWL_hgwr_f_tests_knn_cpp`,
    y, X, Z, as.integer(k), kernel, beta_hat, RSS_full, RSS_ols, DF_ols,
    trS, trStS, tol
  )
}

.hgwr_f_tests_knn_precomputed_cpp <- function(y, X, knn_index, knn_distance,
                                              kernel, beta_hat, RSS_full,
                                              RSS_ols, DF_ols, trS, trStS) {
  .Call(
    `_LHWL_hgwr_f_tests_knn_precomputed_cpp`,
    y, X, knn_index, knn_distance, kernel, beta_hat, RSS_full, RSS_ols,
    DF_ols, trS, trStS
  )
}

.local_collinearity_knn_cpp <- function(X, knn_index, knn_distance, bandwidth,
                                        kernel, center = TRUE, scale = TRUE,
                                        n_threads = 1L) {
  .Call(
    `_LHWL_local_collinearity_knn_cpp`,
    X, knn_index, knn_distance, bandwidth, kernel, center, scale,
    as.integer(n_threads)
  )
}

.bw_hgwr_grid_cpp <- function(y, X, D, bw_grid, approach, kernel, adaptive,
                              bandwidth_type, pilot_k) {
  .Call(
    `_LHWL_bw_hgwr_grid_cpp`,
    y, X, D, bw_grid, approach, kernel, adaptive, bandwidth_type, pilot_k
  )
}

.bw_hgwr_knn_grid_cpp <- function(y, X, Z, bw_grid, approach, kernel,
                                  tol = 1e-6, n_threads = 1L) {
  .Call(
    `_LHWL_bw_hgwr_knn_grid_cpp`,
    y, X, Z, as.numeric(bw_grid), approach, kernel, tol, as.integer(n_threads)
  )
}

.bw_hgwr_knn_golden_cpp <- function(y, X, Z, approach, kernel, lower, upper,
                                    tol, max_iter, dist_tol = 1e-6,
                                    n_threads = 1L) {
  .Call(
    `_LHWL_bw_hgwr_knn_golden_cpp`,
    y, X, Z, approach, kernel, lower, upper, tol, as.integer(max_iter),
    dist_tol, as.integer(n_threads)
  )
}

.bw_hgwr_golden_cpp <- function(y, X, D, approach, kernel, adaptive,
                                bandwidth_type, pilot_k, lower, upper,
                                tol, max_iter) {
  .Call(
    `_LHWL_bw_hgwr_golden_cpp`,
    y, X, D, approach, kernel, adaptive, bandwidth_type, pilot_k,
    lower, upper, tol, max_iter
  )
}

.lhwpr_fit_dense_cpp <- function(y, X, X_focus, D, bw, lambda,
                                 penalty_factor, case_weights, options) {
  .Call(
    `_LHWL_lhwpr_fit_dense_cpp`,
    y, X, X_focus, D, bw, lambda, penalty_factor, case_weights, options
  )
}

.lhwpr_fit_knn_cpp <- function(y, X, X_focus, knn_index, knn_distance,
                               lambda, penalty_factor, case_weights, options) {
  .Call(
    `_LHWL_lhwpr_fit_knn_cpp`,
    y, X, X_focus, knn_index, knn_distance, lambda, penalty_factor,
    case_weights, options
  )
}

.lhwpr_lambda_max_dense_cpp <- function(y, X, D, bw, penalty_factor,
                                        case_weights, options) {
  .Call(
    `_LHWL_lhwpr_lambda_max_dense_cpp`,
    y, X, D, bw, penalty_factor, case_weights, options
  )
}

.lhwpr_lambda_max_knn_cpp <- function(y, X, knn_index, knn_distance,
                                      penalty_factor, case_weights, options) {
  .Call(
    `_LHWL_lhwpr_lambda_max_knn_cpp`,
    y, X, knn_index, knn_distance, penalty_factor, case_weights, options
  )
}

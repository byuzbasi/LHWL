#ifndef LHWL_LHWPR_SOLVER_H
#define LHWL_LHWPR_SOLVER_H

#include "lhwpr_state.h"

namespace lhwl {

inline arma::vec solve_symmetric(const arma::mat& A, const arma::vec& b) {
  arma::vec out;
  bool ok = arma::solve(out, A, b, arma::solve_opts::likely_sympd);
  if (ok && out.is_finite()) return out;
  ok = arma::solve(out, A, b);
  if (ok && out.is_finite()) return out;
  return arma::pinv(A) * b;
}

inline arma::uvec unpenalized_columns(const arma::vec& factors) {
  return arma::find(factors == 0.0);
}

inline void initialize_gaussian_null(const arma::mat& X,
                                     const arma::vec& y,
                                     const arma::vec& weights,
                                     const arma::vec& factors,
                                     const SolverControl& control,
                                     LocalState& state) {
  state.intercept = 0.0;
  state.beta.zeros();
  const arma::uvec free = unpenalized_columns(factors);
  const arma::uword q = free.n_elem + (control.intercept ? 1 : 0);
  if (q == 0) return;

  arma::mat A(X.n_rows, q, arma::fill::none);
  arma::uword col = 0;
  if (control.intercept) {
    A.col(col++).ones();
  }
  for (arma::uword k = 0; k < free.n_elem; ++k) {
    A.col(col++) = X.col(free[k]);
  }
  const arma::mat Aw = A.each_col() % weights;
  arma::mat H = A.t() * Aw;
  H.diag() += 1e-12;
  const arma::vec theta = solve_symmetric(H, A.t() * (weights % y));
  col = 0;
  if (control.intercept) state.intercept = theta[col++];
  for (arma::uword k = 0; k < free.n_elem; ++k) {
    state.beta[free[k]] = theta[col++];
    state.ever_active[free[k]] = 1;
  }
}

inline double binomial_loss(const arma::vec& y,
                            const arma::vec& eta,
                            const arma::vec& weights) {
  double out = 0.0;
  for (arma::uword i = 0; i < y.n_elem; ++i) {
    out += weights[i] * (stable_log1pexp(eta[i]) - y[i] * eta[i]);
  }
  return out;
}

inline void initialize_binomial_null(const arma::mat& X,
                                     const arma::vec& y,
                                     const arma::vec& weights,
                                     const arma::vec& factors,
                                     const SolverControl& control,
                                     LocalState& state) {
  state.intercept = 0.0;
  state.beta.zeros();
  const arma::uvec free = unpenalized_columns(factors);
  const arma::uword q = free.n_elem + (control.intercept ? 1 : 0);
  if (q == 0) return;

  arma::mat A(X.n_rows, q, arma::fill::none);
  arma::uword col = 0;
  if (control.intercept) {
    A.col(col++).ones();
  }
  for (arma::uword k = 0; k < free.n_elem; ++k) {
    A.col(col++) = X.col(free[k]);
  }

  arma::vec theta(q, arma::fill::zeros);
  if (control.intercept) {
    const double mean_y = clamp_probability(arma::dot(weights, y));
    theta[0] = std::log(mean_y / (1.0 - mean_y));
  }

  arma::vec eta = A * theta;
  double old_loss = binomial_loss(y, eta, weights);
  for (int iteration = 0; iteration < control.irls_max_iter; ++iteration) {
    arma::vec probability(y.n_elem, arma::fill::none);
    arma::vec h(y.n_elem, arma::fill::none);
    for (arma::uword i = 0; i < y.n_elem; ++i) {
      probability[i] = clamp_probability(stable_sigmoid(eta[i]));
      h[i] = weights[i] *
        std::max(probability[i] * (1.0 - probability[i]), 1e-8);
    }
    const arma::vec gradient = A.t() * (weights % (y - probability));
    arma::mat H = A.t() * (A.each_col() % h);
    H.diag() += 1e-10;
    arma::vec delta = solve_symmetric(H, gradient);
    if (!delta.is_finite()) break;

    double step = 1.0;
    arma::vec candidate = theta + delta;
    arma::vec eta_candidate = A * candidate;
    double candidate_loss = binomial_loss(y, eta_candidate, weights);
    while (candidate_loss > old_loss && step > 1e-8) {
      step *= 0.5;
      candidate = theta + step * delta;
      eta_candidate = A * candidate;
      candidate_loss = binomial_loss(y, eta_candidate, weights);
    }
    theta = candidate;
    eta = eta_candidate;
    old_loss = candidate_loss;
    if (arma::abs(step * delta).max() <=
        control.irls_tol * (1.0 + arma::abs(theta).max())) {
      break;
    }
  }

  col = 0;
  if (control.intercept) state.intercept = theta[col++];
  for (arma::uword k = 0; k < free.n_elem; ++k) {
    state.beta[free[k]] = theta[col++];
    state.ever_active[free[k]] = 1;
  }
}

inline double model_objective(const arma::mat& X,
                              const arma::vec& y,
                              const arma::vec& weights,
                              const LocalState& state,
                              double lambda,
                              const arma::vec& factors,
                              const SolverControl& control) {
  const arma::vec eta = linear_predictor(X, state.intercept, state.beta);
  double loss = 0.0;
  if (control.family == Family::Gaussian) {
    loss = 0.5 * arma::dot(weights, arma::square(y - eta));
  } else {
    loss = binomial_loss(y, eta, weights);
  }
  return loss + penalty_value(state.beta, lambda, factors, control);
}

inline WlsResult coordinate_descent_wls(const arma::mat& X,
                                        const arma::vec& response,
                                        const arma::vec& weights,
                                        const arma::vec& factors,
                                        double lambda,
                                        const SolverControl& control,
                                        LocalState& state) {
  WlsResult result;
  arma::vec residual = response -
    linear_predictor(X, state.intercept, state.beta);
  const double intercept_curvature = arma::accu(weights);
  arma::mat weighted_square = arma::square(X);
  weighted_square.each_col() %= weights;
  const arma::vec curvature = arma::sum(weighted_square, 0).t();

  const int repair_limit = static_cast<int>(X.n_cols) + 2;
  bool coordinate_converged = false;
  for (int repair = 0; repair < repair_limit; ++repair) {
    coordinate_converged = false;
    for (int iteration = 0; iteration < control.max_iter; ++iteration) {
      double max_change = 0.0;
      if (control.intercept && intercept_curvature > 0.0) {
        const double change = arma::dot(weights, residual) /
          intercept_curvature;
        state.intercept += change;
        residual -= change;
        max_change = std::abs(change) * std::sqrt(intercept_curvature);
      }

      for (arma::uword j = 0; j < X.n_cols; ++j) {
        if (state.active[j] == 0 || curvature[j] <= 0.0) continue;
        const double old_beta = state.beta[j];
        const double partial = arma::dot(
          weights % X.col(j), residual
        ) + curvature[j] * old_beta;
        const double new_beta = coordinate_update(
          partial, curvature[j], lambda, factors[j], old_beta, control
        );
        const double change = new_beta - old_beta;
        if (change != 0.0) {
          state.beta[j] = new_beta;
          residual -= X.col(j) * change;
          state.ever_active[j] = 1;
          max_change = std::max(
            max_change,
            std::abs(change) * std::sqrt(curvature[j])
          );
        }
      }
      ++result.iterations;
      const double scale = 1.0 + std::max(
        std::abs(state.intercept),
        state.beta.n_elem ? arma::abs(state.beta).max() : 0.0
      );
      if (max_change <= control.tol * scale) {
        coordinate_converged = true;
        break;
      }
    }

    const arma::uvec active_before = state.active;
    result.kkt_max = full_kkt_scan(
      state, residual, weights, X, factors, lambda, control, true
    );
    const bool added = arma::accu(state.active != active_before) > 0;
    const double gradient_scale = state.gradient.n_elem ?
      arma::abs(state.gradient).max() : 0.0;
    const double kkt_limit = control.kkt_tol * (1.0 + gradient_scale);
    if (!added && result.kkt_max <= kkt_limit) break;
    if (control.max_active > 0 &&
        static_cast<int>(arma::accu(state.active)) > control.max_active) {
      coordinate_converged = false;
      break;
    }
  }

  const double gradient_scale = state.gradient.n_elem ?
    arma::abs(state.gradient).max() : 0.0;
  result.converged = coordinate_converged &&
    result.kkt_max <= control.kkt_tol * (1.0 + gradient_scale);
  return result;
}

inline WlsResult coordinate_descent_covariance(
    const arma::mat& gram,
    const arma::vec& cross_product,
    const arma::vec& factors,
    double lambda,
    const SolverControl& control,
    LocalState& state) {
  WlsResult result;
  arma::vec gradient = cross_product - gram * state.beta;
  const arma::vec curvature = gram.diag();
  const int repair_limit = static_cast<int>(gram.n_cols) + 2;
  bool coordinate_converged = false;

  for (int repair = 0; repair < repair_limit; ++repair) {
    coordinate_converged = false;
    for (int iteration = 0; iteration < control.max_iter; ++iteration) {
      double max_change = 0.0;
      for (arma::uword j = 0; j < gram.n_cols; ++j) {
        if (state.active[j] == 0 || curvature[j] <= 0.0) continue;
        const double old_beta = state.beta[j];
        const double partial = gradient[j] + curvature[j] * old_beta;
        const double new_beta = coordinate_update(
          partial, curvature[j], lambda, factors[j], old_beta, control
        );
        const double change = new_beta - old_beta;
        if (change != 0.0) {
          state.beta[j] = new_beta;
          gradient -= gram.col(j) * change;
          state.ever_active[j] = 1;
          max_change = std::max(
            max_change,
            std::abs(change) * std::sqrt(curvature[j])
          );
        }
      }
      ++result.iterations;
      const double scale = 1.0 +
        (state.beta.n_elem ? arma::abs(state.beta).max() : 0.0);
      if (max_change <= control.tol * scale) {
        coordinate_converged = true;
        break;
      }
    }

    const arma::uvec active_before = state.active;
    result.kkt_max = full_kkt_gradient_scan(
      state, gradient, factors, lambda, control, true
    );
    const bool added = arma::accu(state.active != active_before) > 0;
    const double gradient_scale = state.gradient.n_elem ?
      arma::abs(state.gradient).max() : 0.0;
    const double kkt_limit = control.kkt_tol * (1.0 + gradient_scale);
    if (!added && result.kkt_max <= kkt_limit) break;
    if (control.max_active > 0 &&
        static_cast<int>(arma::accu(state.active)) > control.max_active) {
      coordinate_converged = false;
      break;
    }
  }
  const double gradient_scale = state.gradient.n_elem ?
    arma::abs(state.gradient).max() : 0.0;
  result.converged = coordinate_converged &&
    result.kkt_max <= control.kkt_tol * (1.0 + gradient_scale);
  return result;
}

inline double lambda_max_one(const arma::mat& X,
                             const arma::vec& y,
                             const arma::vec& weights,
                             const arma::vec& factors,
                             const SolverControl& control) {
  LocalState state(X.n_cols);
  if (control.family == Family::Gaussian) {
    initialize_gaussian_null(X, y, weights, factors, control, state);
  } else {
    initialize_binomial_null(X, y, weights, factors, control, state);
  }
  const arma::vec eta = linear_predictor(X, state.intercept, state.beta);
  arma::vec residual(y.n_elem, arma::fill::none);
  if (control.family == Family::Gaussian) {
    residual = y - eta;
  } else {
    for (arma::uword i = 0; i < y.n_elem; ++i) {
      residual[i] = y[i] - clamp_probability(stable_sigmoid(eta[i]));
    }
  }
  const arma::vec gradient = score_vector(X, residual, weights);
  double maximum = 0.0;
  for (arma::uword j = 0; j < X.n_cols; ++j) {
    if (factors[j] <= 0.0) continue;
    double candidate = std::abs(gradient[j]) / factors[j];
    if (control.alpha > 0.0) {
      candidate /= control.alpha;
    } else {
      candidate /= control.ridge_path_tol;
    }
    if (R_finite(candidate)) maximum = std::max(maximum, candidate);
  }
  if (!R_finite(maximum) || maximum <= 0.0) return 1.0;
  return maximum;
}

inline PathResult allocate_path_result(arma::uword p, arma::uword nlambda) {
  PathResult out;
  out.beta = arma::mat(p, nlambda, arma::fill::zeros);
  out.intercept = arma::vec(nlambda, arma::fill::zeros);
  out.fitted = arma::vec(nlambda, arma::fill::zeros);
  out.linear_predictor = arma::vec(nlambda, arma::fill::zeros);
  out.iterations = arma::ivec(nlambda, arma::fill::zeros);
  out.irls_iterations = arma::ivec(nlambda, arma::fill::zeros);
  out.converged = arma::uvec(nlambda, arma::fill::zeros);
  out.kkt = arma::vec(nlambda, arma::fill::value(R_PosInf));
  out.n_active = arma::ivec(nlambda, arma::fill::zeros);
  return out;
}

inline PathResult fit_gaussian_path(const arma::mat& X,
                                    const arma::rowvec& x_focus,
                                    const arma::vec& y,
                                    const arma::vec& weights,
                                    const arma::vec& lambda,
                                    const arma::vec& factors,
                                    const SolverControl& control) {
  PathResult out = allocate_path_result(X.n_cols, lambda.n_elem);
  LocalState state(X.n_cols);
  initialize_gaussian_null(X, y, weights, factors, control, state);
  arma::vec residual = y -
    linear_predictor(X, state.intercept, state.beta);
  state.gradient = score_vector(X, residual, weights);

  for (arma::uword l = 0; l < lambda.n_elem; ++l) {
    initialize_active_set(
      state,
      factors,
      lambda[l],
      l == 0 ? lambda[l] : lambda[l - 1],
      control,
      l == 0
    );
    const WlsResult fit = coordinate_descent_wls(
      X, y, weights, factors, lambda[l], control, state
    );
    const double eta_focus =
      state.intercept + arma::dot(x_focus, state.beta.t());
    out.beta.col(l) = state.beta;
    out.intercept[l] = state.intercept;
    out.linear_predictor[l] = eta_focus;
    out.fitted[l] = eta_focus;
    out.iterations[l] = fit.iterations;
    out.irls_iterations[l] = 0;
    out.converged[l] = fit.converged ? 1u : 0u;
    out.kkt[l] = fit.kkt_max;
    out.n_active[l] = active_count(state.beta);
  }
  return out;
}

inline PathResult fit_gaussian_covariance_path(
    const arma::mat& X,
    const arma::rowvec& x_focus,
    const arma::vec& y,
    const arma::vec& weights,
    const arma::vec& lambda,
    const arma::vec& factors,
    const SolverControl& control) {
  PathResult out = allocate_path_result(X.n_cols, lambda.n_elem);
  LocalState state(X.n_cols);
  initialize_gaussian_null(X, y, weights, factors, control, state);

  arma::vec x_mean(X.n_cols, arma::fill::zeros);
  double y_mean = 0.0;
  arma::mat weighted_X = X.each_col() % weights;
  arma::mat gram = X.t() * weighted_X;
  arma::vec cross_product = X.t() * (weights % y);
  if (control.intercept) {
    x_mean = X.t() * weights;
    y_mean = arma::dot(weights, y);
    gram -= x_mean * x_mean.t();
    cross_product -= x_mean * y_mean;
    state.intercept = y_mean - arma::dot(x_mean, state.beta);
  }
  state.gradient = cross_product - gram * state.beta;

  for (arma::uword l = 0; l < lambda.n_elem; ++l) {
    initialize_active_set(
      state,
      factors,
      lambda[l],
      l == 0 ? lambda[l] : lambda[l - 1],
      control,
      l == 0
    );
    const WlsResult fit = coordinate_descent_covariance(
      gram, cross_product, factors, lambda[l], control, state
    );
    if (control.intercept) {
      state.intercept = y_mean - arma::dot(x_mean, state.beta);
    }
    const double eta_focus =
      state.intercept + arma::dot(x_focus, state.beta.t());
    out.beta.col(l) = state.beta;
    out.intercept[l] = state.intercept;
    out.linear_predictor[l] = eta_focus;
    out.fitted[l] = eta_focus;
    out.iterations[l] = fit.iterations;
    out.irls_iterations[l] = 0;
    out.converged[l] = fit.converged ? 1u : 0u;
    out.kkt[l] = fit.kkt_max;
    out.n_active[l] = active_count(state.beta);
  }
  return out;
}

inline PathResult fit_binomial_path(const arma::mat& X,
                                    const arma::rowvec& x_focus,
                                    const arma::vec& y,
                                    const arma::vec& weights,
                                    const arma::vec& lambda,
                                    const arma::vec& factors,
                                    const SolverControl& control) {
  PathResult out = allocate_path_result(X.n_cols, lambda.n_elem);
  LocalState state(X.n_cols);
  initialize_binomial_null(X, y, weights, factors, control, state);
  arma::vec eta = linear_predictor(X, state.intercept, state.beta);
  arma::vec probability(y.n_elem, arma::fill::none);
  for (arma::uword i = 0; i < y.n_elem; ++i) {
    probability[i] = clamp_probability(stable_sigmoid(eta[i]));
  }
  state.gradient = score_vector(X, y - probability, weights);

  for (arma::uword l = 0; l < lambda.n_elem; ++l) {
    initialize_active_set(
      state,
      factors,
      lambda[l],
      l == 0 ? lambda[l] : lambda[l - 1],
      control,
      l == 0
    );

    bool outer_converged = false;
    int total_inner = 0;
    int outer_used = 0;
    for (int outer = 0; outer < control.irls_max_iter; ++outer) {
      ++outer_used;
      const double old_intercept = state.intercept;
      const arma::vec old_beta = state.beta;
      eta = linear_predictor(X, state.intercept, state.beta);
      arma::vec h(y.n_elem, arma::fill::none);
      arma::vec working(y.n_elem, arma::fill::none);
      for (arma::uword i = 0; i < y.n_elem; ++i) {
        probability[i] = clamp_probability(stable_sigmoid(eta[i]));
        const double variance = std::max(
          probability[i] * (1.0 - probability[i]), 1e-8
        );
        h[i] = weights[i] * variance;
        working[i] = eta[i] + (y[i] - probability[i]) / variance;
      }

      const double old_objective = model_objective(
        X, y, weights, state, lambda[l], factors, control
      );
      const WlsResult inner = coordinate_descent_wls(
        X, working, h, factors, lambda[l], control, state
      );
      total_inner += inner.iterations;

      double new_objective = model_objective(
        X, y, weights, state, lambda[l], factors, control
      );
      if (!R_finite(new_objective) ||
          new_objective > old_objective + 1e-12) {
        const double candidate_intercept = state.intercept;
        const arma::vec candidate_beta = state.beta;
        double step = 1.0;
        do {
          step *= 0.5;
          state.intercept = old_intercept +
            step * (candidate_intercept - old_intercept);
          state.beta = old_beta + step * (candidate_beta - old_beta);
          new_objective = model_objective(
            X, y, weights, state, lambda[l], factors, control
          );
        } while ((!R_finite(new_objective) ||
                  new_objective > old_objective + 1e-12) &&
                 step > 1e-8);
      }

      const double intercept_change =
        std::abs(state.intercept - old_intercept);
      const double beta_change = state.beta.n_elem ?
        arma::abs(state.beta - old_beta).max() : 0.0;
      const double scale = 1.0 + std::max(
        std::abs(state.intercept),
        state.beta.n_elem ? arma::abs(state.beta).max() : 0.0
      );
      if (std::max(intercept_change, beta_change) <=
          control.irls_tol * scale) {
        outer_converged = true;
        break;
      }
    }

    eta = linear_predictor(X, state.intercept, state.beta);
    for (arma::uword i = 0; i < y.n_elem; ++i) {
      probability[i] = clamp_probability(stable_sigmoid(eta[i]));
    }
    const double kkt = full_kkt_scan(
      state, y - probability, weights, X, factors,
      lambda[l], control, false
    );
    const double gradient_scale = state.gradient.n_elem ?
      arma::abs(state.gradient).max() : 0.0;
    const bool converged = outer_converged &&
      kkt <= control.kkt_tol * (1.0 + gradient_scale);
    const double eta_focus =
      state.intercept + arma::dot(x_focus, state.beta.t());
    out.beta.col(l) = state.beta;
    out.intercept[l] = state.intercept;
    out.linear_predictor[l] = eta_focus;
    out.fitted[l] = stable_sigmoid(eta_focus);
    out.iterations[l] = total_inner;
    out.irls_iterations[l] = outer_used;
    out.converged[l] = converged ? 1u : 0u;
    out.kkt[l] = kkt;
    out.n_active[l] = active_count(state.beta);
  }
  return out;
}

inline PathResult fit_path_one(const arma::mat& X,
                               const arma::rowvec& x_focus,
                               const arma::vec& y,
                               const arma::vec& weights,
                               const arma::vec& lambda,
                               const arma::vec& factors,
                               const SolverControl& control) {
  if (control.family == Family::Gaussian) {
    const bool use_covariance =
      control.solver == "covariance" ||
      (control.solver == "auto" &&
       X.n_cols <= std::max<arma::uword>(20u, X.n_rows / 4u));
    if (use_covariance) {
      return fit_gaussian_covariance_path(
        X, x_focus, y, weights, lambda, factors, control
      );
    }
    return fit_gaussian_path(
      X, x_focus, y, weights, lambda, factors, control
    );
  }
  return fit_binomial_path(
    X, x_focus, y, weights, lambda, factors, control
  );
}

} // namespace lhwl

#endif

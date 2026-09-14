#ifndef LHWL_LHWPR_STATE_H
#define LHWL_LHWPR_STATE_H

#include "lhwpr_penalty.h"

namespace lhwl {

inline arma::vec score_vector(const arma::mat& X,
                              const arma::vec& residual,
                              const arma::vec& weights) {
  return X.t() * (weights % residual);
}

inline void initialize_active_set(LocalState& state,
                                  const arma::vec& factors,
                                  double lambda,
                                  double previous_lambda,
                                  const SolverControl& control,
                                  bool first_lambda) {
  const arma::uword p = state.beta.n_elem;
  if (control.penalty == Penalty::Ridge ||
      control.penalty == Penalty::MCP ||
      control.penalty == Penalty::SCAD ||
      !control.strong_screening ||
      first_lambda) {
    state.active.ones();
    if (first_lambda && control.strong_screening &&
        control.penalty == Penalty::Lasso) {
      state.active.zeros();
      for (arma::uword j = 0; j < p; ++j) {
        if (factors[j] == 0.0 ||
            state.ever_active[j] ||
            std::abs(state.gradient[j]) >=
              control.alpha * factors[j] * lambda) {
          state.active[j] = 1;
        }
      }
    }
    return;
  }

  const double strong_lambda = std::max(2.0 * lambda - previous_lambda, 0.0);
  state.active.zeros();
  for (arma::uword j = 0; j < p; ++j) {
    if (factors[j] == 0.0 ||
        state.ever_active[j] ||
        state.beta[j] != 0.0 ||
        std::abs(state.gradient[j]) >=
          control.alpha * factors[j] * strong_lambda) {
      state.active[j] = 1;
    }
  }
}

inline double kkt_violation_one(double gradient,
                                double beta,
                                double lambda,
                                double factor,
                                const SolverControl& control) {
  const double abs_beta = std::abs(beta);
  if (factor == 0.0) return std::abs(gradient);
  if (abs_beta <= 1e-12) {
    if (control.penalty == Penalty::Ridge) return std::abs(gradient);
    const double threshold = lambda * factor *
      (control.penalty == Penalty::Lasso ? control.alpha : 1.0);
    return std::max(std::abs(gradient) - threshold, 0.0);
  }
  const double sign = beta < 0.0 ? -1.0 : 1.0;
  const double derivative = penalty_derivative_abs(
    abs_beta, lambda, factor, control
  );
  return std::abs(gradient - sign * derivative);
}

inline double full_kkt_scan(LocalState& state,
                            const arma::vec& residual,
                            const arma::vec& weights,
                            const arma::mat& X,
                            const arma::vec& factors,
                            double lambda,
                            const SolverControl& control,
                            bool activate_violations) {
  state.gradient = score_vector(X, residual, weights);
  double max_violation = 0.0;
  for (arma::uword j = 0; j < state.beta.n_elem; ++j) {
    const double violation = kkt_violation_one(
      state.gradient[j], state.beta[j], lambda, factors[j], control
    );
    if (violation > max_violation) max_violation = violation;
    if (activate_violations && violation > control.kkt_tol) {
      state.active[j] = 1;
    }
    if (state.beta[j] != 0.0) state.ever_active[j] = 1;
  }
  return max_violation;
}

inline double full_kkt_gradient_scan(LocalState& state,
                                     const arma::vec& gradient,
                                     const arma::vec& factors,
                                     double lambda,
                                     const SolverControl& control,
                                     bool activate_violations) {
  state.gradient = gradient;
  double max_violation = 0.0;
  for (arma::uword j = 0; j < state.beta.n_elem; ++j) {
    const double violation = kkt_violation_one(
      gradient[j], state.beta[j], lambda, factors[j], control
    );
    if (violation > max_violation) max_violation = violation;
    if (activate_violations && violation > control.kkt_tol) {
      state.active[j] = 1;
    }
    if (state.beta[j] != 0.0) state.ever_active[j] = 1;
  }
  return max_violation;
}

inline int active_count(const arma::vec& beta) {
  return static_cast<int>(arma::accu(arma::abs(beta) > 1e-12));
}

} // namespace lhwl

#endif

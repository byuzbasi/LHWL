#ifndef LHWL_LHWPR_PENALTY_H
#define LHWL_LHWPR_PENALTY_H

#include "lhwpr_types.h"

namespace lhwl {

inline double mcp_value(double t, double lambda, double gamma) {
  if (t <= gamma * lambda) {
    return lambda * t - t * t / (2.0 * gamma);
  }
  return 0.5 * gamma * lambda * lambda;
}

inline double scad_value(double t, double lambda, double gamma) {
  if (t <= lambda) return lambda * t;
  if (t <= gamma * lambda) {
    return (-t * t + 2.0 * gamma * lambda * t - lambda * lambda) /
      (2.0 * (gamma - 1.0));
  }
  return 0.5 * (gamma + 1.0) * lambda * lambda;
}

inline double penalty_value_one(double beta,
                                double lambda,
                                double factor,
                                const SolverControl& control) {
  const double t = std::abs(beta);
  if (control.penalty == Penalty::Ridge) {
    return 0.5 * lambda * factor * beta * beta;
  }
  if (control.penalty == Penalty::Lasso) {
    return lambda * factor *
      (control.alpha * t + 0.5 * (1.0 - control.alpha) * beta * beta);
  }
  if (control.penalty == Penalty::MCP) {
    return factor * mcp_value(t, lambda, control.gamma);
  }
  return factor * scad_value(t, lambda, control.gamma);
}

inline double penalty_value(const arma::vec& beta,
                            double lambda,
                            const arma::vec& factors,
                            const SolverControl& control) {
  double out = 0.0;
  for (arma::uword j = 0; j < beta.n_elem; ++j) {
    out += penalty_value_one(beta[j], lambda, factors[j], control);
  }
  return out;
}

inline void add_candidate(std::vector<double>& candidates,
                          double value,
                          double lower,
                          double upper) {
  if (!R_finite(value)) return;
  const double eps = 1e-12 * (
    1.0 + std::abs(lower) + (R_finite(upper) ? std::abs(upper) : 0.0)
  );
  if (value < lower - eps || value > upper + eps) return;
  if (value < lower) value = lower;
  if (value > upper) value = upper;
  candidates.push_back(value);
}

inline double nonconvex_coordinate(double z,
                                   double curvature,
                                   double lambda,
                                   double factor,
                                   double current,
                                   const SolverControl& control) {
  if (curvature <= 0.0 || factor == 0.0) {
    return curvature > 0.0 ? z / curvature : 0.0;
  }

  const double sign = z < 0.0 ? -1.0 : 1.0;
  const double u = std::abs(z);
  const double zero_threshold = factor * lambda;
  if (std::abs(current) <= 1e-14 &&
      u <= zero_threshold + 1e-12 * (1.0 + zero_threshold)) {
    return 0.0;
  }
  const double boundary = control.gamma * lambda;
  std::vector<double> candidates;
  candidates.reserve(10);
  candidates.push_back(0.0);

  if (control.penalty == Penalty::MCP) {
    candidates.push_back(boundary);
    const double denom = curvature - factor / control.gamma;
    if (std::abs(denom) > 1e-14) {
      add_candidate(
        candidates,
        (u - factor * lambda) / denom,
        0.0,
        boundary
      );
    }
    add_candidate(candidates, u / curvature, boundary,
                  std::numeric_limits<double>::infinity());
  } else {
    candidates.push_back(lambda);
    candidates.push_back(boundary);
    add_candidate(
      candidates,
      (u - factor * lambda) / curvature,
      0.0,
      lambda
    );
    const double denom = curvature - factor / (control.gamma - 1.0);
    if (std::abs(denom) > 1e-14) {
      add_candidate(
        candidates,
        (u - factor * control.gamma * lambda /
          (control.gamma - 1.0)) / denom,
        lambda,
        boundary
      );
    }
    add_candidate(candidates, u / curvature, boundary,
                  std::numeric_limits<double>::infinity());
  }

  double best = 0.0;
  double best_objective = penalty_value_one(
    0.0, lambda, factor, control
  );
  for (std::size_t k = 0; k < candidates.size(); ++k) {
    const double t = candidates[k];
    if (!R_finite(t) || t < 0.0) continue;
    const double objective =
      0.5 * curvature * t * t - u * t +
      penalty_value_one(t, lambda, factor, control);
    if (objective < best_objective) {
      best_objective = objective;
      best = t;
    }
  }
  return sign * best;
}

inline double coordinate_update(double z,
                                double curvature,
                                double lambda,
                                double factor,
                                double current,
                                const SolverControl& control) {
  if (curvature <= 0.0) return 0.0;
  if (factor == 0.0) return z / curvature;
  if (control.penalty == Penalty::Ridge) {
    return z / (curvature + lambda * factor);
  }
  if (control.penalty == Penalty::Lasso) {
    return soft_threshold(z, lambda * control.alpha * factor) /
      (curvature + lambda * (1.0 - control.alpha) * factor);
  }
  return nonconvex_coordinate(
    z, curvature, lambda, factor, current, control
  );
}

inline double penalty_derivative_abs(double t,
                                     double lambda,
                                     double factor,
                                     const SolverControl& control) {
  if (factor == 0.0) return 0.0;
  if (control.penalty == Penalty::Ridge) {
    return lambda * factor * t;
  }
  if (control.penalty == Penalty::Lasso) {
    return lambda * factor *
      (control.alpha + (1.0 - control.alpha) * t);
  }
  if (control.penalty == Penalty::MCP) {
    return factor * std::max(lambda - t / control.gamma, 0.0);
  }
  if (t <= lambda) return factor * lambda;
  if (t <= control.gamma * lambda) {
    return factor * (control.gamma * lambda - t) /
      (control.gamma - 1.0);
  }
  return 0.0;
}

} // namespace lhwl

#endif

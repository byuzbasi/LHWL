#ifndef LHWL_LHWPR_TYPES_H
#define LHWL_LHWPR_TYPES_H

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace lhwl {

enum class Family {
  Gaussian,
  Binomial
};

enum class Penalty {
  Ridge,
  Lasso,
  MCP,
  SCAD
};

struct SolverControl {
  Family family;
  Penalty penalty;
  double alpha;
  double gamma;
  bool intercept;
  std::string solver;
  bool strong_screening;
  double tol;
  double kkt_tol;
  int max_iter;
  double irls_tol;
  int irls_max_iter;
  int max_active;
  bool early_stop;
  double ridge_path_tol;
  int n_threads;
  bool trace;
  bool store_path;
  bool return_inference;
  bool return_path_diagnostics;
  bool return_active_frequency;
  bool compute_lambda_diagnostics;
  int selected_index;
  bool exclude_self;
  int nlambda;
  double lambda_min_ratio;
  std::string kernel;
};

struct LocalState {
  double intercept;
  arma::vec beta;
  arma::uvec active;
  arma::uvec ever_active;
  arma::vec gradient;

  explicit LocalState(arma::uword p)
    : intercept(0.0),
      beta(p, arma::fill::zeros),
      active(p, arma::fill::ones),
      ever_active(p, arma::fill::zeros),
      gradient(p, arma::fill::zeros) {}
};

struct WlsResult {
  int iterations;
  bool converged;
  double kkt_max;

  WlsResult() : iterations(0), converged(false), kkt_max(R_PosInf) {}
};

struct PathResult {
  arma::mat beta;
  arma::vec intercept;
  arma::vec fitted;
  arma::vec linear_predictor;
  arma::ivec iterations;
  arma::ivec irls_iterations;
  arma::uvec converged;
  arma::vec kkt;
  arma::ivec n_active;
};

inline Family parse_family(const std::string& family) {
  if (family == "gaussian") return Family::Gaussian;
  if (family == "binomial") return Family::Binomial;
  Rcpp::stop("Unsupported LHWPR family.");
}

inline Penalty parse_penalty(const std::string& penalty, double alpha) {
  if (alpha == 0.0) return Penalty::Ridge;
  if (penalty == "lasso") return Penalty::Lasso;
  if (penalty == "mcp") return Penalty::MCP;
  if (penalty == "scad") return Penalty::SCAD;
  Rcpp::stop("Unsupported LHWPR penalty.");
}

inline int resolve_threads(int requested) {
  if (requested < 1) requested = 1;
#ifdef _OPENMP
  return requested;
#else
  return 1;
#endif
}

inline double stable_sigmoid(double eta) {
  if (eta >= 0.0) {
    const double e = std::exp(-eta);
    return 1.0 / (1.0 + e);
  }
  const double e = std::exp(eta);
  return e / (1.0 + e);
}

inline double stable_log1pexp(double eta) {
  if (eta > 0.0) return eta + std::log1p(std::exp(-eta));
  return std::log1p(std::exp(eta));
}

inline double clamp_probability(double p) {
  const double eps = 1e-8;
  if (p < eps) return eps;
  if (p > 1.0 - eps) return 1.0 - eps;
  return p;
}

inline double soft_threshold(double z, double threshold) {
  if (z > threshold) return z - threshold;
  if (z < -threshold) return z + threshold;
  return 0.0;
}

inline double weighted_sum(const arma::vec& x, const arma::vec& w) {
  return arma::dot(x, w);
}

inline arma::vec linear_predictor(const arma::mat& X,
                                  double intercept,
                                  const arma::vec& beta) {
  arma::vec eta = X * beta;
  eta += intercept;
  return eta;
}

} // namespace lhwl

#endif

#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "lhwpr_solver.h"

using namespace Rcpp;

namespace {

constexpr double kSurfaceActiveTolerance = 1e-10;

double kernel_weight(double distance,
                     double bandwidth,
                     const std::string& kernel) {
  const double ratio = distance / bandwidth;
  if (kernel == "gaussian") return std::exp(-0.5 * ratio * ratio);
  if (kernel == "exponential") return std::exp(-ratio);
  if (kernel == "tricube") {
    if (ratio >= 1.0) return 0.0;
    const double value = 1.0 - ratio * ratio * ratio;
    return value * value * value;
  }
  if (kernel == "boxcar") return ratio < 1.0 ? 1.0 : 0.0;
  if (ratio >= 1.0) return 0.0;
  const double value = 1.0 - ratio * ratio;
  return value * value;
}

lhwl::SolverControl parse_control(const List& options) {
  lhwl::SolverControl out;
  out.family = lhwl::parse_family(as<std::string>(options["family"]));
  out.alpha = as<double>(options["alpha"]);
  out.penalty = lhwl::parse_penalty(
    as<std::string>(options["penalty"]), out.alpha
  );
  out.gamma = as<double>(options["gamma"]);
  out.intercept = as<bool>(options["intercept"]);
  out.solver = as<std::string>(options["solver"]);
  out.strong_screening =
    as<std::string>(options["screening"]) == "strong";
  out.tol = as<double>(options["tol"]);
  out.kkt_tol = as<double>(options["kkt_tol"]);
  out.max_iter = as<int>(options["max_iter"]);
  out.irls_tol = as<double>(options["irls_tol"]);
  out.irls_max_iter = as<int>(options["irls_max_iter"]);
  out.max_active = as<int>(options["max_active"]);
  out.early_stop = as<bool>(options["early_stop"]);
  out.ridge_path_tol = as<double>(options["ridge_path_tol"]);
  out.n_threads = as<int>(options["n_threads"]);
  out.trace = as<bool>(options["trace"]);
  out.store_path = as<bool>(options["store_path"]);
  out.return_inference =
    options.containsElementNamed("return_inference") ?
    as<bool>(options["return_inference"]) : false;
  out.return_path_diagnostics =
    options.containsElementNamed("return_path_diagnostics") ?
    as<bool>(options["return_path_diagnostics"]) : false;
  out.return_active_frequency =
    options.containsElementNamed("return_active_frequency") ?
    as<bool>(options["return_active_frequency"]) : false;
  out.compute_lambda_diagnostics =
    options.containsElementNamed("compute_lambda_diagnostics") ?
    as<bool>(options["compute_lambda_diagnostics"]) : true;
  out.selected_index = as<int>(options["selected_index"]);
  out.exclude_self = as<bool>(options["exclude_self"]);
  out.nlambda = as<int>(options["nlambda"]);
  out.lambda_min_ratio = as<double>(options["lambda_min_ratio"]);
  out.kernel = as<std::string>(options["kernel"]);
  return out;
}

class DenseWeightProvider {
public:
  DenseWeightProvider(const arma::mat& distance,
                      const arma::vec& bandwidth,
                      const arma::vec& case_weights,
                      const lhwl::SolverControl& control)
    : distance_(distance),
      bandwidth_(bandwidth),
      case_weights_(case_weights),
      control_(control) {}

  arma::uword size() const {
    return distance_.n_cols;
  }

  arma::vec weights(arma::uword focus) const {
    arma::vec out(distance_.n_rows, arma::fill::none);
    const double bandwidth = bandwidth_[focus];
    for (arma::uword i = 0; i < distance_.n_rows; ++i) {
      out[i] = kernel_weight(
        distance_(i, focus), bandwidth, control_.kernel
      ) * case_weights_[i];
    }
    if (control_.exclude_self &&
        distance_.n_cols == distance_.n_rows &&
        focus < out.n_elem) {
      out[focus] = 0.0;
    }
    return out;
  }

private:
  const arma::mat& distance_;
  const arma::vec& bandwidth_;
  const arma::vec& case_weights_;
  const lhwl::SolverControl& control_;
};

class KnnWeightProvider {
public:
  KnnWeightProvider(const arma::imat& index,
                    const arma::mat& distance,
                    const arma::vec& case_weights,
                    const lhwl::SolverControl& control)
    : index_(index),
      distance_(distance),
      case_weights_(case_weights),
      control_(control) {}

  arma::uword size() const {
    return index_.n_rows;
  }

  arma::vec weights(arma::uword focus) const {
    arma::vec out(case_weights_.n_elem, arma::fill::zeros);
    const double bandwidth = distance_(
      focus, distance_.n_cols - 1
    );
    if (!control_.exclude_self && focus < out.n_elem) {
      out[focus] = case_weights_[focus];
    }
    for (arma::uword k = 0; k < index_.n_cols; ++k) {
      const int raw_index = index_(focus, k);
      if (raw_index < 1 ||
          raw_index > static_cast<int>(case_weights_.n_elem)) {
        continue;
      }
      const arma::uword observation =
        static_cast<arma::uword>(raw_index - 1);
      out[observation] = kernel_weight(
        distance_(focus, k), bandwidth, control_.kernel
      ) * case_weights_[observation];
    }
    return out;
  }

private:
  const arma::imat& index_;
  const arma::mat& distance_;
  const arma::vec& case_weights_;
  const lhwl::SolverControl& control_;
};

bool normalize_weights(arma::vec& weights,
                       double& effective_n,
                       int& positive_n) {
  const double total = arma::accu(weights);
  if (!R_finite(total) || total <= 0.0 || !weights.is_finite()) {
    effective_n = 0.0;
    positive_n = 0;
    return false;
  }
  weights /= total;
  const double squared = arma::accu(arma::square(weights));
  effective_n = squared > 0.0 ? 1.0 / squared : 0.0;
  positive_n = static_cast<int>(arma::accu(weights > 0.0));
  return true;
}

arma::mat inverse_symmetric(const arma::mat& matrix) {
  arma::mat out;
  bool ok = arma::inv_sympd(out, matrix);
  if (ok && out.is_finite()) return out;
  ok = arma::inv(out, matrix);
  if (ok && out.is_finite()) return out;
  return arma::pinv(matrix);
}

double matrix_condition_number(const arma::mat& matrix) {
  arma::vec singular_values;
  if (!arma::svd(singular_values, matrix) ||
      singular_values.n_elem == 0) {
    return R_PosInf;
  }
  const double largest = singular_values.max();
  const double smallest = singular_values.min();
  if (!R_finite(largest) || !R_finite(smallest) || smallest <= 0.0) {
    return R_PosInf;
  }
  return largest / smallest;
}

struct GaussianPathDiagnostics {
  arma::vec hat_diag;
  arma::vec post_fitted;
  arma::vec post_hat_diag;
};

arma::uvec conditional_active_set(const arma::vec& beta,
                                  const arma::vec& factors) {
  std::vector<arma::uword> active;
  active.reserve(beta.n_elem);
  for (arma::uword j = 0; j < beta.n_elem; ++j) {
    if (factors[j] == 0.0 ||
        std::abs(beta[j]) > kSurfaceActiveTolerance) {
      active.push_back(j);
    }
  }
  arma::uvec out(active.size());
  for (std::size_t j = 0; j < active.size(); ++j) {
    out[j] = active[j];
  }
  return out;
}

double local_penalty_curvature(double beta,
                               double lambda,
                               double factor,
                               const lhwl::SolverControl& control) {
  if (factor == 0.0) return 0.0;
  if (control.penalty == lhwl::Penalty::Ridge) {
    return lambda * factor;
  }
  if (control.penalty == lhwl::Penalty::Lasso) {
    return lambda * (1.0 - control.alpha) * factor;
  }
  const double magnitude = std::abs(beta);
  if (control.penalty == lhwl::Penalty::MCP) {
    return magnitude < control.gamma * lambda ?
      -factor / control.gamma : 0.0;
  }
  return magnitude > lambda && magnitude < control.gamma * lambda ?
    -factor / (control.gamma - 1.0) : 0.0;
}

GaussianPathDiagnostics gaussian_path_diagnostics(
    const arma::mat& X,
    const arma::rowvec& x_focus,
    const arma::vec& y,
    const arma::vec& weights,
    double self_weight,
    const arma::vec& lambda,
    const arma::vec& factors,
    const lhwl::PathResult& path,
    const lhwl::SolverControl& control) {
  const arma::uword nlambda = lambda.n_elem;
  GaussianPathDiagnostics out;
  out.hat_diag = arma::vec(
    nlambda, arma::fill::value(NA_REAL)
  );
  out.post_fitted = arma::vec(
    nlambda, arma::fill::value(NA_REAL)
  );
  out.post_hat_diag = arma::vec(
    nlambda, arma::fill::value(NA_REAL)
  );

  for (arma::uword l = 0; l < nlambda; ++l) {
    const arma::uvec active = conditional_active_set(
      path.beta.col(l), factors
    );
    const arma::uword q =
      active.n_elem + (control.intercept ? 1u : 0u);
    if (q == 0) {
      out.hat_diag[l] = 0.0;
      out.post_fitted[l] = 0.0;
      out.post_hat_diag[l] = 0.0;
      continue;
    }

    arma::mat design(X.n_rows, q, arma::fill::none);
    arma::vec focus(q, arma::fill::none);
    arma::uword column = 0;
    if (control.intercept) {
      design.col(column).ones();
      focus[column++] = 1.0;
    }
    for (arma::uword k = 0; k < active.n_elem; ++k) {
      design.col(column) = X.col(active[k]);
      focus[column++] = x_focus[active[k]];
    }

    const arma::mat weighted_design = design.each_col() % weights;
    const arma::mat unpenalized_hessian =
      design.t() * weighted_design;
    arma::mat penalized_hessian = unpenalized_hessian;
    column = control.intercept ? 1u : 0u;
    for (arma::uword k = 0; k < active.n_elem; ++k) {
      const arma::uword predictor = active[k];
      penalized_hessian(column, column) += local_penalty_curvature(
        path.beta(predictor, l),
        lambda[l],
        factors[predictor],
        control
      );
      ++column;
    }

    const arma::mat penalized_inverse =
      inverse_symmetric(penalized_hessian);
    const double leverage = self_weight *
      arma::as_scalar(focus.t() * penalized_inverse * focus);
    if (R_finite(leverage)) out.hat_diag[l] = leverage;

    const arma::mat post_inverse =
      inverse_symmetric(unpenalized_hessian);
    const arma::vec theta = post_inverse *
      (design.t() * (weights % y));
    const double post_value = arma::dot(focus, theta);
    const double post_leverage = self_weight *
      arma::as_scalar(focus.t() * post_inverse * focus);
    if (R_finite(post_value)) out.post_fitted[l] = post_value;
    if (R_finite(post_leverage)) {
      out.post_hat_diag[l] = post_leverage;
    }
  }
  return out;
}

template <typename Provider>
double compute_lambda_max(const Provider& provider,
                          const arma::mat& X,
                          const arma::vec& y,
                          const arma::vec& factors,
                          const lhwl::SolverControl& control,
                          arma::vec& local_max,
                          arma::uvec& valid) {
  const int threads = lhwl::resolve_threads(control.n_threads);
#ifndef _OPENMP
  (void) threads;
#endif
  local_max.set_size(provider.size());
  valid.ones(provider.size());

#pragma omp parallel for schedule(static) num_threads(threads) if(threads > 1)
  for (arma::sword ii = 0;
       ii < static_cast<arma::sword>(provider.size()); ++ii) {
    const arma::uword i = static_cast<arma::uword>(ii);
    arma::vec weights = provider.weights(i);
    double effective_n = 0.0;
    int positive_n = 0;
    if (!normalize_weights(weights, effective_n, positive_n)) {
      valid[i] = 0;
      local_max[i] = NA_REAL;
      continue;
    }
    const arma::uvec support = arma::find(weights > 0.0);
    const arma::mat X_local = X.rows(support);
    const arma::vec y_local = y.elem(support);
    const arma::vec weights_local = weights.elem(support);
    local_max[i] = lhwl::lambda_max_one(
      X_local, y_local, weights_local, factors, control
    );
  }
  if (arma::any(valid == 0)) return NA_REAL;
  return local_max.max();
}

arma::vec make_lambda_path(double lambda_max,
                           const lhwl::SolverControl& control) {
  arma::vec lambda(control.nlambda, arma::fill::none);
  if (control.nlambda == 1) {
    lambda[0] = lambda_max;
    return lambda;
  }
  const double lower = lambda_max * control.lambda_min_ratio;
  lambda = arma::exp(arma::linspace<arma::vec>(
    std::log(lambda_max), std::log(lower), control.nlambda
  ));
  return lambda;
}

void validate_common_inputs(const arma::vec& y,
                            const arma::mat& X,
                            const arma::mat& X_focus,
                            const arma::vec& factors,
                            const arma::vec& case_weights,
                            const lhwl::SolverControl& control) {
  if (X.n_rows != y.n_elem ||
      X_focus.n_cols != X.n_cols ||
      factors.n_elem != X.n_cols ||
      case_weights.n_elem != y.n_elem) {
    stop("LHWPR native input dimensions are not compatible.");
  }
  if (!y.is_finite() || !X.is_finite() || !X_focus.is_finite() ||
      !factors.is_finite() || !case_weights.is_finite()) {
    stop("LHWPR native inputs must be finite.");
  }
  if (arma::any(factors < 0.0) || arma::any(case_weights < 0.0)) {
    stop("Penalty factors and case weights must be non-negative.");
  }
  if (arma::accu(case_weights) <= 0.0) {
    stop("At least one case weight must be positive.");
  }
  if (X.n_cols == 0 && !control.intercept) {
    stop("At least one intercept or predictor coefficient is required.");
  }
  if (control.alpha < 0.0 || control.alpha > 1.0) {
    stop("`alpha` must lie in [0, 1].");
  }
  if (control.penalty == lhwl::Penalty::MCP && control.gamma <= 1.0) {
    stop("MCP requires `gamma > 1`.");
  }
  if (control.penalty == lhwl::Penalty::SCAD && control.gamma <= 2.0) {
    stop("SCAD requires `gamma > 2`.");
  }
}

template <typename Provider>
List fit_batch(const Provider& provider,
               const arma::vec& y,
               const arma::mat& X,
               const arma::mat& X_focus,
               arma::vec lambda,
               const arma::vec& factors,
               const lhwl::SolverControl& control) {
  arma::vec local_lambda_max;
  arma::uvec lambda_valid;
  const bool unpenalized_single_fit =
    lambda.n_elem == 1 &&
    lambda[0] == 0.0 &&
    (factors.n_elem == 0 || arma::all(factors == 0.0));
  double common_lambda_max = 1.0;
  if (unpenalized_single_fit) {
    local_lambda_max = arma::vec(provider.size(), arma::fill::ones);
    lambda_valid = arma::uvec(provider.size(), arma::fill::ones);
  } else if (lambda.n_elem == 0 || control.compute_lambda_diagnostics) {
    common_lambda_max = compute_lambda_max(
      provider, X, y, factors, control, local_lambda_max, lambda_valid
    );
  } else {
    common_lambda_max = NA_REAL;
    local_lambda_max = arma::vec(
      provider.size(), arma::fill::value(NA_REAL)
    );
    lambda_valid = arma::uvec(provider.size(), arma::fill::ones);
  }
  if ((lambda.n_elem == 0 || control.compute_lambda_diagnostics) &&
      !R_finite(common_lambda_max)) {
    stop("At least one focal neighborhood has zero total analysis weight.");
  }
  if (lambda.n_elem == 0) {
    lambda = make_lambda_path(common_lambda_max, control);
  }
  if (!lambda.is_finite() || arma::any(lambda < 0.0)) {
    stop("`lambda` must contain finite non-negative values.");
  }

  const arma::uword locations = provider.size();
  const arma::uword p = X.n_cols;
  const arma::uword q = p + (control.intercept ? 1u : 0u);
  const arma::uword nlambda = lambda.n_elem;
  if (X_focus.n_rows != locations) {
    stop("The focal design matrix does not match the number of locations.");
  }
  int selected = control.selected_index;
  if (selected < 0) selected = static_cast<int>(nlambda) - 1;
  if (selected >= static_cast<int>(nlambda)) {
    stop("Selected lambda index is outside the fitted path.");
  }

  arma::mat beta_hat(locations, q, arma::fill::value(NA_REAL));
  arma::cube beta_path;
  if (control.store_path) {
    beta_path = arma::cube(
      locations, q, nlambda, arma::fill::value(NA_REAL)
    );
  }
  arma::mat fitted_path(locations, nlambda, arma::fill::value(NA_REAL));
  arma::mat eta_path(locations, nlambda, arma::fill::value(NA_REAL));
  arma::mat hat_diag_path;
  arma::mat post_fitted_path;
  arma::mat post_hat_diag_path;
  if (control.return_path_diagnostics) {
    hat_diag_path = arma::mat(
      locations, nlambda, arma::fill::value(NA_REAL)
    );
    post_fitted_path = arma::mat(
      locations, nlambda, arma::fill::value(NA_REAL)
    );
    post_hat_diag_path = arma::mat(
      locations, nlambda, arma::fill::value(NA_REAL)
    );
  }
  arma::imat iterations(locations, nlambda, arma::fill::zeros);
  arma::imat irls_iterations(locations, nlambda, arma::fill::zeros);
  arma::imat converged(locations, nlambda, arma::fill::zeros);
  arma::mat kkt(locations, nlambda, arma::fill::value(R_PosInf));
  arma::imat n_active(locations, nlambda, arma::fill::zeros);
  arma::vec effective_n(locations, arma::fill::zeros);
  arma::ivec positive_n(locations, arma::fill::zeros);
  arma::vec weight_sum(locations, arma::fill::zeros);
  arma::mat standard_error(
    locations, q, arma::fill::value(NA_REAL)
  );
  arma::vec local_kappa(locations, arma::fill::value(NA_REAL));
  arma::ivec local_rank(locations);
  local_rank.fill(NA_INTEGER);
  arma::ivec local_class_variation(locations, arma::fill::ones);
  arma::ivec local_events(locations, arma::fill::zeros);
  arma::ivec local_nonevents(locations, arma::fill::zeros);
  arma::uvec valid(locations, arma::fill::ones);

  const int threads = lhwl::resolve_threads(control.n_threads);
  std::vector<arma::imat> active_by_thread;
  if (control.return_active_frequency) {
    active_by_thread.resize(threads);
    for (int thread = 0; thread < threads; ++thread) {
      active_by_thread[thread] = arma::imat(
        p, nlambda, arma::fill::zeros
      );
    }
  }
#ifndef _OPENMP
  (void) threads;
#endif
#pragma omp parallel for schedule(static) num_threads(threads) if(threads > 1)
  for (arma::sword ii = 0;
       ii < static_cast<arma::sword>(locations); ++ii) {
    const arma::uword i = static_cast<arma::uword>(ii);
    arma::vec weights = provider.weights(i);
    weight_sum[i] = arma::accu(weights);
    int positive = 0;
    if (!normalize_weights(weights, effective_n[i], positive)) {
      valid[i] = 0;
      continue;
    }
    const double self_weight =
      i < weights.n_elem ? weights[i] : 0.0;
    positive_n[i] = positive;
    const arma::uvec support = arma::find(weights > 0.0);
    const arma::mat X_local = X.rows(support);
    const arma::vec y_local = y.elem(support);
    const arma::vec weights_local = weights.elem(support);
    if (control.family == lhwl::Family::Binomial) {
      local_events[i] = static_cast<int>(arma::accu(y_local == 1.0));
      local_nonevents[i] = static_cast<int>(arma::accu(y_local == 0.0));
    }
    const lhwl::PathResult path = lhwl::fit_path_one(
      X_local, X_focus.row(i), y_local, weights_local,
      lambda, factors, control
    );
    fitted_path.row(i) = path.fitted.t();
    eta_path.row(i) = path.linear_predictor.t();
    iterations.row(i) = path.iterations.t();
    irls_iterations.row(i) = path.irls_iterations.t();
    converged.row(i) = arma::conv_to<arma::irowvec>::from(
      path.converged.t()
    );
    if (control.family == lhwl::Family::Binomial &&
        y_local.min() == y_local.max()) {
      local_class_variation[i] = 0;
      converged.row(i).zeros();
    }
    kkt.row(i) = path.kkt.t();
    n_active.row(i) = path.n_active.t();

    if (control.return_path_diagnostics) {
      if (control.family != lhwl::Family::Gaussian) {
        valid[i] = 0;
        continue;
      }
      const GaussianPathDiagnostics diagnostics =
        gaussian_path_diagnostics(
          X_local,
          X_focus.row(i),
          y_local,
          weights_local,
          self_weight,
          lambda,
          factors,
          path,
          control
        );
      hat_diag_path.row(i) = diagnostics.hat_diag.t();
      post_fitted_path.row(i) = diagnostics.post_fitted.t();
      post_hat_diag_path.row(i) =
        diagnostics.post_hat_diag.t();

    }

    if (control.return_active_frequency) {
      int thread_id = 0;
#ifdef _OPENMP
      thread_id = omp_get_thread_num();
#endif
      for (arma::uword l = 0; l < nlambda; ++l) {
        for (arma::uword j = 0; j < p; ++j) {
          if (factors[j] > 0.0 &&
              std::abs(path.beta(j, l)) >
                kSurfaceActiveTolerance) {
            ++active_by_thread[thread_id](j, l);
          }
        }
      }
    }

    arma::uword column = 0;
    if (control.intercept) {
      beta_hat(i, column++) = path.intercept[selected];
    }
    if (p > 0) {
      beta_hat.row(i).cols(column, q - 1) =
        path.beta.col(selected).t();
    }

    if (control.return_inference &&
        control.family == lhwl::Family::Binomial) {
      arma::mat design(X_local.n_rows, q, arma::fill::none);
      column = 0;
      if (control.intercept) {
        design.col(column++).ones();
      }
      if (p > 0) {
        design.cols(column, q - 1) = X_local;
      }
      arma::vec eta_local = X_local * path.beta.col(selected);
      eta_local += path.intercept[selected];
      arma::vec information_weight(
        eta_local.n_elem, arma::fill::none
      );
      for (arma::uword observation = 0;
           observation < eta_local.n_elem; ++observation) {
        const double probability = lhwl::clamp_probability(
          lhwl::stable_sigmoid(eta_local[observation])
        );
        information_weight[observation] =
          weights_local[observation] * weight_sum[i] *
          probability * (1.0 - probability);
      }
      const arma::mat weighted_design =
        design.each_col() % information_weight;
      const arma::mat information = design.t() * weighted_design;
      local_kappa[i] = matrix_condition_number(information);
      const arma::uword information_rank = arma::rank(information);
      local_rank[i] = static_cast<int>(information_rank);
      if (local_class_variation[i] != 0 &&
          information_rank == q) {
        const arma::mat covariance = inverse_symmetric(information);
        const arma::vec variance = covariance.diag();
        for (arma::uword coefficient = 0;
             coefficient < q; ++coefficient) {
          standard_error(i, coefficient) =
            R_finite(variance[coefficient]) &&
            variance[coefficient] >= 0.0 ?
            std::sqrt(variance[coefficient]) : NA_REAL;
        }
      }
    }

    if (control.store_path) {
      for (arma::uword l = 0; l < nlambda; ++l) {
        column = 0;
        if (control.intercept) {
          beta_path(i, column++, l) = path.intercept[l];
        }
        if (p > 0) {
          beta_path.slice(l).row(i).cols(column, q - 1) =
            path.beta.col(l).t();
        }
      }
    }
  }
  if (arma::any(valid == 0)) {
    if (control.return_path_diagnostics &&
        control.family != lhwl::Family::Gaussian) {
      stop("Path diagnostics are currently available only for Gaussian LHWPR.");
    }
    stop("At least one focal neighborhood has zero total analysis weight.");
  }

  arma::imat surface_active_count;
  if (control.return_active_frequency) {
    surface_active_count = arma::imat(
      p, nlambda, arma::fill::zeros
    );
    for (int thread = 0; thread < threads; ++thread) {
      surface_active_count += active_by_thread[thread];
    }
  }

  arma::vec lambda_ratio(nlambda, arma::fill::value(NA_REAL));
  if (R_finite(common_lambda_max) && common_lambda_max > 0.0) {
    lambda_ratio = lambda / common_lambda_max;
  }
  return List::create(
    _["beta_hat_standardized"] = beta_hat,
    _["beta_path_standardized"] = control.store_path ?
      static_cast<SEXP>(wrap(beta_path)) : R_NilValue,
    _["fitted_path"] = fitted_path,
    _["linear_predictor_path"] = eta_path,
    _["hat_diag_path"] = control.return_path_diagnostics ?
      static_cast<SEXP>(wrap(hat_diag_path)) : R_NilValue,
    _["post_fitted_path"] = control.return_path_diagnostics ?
      static_cast<SEXP>(wrap(post_fitted_path)) : R_NilValue,
    _["post_hat_diag_path"] = control.return_path_diagnostics ?
      static_cast<SEXP>(wrap(post_hat_diag_path)) : R_NilValue,
    _["surface_active_count"] = control.return_active_frequency ?
      static_cast<SEXP>(wrap(surface_active_count)) : R_NilValue,
    _["lambda"] = lambda,
    _["lambda_ratio"] = lambda_ratio,
    _["lambda_max"] = common_lambda_max,
    _["lambda_local_max"] = local_lambda_max,
    _["selected_index"] = selected + 1,
    _["iterations"] = iterations,
    _["irls_iterations"] = irls_iterations,
    _["converged"] = converged,
    _["kkt"] = kkt,
    _["n_active"] = n_active,
    _["local_neff"] = effective_n,
    _["local_n_positive"] = positive_n,
    _["local_weight_sum"] = weight_sum,
    _["standard_error"] = standard_error,
    _["local_kappa"] = local_kappa,
    _["local_rank"] = local_rank,
    _["local_class_variation"] = local_class_variation,
    _["local_events"] = local_events,
    _["local_nonevents"] = local_nonevents
  );
}

template <typename Provider>
List lambda_max_batch(const Provider& provider,
                      const arma::vec& y,
                      const arma::mat& X,
                      const arma::vec& factors,
                      const lhwl::SolverControl& control) {
  arma::vec local_max;
  arma::uvec valid;
  const double common = compute_lambda_max(
    provider, X, y, factors, control, local_max, valid
  );
  if (!R_finite(common)) {
    stop("At least one focal neighborhood has zero total analysis weight.");
  }
  return List::create(
    _["lambda_max"] = common,
    _["lambda_local_max"] = local_max
  );
}

} // namespace

extern "C" SEXP _LHWL_lhwpr_fit_dense_cpp(SEXP ySEXP,
                                            SEXP XSEXP,
                                            SEXP XFocusSEXP,
                                            SEXP distanceSEXP,
                                            SEXP bandwidthSEXP,
                                            SEXP lambdaSEXP,
                                            SEXP factorsSEXP,
                                            SEXP caseWeightsSEXP,
                                            SEXP optionsSEXP) {
  BEGIN_RCPP
  const arma::vec y = as<arma::vec>(ySEXP);
  const arma::mat X = as<arma::mat>(XSEXP);
  const arma::mat X_focus = as<arma::mat>(XFocusSEXP);
  const arma::mat distance = as<arma::mat>(distanceSEXP);
  arma::vec bandwidth = as<arma::vec>(bandwidthSEXP);
  arma::vec lambda = as<arma::vec>(lambdaSEXP);
  const arma::vec factors = as<arma::vec>(factorsSEXP);
  const arma::vec case_weights = as<arma::vec>(caseWeightsSEXP);
  const List options(optionsSEXP);
  const lhwl::SolverControl control = parse_control(options);
  validate_common_inputs(
    y, X, X_focus, factors, case_weights, control
  );
  if (distance.n_rows != y.n_elem ||
      distance.n_cols != X_focus.n_rows ||
      !distance.is_finite() ||
      arma::accu(distance < 0.0) > 0) {
    stop("Dense LHWPR distance dimensions or values are invalid.");
  }
  if (bandwidth.n_elem == 1) {
    bandwidth = arma::vec(distance.n_cols).fill(bandwidth[0]);
  }
  if (bandwidth.n_elem != distance.n_cols ||
      !bandwidth.is_finite() ||
      arma::any(bandwidth <= 0.0)) {
    stop("Dense LHWPR bandwidth values are invalid.");
  }
  const DenseWeightProvider provider(
    distance, bandwidth, case_weights, control
  );
  return fit_batch(
    provider, y, X, X_focus, lambda, factors, control
  );
  END_RCPP
}

extern "C" SEXP _LHWL_lhwpr_fit_knn_cpp(SEXP ySEXP,
                                          SEXP XSEXP,
                                          SEXP XFocusSEXP,
                                          SEXP indexSEXP,
                                          SEXP distanceSEXP,
                                          SEXP lambdaSEXP,
                                          SEXP factorsSEXP,
                                          SEXP caseWeightsSEXP,
                                          SEXP optionsSEXP) {
  BEGIN_RCPP
  const arma::vec y = as<arma::vec>(ySEXP);
  const arma::mat X = as<arma::mat>(XSEXP);
  const arma::mat X_focus = as<arma::mat>(XFocusSEXP);
  const arma::imat index = as<arma::imat>(indexSEXP);
  const arma::mat distance = as<arma::mat>(distanceSEXP);
  arma::vec lambda = as<arma::vec>(lambdaSEXP);
  const arma::vec factors = as<arma::vec>(factorsSEXP);
  const arma::vec case_weights = as<arma::vec>(caseWeightsSEXP);
  const List options(optionsSEXP);
  const lhwl::SolverControl control = parse_control(options);
  validate_common_inputs(
    y, X, X_focus, factors, case_weights, control
  );
  if (index.n_rows != X_focus.n_rows ||
      distance.n_rows != index.n_rows ||
      distance.n_cols != index.n_cols ||
      index.n_cols == 0 ||
      !distance.is_finite() ||
      arma::accu(distance < 0.0) > 0) {
    stop("kNN LHWPR neighborhood dimensions or values are invalid.");
  }
  const KnnWeightProvider provider(
    index, distance, case_weights, control
  );
  return fit_batch(
    provider, y, X, X_focus, lambda, factors, control
  );
  END_RCPP
}

extern "C" SEXP _LHWL_lhwpr_lambda_max_dense_cpp(SEXP ySEXP,
                                                   SEXP XSEXP,
                                                   SEXP distanceSEXP,
                                                   SEXP bandwidthSEXP,
                                                   SEXP factorsSEXP,
                                                   SEXP caseWeightsSEXP,
                                                   SEXP optionsSEXP) {
  BEGIN_RCPP
  const arma::vec y = as<arma::vec>(ySEXP);
  const arma::mat X = as<arma::mat>(XSEXP);
  const arma::mat distance = as<arma::mat>(distanceSEXP);
  arma::vec bandwidth = as<arma::vec>(bandwidthSEXP);
  const arma::vec factors = as<arma::vec>(factorsSEXP);
  const arma::vec case_weights = as<arma::vec>(caseWeightsSEXP);
  const List options(optionsSEXP);
  const lhwl::SolverControl control = parse_control(options);
  validate_common_inputs(
    y, X, X, factors, case_weights, control
  );
  if (distance.n_rows != y.n_elem || distance.n_cols != y.n_elem) {
    stop("Dense lambda-maximum distances must be square.");
  }
  if (bandwidth.n_elem == 1) {
    bandwidth = arma::vec(distance.n_cols).fill(bandwidth[0]);
  }
  if (bandwidth.n_elem != distance.n_cols) {
    stop("Dense lambda-maximum bandwidth values are invalid.");
  }
  const DenseWeightProvider provider(
    distance, bandwidth, case_weights, control
  );
  return lambda_max_batch(provider, y, X, factors, control);
  END_RCPP
}

extern "C" SEXP _LHWL_lhwpr_lambda_max_knn_cpp(SEXP ySEXP,
                                                 SEXP XSEXP,
                                                 SEXP indexSEXP,
                                                 SEXP distanceSEXP,
                                                 SEXP factorsSEXP,
                                                 SEXP caseWeightsSEXP,
                                                 SEXP optionsSEXP) {
  BEGIN_RCPP
  const arma::vec y = as<arma::vec>(ySEXP);
  const arma::mat X = as<arma::mat>(XSEXP);
  const arma::imat index = as<arma::imat>(indexSEXP);
  const arma::mat distance = as<arma::mat>(distanceSEXP);
  const arma::vec factors = as<arma::vec>(factorsSEXP);
  const arma::vec case_weights = as<arma::vec>(caseWeightsSEXP);
  const List options(optionsSEXP);
  const lhwl::SolverControl control = parse_control(options);
  validate_common_inputs(
    y, X, X, factors, case_weights, control
  );
  if (index.n_rows != y.n_elem ||
      distance.n_rows != index.n_rows ||
      distance.n_cols != index.n_cols ||
      index.n_cols == 0) {
    stop("kNN lambda-maximum neighborhoods are invalid.");
  }
  const KnnWeightProvider provider(
    index, distance, case_weights, control
  );
  return lambda_max_batch(provider, y, X, factors, control);
  END_RCPP
}

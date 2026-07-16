#include <RcppArmadillo.h>
#include <queue>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

namespace {

struct KnnCandidateLess {
  bool operator()(const std::pair<double, int>& a,
                  const std::pair<double, int>& b) const {
    if (a.first != b.first) return a.first < b.first;
    return a.second < b.second;
  }
};

using NeighborHeap = std::priority_queue<
  std::pair<double, int>,
  std::vector<std::pair<double, int> >,
  KnnCandidateLess
>;

int resolve_n_threads(int n_threads) {
  if (n_threads < 1) n_threads = 1;
#ifdef _OPENMP
  return n_threads;
#else
  return 1;
#endif
}

bool candidate_better(const std::pair<double, int>& candidate,
                      const std::pair<double, int>& current_worst) {
  if (candidate.first != current_worst.first) {
    return candidate.first < current_worst.first;
  }
  return candidate.second < current_worst.second;
}

void update_heap(NeighborHeap& heap,
                 const std::pair<double, int>& candidate,
                 int k) {
  if (k <= 0) return;
  if (static_cast<int>(heap.size()) < k) {
    heap.push(candidate);
    return;
  }
  if (candidate_better(candidate, heap.top())) {
    heap.pop();
    heap.push(candidate);
  }
}

double lorentz_distance_rows(const arma::mat& Z,
                             arma::uword i,
                             arma::uword j,
                             double tol) {
  const double inner = -Z(i, 0) * Z(j, 0) + Z(i, 1) * Z(j, 1) + Z(i, 2) * Z(j, 2);
  double arg = -inner;
  if (!R_finite(arg)) {
    stop("Invalid Lorentz distance: non-finite acosh argument.");
  }
  if (arg < 1.0) {
    if (arg >= 1.0 - tol) {
      arg = 1.0;
    } else {
      stop("Invalid Lorentz distance: acosh argument is below 1.");
    }
  }
  return std::acosh(arg);
}

double kernel_weight_one(double d, double bw, const std::string& kernel) {
  if (!R_finite(bw) || bw <= 0.0) {
    stop("Bandwidth values must be positive finite numbers.");
  }
  const double r = d / bw;
  if (kernel == "gaussian") {
    return std::exp(-0.5 * r * r);
  }
  if (kernel == "exponential") {
    return std::exp(-r);
  }
  if (kernel == "tricube") {
    if (r >= 1.0) return 0.0;
    const double u = 1.0 - r * r * r;
    return u * u * u;
  }
  if (kernel == "boxcar") {
    return r < 1.0 ? 1.0 : 0.0;
  }
  if (kernel != "bisquare") {
    stop("Unsupported kernel.");
  }
  if (r >= 1.0) return 0.0;
  const double u = 1.0 - r * r;
  return u * u;
}

std::vector<std::vector<std::pair<double, int> > >
lorentz_neighbor_lists(const arma::mat& Z, int k, double tol) {
  const arma::uword n = Z.n_rows;
  std::vector<NeighborHeap> heaps(n);
  for (arma::uword i = 0; i < n; ++i) {
    for (arma::uword j = i + 1; j < n; ++j) {
      const double dij = lorentz_distance_rows(Z, i, j, tol);
      update_heap(heaps[i], std::make_pair(dij, static_cast<int>(j)), k);
      update_heap(heaps[j], std::make_pair(dij, static_cast<int>(i)), k);
    }
  }

  std::vector<std::vector<std::pair<double, int> > > out(n);
  for (arma::uword i = 0; i < n; ++i) {
    if (static_cast<int>(heaps[i].size()) != k) {
      stop("Could not identify enough Lorentz nearest neighbors.");
    }
    out[i].reserve(k);
    while (!heaps[i].empty()) {
      out[i].push_back(heaps[i].top());
      heaps[i].pop();
    }
    std::sort(out[i].begin(), out[i].end());
  }
  return out;
}

arma::vec kernel_weights(const arma::vec& d, double bw, const std::string& kernel) {
  if (!R_finite(bw) || bw <= 0.0) {
    stop("Bandwidth values must be positive finite numbers.");
  }
  arma::vec r = d / bw;
  if (kernel == "gaussian") {
    return arma::exp(-0.5 * arma::square(r));
  }
  if (kernel == "exponential") {
    return arma::exp(-r);
  }
  if (kernel == "tricube") {
    arma::vec w = arma::pow(1.0 - arma::pow(r, 3), 3);
    w.elem(arma::find(r >= 1.0)).zeros();
    return w;
  }
  if (kernel == "boxcar") {
    arma::vec w(r.n_elem, arma::fill::ones);
    w.elem(arma::find(r >= 1.0)).zeros();
    return w;
  }
  if (kernel != "bisquare") {
    stop("Unsupported kernel.");
  }
  arma::vec w = arma::square(1.0 - arma::square(r));
  w.elem(arma::find(r >= 1.0)).zeros();
  return w;
}

arma::mat solve_xtwx(const arma::mat& XtWX) {
  arma::mat inv;
  bool ok = arma::inv_sympd(inv, XtWX);
  if (ok && inv.is_finite()) return inv;
  
  ok = arma::inv(inv, XtWX);
  if (ok && inv.is_finite()) return inv;
  
  return arma::pinv(XtWX);
}

double condition_number(const arma::mat& A) {
  arma::vec s;
  if (!arma::svd(s, A) || s.n_elem == 0) return R_PosInf;
  const double smax = s.max();
  const double smin = s.min();
  if (!R_finite(smax) || !R_finite(smin) || smin <= 0.0) return R_PosInf;
  return smax / smin;
}

double vec_median(arma::vec x) {
  x = x.elem(arma::find_finite(x));
  if (x.n_elem == 0) return NA_REAL;
  x = arma::sort(x);
  const arma::uword mid = x.n_elem / 2;
  if (x.n_elem % 2 == 1) return x[mid];
  return 0.5 * (x[mid - 1] + x[mid]);
}

double hyperbolic_ball_volume_cpp(double r) {
  return 2.0 * M_PI * (std::cosh(r) - 1.0);
}

double hyperbolic_radius_from_volume_cpp(double volume) {
  return std::acosh(1.0 + volume / (2.0 * M_PI));
}

arma::vec kth_neighbor_radius_cpp(const arma::mat& D, int k) {
  const arma::uword n = D.n_rows;
  arma::vec out(n);
  for (arma::uword j = 0; j < n; ++j) {
    arma::vec col = D.col(j);
    arma::vec positive = col.elem(arma::find(col > 0.0));
    if (positive.n_elem < static_cast<arma::uword>(k)) {
      stop("Adaptive bandwidth must be an integer between 2 and n - 1.");
    }
    positive = arma::sort(positive);
    out[j] = positive[k - 1];
  }
  return out;
}

arma::vec resolve_bandwidth_cpp(const arma::mat& D,
                                double bw_candidate,
                                bool adaptive,
                                const std::string& bandwidth_type,
                                int pilot_k) {
  const int n = static_cast<int>(D.n_rows);
  arma::vec out(n);
  if (adaptive) {
    const int target_mass = static_cast<int>(std::round(bw_candidate));
    if (target_mass < 2 || target_mass >= n) {
      stop("Adaptive bandwidth must be an integer between 2 and n - 1.");
    }
    if (bandwidth_type == "knn") {
      return kth_neighbor_radius_cpp(D, target_mass);
    }
    const int k0 = pilot_k;
    if (k0 < 2 || k0 >= n) {
      stop("`pilot_k` must be an integer between 2 and n - 1.");
    }
    arma::vec pilot_radius = kth_neighbor_radius_cpp(D, k0);
    for (int i = 0; i < n; ++i) {
      const double local_density = static_cast<double>(k0) / hyperbolic_ball_volume_cpp(pilot_radius[i]);
      const double target_volume = static_cast<double>(target_mass) / local_density;
      out[i] = hyperbolic_radius_from_volume_cpp(target_volume);
    }
    return out;
  }
  if (!R_finite(bw_candidate) || bw_candidate <= 0.0) {
    stop("Bandwidth values must be positive finite numbers.");
  }
  const double radius = bandwidth_type == "volume" ?
    hyperbolic_radius_from_volume_cpp(bw_candidate) : bw_candidate;
  out.fill(radius);
  return out;
}

struct FitSummary {
  double score;
  double trS;
  double RSS;
  double CV;
  double kappa_median;
  double kappa_max;
  double neff_min;
  double neff_median;
  double n_positive_min;
  double n_positive_median;
};

FitSummary fit_summary_cpp(const arma::vec& y,
                           const arma::mat& X,
                           const arma::mat& D,
                           const arma::vec& bw,
                           const std::string& kernel,
                           const std::string& approach) {
  const arma::uword n = y.n_elem;
  arma::vec y_hat(n, arma::fill::none);
  arma::vec s_diag(n, arma::fill::zeros);
  arma::vec s2_diag(n, arma::fill::zeros);
  arma::vec local_kappa(n, arma::fill::none);
  arma::vec local_neff(n, arma::fill::none);
  arma::vec local_n_positive(n, arma::fill::none);
  double cv = NA_REAL;
  
  for (arma::uword i = 0; i < n; ++i) {
    arma::vec w = kernel_weights(D.col(i), bw[i], kernel);
    arma::mat Xw = X.each_col() % w;
    arma::mat XtWX = X.t() * Xw;
    arma::vec XtWy = X.t() * (w % y);
    arma::mat XtWX_inv = solve_xtwx(XtWX);
    arma::vec beta = XtWX_inv * XtWy;
    
    y_hat[i] = arma::dot(X.row(i), beta.t());
    
    const double sw = arma::accu(w);
    const double sw2 = arma::accu(arma::square(w));
    local_neff[i] = (sw2 > 0.0 && R_finite(sw) && R_finite(sw2)) ? (sw * sw / sw2) : 0.0;
    local_n_positive[i] = arma::accu(w > 0.0);
    local_kappa[i] = condition_number(XtWX);
    
    arma::mat C = XtWX_inv * Xw.t();
    arma::rowvec s_row = X.row(i) * C;
    s_diag[i] = s_row[i];
    s2_diag[i] = arma::accu(arma::square(s_row));
  }
  
  arma::vec residuals = y - y_hat;
  const double rss = arma::accu(arma::square(residuals));
  const double trS = arma::accu(s_diag);
  const double denom = static_cast<double>(n) - trS - 2.0;
  double aicc = R_PosInf;
  if (R_finite(rss) && rss > 0.0 && R_finite(denom) && denom > 0.0) {
    const double dev = static_cast<double>(n) * std::log(rss / static_cast<double>(n)) +
      static_cast<double>(n) * std::log(2.0 * M_PI) + static_cast<double>(n);
    aicc = dev + 2.0 * static_cast<double>(n) * (trS + 1.0) / denom;
  }
  
  if (approach == "CV") {
    cv = 0.0;
    for (arma::uword i = 0; i < n; ++i) {
      arma::vec w = kernel_weights(D.col(i), bw[i], kernel);
      w[i] = 0.0;
      arma::mat Xw = X.each_col() % w;
      arma::mat XtWX = X.t() * Xw;
      arma::vec XtWy = X.t() * (w % y);
      arma::mat XtWX_inv = solve_xtwx(XtWX);
      arma::vec beta = XtWX_inv * XtWy;
      const double pred = arma::dot(X.row(i), beta.t());
      const double res = y[i] - pred;
      cv += res * res;
    }
  }
  
  FitSummary out;
  out.score = approach == "CV" ? cv : aicc;
  out.trS = trS;
  out.RSS = rss;
  out.CV = approach == "CV" ? cv : NA_REAL;
  out.kappa_median = vec_median(local_kappa);
  out.kappa_max = local_kappa.max();
  out.neff_min = local_neff.min();
  out.neff_median = vec_median(local_neff);
  out.n_positive_min = local_n_positive.min();
  out.n_positive_median = vec_median(local_n_positive);
  return out;
}

FitSummary fit_summary_knn_cpp(const arma::vec& y,
                               const arma::mat& X,
                               const std::vector<std::vector<std::pair<double, int> > >& neighbor_lists,
                               int k,
                               const std::string& kernel,
                               const std::string& approach,
                               int n_threads) {
  const arma::uword n = y.n_elem;
  const arma::uword p = X.n_cols;
  const arma::uword m = static_cast<arma::uword>(k + 1);
  const int threads = resolve_n_threads(n_threads);
#ifndef _OPENMP
  (void) threads;
#endif
  arma::vec y_hat(n, arma::fill::none);
  arma::vec s_diag(n, arma::fill::zeros);
  arma::vec s2_diag(n, arma::fill::zeros);
  arma::vec local_kappa(n, arma::fill::none);
  arma::vec local_neff(n, arma::fill::none);
  arma::vec local_n_positive(n, arma::fill::none);
  double cv = NA_REAL;

#pragma omp parallel for schedule(static) num_threads(threads) if(threads > 1)
  for (arma::sword ii = 0; ii < static_cast<arma::sword>(n); ++ii) {
    const arma::uword i = static_cast<arma::uword>(ii);
    const std::vector<std::pair<double, int> >& neighbors = neighbor_lists[i];
    const double bw_i = neighbors[static_cast<std::size_t>(k - 1)].first;
    if (!R_finite(bw_i) || bw_i <= 0.0) {
      stop("Adaptive kNN radius is not positive and finite.");
    }

    arma::mat X_sub(m, p, arma::fill::none);
    arma::vec y_sub(m, arma::fill::none);
    arma::vec w(m, arma::fill::none);

    X_sub.row(0) = X.row(i);
    y_sub[0] = y[i];
    w[0] = kernel_weight_one(0.0, bw_i, kernel);

    for (int h = 0; h < k; ++h) {
      const arma::uword row = static_cast<arma::uword>(h + 1);
      const int j = neighbors[static_cast<std::size_t>(h)].second;
      X_sub.row(row) = X.row(j);
      y_sub[row] = y[j];
      w[row] = kernel_weight_one(neighbors[static_cast<std::size_t>(h)].first, bw_i, kernel);
    }

    arma::mat Xw = X_sub.each_col() % w;
    arma::mat XtWX = X_sub.t() * Xw;
    arma::vec XtWy = X_sub.t() * (w % y_sub);
    arma::mat XtWX_inv = solve_xtwx(XtWX);
    arma::vec beta = XtWX_inv * XtWy;

    y_hat[i] = arma::dot(X.row(i), beta.t());

    const double sw = arma::accu(w);
    const double sw2 = arma::accu(arma::square(w));
    local_neff[i] = (sw2 > 0.0 && R_finite(sw) && R_finite(sw2)) ? (sw * sw / sw2) : 0.0;
    local_n_positive[i] = arma::accu(w > 0.0);
    local_kappa[i] = condition_number(XtWX);

    arma::mat C = XtWX_inv * Xw.t();
    arma::rowvec s_row = X.row(i) * C;
    s_diag[i] = s_row[0];
    s2_diag[i] = arma::accu(arma::square(s_row));
  }

  arma::vec residuals = y - y_hat;
  const double rss = arma::accu(arma::square(residuals));
  const double trS = arma::accu(s_diag);
  const double denom = static_cast<double>(n) - trS - 2.0;
  double aicc = R_PosInf;
  if (R_finite(rss) && rss > 0.0 && R_finite(denom) && denom > 0.0) {
    const double dev = static_cast<double>(n) * std::log(rss / static_cast<double>(n)) +
      static_cast<double>(n) * std::log(2.0 * M_PI) + static_cast<double>(n);
    aicc = dev + 2.0 * static_cast<double>(n) * (trS + 1.0) / denom;
  }

  if (approach == "CV") {
    cv = 0.0;
#pragma omp parallel for schedule(static) reduction(+:cv) num_threads(threads) if(threads > 1)
    for (arma::sword ii = 0; ii < static_cast<arma::sword>(n); ++ii) {
      const arma::uword i = static_cast<arma::uword>(ii);
      const std::vector<std::pair<double, int> >& neighbors = neighbor_lists[i];
      const double bw_i = neighbors[static_cast<std::size_t>(k - 1)].first;

      arma::mat X_sub(m, p, arma::fill::none);
      arma::vec y_sub(m, arma::fill::none);
      arma::vec w(m, arma::fill::none);

      X_sub.row(0) = X.row(i);
      y_sub[0] = y[i];
      w[0] = 0.0;

      for (int h = 0; h < k; ++h) {
        const arma::uword row = static_cast<arma::uword>(h + 1);
        const int j = neighbors[static_cast<std::size_t>(h)].second;
        X_sub.row(row) = X.row(j);
        y_sub[row] = y[j];
        w[row] = kernel_weight_one(neighbors[static_cast<std::size_t>(h)].first, bw_i, kernel);
      }

      arma::mat Xw = X_sub.each_col() % w;
      arma::mat XtWX = X_sub.t() * Xw;
      arma::vec XtWy = X_sub.t() * (w % y_sub);
      arma::mat XtWX_inv = solve_xtwx(XtWX);
      arma::vec beta = XtWX_inv * XtWy;
      const double pred = arma::dot(X.row(i), beta.t());
      const double res = y[i] - pred;
      cv += res * res;
    }
  }

  FitSummary out;
  out.score = approach == "CV" ? cv : aicc;
  out.trS = trS;
  out.RSS = rss;
  out.CV = approach == "CV" ? cv : NA_REAL;
  out.kappa_median = vec_median(local_kappa);
  out.kappa_max = local_kappa.max();
  out.neff_min = local_neff.min();
  out.neff_median = vec_median(local_neff);
  out.n_positive_min = local_n_positive.min();
  out.n_positive_median = vec_median(local_n_positive);
  return out;
}

int find_candidate(const std::vector<double>& values, double candidate) {
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (values[i] == candidate) return static_cast<int>(i);
  }
  return -1;
}

double normalize_bw_candidate(double x, double lower, double upper, bool adaptive) {
  if (adaptive) {
    x = std::round(x);
  }
  if (x < lower) x = lower;
  if (x > upper) x = upper;
  return x;
}

struct GoldenTrace {
  std::vector<double> bw;
  std::vector<double> score;
  std::vector<double> trS;
  std::vector<double> RSS;
  std::vector<double> CV;
  std::vector<double> kappa_median;
  std::vector<double> kappa_max;
  std::vector<double> neff_min;
  std::vector<double> neff_median;
  std::vector<double> n_positive_min;
  std::vector<double> n_positive_median;
  std::vector<int> iteration;
};

double evaluate_bw_candidate(const arma::vec& y,
                             const arma::mat& X,
                             const arma::mat& D,
                             double candidate,
                             const std::string& approach,
                             const std::string& kernel,
                             bool adaptive,
                             const std::string& bandwidth_type,
                             int pilot_k,
                             GoldenTrace& trace) {
  const int existing = find_candidate(trace.bw, candidate);
  if (existing >= 0) return trace.score[existing];
  
  arma::vec bw_vec = resolve_bandwidth_cpp(D, candidate, adaptive, bandwidth_type, pilot_k);
  FitSummary summary = fit_summary_cpp(y, X, D, bw_vec, kernel, approach);
  
  trace.bw.push_back(candidate);
  trace.score.push_back(summary.score);
  trace.trS.push_back(summary.trS);
  trace.RSS.push_back(summary.RSS);
  trace.CV.push_back(summary.CV);
  trace.kappa_median.push_back(summary.kappa_median);
  trace.kappa_max.push_back(summary.kappa_max);
  trace.neff_min.push_back(summary.neff_min);
  trace.neff_median.push_back(summary.neff_median);
  trace.n_positive_min.push_back(summary.n_positive_min);
  trace.n_positive_median.push_back(summary.n_positive_median);
  trace.iteration.push_back(static_cast<int>(trace.bw.size()));
  
  return summary.score;
}

} // namespace

extern "C" SEXP _LHWR_hgwr_fit_cpp(SEXP ySEXP,
                                           SEXP XSEXP,
                                           SEXP DSEXP,
                                           SEXP bwSEXP,
                                           SEXP kernelSEXP,
                                           SEXP returnHatSEXP) {
  BEGIN_RCPP
  arma::vec y = as<arma::vec>(ySEXP);
  arma::mat X = as<arma::mat>(XSEXP);
  arma::mat D = as<arma::mat>(DSEXP);
  arma::vec bw = as<arma::vec>(bwSEXP);
  const std::string kernel = as<std::string>(kernelSEXP);
  const bool return_hat = as<bool>(returnHatSEXP);
  
  const arma::uword n = y.n_elem;
  const arma::uword p = X.n_cols;
  if (X.n_rows != n || D.n_rows != n || D.n_cols != n) {
    stop("Dimensions of `y`, `X`, and `D` are not compatible.");
  }
  if (bw.n_elem == 1) {
    bw = arma::vec(n).fill(bw[0]);
  }
  if (bw.n_elem != n) {
    stop("`bw` must be scalar or have one value per observation.");
  }
  
  arma::mat beta_hat(n, p, arma::fill::none);
  arma::mat se_unscaled(n, p, arma::fill::value(NA_REAL));
  arma::vec y_hat(n, arma::fill::none);
  arma::vec s_diag(n, arma::fill::zeros);
  arma::vec s2_diag(n, arma::fill::zeros);
  arma::vec local_kappa(n, arma::fill::none);
  arma::vec local_neff(n, arma::fill::none);
  IntegerVector local_n_positive(n);
  
  for (arma::uword i = 0; i < n; ++i) {
    arma::vec w = kernel_weights(D.col(i), bw[i], kernel);
    arma::mat Xw = X.each_col() % w;
    arma::mat XtWX = X.t() * Xw;
    arma::vec XtWy = X.t() * (w % y);
    arma::mat XtWX_inv = solve_xtwx(XtWX);
    arma::vec beta = XtWX_inv * XtWy;
    
    beta_hat.row(i) = beta.t();
    y_hat[i] = arma::dot(X.row(i), beta.t());
    
    const double sw = arma::accu(w);
    const double sw2 = arma::accu(arma::square(w));
    local_neff[i] = (sw2 > 0.0 && R_finite(sw) && R_finite(sw2)) ? (sw * sw / sw2) : 0.0;
    local_n_positive[i] = static_cast<int>(arma::accu(w > 0.0));
    local_kappa[i] = condition_number(XtWX);
    
    if (return_hat) {
      arma::mat C = XtWX_inv * Xw.t();
      arma::rowvec s_row = X.row(i) * C;
      s_diag[i] = s_row[i];
      s2_diag[i] = arma::accu(arma::square(s_row));
      se_unscaled.row(i) = arma::sum(arma::square(C), 1).t();
    }
  }
  
  return List::create(
    _["beta_hat"] = beta_hat,
    _["se_unscaled"] = se_unscaled,
    _["y_hat"] = y_hat,
    _["s_diag"] = s_diag,
    _["s2_diag"] = s2_diag,
    _["local_kappa"] = local_kappa,
    _["local_neff"] = local_neff,
    _["local_n_positive"] = local_n_positive
  );
  END_RCPP
}

extern "C" SEXP _LHWR_hgwr_fit_knn_cpp(SEXP ySEXP,
                                               SEXP XSEXP,
                                               SEXP ZSEXP,
                                               SEXP kSEXP,
                                               SEXP kernelSEXP,
                                               SEXP returnHatSEXP,
                                               SEXP tolSEXP,
                                               SEXP nThreadsSEXP) {
  BEGIN_RCPP
  arma::vec y = as<arma::vec>(ySEXP);
  arma::mat X = as<arma::mat>(XSEXP);
  arma::mat Z = as<arma::mat>(ZSEXP);
  const int k = as<int>(kSEXP);
  const std::string kernel = as<std::string>(kernelSEXP);
  const bool return_hat = as<bool>(returnHatSEXP);
  const double tol = as<double>(tolSEXP);
  const int n_threads = as<int>(nThreadsSEXP);
  const int threads = resolve_n_threads(n_threads);
#ifndef _OPENMP
  (void) threads;
#endif

  const arma::uword n = y.n_elem;
  const arma::uword p = X.n_cols;
  if (X.n_rows != n || Z.n_rows != n || Z.n_cols != 3) {
    stop("Dimensions of `y`, `X`, and `Z` are not compatible.");
  }
  if (k < 2 || k >= static_cast<int>(n)) {
    stop("`k` must be an integer between 2 and n - 1.");
  }
  if (!R_finite(tol) || tol < 0.0) {
    stop("`tol` must be a non-negative finite number.");
  }

  std::vector<std::vector<std::pair<double, int> > > neighbor_lists =
    lorentz_neighbor_lists(Z, k, tol);

  arma::mat beta_hat(n, p, arma::fill::none);
  arma::mat se_unscaled(n, p, arma::fill::value(NA_REAL));
  arma::vec y_hat(n, arma::fill::none);
  arma::vec s_diag(n, arma::fill::zeros);
  arma::vec s2_diag(n, arma::fill::zeros);
  arma::vec local_kappa(n, arma::fill::none);
  arma::vec local_neff(n, arma::fill::none);
  arma::ivec local_n_positive(n, arma::fill::zeros);
  NumericVector bandwidth(n);
  IntegerMatrix knn_index(n, k);
  NumericMatrix knn_distance(n, k);

  for (arma::uword i = 0; i < n; ++i) {
    const std::vector<std::pair<double, int> >& neighbors = neighbor_lists[i];
    for (int h = 0; h < k; ++h) {
      knn_index(i, h) = neighbors[h].second + 1;
      knn_distance(i, h) = neighbors[h].first;
    }
  }

  const arma::uword m = static_cast<arma::uword>(k + 1);
#pragma omp parallel for schedule(static) num_threads(threads) if(threads > 1)
  for (arma::sword ii = 0; ii < static_cast<arma::sword>(n); ++ii) {
    const arma::uword i = static_cast<arma::uword>(ii);
    const std::vector<std::pair<double, int> >& neighbors = neighbor_lists[i];

    const double bw_i = neighbors.back().first;
    if (!R_finite(bw_i) || bw_i <= 0.0) {
      stop("Adaptive kNN radius is not positive and finite.");
    }
    bandwidth[i] = bw_i;

    arma::mat X_sub(m, p, arma::fill::none);
    arma::vec y_sub(m, arma::fill::none);
    arma::vec w(m, arma::fill::none);

    X_sub.row(0) = X.row(i);
    y_sub[0] = y[i];
    w[0] = kernel_weight_one(0.0, bw_i, kernel);

    for (int h = 0; h < k; ++h) {
      const arma::uword row = static_cast<arma::uword>(h + 1);
      const int j = neighbors[h].second;
      X_sub.row(row) = X.row(j);
      y_sub[row] = y[j];
      w[row] = kernel_weight_one(neighbors[h].first, bw_i, kernel);
    }

    arma::mat Xw = X_sub.each_col() % w;
    arma::mat XtWX = X_sub.t() * Xw;
    arma::vec XtWy = X_sub.t() * (w % y_sub);
    arma::mat XtWX_inv = solve_xtwx(XtWX);
    arma::vec beta = XtWX_inv * XtWy;

    beta_hat.row(i) = beta.t();
    y_hat[i] = arma::dot(X.row(i), beta.t());

    const double sw = arma::accu(w);
    const double sw2 = arma::accu(arma::square(w));
    local_neff[i] = (sw2 > 0.0 && R_finite(sw) && R_finite(sw2)) ? (sw * sw / sw2) : 0.0;
    local_n_positive[i] = static_cast<int>(arma::accu(w > 0.0));
    local_kappa[i] = condition_number(XtWX);

    if (return_hat) {
      arma::mat C = XtWX_inv * Xw.t();
      arma::rowvec s_row = X.row(i) * C;
      s_diag[i] = s_row[0];
      s2_diag[i] = arma::accu(arma::square(s_row));
      se_unscaled.row(i) = arma::sum(arma::square(C), 1).t();
    }
  }

  return List::create(
    _["beta_hat"] = beta_hat,
    _["se_unscaled"] = se_unscaled,
    _["y_hat"] = y_hat,
    _["s_diag"] = s_diag,
    _["s2_diag"] = s2_diag,
    _["local_kappa"] = local_kappa,
    _["local_neff"] = local_neff,
    _["local_n_positive"] = local_n_positive,
    _["bandwidth"] = bandwidth,
    _["knn_index"] = knn_index,
    _["knn_distance"] = knn_distance
  );
  END_RCPP
}

extern "C" SEXP _LHWR_hgwr_cv_score_cpp(SEXP ySEXP,
                                                SEXP XSEXP,
                                                SEXP DSEXP,
                                                SEXP bwSEXP,
                                                SEXP kernelSEXP) {
  BEGIN_RCPP
  arma::vec y = as<arma::vec>(ySEXP);
  arma::mat X = as<arma::mat>(XSEXP);
  arma::mat D = as<arma::mat>(DSEXP);
  arma::vec bw = as<arma::vec>(bwSEXP);
  const std::string kernel = as<std::string>(kernelSEXP);
  
  const arma::uword n = y.n_elem;
  if (X.n_rows != n || D.n_rows != n || D.n_cols != n) {
    stop("Dimensions of `y`, `X`, and `D` are not compatible.");
  }
  if (bw.n_elem == 1) {
    bw = arma::vec(n).fill(bw[0]);
  }
  if (bw.n_elem != n) {
    stop("`bw` must be scalar or have one value per observation.");
  }
  
  double cv = 0.0;
  for (arma::uword i = 0; i < n; ++i) {
    arma::vec w = kernel_weights(D.col(i), bw[i], kernel);
    w[i] = 0.0;
    arma::mat Xw = X.each_col() % w;
    arma::mat XtWX = X.t() * Xw;
    arma::vec XtWy = X.t() * (w % y);
    arma::mat XtWX_inv = solve_xtwx(XtWX);
    arma::vec beta = XtWX_inv * XtWy;
    const double pred = arma::dot(X.row(i), beta.t());
    const double res = y[i] - pred;
    cv += res * res;
  }
  
  return wrap(cv);
  END_RCPP
}

extern "C" SEXP _LHWR_hgwr_f_tests_cpp(SEXP ySEXP,
                                                SEXP XSEXP,
                                                SEXP DSEXP,
                                                SEXP bwSEXP,
                                                SEXP kernelSEXP,
                                                SEXP betaHatSEXP,
                                                SEXP rssFullSEXP,
                                                SEXP rssOlsSEXP,
                                                SEXP dfOlsSEXP,
                                                SEXP trSSEXP,
                                                SEXP trStSSEXP) {
  BEGIN_RCPP
  arma::vec y = as<arma::vec>(ySEXP);
  arma::mat X = as<arma::mat>(XSEXP);
  arma::mat D = as<arma::mat>(DSEXP);
  arma::vec bw = as<arma::vec>(bwSEXP);
  const std::string kernel = as<std::string>(kernelSEXP);
  arma::mat beta_hat = as<arma::mat>(betaHatSEXP);
  const double RSS_full = as<double>(rssFullSEXP);
  const double RSS_ols = as<double>(rssOlsSEXP);
  const double DF_ols = as<double>(dfOlsSEXP);
  const double trS = as<double>(trSSEXP);
  const double trStS = as<double>(trStSSEXP);

  const arma::uword n = y.n_elem;
  const arma::uword p = X.n_cols;
  if (X.n_rows != n || D.n_rows != n || D.n_cols != n || beta_hat.n_rows != n || beta_hat.n_cols != p) {
    stop("Dimensions of `y`, `X`, `D`, and `beta_hat` are not compatible.");
  }
  if (bw.n_elem == 1) {
    bw = arma::vec(n).fill(bw[0]);
  }
  if (bw.n_elem != n) {
    stop("`bw` must be scalar or have one value per observation.");
  }

  arma::vec q_diag(n, arma::fill::ones);
  arma::mat B_sum(p, n, arma::fill::zeros);
  arma::mat B_sumsq(p, n, arma::fill::zeros);

  for (arma::uword i = 0; i < n; ++i) {
    arma::vec w = kernel_weights(D.col(i), bw[i], kernel);
    arma::mat Xw = X.each_col() % w;
    arma::mat XtWX = X.t() * Xw;
    arma::mat XtWX_inv = solve_xtwx(XtWX);
    arma::mat C = XtWX_inv * Xw.t();
    arma::rowvec s_row = X.row(i) * C;

    q_diag += arma::square(s_row.t());
    q_diag[i] -= 2.0 * s_row[i];

    B_sum += C;
    B_sumsq += arma::square(C);
  }

  const double L_delta1 = arma::accu(q_diag);
  const double L_delta2 = arma::accu(arma::square(q_diag));
  const double delta1 = static_cast<double>(n) - 2.0 * trS + trStS;
  const double sigma2_delta1 = RSS_full / delta1;

  const double F1 = (RSS_full / L_delta1) / (RSS_ols / DF_ols);
  const double F1_df1 = L_delta1 * L_delta1 / L_delta2;
  const double F1_df2 = DF_ols;
  const double F1_p = R::pf(F1, F1_df1, F1_df2, 1, 0);

  const double F2 = ((RSS_ols - RSS_full) / (DF_ols - L_delta1)) / (RSS_ols / DF_ols);
  const double F2_df1 = (DF_ols - L_delta1) * (DF_ols - L_delta1) /
    (DF_ols - 2.0 * L_delta1 + L_delta2);
  const double F2_df2 = DF_ols;
  const double F2_p = R::pf(F2, F2_df1, F2_df2, 0, 0);

  arma::vec F3(p, arma::fill::value(NA_REAL));
  arma::vec F3_df1(p, arma::fill::value(NA_REAL));
  arma::vec F3_df2(p, arma::fill::value(F1_df1));
  arma::vec F3_p(p, arma::fill::value(NA_REAL));

  for (arma::uword k = 0; k < p; ++k) {
    arma::vec beta_k = beta_hat.col(k);
    const double beta_sum = arma::accu(beta_k);
    const double Vk2 = (arma::accu(arma::square(beta_k)) -
      beta_sum * beta_sum / static_cast<double>(n)) / static_cast<double>(n);

    arma::rowvec diag_BJ = (B_sumsq.row(k) -
      arma::square(B_sum.row(k)) / static_cast<double>(n)) / static_cast<double>(n);
    const double gamma1 = arma::accu(diag_BJ);
    const double gamma2 = arma::accu(arma::square(diag_BJ));
    const double numdf = gamma1 * gamma1 / gamma2;

    F3[k] = (Vk2 / gamma1) / sigma2_delta1;
    F3_df1[k] = numdf;
    F3_p[k] = R::pf(F3[k], F3_df1[k], F3_df2[k], 0, 0);
  }

  return List::create(
    _["F1"] = F1,
    _["F1_df1"] = F1_df1,
    _["F1_df2"] = F1_df2,
    _["F1_p"] = F1_p,
    _["F2"] = F2,
    _["F2_df1"] = F2_df1,
    _["F2_df2"] = F2_df2,
    _["F2_p"] = F2_p,
    _["F3"] = F3,
    _["F3_df1"] = F3_df1,
    _["F3_df2"] = F3_df2,
    _["F3_p"] = F3_p,
    _["q_diag"] = q_diag,
    _["L_delta1"] = L_delta1,
    _["L_delta2"] = L_delta2,
    _["delta1"] = delta1
  );
  END_RCPP
}

extern "C" SEXP _LHWR_hgwr_f_tests_knn_cpp(SEXP ySEXP,
                                                    SEXP XSEXP,
                                                    SEXP ZSEXP,
                                                    SEXP kSEXP,
                                                    SEXP kernelSEXP,
                                                    SEXP betaHatSEXP,
                                                    SEXP rssFullSEXP,
                                                    SEXP rssOlsSEXP,
                                                    SEXP dfOlsSEXP,
                                                    SEXP trSSEXP,
                                                    SEXP trStSSEXP,
                                                    SEXP tolSEXP) {
  BEGIN_RCPP
  arma::vec y = as<arma::vec>(ySEXP);
  arma::mat X = as<arma::mat>(XSEXP);
  arma::mat Z = as<arma::mat>(ZSEXP);
  const int k = as<int>(kSEXP);
  const std::string kernel = as<std::string>(kernelSEXP);
  arma::mat beta_hat = as<arma::mat>(betaHatSEXP);
  const double RSS_full = as<double>(rssFullSEXP);
  const double RSS_ols = as<double>(rssOlsSEXP);
  const double DF_ols = as<double>(dfOlsSEXP);
  const double trS = as<double>(trSSEXP);
  const double trStS = as<double>(trStSSEXP);
  const double tol = as<double>(tolSEXP);

  const arma::uword n = y.n_elem;
  const arma::uword p = X.n_cols;
  if (X.n_rows != n || Z.n_rows != n || Z.n_cols != 3 || beta_hat.n_rows != n || beta_hat.n_cols != p) {
    stop("Dimensions of `y`, `X`, `Z`, and `beta_hat` are not compatible.");
  }
  if (k < 2 || k >= static_cast<int>(n)) {
    stop("`k` must be an integer between 2 and n - 1.");
  }
  if (!R_finite(tol) || tol < 0.0) {
    stop("`tol` must be a non-negative finite number.");
  }

  std::vector<std::vector<std::pair<double, int> > > neighbor_lists =
    lorentz_neighbor_lists(Z, k, tol);

  arma::vec q_diag(n, arma::fill::ones);
  arma::mat B_sum(p, n, arma::fill::zeros);
  arma::mat B_sumsq(p, n, arma::fill::zeros);

  const arma::uword m = static_cast<arma::uword>(k + 1);
  for (arma::uword i = 0; i < n; ++i) {
    const std::vector<std::pair<double, int> >& neighbors = neighbor_lists[i];
    const double bw_i = neighbors.back().first;
    if (!R_finite(bw_i) || bw_i <= 0.0) {
      stop("Adaptive kNN radius is not positive and finite.");
    }

    arma::uvec idx(m);
    arma::mat X_sub(m, p, arma::fill::none);
    arma::vec w(m, arma::fill::none);

    idx[0] = i;
    X_sub.row(0) = X.row(i);
    w[0] = kernel_weight_one(0.0, bw_i, kernel);

    for (int h = 0; h < k; ++h) {
      const arma::uword row = static_cast<arma::uword>(h + 1);
      const arma::uword j = static_cast<arma::uword>(neighbors[h].second);
      idx[row] = j;
      X_sub.row(row) = X.row(j);
      w[row] = kernel_weight_one(neighbors[h].first, bw_i, kernel);
    }

    arma::mat Xw = X_sub.each_col() % w;
    arma::mat XtWX = X_sub.t() * Xw;
    arma::mat XtWX_inv = solve_xtwx(XtWX);
    arma::mat C = XtWX_inv * Xw.t();
    arma::rowvec s_row = X.row(i) * C;

    for (arma::uword h = 0; h < m; ++h) {
      const arma::uword j = idx[h];
      q_diag[j] += s_row[h] * s_row[h];
      B_sum.col(j) += C.col(h);
      B_sumsq.col(j) += arma::square(C.col(h));
    }
    q_diag[i] -= 2.0 * s_row[0];
  }

  const double L_delta1 = arma::accu(q_diag);
  const double L_delta2 = arma::accu(arma::square(q_diag));
  const double delta1 = static_cast<double>(n) - 2.0 * trS + trStS;
  const double sigma2_delta1 = RSS_full / delta1;

  const double F1 = (RSS_full / L_delta1) / (RSS_ols / DF_ols);
  const double F1_df1 = L_delta1 * L_delta1 / L_delta2;
  const double F1_df2 = DF_ols;
  const double F1_p = R::pf(F1, F1_df1, F1_df2, 1, 0);

  const double F2 = ((RSS_ols - RSS_full) / (DF_ols - L_delta1)) / (RSS_ols / DF_ols);
  const double F2_df1 = (DF_ols - L_delta1) * (DF_ols - L_delta1) /
    (DF_ols - 2.0 * L_delta1 + L_delta2);
  const double F2_df2 = DF_ols;
  const double F2_p = R::pf(F2, F2_df1, F2_df2, 0, 0);

  arma::vec F3(p, arma::fill::value(NA_REAL));
  arma::vec F3_df1(p, arma::fill::value(NA_REAL));
  arma::vec F3_df2(p, arma::fill::value(F1_df1));
  arma::vec F3_p(p, arma::fill::value(NA_REAL));

  for (arma::uword kk = 0; kk < p; ++kk) {
    arma::vec beta_k = beta_hat.col(kk);
    const double beta_sum = arma::accu(beta_k);
    const double Vk2 = (arma::accu(arma::square(beta_k)) -
      beta_sum * beta_sum / static_cast<double>(n)) / static_cast<double>(n);

    arma::rowvec diag_BJ = (B_sumsq.row(kk) -
      arma::square(B_sum.row(kk)) / static_cast<double>(n)) / static_cast<double>(n);
    const double gamma1 = arma::accu(diag_BJ);
    const double gamma2 = arma::accu(arma::square(diag_BJ));
    const double numdf = gamma1 * gamma1 / gamma2;

    F3[kk] = (Vk2 / gamma1) / sigma2_delta1;
    F3_df1[kk] = numdf;
    F3_p[kk] = R::pf(F3[kk], F3_df1[kk], F3_df2[kk], 0, 0);
  }

  return List::create(
    _["F1"] = F1,
    _["F1_df1"] = F1_df1,
    _["F1_df2"] = F1_df2,
    _["F1_p"] = F1_p,
    _["F2"] = F2,
    _["F2_df1"] = F2_df1,
    _["F2_df2"] = F2_df2,
    _["F2_p"] = F2_p,
    _["F3"] = F3,
    _["F3_df1"] = F3_df1,
    _["F3_df2"] = F3_df2,
    _["F3_p"] = F3_p,
    _["q_diag"] = q_diag,
    _["L_delta1"] = L_delta1,
    _["L_delta2"] = L_delta2,
    _["delta1"] = delta1
  );
  END_RCPP
}

extern "C" SEXP _LHWR_hgwr_f_tests_knn_precomputed_cpp(SEXP ySEXP,
                                                                SEXP XSEXP,
                                                                SEXP indexSEXP,
                                                                SEXP distanceSEXP,
                                                                SEXP kernelSEXP,
                                                                SEXP betaHatSEXP,
                                                                SEXP rssFullSEXP,
                                                                SEXP rssOlsSEXP,
                                                                SEXP dfOlsSEXP,
                                                                SEXP trSSEXP,
                                                                SEXP trStSSEXP) {
  BEGIN_RCPP
  arma::vec y = as<arma::vec>(ySEXP);
  arma::mat X = as<arma::mat>(XSEXP);
  IntegerMatrix index(indexSEXP);
  NumericMatrix distance(distanceSEXP);
  const std::string kernel = as<std::string>(kernelSEXP);
  arma::mat beta_hat = as<arma::mat>(betaHatSEXP);
  const double RSS_full = as<double>(rssFullSEXP);
  const double RSS_ols = as<double>(rssOlsSEXP);
  const double DF_ols = as<double>(dfOlsSEXP);
  const double trS = as<double>(trSSEXP);
  const double trStS = as<double>(trStSSEXP);

  const arma::uword n = y.n_elem;
  const arma::uword p = X.n_cols;
  const int k = index.ncol();
  if (X.n_rows != n || index.nrow() != static_cast<int>(n) ||
      distance.nrow() != static_cast<int>(n) || distance.ncol() != k ||
      beta_hat.n_rows != n || beta_hat.n_cols != p) {
    stop("Dimensions of `y`, `X`, kNN objects, and `beta_hat` are not compatible.");
  }
  if (k < 2 || k >= static_cast<int>(n)) {
    stop("kNN objects must contain between 2 and n - 1 neighbors.");
  }

  arma::vec q_diag(n, arma::fill::ones);
  arma::mat B_sum(p, n, arma::fill::zeros);
  arma::mat B_sumsq(p, n, arma::fill::zeros);

  const arma::uword m = static_cast<arma::uword>(k + 1);
  for (arma::uword i = 0; i < n; ++i) {
    const double bw_i = distance(i, k - 1);
    if (!R_finite(bw_i) || bw_i <= 0.0) {
      stop("Adaptive kNN radius is not positive and finite.");
    }

    arma::uvec idx(m);
    arma::mat X_sub(m, p, arma::fill::none);
    arma::vec w(m, arma::fill::none);

    idx[0] = i;
    X_sub.row(0) = X.row(i);
    w[0] = kernel_weight_one(0.0, bw_i, kernel);

    for (int h = 0; h < k; ++h) {
      const int j0 = index(i, h) - 1;
      if (j0 < 0 || j0 >= static_cast<int>(n)) {
        stop("kNN index contains an out-of-range neighbor.");
      }
      const arma::uword row = static_cast<arma::uword>(h + 1);
      const arma::uword j = static_cast<arma::uword>(j0);
      idx[row] = j;
      X_sub.row(row) = X.row(j);
      w[row] = kernel_weight_one(distance(i, h), bw_i, kernel);
    }

    arma::mat Xw = X_sub.each_col() % w;
    arma::mat XtWX = X_sub.t() * Xw;
    arma::mat XtWX_inv = solve_xtwx(XtWX);
    arma::mat C = XtWX_inv * Xw.t();
    arma::rowvec s_row = X.row(i) * C;

    for (arma::uword h = 0; h < m; ++h) {
      const arma::uword j = idx[h];
      q_diag[j] += s_row[h] * s_row[h];
      B_sum.col(j) += C.col(h);
      B_sumsq.col(j) += arma::square(C.col(h));
    }
    q_diag[i] -= 2.0 * s_row[0];
  }

  const double L_delta1 = arma::accu(q_diag);
  const double L_delta2 = arma::accu(arma::square(q_diag));
  const double delta1 = static_cast<double>(n) - 2.0 * trS + trStS;
  const double sigma2_delta1 = RSS_full / delta1;

  const double F1 = (RSS_full / L_delta1) / (RSS_ols / DF_ols);
  const double F1_df1 = L_delta1 * L_delta1 / L_delta2;
  const double F1_df2 = DF_ols;
  const double F1_p = R::pf(F1, F1_df1, F1_df2, 1, 0);

  const double F2 = ((RSS_ols - RSS_full) / (DF_ols - L_delta1)) / (RSS_ols / DF_ols);
  const double F2_df1 = (DF_ols - L_delta1) * (DF_ols - L_delta1) /
    (DF_ols - 2.0 * L_delta1 + L_delta2);
  const double F2_df2 = DF_ols;
  const double F2_p = R::pf(F2, F2_df1, F2_df2, 0, 0);

  arma::vec F3(p, arma::fill::value(NA_REAL));
  arma::vec F3_df1(p, arma::fill::value(NA_REAL));
  arma::vec F3_df2(p, arma::fill::value(F1_df1));
  arma::vec F3_p(p, arma::fill::value(NA_REAL));

  for (arma::uword kk = 0; kk < p; ++kk) {
    arma::vec beta_k = beta_hat.col(kk);
    const double beta_sum = arma::accu(beta_k);
    const double Vk2 = (arma::accu(arma::square(beta_k)) -
      beta_sum * beta_sum / static_cast<double>(n)) / static_cast<double>(n);

    arma::rowvec diag_BJ = (B_sumsq.row(kk) -
      arma::square(B_sum.row(kk)) / static_cast<double>(n)) / static_cast<double>(n);
    const double gamma1 = arma::accu(diag_BJ);
    const double gamma2 = arma::accu(arma::square(diag_BJ));
    const double numdf = gamma1 * gamma1 / gamma2;

    F3[kk] = (Vk2 / gamma1) / sigma2_delta1;
    F3_df1[kk] = numdf;
    F3_p[kk] = R::pf(F3[kk], F3_df1[kk], F3_df2[kk], 0, 0);
  }

  return List::create(
    _["F1"] = F1,
    _["F1_df1"] = F1_df1,
    _["F1_df2"] = F1_df2,
    _["F1_p"] = F1_p,
    _["F2"] = F2,
    _["F2_df1"] = F2_df1,
    _["F2_df2"] = F2_df2,
    _["F2_p"] = F2_p,
    _["F3"] = F3,
    _["F3_df1"] = F3_df1,
    _["F3_df2"] = F3_df2,
    _["F3_p"] = F3_p,
    _["q_diag"] = q_diag,
    _["L_delta1"] = L_delta1,
    _["L_delta2"] = L_delta2,
    _["delta1"] = delta1
  );
  END_RCPP
}

extern "C" SEXP _LHWR_local_collinearity_knn_cpp(SEXP XSEXP,
                                                         SEXP indexSEXP,
                                                         SEXP distanceSEXP,
                                                         SEXP bandwidthSEXP,
                                                         SEXP kernelSEXP,
                                                         SEXP centerSEXP,
                                                         SEXP scaleSEXP,
                                                         SEXP nThreadsSEXP) {
  BEGIN_RCPP
  arma::mat X = as<arma::mat>(XSEXP);
  IntegerMatrix index(indexSEXP);
  NumericMatrix distance(distanceSEXP);
  arma::vec bandwidth = as<arma::vec>(bandwidthSEXP);
  const std::string kernel = as<std::string>(kernelSEXP);
  const bool center = as<bool>(centerSEXP);
  const bool scale = as<bool>(scaleSEXP);
  const int n_threads = as<int>(nThreadsSEXP);
  (void) n_threads;

  const arma::uword n = X.n_rows;
  const arma::uword p = X.n_cols;
  const int k = index.ncol();
  if (p < 1) {
    stop("`X` must contain at least one predictor column.");
  }
  if (index.nrow() != static_cast<int>(n) ||
      distance.nrow() != static_cast<int>(n) ||
      distance.ncol() != k ||
      bandwidth.n_elem != n) {
    stop("Dimensions of `X`, kNN objects, and `bandwidth` are not compatible.");
  }
  if (k < 1 || k >= static_cast<int>(n)) {
    stop("kNN objects must contain between 1 and n - 1 neighbors.");
  }

  arma::vec local_CN(n, arma::fill::value(NA_REAL));
  arma::mat local_VIF(n, p, arma::fill::value(NA_REAL));
  arma::vec max_abs_corr(n, arma::fill::value(NA_REAL));
  arma::vec neff(n, arma::fill::zeros);
  arma::ivec n_positive(n, arma::fill::zeros);

  const arma::uword m = static_cast<arma::uword>(k + 1);
  const double eig_tol = 1e-10;

  for (arma::uword i = 0; i < n; ++i) {
    const double bw_i = bandwidth[i];
    if (!R_finite(bw_i) || bw_i <= 0.0) {
      stop("Stored bandwidths must be positive finite numbers.");
    }

    arma::mat X_sub(m, p, arma::fill::none);
    arma::vec w(m, arma::fill::none);
    X_sub.row(0) = X.row(i);
    w[0] = kernel_weight_one(0.0, bw_i, kernel);

    for (int h = 0; h < k; ++h) {
      const int j0 = index(i, h) - 1;
      if (j0 < 0 || j0 >= static_cast<int>(n)) {
        stop("kNN index contains an out-of-range neighbor.");
      }
      const arma::uword row = static_cast<arma::uword>(h + 1);
      X_sub.row(row) = X.row(static_cast<arma::uword>(j0));
      w[row] = kernel_weight_one(distance(i, h), bw_i, kernel);
    }

    const double sw = arma::accu(w);
    const double sw2 = arma::accu(arma::square(w));
    n_positive[i] = static_cast<int>(arma::accu(w > 0.0));
    neff[i] = (sw2 > 0.0 && R_finite(sw) && R_finite(sw2)) ? (sw * sw / sw2) : 0.0;

    arma::mat X_local = X_sub;
    if (center) {
      arma::rowvec mu = (w.t() * X_local) / sw;
      X_local.each_row() -= mu;
    }
    if (scale) {
      arma::rowvec ss = arma::sqrt((w.t() * arma::square(X_local)) / sw);
      for (arma::uword j = 0; j < p; ++j) {
        if (!R_finite(ss[j]) || ss[j] <= std::numeric_limits<double>::epsilon()) {
          ss[j] = 1.0;
        }
      }
      X_local.each_row() /= ss;
    }

    arma::mat Xw = X_local.each_col() % w;
    arma::mat G = X_local.t() * Xw;
    arma::vec eigval;
    bool eig_ok = arma::eig_sym(eigval, G);
    if (!eig_ok || eigval.n_elem == 0 || !eigval.is_finite()) {
      local_CN[i] = R_PosInf;
    } else {
      const double min_ev = eigval.min();
      const double max_ev = eigval.max();
      local_CN[i] = (min_ev <= eig_tol || max_ev <= eig_tol) ?
        R_PosInf : std::sqrt(max_ev / min_ev);
    }

    if (p > 1) {
      arma::mat C = G / sw;
      arma::vec diag_c = C.diag();
      double max_corr = NA_REAL;
      bool any_corr = false;
      for (arma::uword a = 0; a < p; ++a) {
        for (arma::uword b = a + 1; b < p; ++b) {
          const double denom = std::sqrt(diag_c[a] * diag_c[b]);
          if (R_finite(denom) && denom > std::numeric_limits<double>::epsilon()) {
            const double val = std::abs(C(a, b) / denom);
            if (R_finite(val)) {
              max_corr = any_corr ? std::max(max_corr, val) : val;
              any_corr = true;
            }
          }
        }
      }
      max_abs_corr[i] = any_corr ? max_corr : NA_REAL;
    }

    if (p == 1) {
      local_VIF(i, 0) = NA_REAL;
    } else {
      for (arma::uword target = 0; target < p; ++target) {
        const arma::uword q = p;
        if (n_positive[i] <= static_cast<int>(q)) {
          local_VIF(i, target) = R_PosInf;
          continue;
        }

        arma::mat aux_X(m, q, arma::fill::ones);
        arma::uword col = 1;
        for (arma::uword jj = 0; jj < p; ++jj) {
          if (jj == target) continue;
          aux_X.col(col) = X_sub.col(jj);
          col++;
        }
        arma::vec aux_y = X_sub.col(target);
        const double y_bar = arma::dot(w, aux_y) / sw;
        const double tss = arma::accu(w % arma::square(aux_y - y_bar));
        if (!R_finite(tss) || tss <= std::numeric_limits<double>::epsilon()) {
          local_VIF(i, target) = R_PosInf;
          continue;
        }
        arma::mat aux_Xw = aux_X.each_col() % w;
        arma::mat XtWX = aux_X.t() * aux_Xw;
        arma::vec XtWy = aux_X.t() * (w % aux_y);
        arma::vec beta = solve_xtwx(XtWX) * XtWy;
        arma::vec yhat = aux_X * beta;
        const double rss = arma::accu(w % arma::square(aux_y - yhat));
        double r2 = 1.0 - rss / tss;
        if (!R_finite(r2)) {
          local_VIF(i, target) = R_PosInf;
        } else {
          if (r2 < 0.0) r2 = 0.0;
          if (r2 > 1.0) r2 = 1.0;
          local_VIF(i, target) = (r2 >= 1.0 - 1e-10) ? R_PosInf : 1.0 / (1.0 - r2);
        }
      }
    }
  }

  return List::create(
    _["local_CN"] = local_CN,
    _["local_VIF"] = local_VIF,
    _["max_abs_local_correlation"] = max_abs_corr,
    _["neff"] = neff,
    _["n_positive"] = n_positive
  );
  END_RCPP
}

extern "C" SEXP _LHWR_bw_hgwr_grid_cpp(SEXP ySEXP,
                                               SEXP XSEXP,
                                               SEXP DSEXP,
                                               SEXP bwGridSEXP,
                                               SEXP approachSEXP,
                                               SEXP kernelSEXP,
                                               SEXP adaptiveSEXP,
                                               SEXP bandwidthTypeSEXP,
                                               SEXP pilotKSEXP) {
  BEGIN_RCPP
  arma::vec y = as<arma::vec>(ySEXP);
  arma::mat X = as<arma::mat>(XSEXP);
  arma::mat D = as<arma::mat>(DSEXP);
  arma::vec bw_grid = as<arma::vec>(bwGridSEXP);
  const std::string approach = as<std::string>(approachSEXP);
  const std::string kernel = as<std::string>(kernelSEXP);
  const bool adaptive = as<bool>(adaptiveSEXP);
  const std::string bandwidth_type = as<std::string>(bandwidthTypeSEXP);
  const int pilot_k = as<int>(pilotKSEXP);
  
  const arma::uword n = y.n_elem;
  if (X.n_rows != n || D.n_rows != n || D.n_cols != n) {
    stop("Dimensions of `y`, `X`, and `D` are not compatible.");
  }
  
  const int B = static_cast<int>(bw_grid.n_elem);
  NumericVector bw(B), score(B), trS(B), RSS(B), CV(B);
  NumericVector kappa_median(B), kappa_max(B), neff_min(B), neff_median(B);
  NumericVector n_positive_min(B), n_positive_median(B);
  
  for (int b = 0; b < B; ++b) {
    const double bw_candidate = bw_grid[b];
    arma::vec bw_vec = resolve_bandwidth_cpp(
      D,
      bw_candidate,
      adaptive,
      bandwidth_type,
      pilot_k
    );
    FitSummary summary = fit_summary_cpp(y, X, D, bw_vec, kernel, approach);
    
    bw[b] = bw_candidate;
    score[b] = summary.score;
    trS[b] = summary.trS;
    RSS[b] = summary.RSS;
    CV[b] = summary.CV;
    kappa_median[b] = summary.kappa_median;
    kappa_max[b] = summary.kappa_max;
    neff_min[b] = summary.neff_min;
    neff_median[b] = summary.neff_median;
    n_positive_min[b] = summary.n_positive_min;
    n_positive_median[b] = summary.n_positive_median;
  }
  
  return DataFrame::create(
    _["bw"] = bw,
    _["score"] = score,
    _["trS"] = trS,
    _["RSS"] = RSS,
    _["CV"] = CV,
    _["kappa_median"] = kappa_median,
    _["kappa_max"] = kappa_max,
    _["neff_min"] = neff_min,
    _["neff_median"] = neff_median,
    _["n_positive_min"] = n_positive_min,
    _["n_positive_median"] = n_positive_median
  );
  END_RCPP
}

extern "C" SEXP _LHWR_bw_hgwr_knn_grid_cpp(SEXP ySEXP,
                                                   SEXP XSEXP,
                                                   SEXP ZSEXP,
                                                   SEXP bwGridSEXP,
                                                   SEXP approachSEXP,
                                                   SEXP kernelSEXP,
                                                   SEXP tolSEXP,
                                                   SEXP nThreadsSEXP) {
  BEGIN_RCPP
  arma::vec y = as<arma::vec>(ySEXP);
  arma::mat X = as<arma::mat>(XSEXP);
  arma::mat Z = as<arma::mat>(ZSEXP);
  arma::vec bw_grid = as<arma::vec>(bwGridSEXP);
  const std::string approach = as<std::string>(approachSEXP);
  const std::string kernel = as<std::string>(kernelSEXP);
  const double tol = as<double>(tolSEXP);
  const int n_threads = as<int>(nThreadsSEXP);

  const arma::uword n = y.n_elem;
  if (X.n_rows != n || Z.n_rows != n || Z.n_cols != 3) {
    stop("Dimensions of `y`, `X`, and `Z` are not compatible.");
  }
  if (!R_finite(tol) || tol < 0.0) {
    stop("`tol` must be a non-negative finite number.");
  }

  const int B = static_cast<int>(bw_grid.n_elem);
  if (B < 1) {
    stop("`bw_grid` must contain at least one candidate.");
  }
  std::vector<int> candidates(B);
  int max_k = 0;
  for (int b = 0; b < B; ++b) {
    const int k = static_cast<int>(std::round(bw_grid[b]));
    if (k < 2 || k >= static_cast<int>(n)) {
      stop("Adaptive bandwidth candidates must be integers between 2 and n - 1.");
    }
    candidates[b] = k;
    if (k > max_k) max_k = k;
  }

  std::vector<std::vector<std::pair<double, int> > > neighbor_lists =
    lorentz_neighbor_lists(Z, max_k, tol);

  NumericVector bw(B), score(B), trS(B), RSS(B), CV(B);
  NumericVector kappa_median(B), kappa_max(B), neff_min(B), neff_median(B);
  NumericVector n_positive_min(B), n_positive_median(B);

  for (int b = 0; b < B; ++b) {
    const int k = candidates[b];
    FitSummary summary = fit_summary_knn_cpp(
      y,
      X,
      neighbor_lists,
      k,
      kernel,
      approach,
      n_threads
    );

    bw[b] = k;
    score[b] = summary.score;
    trS[b] = summary.trS;
    RSS[b] = summary.RSS;
    CV[b] = summary.CV;
    kappa_median[b] = summary.kappa_median;
    kappa_max[b] = summary.kappa_max;
    neff_min[b] = summary.neff_min;
    neff_median[b] = summary.neff_median;
    n_positive_min[b] = summary.n_positive_min;
    n_positive_median[b] = summary.n_positive_median;
  }

  return DataFrame::create(
    _["bw"] = bw,
    _["score"] = score,
    _["trS"] = trS,
    _["RSS"] = RSS,
    _["CV"] = CV,
    _["kappa_median"] = kappa_median,
    _["kappa_max"] = kappa_max,
    _["neff_min"] = neff_min,
    _["neff_median"] = neff_median,
    _["n_positive_min"] = n_positive_min,
    _["n_positive_median"] = n_positive_median
  );
  END_RCPP
}

double evaluate_knn_bw_candidate(const arma::vec& y,
                                 const arma::mat& X,
                                 const std::vector<std::vector<std::pair<double, int> > >& neighbor_lists,
                                 double candidate,
                                 const std::string& approach,
                                 const std::string& kernel,
                                 GoldenTrace& trace,
                                 int n_threads) {
  candidate = std::round(candidate);
  const int existing = find_candidate(trace.bw, candidate);
  if (existing >= 0) return trace.score[existing];

  const int k = static_cast<int>(candidate);
  FitSummary summary = fit_summary_knn_cpp(y, X, neighbor_lists, k, kernel, approach, n_threads);

  trace.bw.push_back(candidate);
  trace.score.push_back(summary.score);
  trace.trS.push_back(summary.trS);
  trace.RSS.push_back(summary.RSS);
  trace.CV.push_back(summary.CV);
  trace.kappa_median.push_back(summary.kappa_median);
  trace.kappa_max.push_back(summary.kappa_max);
  trace.neff_min.push_back(summary.neff_min);
  trace.neff_median.push_back(summary.neff_median);
  trace.n_positive_min.push_back(summary.n_positive_min);
  trace.n_positive_median.push_back(summary.n_positive_median);
  trace.iteration.push_back(static_cast<int>(trace.bw.size()));

  return summary.score;
}

extern "C" SEXP _LHWR_bw_hgwr_knn_golden_cpp(SEXP ySEXP,
                                                     SEXP XSEXP,
                                                     SEXP ZSEXP,
                                                     SEXP approachSEXP,
                                                     SEXP kernelSEXP,
                                                     SEXP lowerSEXP,
                                                     SEXP upperSEXP,
                                                     SEXP tolSEXP,
                                                     SEXP maxIterSEXP,
                                                     SEXP distTolSEXP,
                                                     SEXP nThreadsSEXP) {
  BEGIN_RCPP
  arma::vec y = as<arma::vec>(ySEXP);
  arma::mat X = as<arma::mat>(XSEXP);
  arma::mat Z = as<arma::mat>(ZSEXP);
  const std::string approach = as<std::string>(approachSEXP);
  const std::string kernel = as<std::string>(kernelSEXP);
  double lower = as<double>(lowerSEXP);
  double upper = as<double>(upperSEXP);
  const double tol = as<double>(tolSEXP);
  const int max_iter = as<int>(maxIterSEXP);
  const double dist_tol = as<double>(distTolSEXP);
  const int n_threads = as<int>(nThreadsSEXP);

  const arma::uword n = y.n_elem;
  if (X.n_rows != n || Z.n_rows != n || Z.n_cols != 3) {
    stop("Dimensions of `y`, `X`, and `Z` are not compatible.");
  }
  if (!R_finite(lower) || !R_finite(upper) || lower >= upper) {
    stop("Invalid bandwidth search interval.");
  }
  if (!R_finite(tol) || tol <= 0.0) {
    stop("`tol` must be a positive finite number.");
  }
  if (max_iter < 1) {
    stop("`max_iter` must be positive.");
  }
  if (!R_finite(dist_tol) || dist_tol < 0.0) {
    stop("`dist_tol` must be a non-negative finite number.");
  }

  lower = std::round(lower);
  upper = std::round(upper);
  if (lower < 2.0) lower = 2.0;
  if (upper > static_cast<double>(n - 1)) upper = static_cast<double>(n - 1);
  if (lower >= upper) {
    stop("Invalid adaptive kNN bandwidth search interval.");
  }

  std::vector<std::vector<std::pair<double, int> > > neighbor_lists =
    lorentz_neighbor_lists(Z, static_cast<int>(upper), dist_tol);

  const double ratio = (std::sqrt(5.0) - 1.0) / 2.0;
  GoldenTrace trace;

  double xL = lower;
  double xU = upper;
  double d = ratio * (xU - xL);
  double x1 = normalize_bw_candidate(xL + d, lower, upper, true);
  double x2 = normalize_bw_candidate(xU - d, lower, upper, true);
  double f1 = evaluate_knn_bw_candidate(y, X, neighbor_lists, x1, approach, kernel, trace, n_threads);
  double f2 = evaluate_knn_bw_candidate(y, X, neighbor_lists, x2, approach, kernel, trace, n_threads);
  double diff_score = f2 - f1;
  int iter = 1;

  while (std::abs(d) > tol && std::abs(diff_score) > tol && iter < max_iter) {
    d = ratio * d;
    if (f1 < f2) {
      xL = x2;
      x2 = x1;
      f2 = f1;
      x1 = normalize_bw_candidate(xL + d, lower, upper, true);
      f1 = evaluate_knn_bw_candidate(y, X, neighbor_lists, x1, approach, kernel, trace, n_threads);
    } else {
      xU = x1;
      x1 = x2;
      f1 = f2;
      x2 = normalize_bw_candidate(xU - d, lower, upper, true);
      f2 = evaluate_knn_bw_candidate(y, X, neighbor_lists, x2, approach, kernel, trace, n_threads);
    }
    diff_score = f2 - f1;
    iter++;
  }

  DataFrame out = DataFrame::create(
    _["bw"] = wrap(trace.bw),
    _["score"] = wrap(trace.score),
    _["trS"] = wrap(trace.trS),
    _["RSS"] = wrap(trace.RSS),
    _["CV"] = wrap(trace.CV),
    _["kappa_median"] = wrap(trace.kappa_median),
    _["kappa_max"] = wrap(trace.kappa_max),
    _["neff_min"] = wrap(trace.neff_min),
    _["neff_median"] = wrap(trace.neff_median),
    _["n_positive_min"] = wrap(trace.n_positive_min),
    _["n_positive_median"] = wrap(trace.n_positive_median),
    _["iteration"] = wrap(trace.iteration)
  );
  out.attr("iterations") = iter;
  return out;
  END_RCPP
}

extern "C" SEXP _LHWR_bw_hgwr_golden_cpp(SEXP ySEXP,
                                                 SEXP XSEXP,
                                                 SEXP DSEXP,
                                                 SEXP approachSEXP,
                                                 SEXP kernelSEXP,
                                                 SEXP adaptiveSEXP,
                                                 SEXP bandwidthTypeSEXP,
                                                 SEXP pilotKSEXP,
                                                 SEXP lowerSEXP,
                                                 SEXP upperSEXP,
                                                 SEXP tolSEXP,
                                                 SEXP maxIterSEXP) {
  BEGIN_RCPP
  arma::vec y = as<arma::vec>(ySEXP);
  arma::mat X = as<arma::mat>(XSEXP);
  arma::mat D = as<arma::mat>(DSEXP);
  const std::string approach = as<std::string>(approachSEXP);
  const std::string kernel = as<std::string>(kernelSEXP);
  const bool adaptive = as<bool>(adaptiveSEXP);
  const std::string bandwidth_type = as<std::string>(bandwidthTypeSEXP);
  const int pilot_k = as<int>(pilotKSEXP);
  double lower = as<double>(lowerSEXP);
  double upper = as<double>(upperSEXP);
  const double tol = as<double>(tolSEXP);
  const int max_iter = as<int>(maxIterSEXP);
  
  const arma::uword n = y.n_elem;
  if (X.n_rows != n || D.n_rows != n || D.n_cols != n) {
    stop("Dimensions of `y`, `X`, and `D` are not compatible.");
  }
  if (!R_finite(lower) || !R_finite(upper) || lower >= upper) {
    stop("Invalid bandwidth search interval.");
  }
  if (!R_finite(tol) || tol <= 0.0) {
    stop("`tol` must be a positive finite number.");
  }
  if (max_iter < 1) {
    stop("`max_iter` must be positive.");
  }
  
  if (adaptive) {
    lower = std::round(lower);
    upper = std::round(upper);
  }
  
  const double ratio = (std::sqrt(5.0) - 1.0) / 2.0;
  GoldenTrace trace;
  
  double xL = lower;
  double xU = upper;
  double d = ratio * (xU - xL);
  double x1 = normalize_bw_candidate(xL + d, lower, upper, adaptive);
  double x2 = normalize_bw_candidate(xU - d, lower, upper, adaptive);
  double f1 = evaluate_bw_candidate(y, X, D, x1, approach, kernel, adaptive,
                                    bandwidth_type, pilot_k, trace);
  double f2 = evaluate_bw_candidate(y, X, D, x2, approach, kernel, adaptive,
                                    bandwidth_type, pilot_k, trace);
  double diff_score = f2 - f1;
  int iter = 1;
  
  while (std::abs(d) > tol && std::abs(diff_score) > tol && iter < max_iter) {
    d = ratio * d;
    if (f1 < f2) {
      xL = x2;
      x2 = x1;
      f2 = f1;
      x1 = normalize_bw_candidate(xL + d, lower, upper, adaptive);
      f1 = evaluate_bw_candidate(y, X, D, x1, approach, kernel, adaptive,
                                  bandwidth_type, pilot_k, trace);
    } else {
      xU = x1;
      x1 = x2;
      f1 = f2;
      x2 = normalize_bw_candidate(xU - d, lower, upper, adaptive);
      f2 = evaluate_bw_candidate(y, X, D, x2, approach, kernel, adaptive,
                                  bandwidth_type, pilot_k, trace);
    }
    diff_score = f2 - f1;
    iter++;
  }
  
  DataFrame out = DataFrame::create(
    _["bw"] = wrap(trace.bw),
    _["score"] = wrap(trace.score),
    _["trS"] = wrap(trace.trS),
    _["RSS"] = wrap(trace.RSS),
    _["CV"] = wrap(trace.CV),
    _["kappa_median"] = wrap(trace.kappa_median),
    _["kappa_max"] = wrap(trace.kappa_max),
    _["neff_min"] = wrap(trace.neff_min),
    _["neff_median"] = wrap(trace.neff_median),
    _["n_positive_min"] = wrap(trace.n_positive_min),
    _["n_positive_median"] = wrap(trace.n_positive_median),
    _["iteration"] = wrap(trace.iteration)
  );
  out.attr("iterations") = iter;
  return out;
  END_RCPP
}

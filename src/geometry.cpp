#include <RcppArmadillo.h>
#include <algorithm>
#include <queue>
#include <vector>

using namespace Rcpp;

namespace {

double lorentz_inner_one(const NumericMatrix& z, int i,
                         const NumericMatrix& w, int j) {
  return -z(i, 0) * w(j, 0) + z(i, 1) * w(j, 1) + z(i, 2) * w(j, 2);
}

double poincare_distance_one(const NumericMatrix& z, int i,
                             const NumericMatrix& w, int j,
                             double tol) {
  const double uz0 = z(i, 1) / (z(i, 0) + 1.0);
  const double uz1 = z(i, 2) / (z(i, 0) + 1.0);
  const double vw0 = w(j, 1) / (w(j, 0) + 1.0);
  const double vw1 = w(j, 2) / (w(j, 0) + 1.0);
  const double du0 = uz0 - vw0;
  const double du1 = uz1 - vw1;
  const double du2 = du0 * du0 + du1 * du1;
  const double nu = uz0 * uz0 + uz1 * uz1;
  const double nv = vw0 * vw0 + vw1 * vw1;
  const double den = (1.0 - nu) * (1.0 - nv);
  
  if (!R_finite(den) || den <= 0.0) return NA_REAL;
  double arg = 1.0 + 2.0 * du2 / den;
  if (arg < 1.0 && arg > 1.0 - std::sqrt(tol)) arg = 1.0;
  if (!R_finite(arg) || arg < 1.0) return NA_REAL;
  return std::acosh(arg);
}

double hyperbolic_distance_one(const NumericMatrix& z, int i,
                               const NumericMatrix& w, int j,
                               double tol) {
  double arg = -lorentz_inner_one(z, i, w, j);
  if (arg >= 1.0) return std::acosh(arg);
  if (arg > 1.0 - std::sqrt(tol)) return 0.0;
  return poincare_distance_one(z, i, w, j, tol);
}

void check_three_columns(const NumericMatrix& x, const char* name) {
  if (x.ncol() != 3) {
    stop("`%s` must have exactly 3 columns.", name);
  }
}

double lorentz_norm_vec(const double x0, const double x1, const double x2) {
  return -x0 * x0 + x1 * x1 + x2 * x2;
}

struct CandidateLess {
  bool operator()(const std::pair<double, int>& a,
                  const std::pair<double, int>& b) const {
    if (a.first < b.first) return true;
    if (a.first > b.first) return false;
    return a.second < b.second;
  }
};

bool candidate_better(const std::pair<double, int>& candidate,
                      const std::pair<double, int>& current_worst) {
  if (candidate.first < current_worst.first) return true;
  if (candidate.first > current_worst.first) return false;
  return candidate.second < current_worst.second;
}

void update_heap(std::priority_queue<
                   std::pair<double, int>,
                   std::vector<std::pair<double, int> >,
                   CandidateLess
                 >& heap,
                 const std::pair<double, int>& candidate,
                 int k) {
  if (static_cast<int>(heap.size()) < k) {
    heap.push(candidate);
    return;
  }
  if (candidate_better(candidate, heap.top())) {
    heap.pop();
    heap.push(candidate);
  }
}

} // namespace

extern "C" SEXP _LHWR_lorentz_inner_cpp(SEXP zSEXP, SEXP wSEXP) {
  BEGIN_RCPP
  NumericMatrix z(zSEXP);
  NumericMatrix w(wSEXP);
  check_three_columns(z, "z");
  check_three_columns(w, "w");
  
  const int nz = z.nrow();
  const int nw = w.nrow();
  if (nz != nw && nz != 1 && nw != 1) {
    stop("`z` and `w` must have compatible row counts.");
  }
  
  const int n = std::max(nz, nw);
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    const int iz = nz == 1 ? 0 : i;
    const int iw = nw == 1 ? 0 : i;
    out[i] = lorentz_inner_one(z, iz, w, iw);
  }
  return out;
  END_RCPP
}

extern "C" SEXP _LHWR_hyperbolic_distance_matrix_cpp(SEXP ZSEXP,
                                                             SEXP WSEXP,
                                                             SEXP symmetricSEXP,
                                                             SEXP tolSEXP) {
  BEGIN_RCPP
  NumericMatrix Z(ZSEXP);
  NumericMatrix W(WSEXP);
  const bool symmetric = as<bool>(symmetricSEXP);
  const double tol = as<double>(tolSEXP);
  check_three_columns(Z, "Z");
  check_three_columns(W, "W");
  
  const int n = Z.nrow();
  const int m = W.nrow();
  NumericMatrix D(n, m);
  
  if (symmetric) {
    if (n != m) stop("Symmetric distance matrix requires equal row counts.");
    for (int i = 0; i < n; ++i) {
      D(i, i) = 0.0;
      for (int j = i + 1; j < n; ++j) {
        const double dij = hyperbolic_distance_one(Z, i, Z, j, tol);
        if (!R_finite(dij)) {
          stop("Invalid Lorentz distance: acosh argument is below 1.");
        }
        D(i, j) = dij;
        D(j, i) = dij;
      }
    }
  } else {
    for (int j = 0; j < m; ++j) {
      for (int i = 0; i < n; ++i) {
        const double dij = hyperbolic_distance_one(Z, i, W, j, tol);
        if (!R_finite(dij)) {
          stop("Invalid Lorentz distance: acosh argument is below 1.");
        }
        D(i, j) = dij;
      }
    }
  }
  
  return D;
  END_RCPP
}

extern "C" SEXP _LHWR_lorentz_knn_cpp(SEXP ZSEXP,
                                              SEXP kSEXP,
                                              SEXP includeSelfSEXP,
                                              SEXP tolSEXP) {
  BEGIN_RCPP
  NumericMatrix Z(ZSEXP);
  const int k = as<int>(kSEXP);
  const bool include_self = as<bool>(includeSelfSEXP);
  const double tol = as<double>(tolSEXP);
  check_three_columns(Z, "Z");
  
  const int n = Z.nrow();
  const int max_k = include_self ? n : n - 1;
  if (k < 1 || k > max_k) {
    stop("`k` must be between 1 and the number of available neighbors.");
  }
  
  IntegerMatrix index(n, k);
  NumericMatrix distance(n, k);
  typedef std::priority_queue<
    std::pair<double, int>,
    std::vector<std::pair<double, int> >,
    CandidateLess
  > NeighborHeap;
  std::vector<NeighborHeap> heaps(n);
  
  for (int i = 0; i < n; ++i) {
    if (include_self) {
      heaps[i].push(std::make_pair(0.0, i));
    }
  }
  
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const double dij = hyperbolic_distance_one(Z, i, Z, j, tol);
      if (!R_finite(dij)) {
        stop("Invalid Lorentz distance: acosh argument is below 1.");
      }
      update_heap(heaps[i], std::make_pair(dij, j), k);
      update_heap(heaps[j], std::make_pair(dij, i), k);
    }
  }
  
  for (int i = 0; i < n; ++i) {
    if (static_cast<int>(heaps[i].size()) != k) {
      stop("Could not compute the requested number of neighbors.");
    }
    std::vector<std::pair<double, int> > neighbors;
    neighbors.reserve(k);
    while (!heaps[i].empty()) {
      neighbors.push_back(heaps[i].top());
      heaps[i].pop();
    }
    std::sort(
      neighbors.begin(),
      neighbors.end(),
      [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
        if (a.first < b.first) return true;
        if (a.first > b.first) return false;
        return a.second < b.second;
      }
    );
    
    for (int ell = 0; ell < k; ++ell) {
      index(i, ell) = neighbors[ell].second + 1;
      distance(i, ell) = neighbors[ell].first;
    }
  }
  
  return List::create(
    _["index"] = index,
    _["distance"] = distance,
    _["k"] = k,
    _["include_self"] = include_self
  );
  END_RCPP
}

extern "C" SEXP _LHWR_tangent_distance_matrix_cpp(SEXP ZSEXP,
                                                          SEXP muSEXP,
                                                          SEXP tolSEXP) {
  BEGIN_RCPP
  NumericMatrix Z(ZSEXP);
  NumericMatrix muMat(muSEXP);
  const double tol = as<double>(tolSEXP);
  check_three_columns(Z, "Z");
  check_three_columns(muMat, "mu");
  if (muMat.nrow() != 1) {
    stop("`mu` must be a single Lorentz point.");
  }
  
  const int n = Z.nrow();
  NumericMatrix V(n, 3);
  NumericMatrix D(n, n);
  
  const double mu0 = muMat(0, 0);
  const double mu1 = muMat(0, 1);
  const double mu2 = muMat(0, 2);
  
  for (int i = 0; i < n; ++i) {
    const double ip = -Z(i, 0) * mu0 + Z(i, 1) * mu1 + Z(i, 2) * mu2;
    double arg = -ip;
    if (arg < 1.0 && arg > 1.0 - std::sqrt(tol)) arg = 1.0;
    if (!R_finite(arg) || arg < 1.0) {
      stop("Invalid Lorentz log-map argument.");
    }
    
    const double d = std::acosh(arg);
    if (d <= tol) {
      V(i, 0) = 0.0;
      V(i, 1) = 0.0;
      V(i, 2) = 0.0;
    } else {
      const double coef = d / std::sinh(d);
      const double cd = std::cosh(d);
      V(i, 0) = coef * (Z(i, 0) - cd * mu0);
      V(i, 1) = coef * (Z(i, 1) - cd * mu1);
      V(i, 2) = coef * (Z(i, 2) - cd * mu2);
    }
  }
  
  for (int i = 0; i < n; ++i) {
    D(i, i) = 0.0;
    for (int j = i + 1; j < n; ++j) {
      const double d0 = V(i, 0) - V(j, 0);
      const double d1 = V(i, 1) - V(j, 1);
      const double d2 = V(i, 2) - V(j, 2);
      double norm = lorentz_norm_vec(d0, d1, d2);
      if (norm < 0.0 && norm > -std::sqrt(tol)) norm = 0.0;
      if (!R_finite(norm) || norm < 0.0) {
        stop("Invalid tangent-space distance under the induced Lorentz metric.");
      }
      const double dij = std::sqrt(norm);
      D(i, j) = dij;
      D(j, i) = dij;
    }
  }
  
  return D;
  END_RCPP
}

#include <RcppArmadillo.h>
#include <algorithm>
#include <queue>
#include <vector>

#include "lorentz_distance.h"

using namespace Rcpp;

namespace {

double lorentz_inner_one(const NumericMatrix& z, int i,
                         const NumericMatrix& w, int j) {
  return -z(i, 0) * w(j, 0) + z(i, 1) * w(j, 1) + z(i, 2) * w(j, 2);
}

double hyperbolic_distance_one(const NumericMatrix& z, int i,
                               const NumericMatrix& w, int j,
                               double tol) {
  static_cast<void>(tol);
  return lhwl_lorentz_distance(
    z(i, 0), z(i, 1), z(i, 2),
    w(j, 0), w(j, 1), w(j, 2)
  );
}

void check_three_columns(const NumericMatrix& x, const char* name) {
  if (x.ncol() != 3) {
    stop("`%s` must have exactly 3 columns.", name);
  }
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

extern "C" SEXP _LHWL_lorentz_inner_cpp(SEXP zSEXP, SEXP wSEXP) {
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

extern "C" SEXP _LHWL_hyperbolic_distance_cpp(SEXP zSEXP,
                                                       SEXP wSEXP,
                                                       SEXP tolSEXP) {
  BEGIN_RCPP
  NumericMatrix z(zSEXP);
  NumericMatrix w(wSEXP);
  const double tol = as<double>(tolSEXP);
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
    out[i] = hyperbolic_distance_one(z, iz, w, iw, tol);
    if (!R_finite(out[i])) {
      stop("Invalid Lorentz-model coordinates for distance calculation.");
    }
  }
  return out;
  END_RCPP
}

extern "C" SEXP _LHWL_hyperbolic_distance_matrix_cpp(SEXP ZSEXP,
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

extern "C" SEXP _LHWL_lorentz_knn_cpp(SEXP ZSEXP,
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

extern "C" SEXP _LHWL_tangent_distance_matrix_cpp(SEXP ZSEXP,
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
  NumericMatrix Q(n, 2);
  NumericMatrix D(n, n);
  
  const double mu0 = muMat(0, 0);
  const double mu1 = muMat(0, 1);
  const double mu2 = muMat(0, 2);
  const double muSpatial = std::hypot(mu1, mu2);
  const double muRadial = std::asinh(muSpatial);
  const double muTime = std::hypot(1.0, muSpatial);
  const double muDirection1 = muSpatial > 0.0 ? mu1 / muSpatial : 1.0;
  const double muDirection2 = muSpatial > 0.0 ? mu2 / muSpatial : 0.0;
  
  for (int i = 0; i < n; ++i) {
    const double d = lhwl_lorentz_distance(
      Z(i, 0), Z(i, 1), Z(i, 2), mu0, mu1, mu2
    );
    if (!R_finite(d)) {
      stop("Invalid Lorentz log-map argument.");
    }

    if (d <= tol) {
      Q(i, 0) = 0.0;
      Q(i, 1) = 0.0;
    } else {
      const double spatial = std::hypot(Z(i, 1), Z(i, 2));
      const double radial = std::asinh(spatial);
      const double direction1 = spatial > 0.0 ? Z(i, 1) / spatial : 0.0;
      const double direction2 = spatial > 0.0 ? Z(i, 2) / spatial : 0.0;
      double radialComponent;
      double angularComponent;
      if (muSpatial > 0.0) {
        const double du1 = direction1 - muDirection1;
        const double du2 = direction2 - muDirection2;
        const double directionGap = 0.5 * (du1 * du1 + du2 * du2);
        radialComponent = std::sinh(radial - muRadial) -
          spatial * muTime * directionGap;
        angularComponent = -Z(i, 1) * muDirection2 +
          Z(i, 2) * muDirection1;
      } else {
        radialComponent = Z(i, 1);
        angularComponent = Z(i, 2);
      }
      const double coef = d < 350.0 ?
        d / std::sinh(d) :
        2.0 * d * std::exp(-d) / (1.0 - std::exp(-2.0 * d));
      Q(i, 0) = coef * radialComponent;
      Q(i, 1) = coef * angularComponent;
    }
  }
  
  for (int i = 0; i < n; ++i) {
    D(i, i) = 0.0;
    for (int j = i + 1; j < n; ++j) {
      const double d1 = Q(i, 0) - Q(j, 0);
      const double d2 = Q(i, 1) - Q(j, 1);
      const double dij = std::hypot(d1, d2);
      D(i, j) = dij;
      D(j, i) = dij;
    }
  }
  
  return D;
  END_RCPP
}

# LHWPR Solver Controls

Construct the numerical controls passed to
[`lhwpr()`](https://byuzbasi.github.io/LHWL/reference/lhwpr.md) and
[`tune.lhwpr()`](https://byuzbasi.github.io/LHWL/reference/tune.lhwpr.md).

## Usage

``` r
lhwpr_control(
  solver = c("auto", "naive", "covariance"),
  screening = c("strong", "none"),
  tol = 1e-07,
  kkt_tol = 1e-07,
  max_iter = 100000L,
  irls_tol = 1e-07,
  irls_max_iter = 100L,
  max_active = NULL,
  early_stop = FALSE,
  ridge_path_tol = 1e-06,
  n_threads = 1L,
  trace = FALSE
)
```

## Arguments

- solver:

  Automatic, residual-updating, or Gaussian covariance-updating
  coordinate descent.

- screening:

  Sequential strong screening with full KKT repair, or no screening.

- tol:

  Coordinate-descent tolerance.

- kkt_tol:

  Full KKT or stationarity tolerance.

- max_iter:

  Maximum coordinate sweeps per weighted least-squares solve.

- irls_tol:

  Binomial outer-loop tolerance.

- irls_max_iter:

  Maximum binomial outer iterations.

- max_active:

  Optional maximum active-set size.

- early_stop:

  Reserved for a later globally synchronized path solver. It must remain
  `FALSE` in this release.

- ridge_path_tol:

  Tolerance defining the ridge path's upper scale.

- n_threads:

  Number of focal-location threads.

- trace:

  Whether numerical tracing is requested.

## Value

A validated object of class `"lhwpr_control"`.

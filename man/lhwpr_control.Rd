\name{lhwpr_control}
\alias{lhwpr_control}
\title{LHWPR Solver Controls}
\usage{
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
}
\arguments{
\item{solver}{Automatic, residual-updating, or Gaussian covariance-updating
coordinate descent.}
\item{screening}{Sequential strong screening with full KKT repair, or no
screening.}
\item{tol}{Coordinate-descent tolerance.}
\item{kkt_tol}{Full KKT or stationarity tolerance.}
\item{max_iter}{Maximum coordinate sweeps per weighted least-squares solve.}
\item{irls_tol}{Binomial outer-loop tolerance.}
\item{irls_max_iter}{Maximum binomial outer iterations.}
\item{max_active}{Optional maximum active-set size.}
\item{early_stop}{Reserved for a later globally synchronized path solver.
It must remain \code{FALSE} in this release.}
\item{ridge_path_tol}{Tolerance defining the ridge path's upper scale.}
\item{n_threads}{Number of focal-location threads.}
\item{trace}{Whether numerical tracing is requested.}
}
\value{A validated object of class \code{"lhwpr_control"}.}
\description{
Construct the numerical controls passed to \code{lhwpr()} and
\code{tune.lhwpr()}.
}

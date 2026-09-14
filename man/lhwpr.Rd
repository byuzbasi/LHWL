\name{lhwpr}
\alias{lhwpr}
\alias{coef.lhwpr}
\alias{fitted.lhwpr}
\alias{plot.lhwpr}
\alias{predict.lhwpr}
\alias{print.lhwpr}
\alias{print.summary.lhwpr}
\alias{residuals.lhwpr}
\alias{summary.lhwpr}
\title{Lorentz-Hyperbolic Weighted Penalized Regression}
\description{
Fit local Gaussian or binomial penalized regressions whose observation
weights are defined by Lorentz-hyperbolic neighborhoods.
}
\usage{
lhwpr(
  formula,
  data,
  Z,
  bw,
  family = c("gaussian", "binomial"),
  penalty = c("lasso", "mcp", "scad"),
  alpha = 1,
  lambda = NULL,
  nlambda = 100L,
  lambda_min_ratio = NULL,
  gamma = NULL,
  penalty_factor = NULL,
  case_weights = NULL,
  standardize = TRUE,
  kernel = hgwr_kernel_names(),
  adaptive = TRUE,
  bandwidth_type = c("knn", "volume"),
  geometry = c("lorentz", "tangent"),
  tangent_ref = NULL,
  dMat = NULL,
  engine = c("dense", "knn"),
  keep_path = c("auto", "all", "selected"),
  control = lhwpr_control()
)

\method{print}{lhwpr}(x, ...)

\method{summary}{lhwpr}(object, ...)

\method{print}{summary.lhwpr}(x, ...)

\method{coef}{lhwpr}(object, s = NULL, ...)

\method{fitted}{lhwpr}(object, ...)

\method{residuals}{lhwpr}(object, ...)

\method{predict}{lhwpr}(
  object,
  newdata = NULL,
  Znew = NULL,
  type = c("response", "link", "coefficients"),
  s = NULL,
  ...
)

\method{plot}{lhwpr}(
  x,
  focus = 1L,
  xvar = c("lambda", "lambda_ratio"),
  ...
)
}
\arguments{
\item{formula}{A model formula with one response and one or more predictors.}
\item{data}{A data frame containing variables in \code{formula}.}
\item{Z}{An \code{n} by 3 matrix of Lorentz-hyperboloid coordinates.}
\item{bw}{Kernel bandwidth, or neighbor count for adaptive kNN fitting.}
\item{family}{Gaussian or binomial response family.}
\item{penalty}{Structural penalty: lasso, MCP, or SCAD. It is ignored when
\code{alpha = 0}.}
\item{alpha}{Mixing value in \code{[0, 1]}. Zero is ridge; a value strictly
between zero and one is elastic net and requires the lasso penalty; one gives
lasso, MCP, or SCAD.}
\item{lambda}{Optional strictly decreasing common absolute lambda path.}
\item{nlambda}{Number of automatically generated lambda values.}
\item{lambda_min_ratio}{Ratio of the final to initial automatic lambda.}
\item{gamma}{MCP or SCAD concavity parameter.}
\item{penalty_factor}{Optional non-negative multiplier for each non-intercept
coefficient.}
\item{case_weights}{Optional non-negative case weights.}
\item{standardize}{Whether predictors are centered and scaled once on the full
analysis sample.}
\item{kernel}{Lorentz-distance kernel.}
\item{adaptive}{Whether \code{bw} is adaptive.}
\item{bandwidth_type}{Neighbor-count or hyperbolic-volume bandwidth scale.}
\item{geometry}{Intrinsic Lorentz or tangent-space distance geometry.}
\item{tangent_ref}{Optional Lorentz point defining the tangent plane.}
\item{dMat}{Optional precomputed distance matrix.}
\item{engine}{Dense or exact compact-kernel kNN computation.}
\item{keep_path}{Whether to retain all local coefficient paths.}
\item{control}{A control object from \code{lhwpr_control()}.}
\item{x, object}{A fitted \code{lhwpr} object, or its summary.}
\item{s}{A lambda value exactly represented in the fitted path.}
\item{newdata}{Optional predictor data for new focal locations.}
\item{Znew}{Lorentz coordinates corresponding to \code{newdata}.}
\item{type}{Prediction scale or local coefficients.}
\item{focus}{One-based focal-location index for a path plot.}
\item{xvar}{Horizontal path scale.}
\item{...}{Additional arguments passed to the plotting method.}
}
\details{
At every focal location, kernel and case weights are multiplied and normalized
before fitting. Predictor standardization is global rather than focal. For a
given bandwidth, every focal model uses the same absolute lambda values.

Lasso and elastic-net paths use warm starts, sequential strong screening, and
a full KKT repair pass. MCP and SCAD use curvature-aware scalar minimization
and full stationarity scans. Binomial models use proximal Newton/IRLS.

The intercept is never penalized. Group penalties are not part of LHWPR.
Classical unpenalized LHWR inference is not attached to penalized fits.

With \code{keep_path = "auto"}, a full dense coefficient cube is retained only
when it contains at most ten million values. Selected coefficients and all
fitted-value paths are retained regardless. Use \code{keep_path = "all"} to
override this memory guard.
}
\value{
\code{lhwpr()} returns an object of class \code{"lhwpr"} with selected local
coefficients and predictions, lambda metadata, optional paths, convergence and
KKT diagnostics, preprocessing metadata, and geometry information. Extractor
methods return the corresponding component.
}
\examples{
set.seed(1)
H <- simulate_H2(40, r_max = 2)
dat <- data.frame(
  y = rnorm(40),
  x1 = rnorm(40),
  x2 = rnorm(40)
)
fit <- lhwpr(
  y ~ x1 + x2,
  data = dat,
  Z = H$Z,
  bw = 15,
  alpha = 0.5,
  nlambda = 10
)
coef(fit)
}

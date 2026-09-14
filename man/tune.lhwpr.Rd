\name{tune.lhwpr}
\alias{tune.lhwpr}
\alias{print.tune.lhwpr}
\title{Joint Bandwidth and Regularization Tuning for LHWPR}
\description{
Select LHWPR bandwidth and regularization parameters with cross-validation or
Gaussian information criteria.
}
\usage{
tune.lhwpr(
  formula,
  data,
  Z,
  bw_grid = NULL,
  family = c("gaussian", "binomial"),
  penalty = c("lasso", "mcp", "scad"),
  alpha_grid = 1,
  gamma_grid = NULL,
  lambda = NULL,
  nlambda = 100L,
  lambda_min_ratio = NULL,
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
  resampling = c("kfold", "loo", "blocked"),
  nfolds = 5L,
  fold_id = NULL,
  rule = c("min", "1se", "sparse_1se"),
  control = lhwpr_control(),
  refit = TRUE,
  bw_lower = NULL,
  bw_upper = NULL,
  n_bw = 15L,
  criterion = c(
    "cv", "aicc", "aic", "bic", "gcv", "ebic", "ebic_refit"
  ),
  ebic_gamma = 0.5,
  refine = NULL,
  n_refine = 7L,
  refine_alpha = FALSE,
  expand_bounds = TRUE,
  compare_lhwr = TRUE,
  eligibility_kkt_tol = 1e-05
)

\method{print}{tune.lhwpr}(x, ...)
}
\arguments{
\item{formula, data, Z}{Model formula, analysis data, and Lorentz coordinates.}
\item{bw_grid}{Positive candidate bandwidths. When \code{NULL}, an adaptive
kNN grid is generated automatically.}
\item{bw_lower, bw_upper}{Optional lower and upper adaptive-kNN bounds used
only when \code{bw_grid = NULL}.}
\item{n_bw}{Number of automatically generated adaptive-kNN candidates.}
\item{family, penalty}{Response family and structural penalty.}
\item{alpha_grid}{Candidate alpha values.}
\item{gamma_grid}{Optional MCP or SCAD gamma values.}
\item{lambda}{Optional common decreasing absolute lambda path.}
\item{nlambda, lambda_min_ratio}{Automatic lambda-path controls.}
\item{penalty_factor, case_weights, standardize}{Penalized preprocessing
controls; see \code{lhwpr()}.}
\item{kernel, adaptive, bandwidth_type, geometry, tangent_ref, dMat, engine}{
Locality and geometry controls; see \code{lhwpr()}.}
\item{resampling}{Leave-one-out, K-fold, or blocked validation.}
\item{nfolds}{Number of automatically constructed deterministic K-folds.}
\item{fold_id}{Optional fold identifier per observation; required for blocked
validation. Without this argument, K-fold allocation cycles through input row
order deterministically; it is not randomly shuffled. Supply a frozen random
or scientifically defined blocked assignment when that is intended.}
\item{rule}{Minimum-score, one-standard-error, or
sparse-one-standard-error selection. The latter two require
\code{criterion = "cv"}.}
\item{control}{A value returned by \code{lhwpr_control()}.}
\item{refit}{Whether to fit the selected model on all observations.}
\item{criterion}{Selection score. Cross-validation is the default. Gaussian
models additionally support AIC, AICc, BIC, GCV, path-based EBIC, and the
refit-based \code{"ebic_refit"} sensitivity criterion.}
\item{ebic_gamma}{EBIC model-space penalty in \eqn{[0,1]}.}
\item{refine}{Whether to refine the bandwidth grid around the provisional
optimum. \code{NULL} refines automatic grids and leaves supplied grids
unchanged.}
\item{n_refine}{Number of candidates in the local bandwidth refinement.}
\item{refine_alpha}{Whether to evaluate alpha midpoints adjacent to the
provisionally selected alpha.}
\item{expand_bounds}{Whether an automatic lower search boundary may expand
toward its model-dimension safety limit. User-supplied bounds are respected.}
\item{compare_lhwr}{Whether to evaluate an unpenalized LHWR reference using
the same folds and bandwidths.}

\item{eligibility_kkt_tol}{Maximum recorded KKT violation allowed for an
eligible tuning row. All relevant focal fits must also converge; CV requires
both validation and full-data path fits to pass.}
\item{x}{A \code{tune.lhwpr} object.}
\item{...}{Unused additional arguments.}
}
\details{
Tuning is conditional on the supplied representation and prepared design.
Predictor standardization, the representation, and its tangent reference are
prepared once for the supplied sample, not relearned within each fold. To
evaluate an entire learned pipeline after tuning, use an outer validation
workflow with transformations and representations fitted on the outer training
sample and applied unchanged to validation observations.

Within each bandwidth and alpha configuration, all validation fits use one
common absolute lambda path. For K-fold and blocked validation, its upper end
is the maximum KKT threshold across folds. The \code{"1se"} rule chooses the
most regularized eligible model, breaking ties toward a larger bandwidth.
The \code{"sparse_1se"} rule first minimizes the number of selected
coefficient surfaces within the same one-standard-error set.

All selection rules use the same numerical eligibility mask. A candidate must
have complete convergence and a maximum recorded KKT violation no larger than
\code{eligibility_kkt_tol}; Gaussian criteria retain their additional effective
degrees-of-freedom restrictions. The final fit reuses the exact decreasing
lambda prefix evaluated for the selected configuration instead of solving the
selected lambda in isolation. The returned \code{selected_lambda_path} and
\code{refit_diagnostics} record this path and its convergence/KKT certificate.

When \code{bw_grid = NULL}, automatic selection is available for adaptive kNN
bandwidths. The lower bound is model-dimension aware and at least 20 neighbors;
the default upper bound is \eqn{n-1}. A cached coarse search is followed by a
local refinement. If a default lower boundary is provisionally selected, the
search may expand downward but never below the model-dimension safety limit.
The returned object records and plots boundary solutions. Supply an explicit
\code{bw_grid} for fixed-radius or volume bandwidths.

For Gaussian paths, ridge uses the exact smoother trace. Lasso and elastic net
use a conditional active-set smoother, while MCP and SCAD additionally use
their local piecewise penalty curvature. Here \eqn{n} is the number of
positive-case-weight observations. If \eqn{\nu} denotes the effective
parameter count including the Gaussian scale parameter and \eqn{L} is minus
twice the maximized Gaussian log likelihood, the reported criteria include
\deqn{\mathrm{AIC}=L+2\nu,}
\deqn{\mathrm{AICc}=\mathrm{AIC}+
2\nu(\nu+1)/(n-\nu-1),}
\deqn{\mathrm{BIC}=L+\log(n)\nu.}
GCV divides the weighted residual mean square by the squared remaining
smoother fraction, giving
\eqn{(\mathrm{RSS}/n)/(1-\mathrm{tr}(S)/n)^2}, where \eqn{S} is the
conditional fitted smoother and \eqn{\mathrm{RSS}} is the residual sum of
squares under the normalized case weights described below.
Gaussian path-information criteria can select only rows with
\code{edf_valid = TRUE}; refit EBIC additionally requires
\code{post_edf_valid = TRUE}.

The primary EBIC augments the pathwise BIC with a model-space term. In this
expression \eqn{\gamma} is \code{ebic_gamma}, not the MCP/SCAD concavity
parameter, and \eqn{p_{search}} counts predictors with positive penalty
factors. The active-surface count \eqn{q_{surface}} is defined below:
\deqn{\mathrm{EBIC}=\mathrm{BIC}_{path}+
2\gamma\log {p_{search} \choose q_{surface}}.}
The optional \code{criterion = "ebic_refit"} sensitivity score replaces
\eqn{\mathrm{BIC}_{path}} with the BIC from an unpenalized local refit on each
selected active set. Refit scores with at most one residual degree of freedom
are recorded as infinite and cannot be selected; \code{post_residual_df} and
\code{post_edf_valid} retain this diagnostic.
Here \eqn{q_{surface}} is the integer count of penalized coefficient surfaces
whose standardized coefficient magnitude exceeds \eqn{10^{-10}} at one or more
focal locations. Predictors with zero penalty factors contribute to smoother
complexity but not to the model-space term. The selection table stores its
logarithmic model-space term in \code{ebic_model_space_log} and flags a
nonempty full-surface model, for which this term is zero, in
\code{ebic_model_space_zero}. Case weights are normalized to sum
to the number of positive-weight observations when computing Gaussian
likelihood criteria. These are path-selection diagnostics, not post-selection
coefficient standard errors. Binomial tuning currently requires
\code{criterion = "cv"}.
}
\value{
A \code{"tune.lhwpr"} object containing the complete selection table,
information criteria, numerical eligibility and support summaries, fold
scores, search history, boundary diagnostics, the same-fold LHWR reference,
selected row and exact lambda prefix, fold identifiers, and optionally the
certified full-data refit.
}
\examples{
set.seed(2)
H <- simulate_H2(35, r_max = 2)
dat <- data.frame(y = rnorm(35), x = rnorm(35))
selection <- tune.lhwpr(
  y ~ x,
  data = dat,
  Z = H$Z,
  bw_lower = 12,
  bw_upper = 28,
  n_bw = 5,
  nlambda = 5,
  nfolds = 3,
  criterion = "cv",
  rule = "sparse_1se"
)
selection
}

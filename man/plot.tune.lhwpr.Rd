\name{plot.tune.lhwpr}
\alias{plot.tune.lhwpr}
\title{Plot an LHWPR Tuning Object}
\description{
Plot bandwidth and lambda profiles, a joint bandwidth--lambda surface, an
exact evaluated-candidate landscape, and an alpha profile for an object
returned by \code{tune.lhwpr()}.
}
\usage{
\method{plot}{tune.lhwpr}(
  x,
  type = c(
    "all", "bandwidth", "lambda", "landscape", "surface", "alpha"
  ),
  criterion = NULL,
  alpha = NULL,
  gamma = NULL,
  include_lhwr = TRUE,
  ...
)
}
\arguments{
\item{x}{A \code{"tune.lhwpr"} object.}
\item{type}{One of \code{"all"}, \code{"bandwidth"}, \code{"lambda"},
\code{"landscape"}, \code{"surface"}, or \code{"alpha"}.}
\item{criterion}{Criterion to display. \code{NULL} uses the criterion that
selected the model.}
\item{alpha, gamma}{Optional alpha and gamma defining the displayed surface
slice. With \code{type = "surface"} or \code{type = "landscape"} and
\code{alpha = NULL}, one panel is drawn for each evaluated alpha.}
\item{include_lhwr}{Whether the same-fold unpenalized LHWR reference is
overlaid on the bandwidth profile.}
\item{...}{Additional graphical arguments passed to profile plots.}
}
\details{
The selected tuning point is marked in red. A selected lower or upper
bandwidth boundary is annotated in the bandwidth and surface panels. The
surface uses relative lambda so configurations whose absolute lambda maxima
differ remain comparable.

The landscape plots only evaluated candidates and therefore does not imply
interpolation between grid points. Point shapes distinguish coarse and refined
stages, grey crosses mark candidates excluded by convergence, KKT, or
criterion-specific validity gates, and the selected eligible point is marked
separately. Profiles and interpolated surfaces are computed only from eligible
candidates.
}
\value{
Invisibly, the data used by the requested plot. For \code{type = "all"}, a
named list containing the four panel data sets is returned; its joint panel is
named \code{landscape}.
}
\examples{
set.seed(3)
H <- simulate_H2(30, r_max = 1.5)
dat <- data.frame(y = rnorm(30), x = rnorm(30))
selection <- tune.lhwpr(
  y ~ x,
  data = dat,
  Z = H$Z,
  bw_grid = c(14, 20),
  nlambda = 3,
  nfolds = 2,
  refit = FALSE
)
plot(selection, type = "bandwidth")
}

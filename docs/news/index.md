# Changelog

## LHWL 0.3.1

- Fixed Moran and LISA diagnostics with Gaussian or exponential kernels:
  adaptive kNN radii still set bandwidths, but stored neighbour lists no
  longer truncate noncompact kernel tails. These diagnostics now use
  full distance matrices. Added cached/full/explicit-weight parity
  tests, including permutation inference. Bisquare diagnostics and model
  estimation are unchanged.

- Added executable LHWR companions for kernels, bandwidth selection,
  F1–F3 and Moran/LISA diagnostics, and manuscript-style
  graphics/export. Documented supported engine combinations, raw versus
  adjusted classes, plotting transformations, and diagnostic-weight
  cache behavior. This documentation revision does not change an
  estimator or historical run.

- Expanded the introduction and Lasso vignettes into linked, executable
  tutorials with quick-start examples, interpreted graphics, numerical
  safeguards and explicit scope limits.

- Added local pkgdown navigation and a shared citation-maintenance
  workflow. Public manuscript identifiers and hosting remain pending; no
  publication or deployment is performed by the documentation tools.

- Synchronized the release candidate with the corrected Tangent-E
  prediction distances. Direct, scaled coordinate differences preserve
  exact zero distances for identical locations without rounding nearby
  distinct points to zero. Added calibration, duplicate, reordered-row,
  and batch tests.

- Made the near-coincident-point test’s absolute error bound explicit;
  this changes the test assertion only, not a distance or estimation
  tolerance.

- Expanded the geometry, reference-center, prediction,
  conditional-validation, tuning, and inference documentation without
  changing estimation settings.

- Matched the C++17 requirement in DESCRIPTION with explicit Unix and
  Windows Makevars declarations; retained portable R-provided
  native-library flags.

- Prepared the public source metadata, package citation, data
  attribution, and self-contained WDI quick start for GitHub and CRAN
  review.

- Recorded the frozen location and scale constants used to construct the
  packaged WDI Lorentz coordinates and clarified that the representation
  variables define locality rather than enter the regression formula.

- Hardened joint tuning, exact selected-model refits, and LHWPR
  graphical diagnostics.

- Required every selected tuning row to pass full-fit convergence and
  KKT gates; cross-validation candidates must also pass the same
  validation-fit gates, and the minimum, one-standard-error, and
  sparse-one-standard-error rules now share exactly the same eligibility
  set.

- Changed the final selected-model refit to reuse the exact decreasing
  lambda prefix evaluated for its tuning configuration. The fitted
  tuning object now records this prefix and a convergence/KKT refit
  certificate.

- Added the non-interpolated `plot.tune.lhwpr(type = "landscape")`
  display, which distinguishes coarse, refined, ineligible, and selected
  candidates.

- Added
  [`plot_lhwpr_surface()`](https://byuzbasi.github.io/LHWL/reference/plot_lhwpr_surface.md)
  for coefficient and focal-selection displays. Numerical zeros are
  hollow and the documentation explicitly distinguishes penalized
  activation from statistical significance.

## LHWL 0.3.0

- Declared the R-provided LAPACK, BLAS, and Fortran runtime libraries in
  the Unix and Windows native-link recipes, preventing unresolved
  Armadillo linear-algebra symbols on clean installations.
- Restored the local R-reference regression tests for the C++ bandwidth,
  fitting, and Lorentz-geometry implementations, and added a
  deterministic Gaussian LHWPR-Lasso test for the `p > n` setting.
- Added a non-overwriting, versioned review audit that records the R
  runtime, compiler/link configuration, user Makevars provenance,
  package versions, test/check logs, source hashes, and SHA-256 artifact
  inventory. The audit now stops before compilation when a required
  compiler or link-library path is unavailable.
- Redesigned
  [`tune.lhwpr()`](https://byuzbasi.github.io/LHWL/reference/tune.lhwpr.md)
  around a cached coarse-to-refined joint search, optional
  alpha-midpoint refinement, automatic adaptive-kNN limits through
  `n - 1`, and explicit boundary diagnostics.
- Added Gaussian pathwise AIC, AICc, BIC, GCV, and EBIC. Effective
  degrees of freedom are computed in C++ from the conditional local
  smoother. Primary EBIC adds the exact coefficient-surface model-space
  penalty to pathwise BIC; the separately named `ebic_refit` sensitivity
  score retains post-selection local refits and excludes saturated
  refits with at most one residual degree of freedom.
- Excluded invalid effective-degree-of-freedom rows from Gaussian
  information-criterion selection and exposed zero EBIC model-space
  terms for nonempty full-surface models.
- Added `rule = "sparse_1se"`, richer convergence/active-support
  selection tables, and a same-fold unpenalized LHWR comparison.
- Added
  [`plot.tune.lhwpr()`](https://byuzbasi.github.io/LHWL/reference/plot.tune.lhwpr.md)
  bandwidth, lambda, joint-surface, and alpha views.
- Added unpenalized binomial LHWR through
  `lhwr(..., family = "binomial")`.
- Added local logistic coefficient, probability, link, convergence,
  effective-neighborhood, class-variation, Fisher-rank, standard-error,
  and z-value outputs.
- Added binomial prediction at calibration and new Lorentz locations.
- Added binomial bandwidth selection by leave-one-out Bernoulli log
  loss.
- Made summaries, reports, diagnostics, and Gaussian-only inference
  checks response-family aware.
- Documented LHWL as one package containing the LHWR and LHWPR modules,
  with later learning modules governed by an explicit staged roadmap.

## LHWL 0.2.0

- Added LHWPR Gaussian and binomial ridge, lasso, elastic-net, MCP, and
  SCAD paths and joint bandwidth-regularization tuning.

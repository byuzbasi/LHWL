# Package index

## Fit and tune local models

LHWR is the unpenalized foundation; LHWPR adds focal regularization.

- [`LHWL`](https://byuzbasi.github.io/LHWL/reference/LHWL-package.md)
  [`LHWL-package`](https://byuzbasi.github.io/LHWL/reference/LHWL-package.md)
  : LHWL Package
- [`lhwr()`](https://byuzbasi.github.io/LHWL/reference/lhwr.md) : Fit a
  Lorentz Hyperbolic Weighted Regression Model
- [`bw.lhwr()`](https://byuzbasi.github.io/LHWL/reference/bw.lhwr.md) :
  Select LHWR Bandwidth
- [`lhwpr()`](https://byuzbasi.github.io/LHWL/reference/lhwpr.md)
  [`print(`*`<lhwpr>`*`)`](https://byuzbasi.github.io/LHWL/reference/lhwpr.md)
  [`summary(`*`<lhwpr>`*`)`](https://byuzbasi.github.io/LHWL/reference/lhwpr.md)
  [`print(`*`<summary.lhwpr>`*`)`](https://byuzbasi.github.io/LHWL/reference/lhwpr.md)
  [`coef(`*`<lhwpr>`*`)`](https://byuzbasi.github.io/LHWL/reference/lhwpr.md)
  [`fitted(`*`<lhwpr>`*`)`](https://byuzbasi.github.io/LHWL/reference/lhwpr.md)
  [`residuals(`*`<lhwpr>`*`)`](https://byuzbasi.github.io/LHWL/reference/lhwpr.md)
  [`predict(`*`<lhwpr>`*`)`](https://byuzbasi.github.io/LHWL/reference/lhwpr.md)
  [`plot(`*`<lhwpr>`*`)`](https://byuzbasi.github.io/LHWL/reference/lhwpr.md)
  : Lorentz-Hyperbolic Weighted Penalized Regression
- [`tune.lhwpr()`](https://byuzbasi.github.io/LHWL/reference/tune.lhwpr.md)
  [`print(`*`<tune.lhwpr>`*`)`](https://byuzbasi.github.io/LHWL/reference/tune.lhwpr.md)
  : Joint Bandwidth and Regularization Tuning for LHWPR
- [`lhwpr_control()`](https://byuzbasi.github.io/LHWL/reference/lhwpr_control.md)
  : LHWPR Solver Controls
- [`plot(`*`<tune.lhwpr>`*`)`](https://byuzbasi.github.io/LHWL/reference/plot.tune.lhwpr.md)
  : Plot an LHWPR Tuning Object

## Coordinates and geometry

- [`check_lhwr_Z()`](https://byuzbasi.github.io/LHWL/reference/make_lhwr_Z.md)
  [`make_lhwr_Z_features()`](https://byuzbasi.github.io/LHWL/reference/make_lhwr_Z.md)
  [`make_lhwr_Z_network()`](https://byuzbasi.github.io/LHWL/reference/make_lhwr_Z.md)
  [`make_lhwr_Z()`](https://byuzbasi.github.io/LHWL/reference/make_lhwr_Z.md)
  [`check_lorentz_Z()`](https://byuzbasi.github.io/LHWL/reference/make_lhwr_Z.md)
  [`make_lorentz_Z_features()`](https://byuzbasi.github.io/LHWL/reference/make_lhwr_Z.md)
  [`make_lorentz_Z_network()`](https://byuzbasi.github.io/LHWL/reference/make_lhwr_Z.md)
  [`make_lorentz_Z()`](https://byuzbasi.github.io/LHWL/reference/make_lhwr_Z.md)
  : Construct Lorentz Hyperboloid Coordinates for LHWR
- [`prepare_lorentz_coordinates()`](https://byuzbasi.github.io/LHWL/reference/prepare_lorentz_coordinates.md)
  : Prepare Coordinates for LHWR
- [`simulate_H2()`](https://byuzbasi.github.io/LHWL/reference/simulate_H2.md)
  : Simulate Points on the Two-Dimensional Hyperboloid
- [`lorentz_inner()`](https://byuzbasi.github.io/LHWL/reference/lorentz_inner.md)
  : Lorentz Inner Product
- [`lorentz_norm()`](https://byuzbasi.github.io/LHWL/reference/lorentz_norm.md)
  : Lorentz Norm
- [`hyperbolic_distance()`](https://byuzbasi.github.io/LHWL/reference/hyperbolic_distance.md)
  : Hyperbolic Distance on the Lorentz Hyperboloid
- [`hyperbolic_distance_matrix()`](https://byuzbasi.github.io/LHWL/reference/hyperbolic_distance_matrix.md)
  : Hyperbolic Distance Matrix
- [`lorentz_knn()`](https://byuzbasi.github.io/LHWL/reference/lorentz_knn.md)
  : Lorentz Hyperbolic k-Nearest Neighbors
- [`check_hyperboloid()`](https://byuzbasi.github.io/LHWL/reference/check_hyperboloid.md)
  : Check Hyperboloid Constraints
- [`normalize_to_H2()`](https://byuzbasi.github.io/LHWL/reference/normalize_to_H2.md)
  : Normalize a Vector to the Lorentz Hyperboloid
- [`project_to_H2()`](https://byuzbasi.github.io/LHWL/reference/project_to_H2.md)
  : Project Spatial Coordinates to the Lorentz Hyperboloid
- [`poincare_to_lorentz()`](https://byuzbasi.github.io/LHWL/reference/poincare_to_lorentz.md)
  : Convert Poincare Disk Coordinates to Lorentz Coordinates
- [`extrinsic_lorentz_mean()`](https://byuzbasi.github.io/LHWL/reference/extrinsic_lorentz_mean.md)
  : Extrinsic Mean Projected to the Lorentz Hyperboloid
- [`lorentz_log_map()`](https://byuzbasi.github.io/LHWL/reference/lorentz_log_map.md)
  : Lorentz Logarithmic Map
- [`tangent_coordinates()`](https://byuzbasi.github.io/LHWL/reference/tangent_coordinates.md)
  : Tangent-Space Coordinates
- [`tangent_distance_matrix()`](https://byuzbasi.github.io/LHWL/reference/tangent_distance_matrix.md)
  : Tangent-Space Distance Matrix

## Extractors, diagnostics and graphics

- [`fitted(`*`<lhwr>`*`)`](https://byuzbasi.github.io/LHWL/reference/fitted.lhwr.md)
  [`residuals(`*`<lhwr>`*`)`](https://byuzbasi.github.io/LHWL/reference/fitted.lhwr.md)
  [`coef(`*`<lhwr>`*`)`](https://byuzbasi.github.io/LHWL/reference/fitted.lhwr.md)
  [`predict(`*`<lhwr>`*`)`](https://byuzbasi.github.io/LHWL/reference/fitted.lhwr.md)
  [`fitted(`*`<hgwr>`*`)`](https://byuzbasi.github.io/LHWL/reference/fitted.lhwr.md)
  [`residuals(`*`<hgwr>`*`)`](https://byuzbasi.github.io/LHWL/reference/fitted.lhwr.md)
  [`coef(`*`<hgwr>`*`)`](https://byuzbasi.github.io/LHWL/reference/fitted.lhwr.md)
  [`predict(`*`<hgwr>`*`)`](https://byuzbasi.github.io/LHWL/reference/fitted.lhwr.md)
  : Extract and Predict from LHWR Models
- [`summary(`*`<lhwr>`*`)`](https://byuzbasi.github.io/LHWL/reference/summary.lhwr.md)
  [`summary(`*`<hgwr>`*`)`](https://byuzbasi.github.io/LHWL/reference/summary.lhwr.md)
  : Summarize an LHWR Model
- [`RSS.lhwr()`](https://byuzbasi.github.io/LHWL/reference/RSS.lhwr.md)
  [`RSS.hgwr()`](https://byuzbasi.github.io/LHWL/reference/RSS.lhwr.md)
  : LHWR Residual Sum of Squares
- [`CV.lhwr()`](https://byuzbasi.github.io/LHWL/reference/CV.lhwr.md)
  [`CV.hgwr()`](https://byuzbasi.github.io/LHWL/reference/CV.lhwr.md) :
  LHWR Cross-Validation Score
- [`AICc.lhwr()`](https://byuzbasi.github.io/LHWL/reference/AICc.lhwr.md)
  [`AICc.hgwr()`](https://byuzbasi.github.io/LHWL/reference/AICc.lhwr.md)
  : LHWR AICc
- [`lhwr_report()`](https://byuzbasi.github.io/LHWL/reference/lhwr_report.md)
  : Professional LHWR Model Report
- [`local_collinearity()`](https://byuzbasi.github.io/LHWL/reference/local_collinearity.md)
  : Local Multicollinearity Diagnostics for LHWR
- [`local_inference()`](https://byuzbasi.github.io/LHWL/reference/local_inference.md)
  [`local_inference.hgwr()`](https://byuzbasi.github.io/LHWL/reference/local_inference.md)
  : Local Inference for LHWR Coefficient Surfaces
- [`diagnostic_tests()`](https://byuzbasi.github.io/LHWL/reference/diagnostic_tests.md)
  : LHWR Diagnostic F-Tests
- [`moran.lhwr()`](https://byuzbasi.github.io/LHWL/reference/moran.lhwr.md)
  [`moran.hgwr()`](https://byuzbasi.github.io/LHWL/reference/moran.lhwr.md)
  : Lorentz-Hyperbolic Moran's I Diagnostic for LHWR Residuals
- [`local_moran.lhwr()`](https://byuzbasi.github.io/LHWL/reference/local_moran.lhwr.md)
  [`local_moran.hgwr()`](https://byuzbasi.github.io/LHWL/reference/local_moran.lhwr.md)
  : Lorentz-Hyperbolic LISA Diagnostic (LH-LISA) for LHWR Residuals
- [`lhwr_plotting_grammar`](https://byuzbasi.github.io/LHWL/reference/lhwr_plotting_grammar.md)
  : LHWR Plotting Grammar
- [`lhwr_plot_control()`](https://byuzbasi.github.io/LHWL/reference/lhwr_plot_control.md)
  : LHWR Plot Control
- [`plot_lhwr_diagnostics()`](https://byuzbasi.github.io/LHWL/reference/plot_lhwr_diagnostics.md)
  : Plot LHWR Diagnostics
- [`plot_lhwr_embedding()`](https://byuzbasi.github.io/LHWL/reference/plot_lhwr_embedding.md)
  : Plot an LHWR Embedding
- [`plot_lhwr_geometry()`](https://byuzbasi.github.io/LHWL/reference/plot_lhwr_geometry.md)
  : Plot Lorentz Geometry
- [`plot_lhwr_hyperboloid()`](https://byuzbasi.github.io/LHWL/reference/plot_lhwr_hyperboloid.md)
  : Plot Lorentz Coordinates on the Hyperboloid
- [`plot_lhwr_inference()`](https://byuzbasi.github.io/LHWL/reference/plot_lhwr_inference.md)
  : Plot LHWR Local Inference
- [`plot_lhwr_model()`](https://byuzbasi.github.io/LHWL/reference/plot_lhwr_model.md)
  : Plot LHWR Model Surfaces
- [`plot_lhwr_network_map()`](https://byuzbasi.github.io/LHWL/reference/plot_lhwr_network_map.md)
  : Plot an LHWR Network Map
- [`plot_lhwr_poincare_map()`](https://byuzbasi.github.io/LHWL/reference/plot_lhwr_poincare_map.md)
  : Plot Lorentz Coordinates in the Poincare Disk
- [`plot_lhwpr_surface()`](https://byuzbasi.github.io/LHWL/reference/plot_lhwpr_surface.md)
  : Plot an LHWPR Coefficient or Focal-Selection Surface
- [`plot(`*`<lhwr>`*`)`](https://byuzbasi.github.io/LHWL/reference/plot.lhwr.md)
  : Plot LHWR Objects
- [`plot(`*`<hgwr_local_moran>`*`)`](https://byuzbasi.github.io/LHWL/reference/plot.hgwr_local_moran.md)
  : Plot Lorentz-Hyperbolic LISA Diagnostics (LH-LISA)
- [`caption_lhwr_map()`](https://byuzbasi.github.io/LHWL/reference/caption_lhwr_map.md)
  : Create a Caption for a Poincare Disk Map
- [`summary(`*`<lhwr_poincare_map>`*`)`](https://byuzbasi.github.io/LHWL/reference/summary.lhwr_poincare_map.md)
  : Summarize a Poincare Disk Map

## Data, kernels and compatibility interfaces

- [`wdi_lhwr`](https://byuzbasi.github.io/LHWL/reference/wdi_lhwr.md) :
  World Bank WDI panel prepared for LHWR examples
- [`hgwr.basic()`](https://byuzbasi.github.io/LHWL/reference/hgwr.basic.md)
  : Deprecated HGWR model-fitting alias
- [`hgwr.fit()`](https://byuzbasi.github.io/LHWL/reference/hgwr.fit.md)
  : Fit the Low-Level LHWR Model
- [`hgwr.fit.knn()`](https://byuzbasi.github.io/LHWL/reference/hgwr.fit.knn.md)
  : Fit the Low-Level LHWR Model from Lorentz kNN Neighborhoods
- [`hgwr_bisquare()`](https://byuzbasi.github.io/LHWL/reference/hgwr_bisquare.md)
  : Bisquare HGWR Kernel
- [`hgwr_boxcar()`](https://byuzbasi.github.io/LHWL/reference/hgwr_boxcar.md)
  : Boxcar HGWR Kernel
- [`hgwr_exponential()`](https://byuzbasi.github.io/LHWL/reference/hgwr_exponential.md)
  : Exponential HGWR Kernel
- [`hgwr_gaussian()`](https://byuzbasi.github.io/LHWL/reference/hgwr_gaussian.md)
  : Gaussian HGWR Kernel
- [`hgwr_kernel()`](https://byuzbasi.github.io/LHWL/reference/hgwr_kernel.md)
  : HGWR Kernel Interface
- [`hgwr_tricube()`](https://byuzbasi.github.io/LHWL/reference/hgwr_tricube.md)
  : Tricube HGWR Kernel
- [`bw.hgwr()`](https://byuzbasi.github.io/LHWL/reference/bw.hgwr.md) :
  Deprecated LHWR bandwidth-selection alias
- [`weighted_cor()`](https://byuzbasi.github.io/LHWL/reference/weighted_cor.md)
  : Weighted Correlation Matrix
- [`weighted_r2()`](https://byuzbasi.github.io/LHWL/reference/weighted_r2.md)
  : Weighted R-squared
- [`safe_condition_number()`](https://byuzbasi.github.io/LHWL/reference/safe_condition_number.md)
  : Safe Condition Number

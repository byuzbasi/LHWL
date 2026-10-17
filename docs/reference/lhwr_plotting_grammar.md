# LHWR Plotting Grammar

Convenience plotting functions that organize LHWR graphics by purpose.
These wrappers keep the package interface predictable for applied users
while still relying on the lower-level plotting functions such as
[`plot.lhwr()`](https://byuzbasi.github.io/LHWL/reference/plot.lhwr.md)
and
[`plot_lhwr_poincare_map()`](https://byuzbasi.github.io/LHWL/reference/plot_lhwr_poincare_map.md).

## Details

Use the functions as a decision guide:
[`plot_lhwr_geometry()`](https://byuzbasi.github.io/LHWL/reference/plot_lhwr_geometry.md)
for the Lorentz geometry itself,
[`plot_lhwr_embedding()`](https://byuzbasi.github.io/LHWL/reference/plot_lhwr_embedding.md)
for embedded observations,
[`plot_lhwr_network_map()`](https://byuzbasi.github.io/LHWL/reference/plot_lhwr_network_map.md)
for graph or network data,
[`plot_lhwr_model()`](https://byuzbasi.github.io/LHWL/reference/plot_lhwr_model.md)
for fitted coefficients and residuals,
[`plot_lhwr_inference()`](https://byuzbasi.github.io/LHWL/reference/plot_lhwr_inference.md)
for local p-values and significance, and
[`plot_lhwr_diagnostics()`](https://byuzbasi.github.io/LHWL/reference/plot_lhwr_diagnostics.md)
for model diagnostic surfaces.

See
[`vignette("lhwr-graphics", package = "LHWL")`](https://byuzbasi.github.io/LHWL/articles/lhwr-graphics.md)
for executable paired coefficient/significance maps, coefficient-colour
plus significance-ring overlays, fitted/residual surfaces, LISA maps,
network overlays and guarded PDF/PNG export. Helpers defined in that
guide are vignette recipes, not new exported package functions. A
default model projection is the Lorentz spatial-coordinate plane, not a
Poincare or geographic map.

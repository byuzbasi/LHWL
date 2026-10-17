# Plot an LHWPR Coefficient or Focal-Selection Surface

Display one fitted LHWPR slope over the Poincare-disk view of the
Lorentz representation. A coefficient surface colors nonzero focal
coefficients and leaves numerical zeros hollow. A selection surface
shows the focal active set for a positively penalized predictor. These
displays summarize penalized estimation and are not significance maps.

## Usage

``` r
plot_lhwpr_surface(
  x,
  predictor,
  type = c("coefficient", "selection"),
  threshold = 1e-08,
  labels = NULL,
  node_names = NULL,
  sizes = NULL,
  title = NULL,
  point_palette = "Blue-Red 3",
  selection_colors = c(`Not selected` = "#D9D9D9", Selected = "#2166AC"),
  show_legend = TRUE,
  radial_rescale = TRUE,
  radial_quantile = 1,
  radial_power = 1,
  output_file = NULL,
  width = 7,
  height = 7,
  res = 220,
  ...
)
```

## Arguments

- x:

  A fitted `"lhwpr"` object.

- predictor:

  Predictor name or one-based slope index. `"Intercept"` is available
  only for `type = "coefficient"`.

- type:

  Draw coefficient values or focal selection.

- threshold:

  Nonnegative numerical threshold used to distinguish stored zeros from
  active coefficients.

- labels, node_names, sizes:

  Optional point labels, longer names, and nonnegative point-size
  values; see
  [`plot_lhwr_poincare_map()`](https://byuzbasi.github.io/LHWL/reference/plot_lhwr_poincare_map.md).

- title:

  Optional plot title. `NULL` creates a descriptive default.

- point_palette:

  HCL palette for nonzero coefficient values.

- selection_colors:

  Two colors for not-selected and selected locations.

- show_legend:

  Logical; show the selection legend or the filled-versus-hollow
  coefficient key.

- radial_rescale, radial_quantile, radial_power:

  Display-only Poincare radial controls; they do not alter fitting
  distances.

- output_file:

  Optional PNG, PDF, EPS, or PS output path.

- width, height, res:

  Output-device dimensions in inches and resolution for raster output.

- ...:

  Additional arguments passed to
  [`plot_lhwr_poincare_map()`](https://byuzbasi.github.io/LHWL/reference/plot_lhwr_poincare_map.md).

## Value

Invisibly, a data frame containing plotted coordinates, the focal
coefficient, the active indicator, predictor, threshold, and plot type.

## Examples

``` r
set.seed(7)
H <- simulate_H2(35, r_max = 2)
dat <- data.frame(y = rnorm(35), x1 = rnorm(35), x2 = rnorm(35))
fit <- lhwpr(
  y ~ x1 + x2, dat, H$Z, bw = 18, lambda = c(0.2, 0.05),
  keep_path = "selected"
)
plot_lhwpr_surface(fit, "x1", type = "coefficient")
plot_lhwpr_surface(fit, "x1", type = "selection")
```

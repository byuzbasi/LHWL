# Plot Lorentz Coordinates on the Hyperboloid

Visualize Lorentz hyperboloid coordinates on the upper sheet
`z0 = sqrt(1 + z1^2 + z2^2)`. This plot emphasizes the geometry used by
LHWR estimation: observations live on the Lorentz hyperboloid, and local
neighborhoods are defined by Lorentz-hyperbolic distances. For readable
two-dimensional displays, use
[`plot_lhwr_poincare_map()`](https://byuzbasi.github.io/LHWL/reference/plot_lhwr_poincare_map.md).

## Usage

``` r
plot_lhwr_hyperboloid(
  Z,
  values = NULL,
  labels = NULL,
  sizes = NULL,
  title = "LHWR Lorentz hyperboloid view",
  value_label = NULL,
  point_palette = "Blue-Red 3",
  surface_col = "#D7E8F0",
  surface_alpha = 0.55,
  surface_border = "#FFFFFF",
  theta = 38,
  phi = 24,
  expand = 0.8,
  grid_n = 55L,
  surface_quantile = "auto",
  point_size = 1.1,
  show_labels = FALSE,
  n_label = 30L,
  mar = c(1.4, 1.4, 2.8, 0.8),
  label_cex = NULL,
  colorbar_cex = NULL,
  title_cex = NULL,
  point_min = NULL,
  point_max = NULL,
  plot_control = NULL,
  show_colorbar = TRUE,
  output_file = NULL,
  device = c("auto", "png", "eps", "pdf"),
  width = 7,
  height = 6,
  res = 220,
  ...
)
```

## Arguments

- Z:

  An `n x 3` Lorentz hyperboloid coordinate matrix with columns
  `(z0, z1, z2)`.

- values:

  Optional numeric vector used to color points. If omitted, the
  hyperbolic radius `acosh(z0)` is used.

- labels:

  Optional point labels.

- sizes:

  Optional numeric vector used to scale point sizes.

- title:

  Plot title. Use `NULL` for manuscript-style panels.

- value_label:

  Label for the colorbar. If omitted, a sensible default is chosen from
  `values`.

- point_palette:

  HCL palette for numeric `values`.

- surface_col:

  Color of the hyperboloid surface.

- surface_alpha:

  Transparency of the hyperboloid surface.

- surface_border:

  Grid-line color for the surface.

- theta, phi:

  Perspective angles passed to
  [`graphics::persp()`](https://rdrr.io/r/graphics/persp.html).

- expand:

  Expansion factor passed to
  [`graphics::persp()`](https://rdrr.io/r/graphics/persp.html).

- grid_n:

  Number of grid points per axis used to draw the surface.

- surface_quantile:

  Quantile in `(0, 1]` used to set the displayed surface extent from
  `abs(z1)` and `abs(z2)`, or `"auto"`. The automatic setting keeps the
  full range unless a few extreme coordinates would make the hyperboloid
  surface visually too large. Point coordinates are not changed.

- point_size:

  Fixed point size used when `sizes` is omitted.

- show_labels:

  Logical; draw labels for selected points.

- n_label:

  Number of labels to draw when `show_labels = TRUE`.

- mar:

  Base graphics margins passed to
  [`graphics::par()`](https://rdrr.io/r/graphics/par.html) for the
  hyperboloid panel.

- label_cex, colorbar_cex, title_cex:

  Character sizes. These can also be supplied through `plot_control`.
  For `label_cex`, `"auto"` chooses a small value from the number of
  observations when labels are drawn.

- point_min, point_max:

  Minimum and maximum point sizes when `sizes` is supplied.

- plot_control:

  Optional list created by
  [`lhwr_plot_control()`](https://byuzbasi.github.io/LHWL/reference/lhwr_plot_control.md).
  Direct arguments such as `label_cex` or `colorbar_cex` override values
  in this list.

- show_colorbar:

  Logical; draw a numeric colorbar.

- output_file:

  Optional path. If supplied, the plot is written to this file. The
  device is inferred from the file extension unless `device` is
  supplied.

- device:

  Output device: `"auto"`, `"png"`, `"eps"`, or `"pdf"`.

- width, height:

  Output size in inches when `output_file` is supplied.

- res:

  Raster resolution for PNG output and fallback resolution for EPS.

- ...:

  Additional arguments passed to
  [`graphics::persp()`](https://rdrr.io/r/graphics/persp.html).

## Value

Invisibly returns a data frame containing the plotted Lorentz
coordinates, hyperbolic radii, values, and point sizes.

## Examples

``` r
H <- simulate_H2(30, r_max = 2)
plot_lhwr_hyperboloid(H$Z, values = H$r, title = NULL)
```

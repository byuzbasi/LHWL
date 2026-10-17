# Plot Lorentz Coordinates in the Poincare Disk

Visualize Lorentz hyperboloid coordinates in the Poincare disk. This is
a display tool: distances and model estimation in LHWR are still based
on the original Lorentz-hyperbolic geometry unless the user explicitly
uses another geometry. The function is useful for network, similarity,
and trade-map applications where nodes, edge backbones, labels, node
sizes, and value-based colors should be shown in a hyperbolic disk view.

## Usage

``` r
plot_lhwr_poincare_map(
  Z,
  values = NULL,
  groups = NULL,
  significant = NULL,
  significance_adjusted = NULL,
  significance_direction = NULL,
  labels = NULL,
  node_names = NULL,
  sizes = NULL,
  edges = NULL,
  edge_weight = NULL,
  edge_alpha_range = c(0.08, 0.3),
  edge_lwd_range = c(0.25, 1.6),
  node_id = NULL,
  title = "LHWR Poincare disk map",
  value_label = "Value",
  group_label = "groups",
  size_label = "Node size",
  edge_label = "Edges",
  disk_radius = 0.94,
  radial_rescale = TRUE,
  radial_quantile = 1,
  radial_power = 1,
  show_all_labels = TRUE,
  n_label = 30L,
  label_cex = NULL,
  degree_cex = NULL,
  legend_cex = NULL,
  colorbar_cex = NULL,
  title_cex = NULL,
  info_cex = NULL,
  point_min = NULL,
  point_max = NULL,
  plot_control = NULL,
  point_palette = "Blue-Red 3",
  group_palette = NULL,
  significance_colors = c(`Not significant` = "#BDBDBD", Significant = "#2E7D32"),
  signed_significance_colors = c(`Not significant` = "#C9C9C9", `Negative significant` =
    "#7B3294", `Positive significant` = "#2E7D32"),
  adjusted_significance_colors = c(`Not significant` = "#C9C9C9", `Raw only` = "#E69F00",
    `Raw and BH` = "#1B7837"),
  significance_label = "Local significance",
  background_col = "#F6F1E7",
  boundary_col = "#0A9396",
  edge_col = "#607D8B",
  show_colorbar = TRUE,
  show_group_legend = TRUE,
  group_legend_position = c("right", "bottom"),
  show_size_legend = TRUE,
  show_info = TRUE,
  output_file = NULL,
  width = 1500,
  height = 1500,
  res = 170,
  ...
)
```

## Arguments

- Z:

  An `n x 3` Lorentz hyperboloid coordinate matrix with columns
  `(z0, z1, z2)`.

- values:

  Optional numeric vector used to color points. If omitted and `groups`
  is also omitted, points use a single color.

- groups:

  Optional grouping variable used for categorical point colors. Ignored
  when `values` or `significant` is supplied.

- significant:

  Optional logical vector. If supplied, points are colored as
  significant/non-significant using `significance_colors`; useful for
  local p-value maps.

- significance_adjusted:

  Optional logical vector used together with `significant`. If supplied,
  local significance is classified as not significant, raw-only
  significant, or significant after adjustment. This is useful for
  comparing raw local p-values with multiplicity-adjusted p-values
  without drawing a second map.

- significance_direction:

  Optional numeric vector used together with `significant`. Negative
  values are shown as negative significant, positive values as positive
  significant, and non-significant locations remain grey. This is useful
  for local coefficient maps where both significance and sign should be
  visible.

- labels:

  Optional labels for points, for example country ISO3 codes.

- node_names:

  Optional longer node names used in map summaries.

- sizes:

  Optional nonnegative numeric vector for point sizes. A common choice
  is original total trade or a node degree/strength measure. The map
  already applies `log1p` before rescaling sizes; do not log-transform
  twice unless that additional visual transformation is explicitly
  intended.

- edges:

  Optional edge list as a data frame or matrix with two columns. Entries
  may be integer row indices or labels matching `node_id`.

- edge_weight:

  Optional numeric vector of edge weights used for line widths and
  transparency.

- edge_alpha_range:

  Length-2 numeric vector giving the minimum and maximum edge
  transparency. Larger values make network links more visible.

- edge_lwd_range:

  Length-2 numeric vector giving the minimum and maximum edge line
  width.

- node_id:

  Optional node identifiers used to match `edges` when edge endpoints
  are not integer row indices.

- title:

  Plot title. Use `NULL` to suppress the title, which is useful when the
  title is supplied in a manuscript caption.

- value_label:

  Label for the colorbar when `values` is supplied.

- group_label:

  Label used in summaries and captions when `groups` is supplied.

- size_label:

  Label explaining point size.

- edge_label:

  Label explaining edges.

- disk_radius:

  Visual radius used by the largest point after optional radial scaling.
  This changes only the display, not `Z` or distances.

- radial_rescale:

  If `TRUE`, rescales Poincare coordinates for display. This changes
  only the graphic, not `Z` or the distances used by LHWR.

- radial_quantile:

  Quantile of the Poincare radii used as the visual reference radius
  when `radial_rescale = TRUE`. The default `1` preserves the previous
  behavior. Values below 1 are useful when a few peripheral points would
  otherwise compress the rest of the map.

- radial_power:

  Positive power applied to rescaled radii for display. Values below 1
  spread central points outward while preserving angles.

- show_all_labels:

  If `TRUE`, all labels are drawn. If `FALSE`, only the largest
  `n_label` points according to `sizes` are labeled.

- n_label:

  Number of labels to draw when `show_all_labels = FALSE`.

- label_cex:

  Label size. If `NULL`, the value in `plot_control` is used.

- degree_cex:

  Character size for degree labels on the disk boundary.

- legend_cex:

  Character size for categorical, significance, and size legends.

- colorbar_cex:

  Character size for numeric colorbar labels and title.

- title_cex:

  Character size for the plot title.

- info_cex:

  Character size for the optional information box.

- point_min, point_max:

  Minimum and maximum point sizes after scaling.

- plot_control:

  Optional list created by
  [`lhwr_plot_control()`](https://byuzbasi.github.io/LHWL/reference/lhwr_plot_control.md).
  Direct arguments such as `label_cex` or `legend_cex` override values
  in this list.

- point_palette:

  HCL palette for numeric `values`.

- group_palette:

  Optional vector of colors for `groups`.

- significance_colors:

  Named or unnamed length-2 color vector for `significant = FALSE` and
  `significant = TRUE`, respectively.

- signed_significance_colors:

  Named or unnamed length-3 color vector for not significant, negative
  significant, and positive significant locations when both
  `significant` and `significance_direction` are supplied.

- adjusted_significance_colors:

  Named or unnamed length-3 color vector for not significant, raw-only
  significant, and raw-and-adjusted significant locations when both
  `significant` and `significance_adjusted` are supplied.

- significance_label:

  Legend title used when `significant` is supplied.

- background_col:

  Disk background color.

- boundary_col:

  Disk boundary color.

- edge_col:

  Edge color.

- show_colorbar:

  Logical; draw a colorbar for numeric `values`.

- show_group_legend:

  Logical; draw a legend for categorical `groups`.

- group_legend_position:

  Position for the categorical/significance legend. Use `"right"` for a
  side legend or `"bottom"` to keep the legend below the disk when side
  placement overlaps the plot.

- show_size_legend:

  Logical; draw point-size legend.

- show_info:

  Logical; draw an explanatory annotation box.

- output_file:

  Optional output path. PDF and EPS/PS extensions select the
  corresponding device; other extensions use PNG. Existing files may be
  overwritten: check the target first if preservation is required.

- width, height, res:

  Width and height in pixels and resolution in pixels per inch. PDF/EPS
  dimensions in inches are `width/res` and `height/res`. These differ
  from the inch-based hyperboloid output settings.

- ...:

  Additional graphical arguments currently ignored.

## Value

Invisibly returns a data frame containing the Poincare coordinates used
for plotting.

## Examples

``` r
set.seed(1)
H <- simulate_H2(40, r_max = 3)
Z <- H$Z
plot_lhwr_poincare_map(
  Z,
  values = H$r,
  labels = paste0("n", seq_len(nrow(Z))),
  title = "Simulated points in the Poincare disk",
  value_label = "radius"
)
```

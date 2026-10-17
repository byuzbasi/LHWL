# Plot an LHWR Network Map

Display graph or network data whose nodes have Lorentz hyperboloid
coordinates. This is the standard plotting interface for trade,
citation, social, transportation, similarity, and general graph
embeddings.

## Usage

``` r
plot_lhwr_network_map(
  Z,
  edges,
  edge_weight = NULL,
  values = NULL,
  groups = NULL,
  significant = NULL,
  significance_adjusted = NULL,
  significance_direction = NULL,
  labels = NULL,
  node_names = NULL,
  sizes = NULL,
  node_id = NULL,
  title = NULL,
  value_label = "Value",
  group_label = "network communities",
  size_label = "Node size",
  edge_label = "Network edges",
  edge_alpha_range = c(0.008, 0.055),
  edge_lwd_range = c(0.025, 0.3),
  show_info = FALSE,
  ...
)
```

## Arguments

- Z:

  An `n x 3` Lorentz hyperboloid coordinate matrix with columns
  `(z0, z1, z2)`.

- edges:

  Optional edge list as a data frame or matrix with two columns. Entries
  may be integer row indices or labels matching `node_id`.

- edge_weight:

  Optional numeric vector of edge weights used for line widths and
  transparency.

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

- edge_alpha_range:

  Length-2 numeric vector giving the minimum and maximum edge
  transparency. Larger values make network links more visible.

- edge_lwd_range:

  Length-2 numeric vector giving the minimum and maximum edge line
  width.

- show_info:

  Logical; draw an explanatory annotation box.

- ...:

  Additional graphical arguments currently ignored.

## Value

Invisibly returns an object of class `lhwr_poincare_map`.

## Examples

``` r
H <- simulate_H2(20, r_max = 2)
edges <- data.frame(from = 1:5, to = 6:10)
plot_lhwr_network_map(H$Z, edges = edges, title = NULL)
```

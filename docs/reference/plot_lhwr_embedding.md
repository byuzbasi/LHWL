# Plot an LHWR Embedding

Display embedded observations in the Poincare disk. This is the standard
plot for non-network feature maps, economic indicators, clinical
variables, functional scores, shape scores, and other two-dimensional
representations converted to Lorentz coordinates.

## Usage

``` r
plot_lhwr_embedding(
  Z,
  values = NULL,
  groups = NULL,
  significant = NULL,
  significance_adjusted = NULL,
  significance_direction = NULL,
  labels = NULL,
  node_names = NULL,
  sizes = NULL,
  title = NULL,
  value_label = "Value",
  group_label = "groups",
  show_info = FALSE,
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

- title:

  Plot title. Use `NULL` to suppress the title, which is useful when the
  title is supplied in a manuscript caption.

- value_label:

  Label for the colorbar when `values` is supplied.

- group_label:

  Label used in summaries and captions when `groups` is supplied.

- show_info:

  Logical; draw an explanatory annotation box.

- ...:

  Additional graphical arguments currently ignored.

## Value

Invisibly returns an object of class `lhwr_poincare_map`.

## Examples

``` r
H <- simulate_H2(30, r_max = 2)
plot_lhwr_embedding(H$Z, values = H$r, title = NULL)
```

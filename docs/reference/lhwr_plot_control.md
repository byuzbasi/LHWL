# LHWR Plot Control

Create a small list of graphical defaults used by the high-level LHWR
Poincare disk plotting functions. Most users can ignore this helper and
use the plotting functions directly. It is useful when preparing
publication figures or when labels and legends need to be larger or
smaller.

## Usage

``` r
lhwr_plot_control(
  label_cex = 0.55,
  degree_cex = 0.5,
  legend_cex = 0.55,
  colorbar_cex = 0.55,
  title_cex = 1.2,
  info_cex = 0.55,
  point_min = 0.45,
  point_max = 2.8
)
```

## Arguments

- label_cex:

  Character expansion for observation labels.

- degree_cex:

  Character expansion for degree labels on the disk boundary.

- legend_cex:

  Character expansion for categorical and size legends.

- colorbar_cex:

  Character expansion for colorbar text.

- title_cex:

  Character expansion for plot titles.

- info_cex:

  Character expansion for the optional information box.

- point_min, point_max:

  Minimum and maximum point sizes after scaling.

## Value

A list that can be passed to `plot_control` in
[`plot_lhwr_poincare_map()`](https://byuzbasi.github.io/LHWL/reference/plot_lhwr_poincare_map.md),
[`plot_lhwr_embedding()`](https://byuzbasi.github.io/LHWL/reference/plot_lhwr_embedding.md),
or
[`plot_lhwr_network_map()`](https://byuzbasi.github.io/LHWL/reference/plot_lhwr_network_map.md).

## Examples

``` r
H <- simulate_H2(30, r_max = 2)
ctl <- lhwr_plot_control(label_cex = 0.65, legend_cex = 0.65)
plot_lhwr_embedding(H$Z, values = H$r, plot_control = ctl)
```

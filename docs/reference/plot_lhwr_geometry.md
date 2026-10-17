# Plot Lorentz Geometry

Inspect the Lorentz hyperboloid coordinates used by LHWR. The
hyperboloid view emphasizes the theoretical model, the Poincare view
gives a disk display, and the radial view shows the radial/angular
distribution of observations.

## Usage

``` r
plot_lhwr_geometry(
  Z,
  type = c("poincare", "hyperboloid", "radial"),
  values = NULL,
  labels = NULL,
  title = NULL,
  palette = "Viridis",
  show_size_legend = FALSE,
  show_info = FALSE,
  ...
)
```

## Arguments

- Z:

  An `n x 3` Lorentz hyperboloid coordinate matrix.

- type:

  Plot type: `"poincare"`, `"hyperboloid"`, or `"radial"`.

- values:

  Optional numeric vector for point colors. Defaults to the hyperbolic
  radius.

- labels:

  Optional point labels.

- title:

  Optional title. `NULL` suppresses the Poincare title but selects a
  default title for hyperboloid/radial views; use `""` for those views
  when the title is provided by a manuscript caption.

- palette:

  HCL palette for numeric values.

- show_size_legend:

  Logical; draw the point-size legend in the Poincare view.

- show_info:

  Logical; draw the explanatory information box in the Poincare view.

- ...:

  Additional arguments passed to the underlying plotting function.

## Value

Invisibly returns plotting coordinates or the input `Z`.

## Examples

``` r
H <- simulate_H2(40, r_max = 3)
plot_lhwr_geometry(H$Z, type = "poincare", title = NULL)
plot_lhwr_geometry(H$Z, type = "radial")
```

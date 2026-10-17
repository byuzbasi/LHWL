# Summarize a Poincare Disk Map

Print a compact interpretation summary for an object returned by
[`plot_lhwr_poincare_map()`](https://byuzbasi.github.io/LHWL/reference/plot_lhwr_poincare_map.md).
The summary keeps the graphic clean while making the figure metadata
reproducible, including the number of nodes and edges, radial-size
correlations, and the most central observations.

## Usage

``` r
# S3 method for class 'lhwr_poincare_map'
summary(object, top_n = 10L, ...)
```

## Arguments

- object:

  An object returned invisibly by
  [`plot_lhwr_poincare_map()`](https://byuzbasi.github.io/LHWL/reference/plot_lhwr_poincare_map.md).

- top_n:

  Number of most central observations to print.

- ...:

  Unused.

## Value

Invisibly returns a list containing map metadata and the central-node
table.

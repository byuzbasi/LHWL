# Simulate Points on the Two-Dimensional Hyperboloid

Simulate Points on the Two-Dimensional Hyperboloid

## Usage

``` r
simulate_H2(n, r_max = 8, uniform_area = TRUE)
```

## Arguments

- n:

  Number of points.

- r_max:

  Maximum hyperbolic radius.

- uniform_area:

  If `TRUE`, sample uniformly with respect to hyperbolic area inside the
  radius `r_max`; otherwise sample radius uniformly.

## Value

A list with Lorentz coordinates `Z`, spatial projection `S`, radius `r`,
and angle `theta`.

## Examples

``` r
set.seed(1)
H <- simulate_H2(10, r_max = 2)
head(H$Z)
check_hyperboloid(H$Z)
```

# Plot an LHWPR Tuning Object

Plot bandwidth and lambda profiles, a joint bandwidth–lambda surface, an
exact evaluated-candidate landscape, and an alpha profile for an object
returned by
[`tune.lhwpr()`](https://byuzbasi.github.io/LHWL/reference/tune.lhwpr.md).

## Usage

``` r
# S3 method for class 'tune.lhwpr'
plot(
  x,
  type = c(
    "all", "bandwidth", "lambda", "landscape", "surface", "alpha"
  ),
  criterion = NULL,
  alpha = NULL,
  gamma = NULL,
  include_lhwr = TRUE,
  ...
)
```

## Arguments

- x:

  A `"tune.lhwpr"` object.

- type:

  One of `"all"`, `"bandwidth"`, `"lambda"`, `"landscape"`, `"surface"`,
  or `"alpha"`.

- criterion:

  Criterion to display. `NULL` uses the criterion that selected the
  model.

- alpha, gamma:

  Optional alpha and gamma defining the displayed surface slice. With
  `type = "surface"` or `type = "landscape"` and `alpha = NULL`, one
  panel is drawn for each evaluated alpha.

- include_lhwr:

  Whether the same-fold unpenalized LHWR reference is overlaid on the
  bandwidth profile.

- ...:

  Additional graphical arguments passed to profile plots.

## Details

The selected tuning point is marked in red. A selected lower or upper
bandwidth boundary is annotated in the bandwidth and surface panels. The
surface uses relative lambda so configurations whose absolute lambda
maxima differ remain comparable.

The landscape plots only evaluated candidates and therefore does not
imply interpolation between grid points. Point shapes distinguish coarse
and refined stages, grey crosses mark candidates excluded by
convergence, KKT, or criterion-specific validity gates, and the selected
eligible point is marked separately. Profiles and interpolated surfaces
are computed only from eligible candidates.

## Value

Invisibly, the data used by the requested plot. For `type = "all"`, a
named list containing the four panel data sets is returned; its joint
panel is named `landscape`.

## Examples

``` r
set.seed(3)
H <- simulate_H2(30, r_max = 1.5)
dat <- data.frame(y = rnorm(30), x = rnorm(30))
selection <- tune.lhwpr(
  y ~ x,
  data = dat,
  Z = H$Z,
  bw_grid = c(14, 20),
  nlambda = 3,
  nfolds = 2,
  refit = FALSE
)
plot(selection, type = "bandwidth")
```

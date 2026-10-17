# Construct Lorentz Hyperboloid Coordinates for LHWR

These helpers turn raw two-component representations or relation
matrices into the `n x 3` Lorentz hyperboloid coordinate matrix `Z` used
by [`lhwr()`](https://byuzbasi.github.io/LHWL/reference/lhwr.md) and
[`bw.lhwr()`](https://byuzbasi.github.io/LHWL/reference/bw.lhwr.md). The
coordinates are not required to be longitude and latitude. They may come
from geographic coordinates, economic scores, clinical variables,
functional or shape summaries, similarity matrices, trade networks,
social networks, transportation networks, or graph embeddings.

## Usage

``` r
check_lhwr_Z(Z, tol = 1e-06)

make_lhwr_Z_features(X, scale = TRUE, tol = 1e-06)

make_lhwr_Z_network(
  A,
  method = c("shortest_path_mds", "spectral"),
  scale = TRUE,
  disconnected = c("stop", "impute"),
  symmetrize = TRUE,
  tol = 1e-06
)

make_lhwr_Z(x, method = c("features", "network"), ...)

check_lorentz_Z(Z, tol = 1e-06)

make_lorentz_Z_features(X, scale = TRUE, tol = 1e-06)

make_lorentz_Z_network(
  A,
  method = c("shortest_path_mds", "spectral"),
  scale = TRUE,
  disconnected = c("stop", "impute"),
  symmetrize = TRUE,
  tol = 1e-06
)

make_lorentz_Z(x, method = c("features", "network"), ...)
```

## Arguments

- Z:

  An `n x 3` Lorentz hyperboloid coordinate matrix to validate.

- tol:

  Numerical tolerance for checking the hyperboloid constraint.

- X:

  An `n x 2` numeric matrix or data frame containing a two-component
  representation of the observations.

- scale:

  If `TRUE`, standardize the two input components before constructing
  `Z`. A constant component has undefined sample standard deviation; it
  is mapped to zero with an explicit warning.

- A:

  An `n x n` numeric adjacency, similarity, or relation matrix.

- method:

  Construction or network embedding method. In `make_lhwr_Z()`, use
  `"features"` for a two-component feature map or `"network"` for a
  relation matrix. In `make_lhwr_Z_network()`, use `"shortest_path_mds"`
  or `"spectral"`.

- disconnected:

  Handling of disconnected graphs for `"shortest_path_mds"`. `"stop"`
  reports an error; `"impute"` replaces infinite graph distances by a
  large finite distance before MDS.

- symmetrize:

  If `TRUE`, non-symmetric relation matrices are replaced by
  `(A + t(A)) / 2` with a warning.

- x:

  Input object passed to the generic constructor. Use an `n x 2` feature
  matrix when `method = "features"`, or an `n x n` relation matrix when
  `method = "network"`.

- ...:

  Arguments passed from `make_lhwr_Z()` to the selected construction
  helper.

## Value

An `n x 3` numeric matrix with columns `z0`, `z1`, and `z2`.

## Details

Let `x_i` denote observation `i`. Suppose a two-component representation
maps each observation to `(z_i1, z_i2)`. LHWR embeds this pair into the
upper Lorentz hyperboloid by

\$\$ z_i = \left(\sqrt{1 + z\_{i1}^2 + z\_{i2}^2}, z\_{i1},
z\_{i2}\right). \$\$

The resulting point satisfies

\$\$-z\_{i0}^2 + z\_{i1}^2 + z\_{i2}^2 = -1,\quad z\_{i0} \> 0.\$\$

For non-geographic data, LHWR should be interpreted as locally weighted
regression over a latent hyperbolic geometry. The validity of the model
depends on whether the chosen representation captures meaningful
neighborhood relations among observations.

## Examples

``` r
X <- cbind(score1 = c(-1, 0, 1), score2 = c(0.5, 0, -0.5))
Z <- make_lhwr_Z_features(X)
check_lhwr_Z(Z)

A <- matrix(c(
  0, 1, 1, 0,
  1, 0, 1, 1,
  1, 1, 0, 1,
  0, 1, 1, 0
), nrow = 4, byrow = TRUE)
Z_net <- make_lhwr_Z_network(A, method = "spectral")
check_lhwr_Z(Z_net)
```

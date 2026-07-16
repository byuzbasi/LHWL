#' Construct Lorentz Hyperboloid Coordinates for LHWR
#'
#' These helpers turn raw two-component representations or relation matrices
#' into the `n x 3` Lorentz hyperboloid coordinate matrix `Z` used by
#' [lhwr()] and [bw.lhwr()]. The coordinates are not required to be longitude
#' and latitude. They may come from geographic coordinates, economic scores,
#' clinical variables, functional or shape summaries, similarity matrices,
#' trade networks, social networks, transportation networks, or graph
#' embeddings.
#'
#' Let `x_i` denote observation `i`. Suppose a two-component representation
#' maps each observation to `(z_i1, z_i2)`. LHWR embeds this pair into the
#' upper Lorentz hyperboloid by
#'
#' \deqn{
#' z_i = \left(\sqrt{1 + z_{i1}^2 + z_{i2}^2}, z_{i1}, z_{i2}\right).
#' }
#'
#' The resulting point satisfies
#'
#' \deqn{-z_{i0}^2 + z_{i1}^2 + z_{i2}^2 = -1,\quad z_{i0} > 0.}
#'
#' For non-geographic data, LHWR should be interpreted as locally weighted
#' regression over a latent hyperbolic geometry. The validity of the model
#' depends on whether the chosen representation captures meaningful
#' neighborhood relations among observations.
#'
#' @param X An `n x 2` numeric matrix or data frame containing a two-component
#'   representation of the observations.
#' @param x Input object passed to the generic constructor. Use an `n x 2`
#'   feature matrix when `method = "features"`, or an `n x n` relation matrix
#'   when `method = "network"`.
#' @param Z An `n x 3` Lorentz hyperboloid coordinate matrix to validate.
#' @param scale If `TRUE`, standardize the two input components before
#'   constructing `Z`.
#' @param A An `n x n` numeric adjacency, similarity, or relation matrix.
#' @param method Construction or network embedding method. In `make_lhwr_Z()`,
#'   use `"features"` for a two-component feature map or `"network"` for a
#'   relation matrix. In `make_lhwr_Z_network()`, use `"shortest_path_mds"` or
#'   `"spectral"`.
#' @param disconnected Handling of disconnected graphs for
#'   `"shortest_path_mds"`. `"stop"` reports an error; `"impute"` replaces
#'   infinite graph distances by a large finite distance before MDS.
#' @param symmetrize If `TRUE`, non-symmetric relation matrices are replaced by
#'   `(A + t(A)) / 2` with a warning.
#' @param tol Numerical tolerance for checking the hyperboloid constraint.
#' @param ... Arguments passed from `make_lhwr_Z()` to the selected
#'   construction helper.
#' @return An `n x 3` numeric matrix with columns `z0`, `z1`, and `z2`.
#' @examples
#' X <- cbind(score1 = c(-1, 0, 1), score2 = c(0.5, 0, -0.5))
#' Z <- make_lhwr_Z_features(X)
#' check_lhwr_Z(Z)
#'
#' A <- matrix(c(
#'   0, 1, 1, 0,
#'   1, 0, 1, 1,
#'   1, 1, 0, 1,
#'   0, 1, 1, 0
#' ), nrow = 4, byrow = TRUE)
#' Z_net <- make_lhwr_Z_network(A, method = "spectral")
#' check_lhwr_Z(Z_net)
#' @name make_lhwr_Z
NULL

#' @rdname make_lhwr_Z
#' @export
check_lhwr_Z <- function(Z, tol = 1e-6) {
  if (is.data.frame(Z)) {
    Z <- as.matrix(Z)
  }
  if (!is.matrix(Z) || !is.numeric(Z)) {
    stop("`Z` must be a numeric matrix or data frame.", call. = FALSE)
  }
  if (ncol(Z) != 3L) {
    stop("`Z` must have exactly 3 columns: z0, z1, and z2.", call. = FALSE)
  }
  if (!all(is.finite(Z))) {
    stop("`Z` must not contain missing, NaN, or infinite values.", call. = FALSE)
  }
  ok <- check_hyperboloid(Z, tol = tol, stop_on_fail = FALSE)
  if (!all(ok)) {
    bad <- which(!ok)
    stop(
      "`Z` is not on the upper Lorentz hyperboloid for row(s): ",
      paste(utils::head(bad, 5L), collapse = ", "),
      if (length(bad) > 5L) ", ...",
      ". Required: -z0^2 + z1^2 + z2^2 = -1 and z0 > 0.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' @rdname make_lhwr_Z
#' @export
make_lhwr_Z_features <- function(X, scale = TRUE, tol = 1e-6) {
  X <- lhwr_validate_two_component_input(X, name = "X")
  feature_names <- colnames(X)
  if (isTRUE(scale)) {
    X <- scale(X)
    X[!is.finite(X)] <- 0
  } else if (!identical(scale, FALSE)) {
    stop("`scale` must be TRUE or FALSE.", call. = FALSE)
  }

  Z <- cbind(
    z0 = sqrt(1 + X[, 1]^2 + X[, 2]^2),
    z1 = X[, 1],
    z2 = X[, 2]
  )
  storage.mode(Z) <- "double"
  check_lhwr_Z(Z, tol = tol)
  attr(Z, "lhwr_Z_source") <- "features"
  attr(Z, "feature_names") <- feature_names
  Z
}

#' @rdname make_lhwr_Z
#' @export
make_lhwr_Z_network <- function(A,
                                method = c("shortest_path_mds", "spectral"),
                                scale = TRUE,
                                disconnected = c("stop", "impute"),
                                symmetrize = TRUE,
                                tol = 1e-6) {
  method <- match.arg(method)
  disconnected <- match.arg(disconnected)
  A <- lhwr_validate_relation_matrix(A, symmetrize = symmetrize)

  embedding <- switch(
    method,
    shortest_path_mds = lhwr_network_shortest_path_mds(A, disconnected = disconnected),
    spectral = lhwr_network_spectral_embedding(A)
  )
  colnames(embedding) <- c("component1", "component2")

  Z <- make_lhwr_Z_features(embedding, scale = scale, tol = tol)
  attr(Z, "lhwr_Z_source") <- "network"
  attr(Z, "network_method") <- method
  Z
}

#' @rdname make_lhwr_Z
#' @export
make_lhwr_Z <- function(x, method = c("features", "network"), ...) {
  method <- match.arg(method)
  switch(
    method,
    features = make_lhwr_Z_features(x, ...),
    network = make_lhwr_Z_network(x, ...)
  )
}

#' @rdname make_lhwr_Z
#' @export
check_lorentz_Z <- function(Z, tol = 1e-6) {
  check_lhwr_Z(Z, tol = tol)
}

#' @rdname make_lhwr_Z
#' @export
make_lorentz_Z_features <- function(X, scale = TRUE, tol = 1e-6) {
  make_lhwr_Z_features(X = X, scale = scale, tol = tol)
}

#' @rdname make_lhwr_Z
#' @export
make_lorentz_Z_network <- function(A,
                                   method = c("shortest_path_mds", "spectral"),
                                   scale = TRUE,
                                   disconnected = c("stop", "impute"),
                                   symmetrize = TRUE,
                                   tol = 1e-6) {
  make_lhwr_Z_network(
    A = A,
    method = method,
    scale = scale,
    disconnected = disconnected,
    symmetrize = symmetrize,
    tol = tol
  )
}

#' @rdname make_lhwr_Z
#' @export
make_lorentz_Z <- function(x, method = c("features", "network"), ...) {
  make_lhwr_Z(x = x, method = method, ...)
}

lhwr_validate_two_component_input <- function(X, name = "X") {
  if (is.data.frame(X)) {
    numeric_cols <- vapply(X, is.numeric, logical(1))
    if (!all(numeric_cols)) {
      stop("`", name, "` must contain only numeric columns.", call. = FALSE)
    }
    X <- as.matrix(X)
  }
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("`", name, "` must be a numeric matrix or data frame.", call. = FALSE)
  }
  if (ncol(X) != 2L) {
    stop("`", name, "` must have exactly 2 columns.", call. = FALSE)
  }
  if (!all(is.finite(X))) {
    stop("`", name, "` must not contain missing, NaN, or infinite values.", call. = FALSE)
  }
  storage.mode(X) <- "double"
  X
}

lhwr_validate_relation_matrix <- function(A, symmetrize = TRUE) {
  if (is.data.frame(A)) {
    A <- as.matrix(A)
  }
  if (!is.matrix(A) || !is.numeric(A)) {
    stop("`A` must be a numeric square matrix.", call. = FALSE)
  }
  if (nrow(A) != ncol(A)) {
    stop("`A` must be a square matrix.", call. = FALSE)
  }
  if (nrow(A) < 3L) {
    stop("`A` must contain at least 3 observations.", call. = FALSE)
  }
  if (!all(is.finite(A))) {
    stop("`A` must not contain missing, NaN, or infinite values.", call. = FALSE)
  }
  if (any(A < 0)) {
    stop("`A` must contain nonnegative relation or similarity values.", call. = FALSE)
  }
  if (!isTRUE(all.equal(A, t(A), tolerance = 1e-10))) {
    if (!isTRUE(symmetrize)) {
      stop("`A` must be symmetric, or use `symmetrize = TRUE`.", call. = FALSE)
    }
    warning("`A` is not symmetric; using (A + t(A)) / 2.", call. = FALSE)
    A <- (A + t(A)) / 2
  }
  diag(A) <- 0
  storage.mode(A) <- "double"
  A
}

lhwr_network_spectral_embedding <- function(A) {
  deg <- rowSums(A)
  if (sum(deg > 0) < 3L) {
    stop("`A` does not contain enough connected information for spectral embedding.", call. = FALSE)
  }
  L <- diag(deg, nrow(A)) - A
  eig <- eigen(L, symmetric = TRUE)
  ord <- order(eig$values, decreasing = FALSE)
  if (length(ord) < 3L) {
    stop("Spectral embedding requires at least 3 observations.", call. = FALSE)
  }
  eig$vectors[, ord[2:3], drop = FALSE]
}

lhwr_network_shortest_path_mds <- function(A, disconnected = c("stop", "impute")) {
  disconnected <- match.arg(disconnected)
  n <- nrow(A)
  D <- matrix(Inf, n, n)
  diag(D) <- 0
  positive <- A > 0
  D[positive] <- 1 / A[positive]

  for (k in seq_len(n)) {
    candidate <- outer(D[, k], D[k, ], "+")
    D <- pmin(D, candidate)
  }

  finite_offdiag <- is.finite(D[row(D) != col(D)])
  if (!all(finite_offdiag)) {
    if (identical(disconnected, "stop")) {
      stop(
        "`A` is disconnected; shortest-path MDS requires a connected graph. ",
        "Use `disconnected = \"impute\"` or `method = \"spectral\"`.",
        call. = FALSE
      )
    }
    finite_values <- D[is.finite(D) & D > 0]
    if (!length(finite_values)) {
      stop("`A` has no positive relations to embed.", call. = FALSE)
    }
    replacement <- 1.25 * max(finite_values)
    warning("`A` is disconnected; replacing infinite shortest-path distances by ", signif(replacement, 4), ".", call. = FALSE)
    D[!is.finite(D)] <- replacement
  }

  fit <- stats::cmdscale(stats::as.dist(D), k = 2, eig = TRUE)
  points <- fit$points
  if (is.null(points) || ncol(points) == 0L) {
    stop("Classical MDS did not produce an embedding.", call. = FALSE)
  }
  if (ncol(points) == 1L) {
    points <- cbind(points[, 1], 0)
  }
  points[, 1:2, drop = FALSE]
}

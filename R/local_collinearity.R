#' Local Multicollinearity Diagnostics for LHWR
#'
#' Diagnose local multicollinearity in a fitted LHWR model. The
#' diagnostics follow the classical GWR weighted linear algebra discussed for
#' local collinearity diagnostics, including Wheeler and Tiefelsdorf (2005), but
#' use the LHWR kernel weights induced by the fitted distance geometry.
#'
#' @param object A fitted model object.
#' @param ... Additional arguments passed to methods.
#' @return A local collinearity diagnostic object.
#' @references
#' Wheeler, D. C. and Tiefelsdorf, M. (2005). Multicollinearity and correlation
#' among local regression coefficients in LHWR and related local weighted regression models.
#' \emph{Journal of Geographical Systems}, 7, 161--187.
#' @examples
#' set.seed(31)
#' H <- simulate_H2(60, r_max = 3)
#' x1 <- rnorm(60)
#' x2 <- 0.7 * x1 + rnorm(60, sd = 0.4)
#' y <- 1 + 2 * x1 - 0.5 * x2 + rnorm(60, sd = 0.3)
#' fit <- lhwr(y ~ x1 + x2, data.frame(y, x1, x2), Z = H$Z, bw = 25)
#' lc <- local_collinearity(fit)
#' lc
#' @export
local_collinearity <- function(object, ...) {
  UseMethod("local_collinearity")
}

#' @rdname local_collinearity
#' @param threshold_cn Threshold used to flag large local condition numbers.
#' @param threshold_vif Threshold used to flag large local VIF values.
#' @param center Logical. If `TRUE`, locally center predictor columns before
#'   computing condition numbers and local correlations.
#' @param scale Logical. If `TRUE`, locally scale predictor columns before
#'   computing condition numbers and local correlations.
#' @param return_pairwise Logical. If `TRUE`, return all pairwise local weighted
#'   correlations in long format.
#' @param engine Computational engine. `"cpp"` uses a compiled kNN
#'   implementation when cached kNN neighborhoods are stored in the fitted
#'   object. `"R"` uses the reference R implementation. `"auto"` selects
#'   `"cpp"` when possible and otherwise falls back to `"R"`.
#' @param digits Number of digits used by the print method.
#' @export
local_collinearity.lhwr <- function(object,
                                    threshold_cn = 30,
                                    threshold_vif = 10,
                                    center = TRUE,
                                    scale = TRUE,
                                    return_pairwise = FALSE,
                                    engine = c("auto", "cpp", "R"),
                                    digits = 3,
                                    ...) {
  engine <- match.arg(engine)
  if (!is_lhwr_object(object)) {
    stop("`object` must inherit from class 'lhwr'.", call. = FALSE)
  }
  X <- local_collinearity_model_matrix(object)
  colnames(X) <- clean_coef_names(colnames(X))
  non_intercept <- colnames(X) != "Intercept"
  Xp <- X[, non_intercept, drop = FALSE]
  variable_names <- colnames(Xp)
  if (!ncol(Xp)) {
    stop("Local collinearity diagnostics require at least one non-intercept predictor.", call. = FALSE)
  }

  n <- nrow(Xp)
  p <- ncol(Xp)
  use_cpp <- local_collinearity_use_cpp(object, engine = engine, return_pairwise = return_pairwise)
  components <- if (use_cpp) {
    local_collinearity_cpp_components(object, Xp, center = center, scale = scale)
  } else {
    local_collinearity_R_components(object, Xp, center = center, scale = scale, return_pairwise = return_pairwise)
  }
  cn_values <- components$local_CN
  neff <- components$neff
  n_positive <- components$n_positive
  vif_matrix <- components$local_VIF
  dimnames(vif_matrix) <- list(NULL, variable_names)
  max_abs_corr_values <- components$max_abs_local_correlation

  local_vif <- data.frame(
    observation_id = rep(seq_len(n), each = p),
    variable = rep(variable_names, times = n),
    local_VIF = as.vector(t(vif_matrix)),
    flag_VIF = as.vector(t(vif_matrix)) > threshold_vif,
    check.names = FALSE
  )
  local_vif$flag_VIF[is.na(local_vif$flag_VIF)] <- FALSE

  out <- list(
    condition_number = data.frame(
      observation_id = seq_len(n),
      local_CN = cn_values,
      flag_CN = cn_values > threshold_cn,
      check.names = FALSE
    ),
    local_vif = local_vif,
    local_correlation_summary = data.frame(
      observation_id = seq_len(n),
      max_abs_local_correlation = max_abs_corr_values,
      check.names = FALSE
    ),
    local_correlation_pairs = if (isTRUE(return_pairwise) && p > 1L) {
      components$local_correlation_pairs
    } else {
      NULL
    },
    thresholds = list(
      condition_number = threshold_cn,
      vif = threshold_vif
    ),
    diagnostics = list(
      n = n,
      p = p,
      variables = variable_names,
      bandwidth = object$bandwidth,
      kernel = object$LHWR.arguments$kernel,
      distance = object$distance_label %||% object$LHWR.arguments$distance,
      geometry = object$LHWR.arguments$geometry,
      bandwidth_type = object$LHWR.arguments$bandwidth_type,
      adaptive = object$LHWR.arguments$adaptive,
      neighbor_count = object$LHWR.arguments$bw,
      engine = object$LHWR.arguments$engine,
      collinearity_engine = if (use_cpp) "cpp" else "R",
      center = center,
      scale = scale,
      neff = neff,
      n_positive = n_positive
    ),
    digits = digits
  )
  class(out) <- "local_collinearity_lhwr"
  out
}

#' @export
local_collinearity.hgwr <- local_collinearity.lhwr
#' @export
local_collinearity.lorentzhgwr <- local_collinearity.lhwr

#' @export
print.local_collinearity_lhwr <- function(x, digits = x$digits %||% 3, ...) {
  cn <- x$condition_number$local_CN
  cn_flags <- x$condition_number$flag_CN
  vif <- x$local_vif
  corr <- x$local_correlation_summary$max_abs_local_correlation
  diagnostics <- x$diagnostics

  hgwr_rule()
  cat(" LHWR Local Multicollinearity Diagnostics\n")
  hgwr_rule()
  cat(" Distance metric     : ", diagnostics$distance, "\n", sep = "")
  cat(" Kernel              : ", diagnostics$kernel, "\n", sep = "")
  cat(" Bandwidth type      : ", local_collinearity_bandwidth_type_label(diagnostics), "\n", sep = "")
  if (isTRUE(diagnostics$adaptive) && identical(diagnostics$bandwidth_type, "knn")) {
    cat(" Neighbor count      : ", local_collinearity_neighbor_count(diagnostics), "\n", sep = "")
    cat(" Median radius       : ", local_collinearity_bandwidth_label(diagnostics$bandwidth), "\n", sep = "")
  } else {
    cat(" Bandwidth           : ", local_collinearity_bandwidth_label(diagnostics$bandwidth), "\n", sep = "")
  }

  cat("\n Condition number summary\n")
  hgwr_rule()
  cn_out <- data.frame(
    Metric = c("Min", "Median", "Mean", "Max", paste0("Locations CN > ", x$thresholds$condition_number)),
    Value = c(
      format_diag_number(min(cn, na.rm = TRUE), digits),
      format_diag_number(stats::median(cn, na.rm = TRUE), digits),
      format_diag_number(mean(cn, na.rm = TRUE), digits),
      format_diag_number(max(cn, na.rm = TRUE), digits),
      sum(cn_flags, na.rm = TRUE)
    ),
    check.names = FALSE
  )
  print(cn_out, row.names = FALSE, right = FALSE)

  cat("\n Local VIF summary\n")
  hgwr_rule()
  vif_summary <- do.call(rbind, lapply(split(vif, vif$variable), function(d) {
    data.frame(
      Variable = unique(d$variable),
      Min = format_diag_number(min(d$local_VIF, na.rm = TRUE), digits),
      Median = format_diag_number(stats::median(d$local_VIF, na.rm = TRUE), digits),
      Mean = format_diag_number(mean(d$local_VIF, na.rm = TRUE), digits),
      Max = format_diag_number(max(d$local_VIF, na.rm = TRUE), digits),
      `# Locations VIF > threshold` = sum(d$flag_VIF, na.rm = TRUE),
      check.names = FALSE
    )
  }))
  names(vif_summary)[ncol(vif_summary)] <- paste0("# Locations VIF > ", x$thresholds$vif)
  print(vif_summary, row.names = FALSE, right = FALSE)

  cat("\n Local weighted correlation summary\n")
  hgwr_rule()
  corr_out <- data.frame(
    Metric = c("Median max |corr|", "Max max |corr|"),
    Value = c(
      format_diag_number(stats::median(corr, na.rm = TRUE), digits),
      format_diag_number(max(corr, na.rm = TRUE), digits)
    ),
    check.names = FALSE
  )
  print(corr_out, row.names = FALSE, right = FALSE)

  cat("\n Interpretation\n")
  hgwr_rule()
  cat(" CN > ", x$thresholds$condition_number, " suggests serious local multicollinearity.\n", sep = "")
  cat(" VIF > ", x$thresholds$vif, " suggests coefficient-specific local collinearity.\n", sep = "")
  hgwr_rule()
  invisible(x)
}

local_collinearity_model_matrix <- function(object) {
  if (!is.null(object$lm) && inherits(object$lm, "lm")) {
    return(stats::model.matrix(object$lm))
  }
  stop("Could not reconstruct the model matrix from `object`.", call. = FALSE)
}

local_collinearity_use_cpp <- function(object, engine, return_pairwise) {
  if (engine == "R") {
    return(FALSE)
  }
  if (isTRUE(return_pairwise)) {
    if (engine == "cpp") {
      warning("`return_pairwise = TRUE` uses the R reference path for pairwise correlations.", call. = FALSE)
    }
    return(FALSE)
  }
  ok <- identical(object$LHWR.arguments$engine, "knn") &&
    !is.null(object$knn_index) &&
    !is.null(object$knn_distance)
  if (engine == "cpp" && !ok) {
    stop("`engine = \"cpp\"` for local collinearity requires cached kNN neighborhoods from `lhwr(..., engine = \"knn\")`.", call. = FALSE)
  }
  ok
}

local_collinearity_cpp_components <- function(object, Xp, center, scale) {
  out <- .local_collinearity_knn_cpp(
    X = Xp,
    knn_index = object$knn_index,
    knn_distance = object$knn_distance,
    bandwidth = object$bandwidth,
    kernel = object$LHWR.arguments$kernel,
    center = isTRUE(center),
    scale = isTRUE(scale),
    n_threads = object$LHWR.arguments$n_threads %||% 1L
  )
  list(
    local_CN = as.numeric(out$local_CN),
    local_VIF = as.matrix(out$local_VIF),
    max_abs_local_correlation = as.numeric(out$max_abs_local_correlation),
    neff = as.numeric(out$neff),
    n_positive = as.integer(out$n_positive),
    local_correlation_pairs = NULL
  )
}

local_collinearity_R_components <- function(object, Xp, center, scale, return_pairwise) {
  weights <- local_collinearity_weight_source(object)
  n <- nrow(Xp)
  p <- ncol(Xp)
  cn_values <- rep(NA_real_, n)
  neff <- rep(NA_real_, n)
  n_positive <- rep(NA_integer_, n)
  vif_matrix <- matrix(NA_real_, nrow = n, ncol = p)
  max_abs_corr_values <- rep(NA_real_, n)
  pairwise_rows <- if (isTRUE(return_pairwise) && p > 1L) vector("list", n) else NULL

  for (i in seq_len(n)) {
    w <- weights$get(i)
    positive <- w > 0 & is.finite(w)
    n_positive[i] <- sum(positive)
    sw <- sum(w)
    sw2 <- sum(w^2)
    neff[i] <- if (is.finite(sw) && is.finite(sw2) && sw2 > 0) sw^2 / sw2 else 0

    X_local <- local_standardize_matrix(Xp, w, center = center, scale = scale)
    G <- crossprod(X_local, X_local * w)
    cn_values[i] <- safe_condition_number(G)

    if (p > 1L) {
      C <- weighted_cor(Xp, w, center = center, scale = scale)
      off_diag <- abs(C[upper.tri(C)])
      max_abs_corr_values[i] <- if (length(off_diag)) max(off_diag, na.rm = TRUE) else NA_real_
      if (isTRUE(return_pairwise)) {
        pairwise_rows[[i]] <- local_correlation_long(C, i)
      }
    } else {
      max_abs_corr_values[i] <- NA_real_
    }

    for (j in seq_len(p)) {
      if (p == 1L) {
        vif_matrix[i, j] <- NA_real_
      } else {
        aux_y <- Xp[, j]
        aux_X <- cbind(Intercept = 1, Xp[, -j, drop = FALSE])
        r2 <- weighted_r2(aux_y, aux_X, w)
        vif_matrix[i, j] <- if (!is.finite(r2) || r2 >= 1 - 1e-10) {
          Inf
        } else {
          1 / (1 - r2)
        }
      }
    }
  }

  list(
    local_CN = cn_values,
    local_VIF = vif_matrix,
    max_abs_local_correlation = max_abs_corr_values,
    neff = neff,
    n_positive = n_positive,
    local_correlation_pairs = if (isTRUE(return_pairwise) && p > 1L) {
      do.call(rbind, pairwise_rows)
    } else {
      NULL
    }
  )
}

local_collinearity_weight_source <- function(object) {
  if (
    identical(object$LHWR.arguments$engine, "knn") &&
      !is.null(object$knn_index) &&
      !is.null(object$knn_distance)
  ) {
    index <- as.matrix(object$knn_index)
    distance <- as.matrix(object$knn_distance)
    bw <- as.numeric(object$bandwidth)
    kernel <- object$LHWR.arguments$kernel
    n <- nrow(index)
    if (nrow(distance) == n && ncol(distance) == ncol(index) && length(bw) == n) {
      return(list(
        get = function(i) {
          w <- numeric(n)
          w[i] <- hgwr_kernel(0, bw[i], kernel = kernel)
          idx <- index[i, ]
          dij <- distance[i, ]
          ok <- is.finite(idx) & idx >= 1L & idx <= n & is.finite(dij)
          if (any(ok)) {
            w[idx[ok]] <- hgwr_kernel(dij[ok], bw[i], kernel = kernel)
          }
          w
        }
      ))
    }
  }

  D <- local_collinearity_distance_matrix(object)
  bw <- as.numeric(object$bandwidth)
  if (length(bw) == 1L) {
    bw <- rep(bw, nrow(D))
  }
  if (length(bw) != nrow(D)) {
    stop("Stored bandwidths are not compatible with the model size.", call. = FALSE)
  }
  kernel <- object$LHWR.arguments$kernel
  list(
    get = function(i) {
      hgwr_kernel(D[i, ], bw[i], kernel = kernel)
    }
  )
}

local_collinearity_distance_matrix <- function(object) {
  if (
    !is.null(object$distance) &&
      is.matrix(object$distance) &&
      nrow(object$distance) == nrow(object$Z) &&
      ncol(object$distance) == nrow(object$Z)
  ) {
    return(object$distance)
  }
  if (identical(object$LHWR.arguments$geometry, "tangent")) {
    return(tangent_distance_matrix(object$Z, mu = object$tangent_ref))
  }
  hyperbolic_distance_matrix(object$Z, check = FALSE)
}

local_standardize_matrix <- function(X, w, center = TRUE, scale = TRUE) {
  X <- as.matrix(X)
  if (isTRUE(center)) {
    mu <- weighted_col_means(X, w)
    X <- sweep(X, 2L, mu, "-")
  }
  if (isTRUE(scale)) {
    ss <- sqrt(weighted_col_means(X^2, w))
    ss[!is.finite(ss) | ss <= .Machine$double.eps] <- 1
    X <- sweep(X, 2L, ss, "/")
  }
  X
}

weighted_col_means <- function(X, w) {
  sw <- sum(w)
  if (!is.finite(sw) || sw <= .Machine$double.eps) {
    return(rep(NA_real_, ncol(X)))
  }
  as.numeric(crossprod(w, X) / sw)
}

#' Weighted Correlation Matrix
#'
#' Compute a weighted correlation matrix for numeric columns.
#'
#' @param X Numeric matrix.
#' @param w Nonnegative weight vector.
#' @param center Logical. If `TRUE`, center columns by weighted means.
#' @param scale Logical. If `TRUE`, scale columns by weighted standard
#'   deviations.
#' @return Weighted correlation matrix.
weighted_cor <- function(X, w, center = TRUE, scale = TRUE) {
  Xs <- local_standardize_matrix(X, w, center = center, scale = scale)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= .Machine$double.eps) {
    return(matrix(NA_real_, ncol(X), ncol(X), dimnames = list(colnames(X), colnames(X))))
  }
  C <- crossprod(Xs, Xs * w) / sw
  d <- sqrt(diag(C))
  d[!is.finite(d) | d <= .Machine$double.eps] <- NA_real_
  C <- C / tcrossprod(d)
  diag(C) <- 1
  dimnames(C) <- list(colnames(X), colnames(X))
  C
}

#' Weighted R-squared
#'
#' Compute weighted R-squared for a weighted least-squares auxiliary regression.
#'
#' @param y Numeric response vector.
#' @param X Numeric design matrix, usually including an intercept.
#' @param w Nonnegative weight vector.
#' @return Weighted R-squared, or `NA_real_` if it cannot be computed.
weighted_r2 <- function(y, X, w) {
  y <- as.numeric(y)
  X <- as.matrix(X)
  w <- as.numeric(w)
  positive <- w > 0 & is.finite(w)
  if (sum(positive) <= ncol(X)) {
    return(NA_real_)
  }
  sw <- sum(w)
  if (!is.finite(sw) || sw <= .Machine$double.eps) {
    return(NA_real_)
  }
  y_bar <- sum(w * y) / sw
  tss <- sum(w * (y - y_bar)^2)
  if (!is.finite(tss) || tss <= .Machine$double.eps) {
    return(NA_real_)
  }
  Xw <- X * w
  XtWX <- crossprod(X, Xw)
  XtWy <- crossprod(X, w * y)
  beta <- tryCatch(
    solve(XtWX, XtWy),
    error = function(e) MASS::ginv(XtWX) %*% XtWy
  )
  yhat <- as.numeric(X %*% beta)
  rss <- sum(w * (y - yhat)^2)
  r2 <- 1 - rss / tss
  if (!is.finite(r2)) NA_real_ else max(0, min(1, r2))
}

#' Safe Condition Number
#'
#' Compute a condition number from a symmetric cross-product matrix.
#'
#' @param G Symmetric positive semidefinite matrix.
#' @param tol Eigenvalue tolerance.
#' @return Numeric condition number, or `Inf` if the matrix is singular.
safe_condition_number <- function(G, tol = 1e-10) {
  ev <- tryCatch(
    eigen(G, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) NA_real_
  )
  if (any(!is.finite(ev))) return(Inf)
  max_ev <- max(ev)
  min_ev <- min(ev)
  if (min_ev <= tol || max_ev <= tol) return(Inf)
  sqrt(max_ev / min_ev)
}

local_correlation_long <- function(C, observation_id) {
  if (ncol(C) < 2L) return(NULL)
  ij <- which(upper.tri(C), arr.ind = TRUE)
  data.frame(
    observation_id = observation_id,
    variable_1 = rownames(C)[ij[, 1]],
    variable_2 = colnames(C)[ij[, 2]],
    local_correlation = C[ij],
    check.names = FALSE
  )
}

local_collinearity_bandwidth_type_label <- function(diagnostics) {
  if (isTRUE(diagnostics$adaptive)) {
    if (identical(diagnostics$bandwidth_type, "knn")) {
      return("adaptive kNN")
    }
    return("adaptive volume")
  }
  "fixed"
}

local_collinearity_bandwidth_label <- function(bw) {
  if (length(bw) == 1L) {
    return(format(signif(bw, 6)))
  }
  format(signif(stats::median(bw, na.rm = TRUE), 6))
}

local_collinearity_neighbor_count <- function(diagnostics) {
  k <- diagnostics$neighbor_count
  if (is.null(k) || length(k) != 1L || !is.finite(k)) {
    return("not stored")
  }
  as.character(as.integer(round(k)))
}

format_diag_number <- function(x, digits) {
  ifelse(is.na(x), "NA", ifelse(is.infinite(x), "Inf", formatC(x, digits = digits, format = "f")))
}

#' @export
print.local_collinearity_hgwr <- print.local_collinearity_lhwr

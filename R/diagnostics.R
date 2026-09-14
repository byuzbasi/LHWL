#' LHWR AICc
#'
#' @param object An LHWR model object returned by [lhwr()] or a low-level
#'   fit returned by [hgwr.fit()].
#' @return A single numeric value containing the corrected Akaike information
#'   criterion. Smaller values indicate a better fitted model under this
#'   criterion.
#' @details
#' This extractor does not refit or select a model. AICc requires hat-matrix
#' diagnostics; a nonpositive n - tr(S) - 2 denominator produces infinity.
#' Gaussian AIC and BIC are also stored in `object$LHWR.diagnostic`, but
#' `bw.lhwr()` directly selects only by AICc or CV. See
#' `vignette("lhwr-diagnostics-and-tuning", package = "LHWL")` for formulas.
#' @export
AICc.lhwr <- function(object) {
  if (identical(object$family %||% object$LHWR.arguments$family, "binomial")) {
    stop(
      "AICc is not currently defined for binomial LHWR; use binomial leave-one-out deviance through `bw.lhwr(..., approach = \"CV\")`.",
      call. = FALSE
    )
  }
  if (!is.null(object$LHWR.diagnostic$AICc)) {
    return(object$LHWR.diagnostic$AICc)
  }
  if (!is.null(object$AICc)) object$AICc else NA_real_
}

#' LHWR Cross-Validation Score
#'
#' @param object An LHWR bandwidth or model object.
#' @return A single numeric leave-one-out cross-validation score if it was
#'   computed for `object`; otherwise `NA_real_`.
#' @details
#' Gaussian CV is a sum of squared leave-one-out errors; binomial tuning uses
#' Bernoulli log loss. The representation and candidate radii are held fixed
#' for internal bandwidth scoring. This extractor does not run cross-validation
#' and a stored score is not an independent outer evaluation after tuning.
#' @export
CV.lhwr <- function(object) {
  if (!is.null(object$CV)) return(object$CV)
  if (!is.null(object$LHWR.diagnostic$CV)) return(object$LHWR.diagnostic$CV)
  NA_real_
}

#' LHWR Residual Sum of Squares
#'
#' @param object An LHWR model object.
#' @return A single numeric residual sum of squares computed from the fitted
#'   values and observed response.
#' @export
RSS.lhwr <- function(object) {
  if (identical(object$family %||% object$LHWR.arguments$family, "binomial")) {
    stop(
      "Residual sum of squares is a Gaussian LHWR diagnostic; use binomial deviance or Brier score instead.",
      call. = FALSE
    )
  }
  if (!is.null(object$LHWR.diagnostic$RSS.gw)) {
    return(object$LHWR.diagnostic$RSS.gw)
  }
  object$RSS
}

hgwr_diagnostics <- function(y, X, beta, y_hat, residuals, s_diag,
                             s2_diag = NULL,
                             local_kappa = NULL,
                             local_neff = NULL,
                             local_n_positive = NULL) {
  n <- length(y)
  p <- ncol(X)
  rss <- sum(residuals^2)
  tss <- sum((y - mean(y))^2)
  trS <- sum(s_diag)
  trStS <- if (is.null(s2_diag)) sum(s_diag^2) else sum(s2_diag)
  enp <- 2 * trS - trStS
  edf <- n - enp
  sigma2 <- rss / max(edf, .Machine$double.eps)
  dev <- n * log(rss / n) + n * log(2 * pi) + n
  AIC <- dev + 2 * (trS + 1)
  AICc_denom <- n - trS - 2
  AICc <- if (AICc_denom > 0) {
    dev + 2 * n * (trS + 1) / AICc_denom
  } else {
    Inf
  }
  BIC <- dev + log(n) * (trS + 1)
  R2 <- if (tss > 0) 1 - rss / tss else NA_real_
  adjR2 <- if (is.finite(R2)) 1 - (1 - R2) * (n - 1) / max(edf, 1) else NA_real_
  
  out <- list(
    RSS.gw = rss,
    AIC = AIC,
    AICc = AICc,
    BIC = BIC,
    tr.S = trS,
    tr.StS = trStS,
    enp = enp,
    edf = edf,
    sigma2 = sigma2,
    gw.R2 = R2,
    gwR2.adj = adjR2
  )
  
  if (!is.null(local_kappa)) {
    out$kappa.median <- stats::median(local_kappa, na.rm = TRUE)
    out$kappa.max <- max(local_kappa, na.rm = TRUE)
  }
  if (!is.null(local_neff)) {
    out$neff.min <- min(local_neff, na.rm = TRUE)
    out$neff.median <- stats::median(local_neff, na.rm = TRUE)
  }
  if (!is.null(local_n_positive)) {
    out$n.positive.min <- min(local_n_positive, na.rm = TRUE)
    out$n.positive.median <- stats::median(local_n_positive, na.rm = TRUE)
  }
  
  out
}

hgwr_minimal_diagnostics <- function(y, y_hat, residuals) {
  rss <- sum(residuals^2)
  tss <- sum((y - mean(y))^2)
  R2 <- if (tss > 0) 1 - rss / tss else NA_real_
  list(
    RSS.gw = rss,
    gw.R2 = R2
  )
}

local_r2 <- function(y, y_hat, D, bw, kernel) {
  W <- hgwr_kernel(D, bw, kernel = kernel)
  y_bar <- as.numeric(crossprod(W, y) / pmax(colSums(W), .Machine$double.eps))
  tss <- colSums(W * (y - matrix(y_bar, nrow(D), ncol(D), byrow = TRUE))^2)
  rss <- colSums(W * (y - y_hat)^2)
  out <- 1 - rss / tss
  out[!is.finite(out)] <- NA_real_
  out
}

local_r2_knn <- function(y, y_hat, knn_index, knn_distance, bandwidth, kernel) {
  y <- as.numeric(y)
  y_hat <- as.numeric(y_hat)
  knn_index <- as.matrix(knn_index)
  knn_distance <- as.matrix(knn_distance)
  n <- length(y)
  if (nrow(knn_index) != n || nrow(knn_distance) != n || ncol(knn_index) != ncol(knn_distance)) {
    return(rep(NA_real_, n))
  }
  if (length(bandwidth) == 1L) bandwidth <- rep(bandwidth, n)
  out <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    idx <- as.integer(knn_index[i, ])
    dij <- as.numeric(knn_distance[i, ])
    ok <- is.finite(idx) & idx >= 1L & idx <= n & is.finite(dij)
    idx <- idx[ok]
    dij <- dij[ok]
    bw_i <- bandwidth[i]
    if (!is.finite(bw_i) || bw_i <= 0) next
    w <- c(1, hgwr_kernel(dij, bw_i, kernel = kernel))
    yy <- c(y[i], y[idx])
    yy_hat <- c(y_hat[i], y_hat[idx])
    good <- is.finite(w) & w > 0 & is.finite(yy) & is.finite(yy_hat)
    if (!any(good)) next
    w <- w[good]
    yy <- yy[good]
    yy_hat <- yy_hat[good]
    y_bar <- sum(w * yy) / sum(w)
    tss <- sum(w * (yy - y_bar)^2)
    rss <- sum(w * (yy - yy_hat)^2)
    out[i] <- if (tss > .Machine$double.eps) 1 - rss / tss else NA_real_
  }
  out
}

#' Summarize an LHWR Model
#'
#' @param object An LHWR model object.
#' @param ... Unused.
#' @return An object of class `"summary.hgwr"` containing the original call,
#'   model arguments, diagnostic values, and summary statistics for each local
#'   coefficient surface.
#' @export
summary.lhwr <- function(object, ...) {
  if (!is_lhwr_object(object)) {
    stop("`object` must inherit from class 'lhwr'.", call. = FALSE)
  }
  
  coef_cols <- object$coef_names
  coef_summary <- t(vapply(
    object$SDF[, coef_cols, drop = FALSE],
    summary,
    numeric(6L)
  ))
  
  ans <- list(
    call = object$this.call,
    arguments = object$LHWR.arguments,
    family = object$family %||% object$LHWR.arguments$family %||% "gaussian",
    diagnostic = object$LHWR.diagnostic,
    coefficients = coef_summary
  )
  class(ans) <- "summary.lhwr"
  ans
}

#' @export
print.lhwr <- function(x, ...) {
  if (!is_lhwr_object(x)) {
    stop("`x` must inherit from class 'lhwr'.", call. = FALSE)
  }
  args <- x$LHWR.arguments
  diag <- x$LHWR.diagnostic
  n <- length(x$residuals)
  p <- length(x$coef_names)
  
  hgwr_rule()
  cat(" LHWR model fit\n")
  family <- x$family %||% args$family %||% "gaussian"
  cat(
    if (family == "binomial") {
      " Local logistic regression with geometry-aware kernel weights\n"
    } else {
      " Local linear regression with geometry-aware kernel weights\n"
    }
  )
  hgwr_rule()
  cat("Call:\n")
  print(x$this.call)
  
  info <- data.frame(
    Setting = c(
      "Model",
      "Family",
      "Geometry",
      "Distance",
      "Kernel",
      "Engine",
      "Diagnostics",
      "Bandwidth",
      "Bandwidth type",
      "Observations",
      "Coefficients"
    ),
    Value = c(
      "Lorentz Hyperbolic Weighted Regression",
      family,
      hgwr_geometry_label(args$geometry),
      args$distance,
      hgwr_kernel_label(args$kernel),
      hgwr_engine_label(args$engine),
      hgwr_diagnostics_label(args$diagnostics),
      hgwr_bandwidth_label(args),
      hgwr_bandwidth_type_label(args$bandwidth_type),
      n,
      p
    ),
    check.names = FALSE
  )
  cat("\nCalibration\n")
  print(info, row.names = FALSE, right = FALSE)
  
  cat("\nDiagnostics\n")
  print(hgwr_diagnostic_table(diag), row.names = FALSE, right = FALSE)
  
  cat("\nLocal coefficient summaries\n")
  cat("Rows summarize the estimated coefficient surface over all locations.\n")
  printCoefmat(hgwr_coef_summary(x), digits = 4, signif.stars = FALSE)

  if (!is.null(x$Ftests)) {
    cat("\n")
    print(x$Ftests)
  }
  
  if (!is.null(x$timings$stop) && !is.null(x$timings$start)) {
    elapsed <- as.numeric(difftime(x$timings$stop, x$timings$start, units = "secs"))
    cat("\nElapsed:", hgwr_format_number(elapsed), "seconds\n")
  }
  hgwr_rule()
  invisible(x)
}

#' @export
print.summary.lhwr <- function(x, ...) {
  cat("Summary of LHWR\n")
  family <- x$family %||% x$arguments$family %||% "gaussian"
  cat(
    if (family == "binomial") {
      "Geometry-aware local logistic regression on hyperboloid locations\n"
    } else {
      "Geometry-aware local linear regression on hyperboloid locations\n"
    }
  )
  hgwr_rule()
  cat("Call:\n")
  print(x$call)
  cat("\nCalibration\n")
  info <- data.frame(
    Setting = c("Family", "Geometry", "Distance", "Kernel", "Diagnostics", "Bandwidth", "Bandwidth type"),
    Value = c(
      family,
      hgwr_geometry_label(x$arguments$geometry),
      x$arguments$distance,
      hgwr_kernel_label(x$arguments$kernel),
      hgwr_diagnostics_label(x$arguments$diagnostics),
      hgwr_bandwidth_label(x$arguments),
      hgwr_bandwidth_type_label(x$arguments$bandwidth_type)
    ),
    check.names = FALSE
  )
  print(info, row.names = FALSE, right = FALSE)
  cat("\nDiagnostics\n")
  print(hgwr_diagnostic_table(x$diagnostic), row.names = FALSE, right = FALSE)
  cat("\nLocal coefficient summaries:\n")
  printCoefmat(x$coefficients, digits = 4, signif.stars = FALSE)
  hgwr_rule()
  invisible(x)
}

# Backward-compatible diagnostic aliases.
#' @rdname AICc.lhwr
#' @export
AICc.hgwr <- AICc.lhwr
#' @rdname CV.lhwr
#' @export
CV.hgwr <- CV.lhwr
#' @rdname RSS.lhwr
#' @export
RSS.hgwr <- RSS.lhwr
#' @rdname summary.lhwr
#' @export
summary.hgwr <- summary.lhwr
#' @export
print.hgwr <- print.lhwr
#' @export
print.summary.hgwr <- print.summary.lhwr

is_lhwr_object <- function(object) {
  inherits(object, "lhwr") || inherits(object, "hgwr")
}

hgwr_rule <- function() {
  cat(strrep("-", 72), "\n", sep = "")
}

hgwr_bandwidth_label <- function(args) {
  bw <- if (length(args$bw) == 1L) {
    hgwr_format_number(args$bw)
  } else {
    paste0("vector[n=", length(args$bw), "]")
  }
  if (isTRUE(args$adaptive)) {
    paste0(bw, " nearest neighbors")
  } else {
    paste0(bw, " fixed radius")
  }
}

hgwr_format_number <- function(x, digits = 6) {
  ifelse(
    is.na(x),
    NA_character_,
    format(signif(as.numeric(x), digits), scientific = FALSE, trim = TRUE)
  )
}

hgwr_diagnostic_table <- function(diag) {
  labels <- c(
    RSS.gw = "Residual sum of squares",
    deviance = "Binomial deviance",
    null.deviance = "Null binomial deviance",
    pseudo.R2 = "McFadden pseudo-R-squared",
    logLik = "Binomial log-likelihood",
    brier = "Brier score",
    converged = "Converged local fits",
    not.converged = "Non-converged local fits",
    class.variation.failures = "Single-class neighborhoods",
    rank.deficient = "Rank-deficient local Fisher matrices",
    AIC = "AIC",
    AICc = "AICc",
    BIC = "BIC",
    tr.S = "Trace of hat matrix",
    tr.StS = "Trace of S'S",
    enp = "Effective parameters",
    edf = "Effective degrees of freedom",
    sigma2 = "Residual variance",
    gw.R2 = "Local-model R-squared",
    gwR2.adj = "Adjusted local R-squared",
    kappa.median = "Median local condition number",
    kappa.max = "Maximum local condition number",
    neff.min = "Minimum effective sample size",
    neff.median = "Median effective sample size",
    n.positive.min = "Minimum positive kernel weights",
    n.positive.median = "Median positive kernel weights"
  )
  present <- intersect(names(labels), names(diag))
  data.frame(
    Metric = unname(labels[present]),
    Value = hgwr_format_number(unlist(diag[present])),
    check.names = FALSE
  )
}

hgwr_geometry_label <- function(geometry) {
  switch(
    geometry,
    lorentz = "Lorentz hyperboloid (intrinsic)",
    tangent = "Tangent-space approximation",
    geometry
  )
}

hgwr_kernel_label <- function(kernel) {
  paste0(kernel, " kernel")
}

hgwr_bandwidth_type_label <- function(bandwidth_type) {
  switch(
    bandwidth_type,
    knn = "Adaptive nearest-neighbor radius",
    volume = "Hyperbolic volume-adaptive radius",
    bandwidth_type
  )
}

hgwr_criterion_label <- function(approach) {
  switch(
    approach,
    AICc = "AICc, lower is better",
    CV = "Leave-one-out CV, lower is better",
    approach
  )
}

hgwr_search_label <- function(search) {
  switch(
    search,
    golden = "Golden-section search",
    grid = "Candidate grid search",
    search
  )
}

hgwr_engine_label <- function(engine) {
  if (is.null(engine)) {
    return("Dense distance matrix")
  }
  switch(
    engine,
    dense = "Dense distance matrix",
    knn = "Exact Lorentz kNN",
    engine
  )
}

hgwr_diagnostics_label <- function(diagnostics) {
  if (is.null(diagnostics)) {
    return("Full")
  }
  switch(
    diagnostics,
    full = "Full",
    aicc = "AICc only",
    none = "None",
    diagnostics
  )
}

hgwr_coef_summary <- function(x) {
  t(vapply(
    x$SDF[, x$coef_names, drop = FALSE],
    summary,
    numeric(6L)
  ))
}

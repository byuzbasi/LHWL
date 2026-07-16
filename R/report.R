#' Professional LHWR Model Report
#'
#' Print a structured report for an LHWR fitted model. The report is designed
#' for applied analysis: it combines the model call, calibration settings,
#' global OLS baseline, LHWR diagnostics, local coefficient summaries, optional
#' local inference, optional diagnostic F-tests, and optional residual or
#' collinearity diagnostics.
#'
#' The function does not refit the model. Local inference is computed only from
#' quantities already stored in the fitted object, and therefore requires a fit
#' created with `diagnostics = "full"`. Moran and local collinearity diagnostics
#' can be passed through the `moran` and `collinearity` arguments when they have
#' already been computed separately.
#'
#' @param object An LHWR model object returned by [lhwr()].
#' @param include_global Logical. If `TRUE`, include the global OLS baseline
#'   stored in `object$lm` when available.
#' @param include_inference Logical. If `TRUE`, include exploratory local
#'   coefficient inference when available or computable from the fit.
#' @param include_f_tests Logical. If `TRUE`, include diagnostic F-tests when
#'   stored in the fit or computable from the fit.
#' @param include_moran Logical. If `TRUE`, include the object supplied through
#'   `moran`. The function does not compute Moran's I automatically.
#' @param include_collinearity Logical. If `TRUE`, include the object supplied
#'   through `collinearity`. The function does not compute local collinearity
#'   diagnostics automatically.
#' @param local_inference Optional object returned by [local_inference()]. If
#'   `NULL` and `include_inference = TRUE`, the function attempts to compute it
#'   from `object`.
#' @param moran Optional object returned by [moran.lhwr()] or [moran.hgwr()].
#' @param collinearity Optional object returned by [local_collinearity()].
#' @param alpha Significance level used when local inference must be computed.
#' @param adjust Multiplicity adjustment used when local inference must be
#'   computed. Passed to [local_inference()].
#' @param digits Number of significant digits used in printed numeric values.
#' @param ... Unused.
#' @return Invisibly returns an object of class `"lhwr_report"` containing the
#'   tables printed in the report.
#' @examples
#' set.seed(7)
#' H <- simulate_H2(80, r_max = 3)
#' Z <- H$Z
#' x1 <- rnorm(80)
#' x2 <- rnorm(80)
#' y <- 1 + 2 * x1 - 0.5 * x2 + rnorm(80, sd = 0.3)
#' dat <- data.frame(y, x1, x2)
#' fit <- lhwr(y ~ x1 + x2, data = dat, Z = Z, bw = 35,
#'             adaptive = TRUE, diagnostics = "full", f_tests = TRUE)
#' lhwr_report(fit, include_moran = FALSE, include_collinearity = FALSE)
#' @export
lhwr_report <- function(object,
                        include_global = TRUE,
                        include_inference = TRUE,
                        include_f_tests = TRUE,
                        include_moran = !is.null(moran),
                        include_collinearity = !is.null(collinearity),
                        local_inference = NULL,
                        moran = NULL,
                        collinearity = NULL,
                        alpha = 0.05,
                        adjust = "BH",
                        digits = 4,
                        ...) {
  if (!is_lhwr_object(object)) {
    stop("`object` must be an LHWR model object.", call. = FALSE)
  }
  if (!is.numeric(digits) || length(digits) != 1L || !is.finite(digits) || digits < 1) {
    stop("`digits` must be a positive number.", call. = FALSE)
  }

  args <- object$LHWR.arguments
  out <- list(
    call = object$this.call,
    calibration = lhwr_report_calibration_table(object),
    global = if (isTRUE(include_global)) lhwr_report_global_table(object) else NULL,
    diagnostics = hgwr_diagnostic_table(object$LHWR.diagnostic),
    coefficients = hgwr_coef_summary(object),
    local_inference = NULL,
    f_tests = NULL,
    moran = if (isTRUE(include_moran)) moran else NULL,
    collinearity = if (isTRUE(include_collinearity)) collinearity else NULL,
    elapsed = lhwr_report_fit_elapsed(object),
    notes = character()
  )

  if (isTRUE(include_inference)) {
    if (is.null(local_inference)) {
      local_inference <- tryCatch(
        local_inference(object, alpha = alpha, adjust = adjust),
        error = function(e) e
      )
    }
    if (inherits(local_inference, "error")) {
      out$notes <- c(out$notes, paste0("Local inference not shown: ", conditionMessage(local_inference)))
    } else {
      out$local_inference <- local_inference
    }
  }

  if (isTRUE(include_f_tests)) {
    f_obj <- object$Ftests
    if (is.null(f_obj)) {
      f_obj <- tryCatch(diagnostic_tests(object), error = function(e) e)
    }
    if (inherits(f_obj, "error")) {
      out$notes <- c(out$notes, paste0("Diagnostic F-tests not shown: ", conditionMessage(f_obj)))
    } else {
      out$f_tests <- f_obj
    }
  }

  if (isTRUE(include_moran) && is.null(moran)) {
    out$notes <- c(out$notes, "Moran's I not shown: pass a precomputed `moran` object to lhwr_report().")
  }
  if (isTRUE(include_collinearity) && is.null(collinearity)) {
    out$notes <- c(out$notes, "Local collinearity not shown: pass a precomputed `collinearity` object to lhwr_report().")
  }

  out$settings <- list(
    include_global = include_global,
    include_inference = include_inference,
    include_f_tests = include_f_tests,
    include_moran = include_moran,
    include_collinearity = include_collinearity,
    digits = digits,
    geometry = args$geometry,
    kernel = args$kernel
  )
  class(out) <- "lhwr_report"
  print(out, digits = digits)
  invisible(out)
}

#' @export
print.lhwr_report <- function(x, digits = NULL, ...) {
  digits <- digits %||% x$settings$digits %||% 4
  hgwr_rule()
  cat(" LHWR Model Report\n")
  cat(" Lorentz Hyperbolic Weighted Regression\n")
  hgwr_rule()

  cat("Call:\n")
  print(x$call)

  cat("\nData and calibration\n")
  hgwr_rule()
  print(x$calibration, row.names = FALSE, right = FALSE)

  if (!is.null(x$global)) {
    cat("\nGlobal OLS baseline\n")
    hgwr_rule()
    print(x$global, row.names = FALSE, right = FALSE)
  }

  cat("\nLHWR diagnostics\n")
  hgwr_rule()
  print(x$diagnostics, row.names = FALSE, right = FALSE)

  cat("\nLocal coefficient summaries\n")
  hgwr_rule()
  cat("Rows summarize the estimated coefficient surface over all locations.\n")
  printCoefmat(x$coefficients, digits = digits, signif.stars = FALSE)

  if (!is.null(x$local_inference)) {
    cat("\nLocal coefficient inference\n")
    hgwr_rule()
    lhwr_report_print_inference(x$local_inference)
  }

  if (!is.null(x$f_tests)) {
    cat("\n")
    print(x$f_tests)
  }

  if (!is.null(x$moran)) {
    cat("\nResidual hyperbolic autocorrelation\n")
    hgwr_rule()
    print(x$moran)
  }

  if (!is.null(x$collinearity)) {
    cat("\n")
    print(x$collinearity)
  }

  if (length(x$notes)) {
    cat("\nNotes\n")
    hgwr_rule()
    for (note in x$notes) cat(" - ", note, "\n", sep = "")
  }

  elapsed <- x$elapsed
  if (!is.na(elapsed)) {
    cat("\nElapsed time\n")
    hgwr_rule()
    cat(" ", hgwr_format_number(elapsed), " seconds\n", sep = "")
  }
  hgwr_rule()
  invisible(x)
}

lhwr_report_calibration_table <- function(object) {
  args <- object$LHWR.arguments
  n <- length(object$residuals)
  response <- lhwr_report_response_name(object)
  predictors <- setdiff(object$coef_names, "Intercept")
  if (!length(predictors)) predictors <- "Intercept only"
  data.frame(
    Setting = c(
      "Response",
      "Predictors",
      "Observations",
      "Local coefficients",
      "Geometry",
      "Distance",
      "Kernel",
      "Engine",
      "Diagnostics",
      "Bandwidth",
      "Bandwidth type"
    ),
    Value = c(
      response,
      paste(predictors, collapse = ", "),
      n,
      length(object$coef_names),
      hgwr_geometry_label(args$geometry),
      args$distance,
      hgwr_kernel_label(args$kernel),
      hgwr_engine_label(args$engine),
      hgwr_diagnostics_label(args$diagnostics),
      hgwr_bandwidth_label(args),
      hgwr_bandwidth_type_label(args$bandwidth_type)
    ),
    check.names = FALSE
  )
}

lhwr_report_global_table <- function(object) {
  lm_fit <- object$lm
  if (is.null(lm_fit) || !inherits(lm_fit, "lm")) {
    return(data.frame(Metric = "Global OLS baseline", Value = "Not available", check.names = FALSE))
  }
  sm <- summary(lm_fit)
  rss <- sum(stats::residuals(lm_fit)^2)
  vals <- c(
    "Residual sum of squares" = rss,
    "AIC" = stats::AIC(lm_fit),
    "BIC" = stats::BIC(lm_fit),
    "R-squared" = unname(sm$r.squared),
    "Adjusted R-squared" = unname(sm$adj.r.squared),
    "Residual df" = stats::df.residual(lm_fit)
  )
  data.frame(
    Metric = names(vals),
    Value = hgwr_format_number(vals),
    check.names = FALSE
  )
}

lhwr_report_print_inference <- function(inf) {
  settings <- inf$settings
  cat("Reference distribution : ", if (is.infinite(settings$reference_df)) "normal" else paste0("t(df = ", signif(settings$reference_df, 5), ")"), "\n", sep = "")
  cat("Alpha                  : ", settings$alpha, "\n", sep = "")
  cat("P-value adjustment     : ", settings$adjust, " (", settings$adjust_scope, ")\n", sep = "")
  cat("Number of local tests  : ", settings$n_tests, "\n", sep = "")
  tab <- inf$summary
  tab$significant_percent <- sprintf("%.1f", tab$significant_percent)
  tab$median_abs_t <- sprintf("%.3f", tab$median_abs_t)
  tab$min_p_value <- format.pval(tab$min_p_value, digits = 3, eps = 1e-4)
  tab$min_p_adjusted <- format.pval(tab$min_p_adjusted, digits = 3, eps = 1e-4)
  tab$max_abs_t <- sprintf("%.3f", tab$max_abs_t)
  names(tab) <- c("Coefficient", "Significant %", "Median |t|", "Min p", "Min adjusted p", "Max |t|")
  cat("\nCoefficient-level summary\n")
  print(tab, row.names = FALSE, right = FALSE)
  cat("Note: local p-values are exploratory and should be interpreted with\n")
  cat("multiplicity and neighborhood dependence in mind.\n")
}

lhwr_report_response_name <- function(object) {
  formula <- object$LHWR.arguments$formula
  if (is.null(formula)) return("Unknown")
  all.vars(formula)[1L] %||% "Unknown"
}

lhwr_report_fit_elapsed <- function(object) {
  if (!is.null(object$timings$stop) && !is.null(object$timings$start)) {
    return(as.numeric(difftime(object$timings$stop, object$timings$start, units = "secs")))
  }
  NA_real_
}

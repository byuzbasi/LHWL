#' Local Inference for LHWR Coefficient Surfaces
#'
#' Compute exploratory local t-tests and p-values for the coefficient surfaces
#' estimated by an LHWR model. The function uses the local coefficient estimates,
#' local standard errors, and local t-values stored in the fitted object; it does
#' not refit the model.
#'
#' These p-values should be interpreted as local exploratory diagnostics. Because
#' LHWR produces many local tests over dependent neighborhoods, multiplicity and
#' local dependence matter. Use `adjust = "BH"`, `"holm"`, or `"bonferroni"`
#' when the goal is to identify robust local patterns.
#'
#' Inference is conditional on the Lorentz-model coordinates and neighborhood
#' weights stored in `object`. If the representation was estimated, the reported
#' standard errors and p-values do not include first-stage representation
#' uncertainty or uncertainty from bandwidth selection.
#'
#' @param object An LHWR model object returned by [lhwr()].
#' @param alpha Significance level used to create the `significant` flag.
#' @param df Degrees-of-freedom rule. `"edf"` uses the model effective degrees of
#'   freedom stored in `object$LHWR.diagnostic$edf`; `"residual"` uses
#'   `n - tr(S)`; `"normal"` uses a standard normal reference distribution.
#' @param adjust Multiplicity adjustment passed to [stats::p.adjust()]. Use
#'   `"none"` for raw p-values.
#' @param adjust_scope Adjustment scope. `"coefficient"` adjusts p-values
#'   separately within each coefficient surface; `"global"` adjusts over all
#'   local coefficient tests at once.
#' @param ... Unused.
#' @return An object of class `"lhwr_local_inference"` with components:
#'   `local`, a long data frame of local estimates, standard errors, t-values,
#'   raw p-values, adjusted p-values, and significance flags; `summary`, a
#'   coefficient-level summary table; and `settings`.
#' @examples
#' set.seed(5)
#' H <- simulate_H2(80, r_max = 3)
#' Z <- H$Z
#' x1 <- rnorm(80)
#' y <- 1 + 2 * x1 + rnorm(80, sd = 0.3)
#' dat <- data.frame(y, x1)
#' fit <- lhwr(y ~ x1, dat, Z = Z, bw = 35, adaptive = TRUE, diagnostics = "full")
#' inf <- local_inference(fit, adjust = "BH")
#' inf
#' head(inf$local)
#' @export
local_inference <- function(object,
                            alpha = 0.05,
                            df = c("edf", "residual", "normal"),
                            adjust = c("BH", "none", "holm", "bonferroni", "BY", "fdr"),
                            adjust_scope = c("coefficient", "global"),
                            ...) {
  if (!is_lhwr_object(object)) {
    stop("`object` must be an LHWR model object.", call. = FALSE)
  }
  if (identical(object$family %||% object$LHWR.arguments$family, "binomial")) {
    stop(
      "Classical LHWR local t-inference is defined only for Gaussian models. Binomial fits retain Fisher-information standard errors and z-values in `SDF` for exploratory use.",
      call. = FALSE
    )
  }
  df <- match.arg(df)
  adjust <- match.arg(adjust)
  adjust_scope <- match.arg(adjust_scope)
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number between 0 and 1.", call. = FALSE)
  }

  coef_names <- object$coef_names
  sdf <- object$SDF
  missing_se <- setdiff(paste0(coef_names, "_SE"), names(sdf))
  missing_tv <- setdiff(paste0(coef_names, "_TV"), names(sdf))
  if (length(missing_se) || length(missing_tv)) {
    stop("The fitted object does not contain local standard errors and t-values. Refit with diagnostics = 'full' before calling local_inference().", call. = FALSE)
  }

  ref_df <- lhwr_inference_df(object, df)
  rows <- vector("list", length(coef_names))
  for (j in seq_along(coef_names)) {
    nm <- coef_names[j]
    estimate <- as.numeric(sdf[[nm]])
    se <- as.numeric(sdf[[paste0(nm, "_SE")]])
    t_value <- as.numeric(sdf[[paste0(nm, "_TV")]])
    p_value <- lhwr_t_pvalue(t_value, ref_df)
    rows[[j]] <- data.frame(
      observation_id = seq_along(estimate),
      coefficient = nm,
      estimate = estimate,
      se = se,
      t_value = t_value,
      p_value = p_value,
      stringsAsFactors = FALSE
    )
  }
  local <- do.call(rbind, rows)

  if (adjust_scope == "global") {
    local$p_adjusted <- stats::p.adjust(local$p_value, method = adjust)
  } else {
    local$p_adjusted <- stats::ave(
      local$p_value,
      local$coefficient,
      FUN = function(z) stats::p.adjust(z, method = adjust)
    )
  }
  local$significant <- is.finite(local$p_adjusted) & local$p_adjusted < alpha

  summary <- lhwr_local_inference_summary(local)
  settings <- list(
    alpha = alpha,
    df_rule = df,
    reference_df = ref_df,
    adjust = adjust,
    adjust_scope = adjust_scope,
    n_tests = nrow(local)
  )
  out <- list(local = local, summary = summary, settings = settings)
  class(out) <- "lhwr_local_inference"
  out
}

lhwr_inference_df <- function(object, df_rule) {
  n <- nrow(object$SDF)
  diag <- object$LHWR.diagnostic
  out <- switch(
    df_rule,
    edf = if (!is.null(diag$edf)) diag$edf else NA_real_,
    residual = if (!is.null(diag$tr.S)) n - diag$tr.S else NA_real_,
    normal = Inf
  )
  if (df_rule != "normal" && (!is.finite(out) || out <= 0)) {
    warning("Could not compute a positive reference degrees of freedom; using normal approximation.", call. = FALSE)
    out <- Inf
  }
  out
}

lhwr_t_pvalue <- function(t_value, df) {
  t_abs <- abs(as.numeric(t_value))
  if (is.infinite(df)) {
    out <- 2 * stats::pnorm(t_abs, lower.tail = FALSE)
  } else {
    out <- 2 * stats::pt(t_abs, df = df, lower.tail = FALSE)
  }
  out[!is.finite(t_abs)] <- NA_real_
  out
}

lhwr_local_inference_summary <- function(local) {
  rows <- lapply(split(local, local$coefficient), function(z) {
    data.frame(
      coefficient = z$coefficient[1L],
      significant_percent = 100 * mean(z$significant, na.rm = TRUE),
      median_abs_t = stats::median(abs(z$t_value), na.rm = TRUE),
      min_p_value = suppressWarnings(min(z$p_value, na.rm = TRUE)),
      min_p_adjusted = suppressWarnings(min(z$p_adjusted, na.rm = TRUE)),
      max_abs_t = suppressWarnings(max(abs(z$t_value), na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' @export
print.lhwr_local_inference <- function(x, ...) {
  cat("------------------------------------------------------------------------\n")
  cat(" LHWR Local Coefficient Inference\n")
  cat("------------------------------------------------------------------------\n")
  cat("Reference distribution : ", if (is.infinite(x$settings$reference_df)) "normal" else paste0("t(df = ", signif(x$settings$reference_df, 5), ")"), "\n", sep = "")
  cat("Alpha                  : ", x$settings$alpha, "\n", sep = "")
  cat("P-value adjustment     : ", x$settings$adjust, " (", x$settings$adjust_scope, ")\n", sep = "")
  cat("Number of local tests  : ", x$settings$n_tests, "\n", sep = "")
  cat("\nCoefficient-level summary\n")
  cat("------------------------------------------------------------------------\n")
  tab <- x$summary
  tab$significant_percent <- sprintf("%.1f", tab$significant_percent)
  tab$median_abs_t <- sprintf("%.3f", tab$median_abs_t)
  tab$min_p_value <- format.pval(tab$min_p_value, digits = 3, eps = 1e-4)
  tab$min_p_adjusted <- format.pval(tab$min_p_adjusted, digits = 3, eps = 1e-4)
  tab$max_abs_t <- sprintf("%.3f", tab$max_abs_t)
  names(tab) <- c("Coefficient", "Significant %", "Median |t|", "Min p", "Min adjusted p", "Max |t|")
  print(tab, row.names = FALSE, right = FALSE)
  cat("------------------------------------------------------------------------\n")
  cat("Note: local p-values are exploratory and should be interpreted with\n")
  cat("multiplicity and neighborhood dependence in mind.\n")
  cat("------------------------------------------------------------------------\n")
  invisible(x)
}

# Backward-compatible alias for objects that still inherit from hgwr.
#' @rdname local_inference
#' @export
local_inference.hgwr <- local_inference

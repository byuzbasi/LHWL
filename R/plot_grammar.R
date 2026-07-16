#' LHWR Plotting Grammar
#'
#' Convenience plotting functions that organize LHWR graphics by purpose. These
#' wrappers keep the package interface predictable for applied users while still
#' relying on the lower-level plotting functions such as [plot.lhwr()] and
#' [plot_lhwr_poincare_map()].
#'
#' Use the functions as a decision guide:
#' `plot_lhwr_geometry()` for the Lorentz geometry itself,
#' `plot_lhwr_embedding()` for embedded observations,
#' `plot_lhwr_network_map()` for graph or network data,
#' `plot_lhwr_model()` for fitted coefficients and residuals,
#' `plot_lhwr_inference()` for local p-values and significance,
#' and `plot_lhwr_diagnostics()` for model diagnostic surfaces.
#'
#' @name lhwr_plotting_grammar
NULL

#' Plot Lorentz Geometry
#'
#' Inspect the Lorentz hyperboloid coordinates used by LHWR. The hyperboloid view
#' emphasizes the theoretical model, the Poincare view gives a disk display, and
#' the radial view shows the radial/angular distribution of observations.
#'
#' @param Z An `n x 3` Lorentz hyperboloid coordinate matrix.
#' @param type Plot type: `"poincare"`, `"hyperboloid"`, or `"radial"`.
#' @param values Optional numeric vector for point colors. Defaults to the
#'   hyperbolic radius.
#' @param labels Optional point labels.
#' @param title Optional title. Use `NULL` for manuscript-style figures.
#' @param palette HCL palette for numeric values.
#' @param ... Additional arguments passed to the underlying plotting function.
#' @return Invisibly returns plotting coordinates or the input `Z`.
#' @examples
#' H <- simulate_H2(40, r_max = 3)
#' plot_lhwr_geometry(H$Z, type = "poincare", title = NULL)
#' plot_lhwr_geometry(H$Z, type = "radial")
#' @export
plot_lhwr_geometry <- function(Z,
                               type = c("poincare", "hyperboloid", "radial"),
                               values = NULL,
                               labels = NULL,
                               title = NULL,
                               palette = "Viridis",
                               show_size_legend = FALSE,
                               show_info = FALSE,
                               ...) {
  Z <- as_lorentz_matrix(Z, "Z")
  check_hyperboloid(Z, stop_on_fail = TRUE)
  type <- match.arg(type)
  r <- acosh(pmax(Z[, 1], 1))
  theta <- atan2(Z[, 3], Z[, 2])
  if (is.null(values)) values <- r

  if (type == "poincare") {
    return(plot_lhwr_poincare_map(
      Z,
      values = values,
      labels = labels,
      title = title,
      value_label = "Hyperbolic radius",
      point_palette = palette,
      show_size_legend = show_size_legend,
      show_info = show_info,
      ...
    ))
  }
  if (type == "hyperboloid") {
    return(plot_lhwr_hyperboloid(
      Z = Z,
      values = values,
      labels = labels,
      title = if (is.null(title)) "Lorentz hyperboloid" else title,
      point_palette = palette,
      ...
    ))
  }

  cols <- plot_value_colors(values, palette = palette)
  graphics::plot(
    theta,
    r,
    pch = 19,
    col = cols,
    xlab = "Angular coordinate theta",
    ylab = "Hyperbolic radius r",
    main = if (is.null(title)) "Lorentz radial coordinates" else title,
    ...
  )
  plot_value_legend(values, palette = palette)
  invisible(data.frame(theta = theta, radius = r, value = values))
}

#' Plot an LHWR Embedding
#'
#' Display embedded observations in the Poincare disk. This is the standard plot
#' for non-network feature maps, economic indicators, clinical variables,
#' functional scores, shape scores, and other two-dimensional representations
#' converted to Lorentz coordinates.
#'
#' @inheritParams plot_lhwr_poincare_map
#' @return Invisibly returns an object of class `lhwr_poincare_map`.
#' @examples
#' H <- simulate_H2(30, r_max = 2)
#' plot_lhwr_embedding(H$Z, values = H$r, title = NULL)
#' @export
plot_lhwr_embedding <- function(Z,
                                values = NULL,
                                groups = NULL,
                                significant = NULL,
                                labels = NULL,
                                node_names = NULL,
                                sizes = NULL,
                                title = NULL,
                                value_label = "Value",
                                group_label = "groups",
                                show_info = FALSE,
                                ...) {
  plot_lhwr_poincare_map(
    Z,
    values = values,
    groups = groups,
    significant = significant,
    labels = labels,
    node_names = node_names,
    sizes = sizes,
    title = title,
    value_label = value_label,
    group_label = group_label,
    show_info = show_info,
    ...
  )
}

#' Plot an LHWR Network Map
#'
#' Display graph or network data whose nodes have Lorentz hyperboloid
#' coordinates. This is the standard plotting interface for trade, citation,
#' social, transportation, similarity, and general graph embeddings.
#'
#' @inheritParams plot_lhwr_poincare_map
#' @return Invisibly returns an object of class `lhwr_poincare_map`.
#' @examples
#' H <- simulate_H2(20, r_max = 2)
#' edges <- data.frame(from = 1:5, to = 6:10)
#' plot_lhwr_network_map(H$Z, edges = edges, title = NULL)
#' @export
plot_lhwr_network_map <- function(Z,
                                  edges,
                                  edge_weight = NULL,
                                  values = NULL,
                                  groups = NULL,
                                  significant = NULL,
                                  labels = NULL,
                                  node_names = NULL,
                                  sizes = NULL,
                                  node_id = NULL,
                                  title = NULL,
                                  value_label = "Value",
                                  group_label = "network communities",
                                  size_label = "Node size",
                                  edge_label = "Network edges",
                                  edge_alpha_range = c(0.008, 0.055),
                                  edge_lwd_range = c(0.025, 0.30),
                                  show_info = FALSE,
                                  ...) {
  plot_lhwr_poincare_map(
    Z,
    values = values,
    groups = groups,
    significant = significant,
    labels = labels,
    node_names = node_names,
    sizes = sizes,
    edges = edges,
    edge_weight = edge_weight,
    node_id = node_id,
    title = title,
    value_label = value_label,
    group_label = group_label,
    size_label = size_label,
    edge_label = edge_label,
    edge_alpha_range = edge_alpha_range,
    edge_lwd_range = edge_lwd_range,
    show_info = show_info,
    ...
  )
}

#' Plot LHWR Model Surfaces
#'
#' User-friendly wrapper around [plot.lhwr()] for coefficient, fitted-value,
#' residual, and local R-squared maps.
#'
#' @param object An LHWR model object.
#' @param type Model surface type: `"coefficient"`, `"se"`, `"t_value"`,
#'   `"fitted"`, `"residual"`, or `"local_r2"`.
#' @param coefficient Coefficient name for coefficient-level surfaces.
#' @param ... Additional arguments passed to [plot.lhwr()].
#' @return Invisibly returns `object`.
#' @examples
#' H <- simulate_H2(35, r_max = 2)
#' x <- rnorm(35)
#' y <- 1 + x + rnorm(35, sd = 0.2)
#' fit <- lhwr(y ~ x, data.frame(y, x), Z = H$Z, bw = 15, adaptive = TRUE)
#' plot_lhwr_model(fit, type = "coefficient", coefficient = "x")
#' @export
plot_lhwr_model <- function(object,
                            type = c("coefficient", "se", "t_value", "fitted", "residual", "local_r2"),
                            coefficient = NULL,
                            ...) {
  type <- match.arg(type)
  plot(object, type = type, coefficient = coefficient, ...)
  invisible(object)
}

#' Plot LHWR Local Inference
#'
#' Plot local p-values, adjusted p-values, or local significance for a fitted
#' LHWR model. The Poincare geometry gives a clean disk map; projection and
#' hyperboloid modes reuse [plot.lhwr()].
#'
#' @param object An LHWR model object fitted with `diagnostics = "full"`.
#' @param coefficient Coefficient name.
#' @param type Inference surface: `"p_value"`, `"p_adjusted"`, or
#'   `"significance"`.
#' @param geometry Plot geometry: `"poincare"`, `"projection"`, or
#'   `"hyperboloid"`.
#' @param alpha Significance level.
#' @param p_adjust Multiplicity adjustment.
#' @param inference_df Degrees-of-freedom rule.
#' @param ... Additional arguments passed to the plotting function.
#' @return Invisibly returns the map object for Poincare plots or `object`.
#' @examples
#' H <- simulate_H2(45, r_max = 2)
#' x <- rnorm(45)
#' y <- 1 + x + rnorm(45, sd = 0.2)
#' fit <- lhwr(y ~ x, data.frame(y, x), Z = H$Z, bw = 18, adaptive = TRUE,
#'             diagnostics = "full")
#' plot_lhwr_inference(fit, coefficient = "x", type = "significance")
#' @export
plot_lhwr_inference <- function(object,
                                coefficient = NULL,
                                type = c("p_value", "p_adjusted", "significance"),
                                geometry = c("poincare", "projection", "hyperboloid"),
                                alpha = 0.05,
                                p_adjust = c("BH", "none", "holm", "bonferroni", "BY", "fdr"),
                                inference_df = c("edf", "residual", "normal"),
                                show_info = FALSE,
                                ...) {
  if (!is_lhwr_object(object)) stop("`object` must be an LHWR model object.", call. = FALSE)
  type <- match.arg(type)
  geometry <- match.arg(geometry)
  p_adjust <- match.arg(p_adjust)
  inference_df <- match.arg(inference_df)
  if (is.null(coefficient)) coefficient <- object$coef_names[1L]
  if (!coefficient %in% object$coef_names) stop("Unknown coefficient: ", coefficient, call. = FALSE)

  if (geometry != "poincare") {
    plot(object, type = type, coefficient = coefficient, geometry = geometry,
         alpha = alpha, p_adjust = p_adjust, inference_df = inference_df, ...)
    return(invisible(object))
  }

  inf <- local_inference(object, alpha = alpha, df = inference_df, adjust = p_adjust)
  z <- inf$local[inf$local$coefficient == coefficient, , drop = FALSE]
  Z <- resolve_plot_Z(object, NULL)
  if (type == "significance") {
    return(plot_lhwr_poincare_map(
      Z,
      significant = z$significant,
      labels = seq_len(nrow(Z)),
      title = NULL,
      value_label = paste0("Local significance: ", coefficient),
      show_size_legend = FALSE,
      show_info = show_info,
      ...
    ))
  }
  vals <- if (type == "p_value") z$p_value else z$p_adjusted
  label <- if (type == "p_value") {
    paste0("Local p-value: ", coefficient)
  } else {
    paste0("Local adjusted p-value: ", coefficient)
  }
  plot_lhwr_poincare_map(
    Z,
    values = vals,
    labels = seq_len(nrow(Z)),
    title = NULL,
    value_label = label,
    point_palette = "Inferno",
    show_size_legend = FALSE,
    show_info = show_info,
    ...
  )
}

#' Plot LHWR Diagnostics
#'
#' Plot common diagnostics from bandwidth and fitted LHWR objects.
#'
#' @param object An LHWR model object or a bandwidth object returned by
#'   [bw.lhwr()].
#' @param type Diagnostic type. For bandwidth objects use `"bandwidth"`. For
#'   model objects use `"residual"`, `"local_r2"`, `"condition_number"`, or
#'   `"effective_sample_size"`.
#' @param ... Additional plotting arguments.
#' @return Invisibly returns `object`.
#' @examples
#' H <- simulate_H2(35, r_max = 2)
#' x <- rnorm(35)
#' y <- 1 + x + rnorm(35, sd = 0.2)
#' fit <- lhwr(y ~ x, data.frame(y, x), Z = H$Z, bw = 15, adaptive = TRUE)
#' plot_lhwr_diagnostics(fit, type = "residual")
#' @export
plot_lhwr_diagnostics <- function(object,
                                  type = c("bandwidth", "residual", "local_r2", "condition_number", "effective_sample_size"),
                                  ...) {
  type <- match.arg(type)
  if (inherits(object, "bw.lhwr") || inherits(object, "bw.hgwr")) {
    if (type != "bandwidth") stop("Bandwidth objects only support type = 'bandwidth'.", call. = FALSE)
    plot(object, ...)
    return(invisible(object))
  }
  if (!is_lhwr_object(object)) stop("`object` must be an LHWR model or bandwidth object.", call. = FALSE)
  if (type %in% c("residual", "local_r2")) {
    plot(object, type = type, ...)
    return(invisible(object))
  }
  if (type == "condition_number") {
    lc <- local_collinearity(object)
    Z <- resolve_plot_Z(object, NULL)
    plot_lorentz_projection(Z, lc$condition_number$local_CN,
                            main = "Local condition number", palette = "Inferno", ...)
    return(invisible(object))
  }
  vals <- object$local_neff
  if (is.null(vals) && !is.null(object$fit.diagnostics$local_neff)) vals <- object$fit.diagnostics$local_neff
  if (is.null(vals)) stop("Effective sample size values are not available in this fitted object.", call. = FALSE)
  Z <- resolve_plot_Z(object, NULL)
  plot_lorentz_projection(Z, vals, main = "Local effective sample size", palette = "Viridis", ...)
  invisible(object)
}

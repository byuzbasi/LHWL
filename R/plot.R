
#' LHWR Plot Control
#'
#' Create a small list of graphical defaults used by the high-level LHWR
#' Poincare disk plotting functions. Most users can ignore this helper and use
#' the plotting functions directly. It is useful when preparing publication
#' figures or when labels and legends need to be larger or smaller.
#'
#' @param label_cex Character expansion for observation labels.
#' @param degree_cex Character expansion for degree labels on the disk boundary.
#' @param legend_cex Character expansion for categorical and size legends.
#' @param colorbar_cex Character expansion for colorbar text.
#' @param title_cex Character expansion for plot titles.
#' @param info_cex Character expansion for the optional information box.
#' @param point_min,point_max Minimum and maximum point sizes after scaling.
#' @return A list that can be passed to `plot_control` in
#'   [plot_lhwr_poincare_map()], [plot_lhwr_embedding()], or
#'   [plot_lhwr_network_map()].
#' @examples
#' H <- simulate_H2(30, r_max = 2)
#' ctl <- lhwr_plot_control(label_cex = 0.65, legend_cex = 0.65)
#' plot_lhwr_embedding(H$Z, values = H$r, plot_control = ctl)
#' @export
lhwr_plot_control <- function(label_cex = 0.55,
                              degree_cex = 0.50,
                              legend_cex = 0.55,
                              colorbar_cex = 0.55,
                              title_cex = 1.20,
                              info_cex = 0.55,
                              point_min = 0.45,
                              point_max = 2.80) {
  out <- list(
    label_cex = label_cex,
    degree_cex = degree_cex,
    legend_cex = legend_cex,
    colorbar_cex = colorbar_cex,
    title_cex = title_cex,
    info_cex = info_cex,
    point_min = point_min,
    point_max = point_max
  )
  for (nm in names(out)) {
    if (!is.numeric(out[[nm]]) || length(out[[nm]]) != 1L || !is.finite(out[[nm]]) || out[[nm]] <= 0) {
      stop("`", nm, "` must be a positive finite number.", call. = FALSE)
    }
  }
  out
}

lhwr_resolve_plot_control <- function(plot_control = NULL,
                                      label_cex = NULL,
                                      degree_cex = NULL,
                                      legend_cex = NULL,
                                      colorbar_cex = NULL,
                                      title_cex = NULL,
                                      info_cex = NULL,
                                      point_min = NULL,
                                      point_max = NULL) {
  ctl <- lhwr_plot_control()
  if (!is.null(plot_control)) {
    if (!is.list(plot_control)) stop("`plot_control` must be a list created by `lhwr_plot_control()`.", call. = FALSE)
    unknown <- setdiff(names(plot_control), names(ctl))
    if (length(unknown)) stop("Unknown `plot_control` field: ", paste(unknown, collapse = ", "), call. = FALSE)
    ctl[names(plot_control)] <- plot_control
  }
  overrides <- list(
    label_cex = label_cex,
    degree_cex = degree_cex,
    legend_cex = legend_cex,
    colorbar_cex = colorbar_cex,
    title_cex = title_cex,
    info_cex = info_cex,
    point_min = point_min,
    point_max = point_max
  )
  for (nm in names(overrides)) {
    if (!is.null(overrides[[nm]])) ctl[[nm]] <- overrides[[nm]]
  }
  for (nm in names(ctl)) {
    if (!is.numeric(ctl[[nm]]) || length(ctl[[nm]]) != 1L || !is.finite(ctl[[nm]]) || ctl[[nm]] <= 0) {
      stop("`", nm, "` must be a positive finite number.", call. = FALSE)
    }
  }
  ctl
}

#' Plot LHWR Objects
#'
#' Plot local model outputs from a fitted LHWR model or plot a bandwidth
#' selection curve. Model plots can show local coefficients, standard errors,
#' local t-values, local p-values, local significance flags, fitted values,
#' residuals, or local R-squared values. The same
#' quantity can be drawn either on the projected `(z1, z2)` coordinates or on a
#' perspective view of the Lorentz hyperboloid.
#'
#' @param x An object returned by [lhwr()] or [bw.lhwr()].
#' @param type Quantity to plot for model objects. `"coefficient"` plots a local
#'   coefficient estimate. `"se"` plots its local standard error. `"t_value"`
#'   plots its local t-value. `"p_value"` and `"p_adjusted"` plot local
#'   inference p-values. `"significance"` plots a 0/1 local significance flag
#'   based on adjusted p-values. `"residual"`, `"fitted"`, and `"local_r2"`
#'   plot model diagnostics stored in the fitted object.
#' @param coefficient Coefficient name to plot for coefficient-level plot types.
#'   Defaults to the first local coefficient.
#' @param Z Optional `n x 3` Lorentz coordinate matrix. If omitted for model
#'   objects, the coordinates stored in the fitted object are used.
#' @param geometry Plotting geometry. `"projection"` draws a fast two-dimensional
#'   map using `(z1, z2)`. `"hyperboloid"` draws the upper Lorentz hyperboloid
#'   surface `x0 = sqrt(1 + x1^2 + x2^2)` and overlays the observations.
#' @param palette Color palette passed to [grDevices::hcl.colors()].
#' @param point_size Point size.
#' @param main Optional plot title.
#' @param inference_df Degrees-of-freedom rule used when plotting local p-values;
#'   passed to [local_inference()].
#' @param p_adjust Multiplicity adjustment used when plotting adjusted p-values
#'   or significance; passed to [local_inference()].
#' @param alpha Significance level used when `type = "significance"`.
#' @param ... Additional arguments passed to base plotting functions.
#' @return Invisibly returns `x`.
#' @examples
#' set.seed(4)
#' H <- simulate_H2(60, r_max = 3)
#' Z <- H$Z
#' x1 <- rnorm(60)
#' y <- 1 + 2 * x1 + rnorm(60, sd = 0.3)
#' dat <- data.frame(y = y, x1 = x1)
#' fit <- lhwr(y ~ x1, dat, Z = Z, bw = 30, adaptive = TRUE, diagnostics = "full")
#' plot(fit, coefficient = "x1")
#' plot(fit, type = "p_value", coefficient = "x1")
#' plot(fit, type = "significance", coefficient = "x1")
#' plot(fit, type = "residual")
#' plot(fit, coefficient = "x1", geometry = "hyperboloid")
#' @export
plot.lhwr <- function(x,
                      type = c("coefficient", "se", "t_value", "p_value", "p_adjusted", "significance", "residual", "fitted", "local_r2"),
                      coefficient = NULL,
                      Z = NULL,
                      geometry = c("projection", "hyperboloid"),
                      palette = "Viridis",
                      point_size = 1.1,
                      main = NULL,
                      inference_df = c("edf", "residual", "normal"),
                      p_adjust = c("BH", "none", "holm", "bonferroni", "BY", "fdr"),
                      alpha = 0.05,
                      ...) {
  if (!is_lhwr_object(x)) {
    stop("`x` must be an LHWR model object.", call. = FALSE)
  }
  type <- match.arg(type)
  geometry <- match.arg(geometry)
  inference_df <- match.arg(inference_df)
  p_adjust <- match.arg(p_adjust)
  Z <- resolve_plot_Z(x, Z)
  value_info <- hgwr_plot_values(
    x,
    type = type,
    coefficient = coefficient,
    inference_df = inference_df,
    p_adjust = p_adjust,
    alpha = alpha
  )
  
  if (is.null(main)) {
    main <- value_info$title
  }
  if (geometry == "projection") {
    plot_lorentz_projection(
      Z = Z,
      values = value_info$values,
      main = main,
      palette = palette,
      point_size = point_size,
      ...
    )
  } else {
    plot_lhwr_hyperboloid(
      Z = Z,
      values = value_info$values,
      title = main,
      point_palette = palette,
      point_size = point_size,
      ...
    )
  }
  invisible(x)
}

#' @export
plot.bw.lhwr <- function(x, main = NULL, point_size = 1.1, ...) {
  if (!inherits(x, "bw.lhwr")) {
    stop("`x` must inherit from class 'bw.lhwr'.", call. = FALSE)
  }
  if (is.null(main)) {
    main <- paste("Bandwidth selection:", x$approach)
  }
  graphics::plot(
    x$selection.table$bw,
    x$selection.table$score,
    type = "b",
    pch = 19,
    cex = point_size,
    xlab = "Bandwidth",
    ylab = x$approach,
    main = main,
    ...
  )
  graphics::abline(v = x$bw, col = "firebrick", lwd = 2, lty = 2)
  graphics::legend(
    "topright",
    legend = paste("selected =", signif(x$bw, 5)),
    col = "firebrick",
    lwd = 2,
    lty = 2,
    bty = "n"
  )
  invisible(x)
}

resolve_plot_Z <- function(x, Z = NULL) {
  if (is.null(Z)) {
    if (!is.null(x$Z)) {
      Z <- x$Z
    } else {
      stop("`Z` is required for model plots.", call. = FALSE)
    }
  }
  Z <- as_lorentz_matrix(Z, "Z")
  check_hyperboloid(Z, stop_on_fail = TRUE)
  if (nrow(Z) != nrow(x$SDF)) {
    stop("`Z` must have one row per fitted observation.", call. = FALSE)
  }
  Z
}

hgwr_plot_values <- function(x, type, coefficient = NULL,
                             inference_df = "edf", p_adjust = "BH", alpha = 0.05) {
  if (type %in% c("coefficient", "se", "t_value", "p_value", "p_adjusted", "significance")) {
    if (is.null(coefficient)) coefficient <- x$coef_names[1L]
    if (!coefficient %in% x$coef_names) {
      stop("Unknown coefficient: ", coefficient, call. = FALSE)
    }
    if (type %in% c("p_value", "p_adjusted", "significance")) {
      inf <- local_inference(x, alpha = alpha, df = inference_df, adjust = p_adjust)
      z <- inf$local[inf$local$coefficient == coefficient, , drop = FALSE]
      values <- switch(
        type,
        p_value = z$p_value,
        p_adjusted = z$p_adjusted,
        significance = as.numeric(z$significant)
      )
      label <- switch(
        type,
        p_value = "Local p-value",
        p_adjusted = paste0("Local adjusted p-value (", p_adjust, ")"),
        significance = paste0("Local significance (alpha = ", alpha, ")")
      )
      return(list(values = values, title = paste(label, coefficient)))
    }
    suffix <- switch(
      type,
      coefficient = "",
      se = "_SE",
      t_value = "_TV"
    )
    field <- paste0(coefficient, suffix)
    if (!field %in% names(x$SDF)) {
      stop("Field not available in fitted object: ", field, call. = FALSE)
    }
    label <- switch(
      type,
      coefficient = "Local coefficient",
      se = "Local standard error",
      t_value = "Local t-value"
    )
    return(list(
      values = x$SDF[[field]],
      title = paste(label, coefficient)
    ))
  }
  
  field <- switch(
    type,
    residual = "residual",
    fitted = "yhat",
    local_r2 = "Local_R2"
  )
  label <- switch(
    type,
    residual = "Residuals",
    fitted = "Fitted values",
    local_r2 = "Local R-squared"
  )
  if (!field %in% names(x$SDF)) {
    stop("Field not available in fitted object: ", field, call. = FALSE)
  }
  list(values = x$SDF[[field]], title = label)
}

plot_value_colors <- function(values, palette = "Viridis", n = 100) {
  pal <- grDevices::hcl.colors(n, palette)
  if (all(!is.finite(values))) {
    return(rep(pal[1L], length(values)))
  }
  rng <- range(values, finite = TRUE)
  if (diff(rng) < .Machine$double.eps) {
    return(rep(pal[ceiling(n / 2)], length(values)))
  }
  idx <- cut(values, breaks = n, include.lowest = TRUE, labels = FALSE)
  pal[pmax(1L, pmin(n, idx))]
}

plot_lorentz_projection <- function(Z, values, main, palette = "Viridis",
                                    point_size = 1.1, ...) {
  cols <- plot_value_colors(values, palette = palette)
  graphics::plot(
    Z[, 2],
    Z[, 3],
    col = cols,
    pch = 19,
    cex = point_size,
    xlab = "x1",
    ylab = "x2",
    main = main,
    asp = 1,
    ...
  )
  plot_value_legend(values, palette = palette)
  invisible(NULL)
}

plot_lorentz_hyperboloid_base <- function(Z, values, main, palette = "Viridis",
                                          point_size = 1.1, grid_n = 60, ...) {
  lim <- max(abs(Z[, 2:3]), finite = TRUE) * 1.08
  x1 <- seq(-lim, lim, length.out = grid_n)
  x2 <- seq(-lim, lim, length.out = grid_n)
  x0 <- outer(x1, x2, function(a, b) sqrt(1 + a^2 + b^2))
  cols <- plot_value_colors(values, palette = palette)
  zlim <- range(c(x0, Z[, 1]), finite = TRUE)
  
  persp_obj <- graphics::persp(
    x1,
    x2,
    x0,
    theta = 35,
    phi = 24,
    expand = 0.78,
    col = grDevices::adjustcolor("grey85", alpha.f = 0.38),
    border = grDevices::adjustcolor("grey65", alpha.f = 0.18),
    xlab = "x1",
    ylab = "x2",
    zlab = "x0",
    zlim = zlim,
    ticktype = "detailed",
    main = main,
    ...
  )
  pts <- grDevices::trans3d(Z[, 2], Z[, 3], Z[, 1], pmat = persp_obj)
  graphics::points(pts, pch = 19, col = grDevices::adjustcolor(cols, alpha.f = 0.9), cex = point_size)
  plot_value_legend(values, palette = palette)
  invisible(NULL)
}

plot_value_legend <- function(values, palette = "Viridis") {
  if (all(!is.finite(values))) return(invisible(NULL))
  rng <- range(values, finite = TRUE)
  graphics::legend(
    "topleft",
    legend = c(sprintf("min %.3f", rng[1L]), sprintf("max %.3f", rng[2L])),
    pch = 19,
    col = plot_value_colors(rng, palette = palette),
    bty = "n",
    cex = 0.85
  )
  invisible(NULL)
}



#' Plot Lorentz Coordinates on the Hyperboloid
#'
#' Visualize Lorentz hyperboloid coordinates on the upper sheet
#' `z0 = sqrt(1 + z1^2 + z2^2)`. This plot emphasizes the geometry used by
#' LHWR estimation: observations live on the Lorentz hyperboloid, and local
#' neighborhoods are defined by Lorentz-hyperbolic distances. For readable
#' two-dimensional displays, use [plot_lhwr_poincare_map()].
#'
#' @param Z An `n x 3` Lorentz hyperboloid coordinate matrix with columns
#'   `(z0, z1, z2)`.
#' @param values Optional numeric vector used to color points. If omitted, the
#'   hyperbolic radius `acosh(z0)` is used.
#' @param labels Optional point labels.
#' @param sizes Optional numeric vector used to scale point sizes.
#' @param title Plot title. Use `NULL` for manuscript-style panels.
#' @param value_label Label for the colorbar. If omitted, a sensible default is
#'   chosen from `values`.
#' @param point_palette HCL palette for numeric `values`.
#' @param surface_col Color of the hyperboloid surface.
#' @param surface_alpha Transparency of the hyperboloid surface.
#' @param surface_border Grid-line color for the surface.
#' @param theta,phi Perspective angles passed to [graphics::persp()].
#' @param expand Expansion factor passed to [graphics::persp()].
#' @param grid_n Number of grid points per axis used to draw the surface.
#' @param surface_quantile Quantile in `(0, 1]` used to set the displayed
#'   surface extent from `abs(z1)` and `abs(z2)`, or `"auto"`. The automatic
#'   setting keeps the full range unless a few extreme coordinates would make the
#'   hyperboloid surface visually too large. Point coordinates are not changed.
#' @param point_size Fixed point size used when `sizes` is omitted.
#' @param show_labels Logical; draw labels for selected points.
#' @param n_label Number of labels to draw when `show_labels = TRUE`.
#' @param mar Base graphics margins passed to [graphics::par()] for the
#'   hyperboloid panel.
#' @param label_cex,colorbar_cex,title_cex Character sizes. These can also be
#'   supplied through `plot_control`. For `label_cex`, `"auto"` chooses a small
#'   value from the number of observations when labels are drawn.
#' @param point_min,point_max Minimum and maximum point sizes when `sizes` is
#'   supplied.
#' @param plot_control Optional list created by [lhwr_plot_control()]. Direct
#'   arguments such as `label_cex` or `colorbar_cex` override values in this
#'   list.
#' @param show_colorbar Logical; draw a numeric colorbar.
#' @param output_file Optional path. If supplied, the plot is written to this
#'   file. The device is inferred from the file extension unless `device` is
#'   supplied.
#' @param device Output device: `"auto"`, `"png"`, `"eps"`, or `"pdf"`.
#' @param width,height Output size in inches when `output_file` is supplied.
#' @param res Raster resolution for PNG output and fallback resolution for EPS.
#' @param ... Additional arguments passed to [graphics::persp()].
#' @return Invisibly returns a data frame containing the plotted Lorentz
#'   coordinates, hyperbolic radii, values, and point sizes.
#' @examples
#' H <- simulate_H2(30, r_max = 2)
#' plot_lhwr_hyperboloid(H$Z, values = H$r, title = NULL)
#' @export
plot_lhwr_hyperboloid <- function(Z,
                                  values = NULL,
                                  labels = NULL,
                                  sizes = NULL,
                                  title = "LHWR Lorentz hyperboloid view",
                                  value_label = NULL,
                                  point_palette = "Blue-Red 3",
                                  surface_col = "#D7E8F0",
                                  surface_alpha = 0.55,
                                  surface_border = "#FFFFFF",
                                  theta = 38,
                                  phi = 24,
                                  expand = 0.80,
                                  grid_n = 55L,
                                  surface_quantile = "auto",
                                  point_size = 1.1,
                                  show_labels = FALSE,
                                  n_label = 30L,
                                  mar = c(1.4, 1.4, 2.8, 0.8),
                                  label_cex = NULL,
                                  colorbar_cex = NULL,
                                  title_cex = NULL,
                                  point_min = NULL,
                                  point_max = NULL,
                                  plot_control = NULL,
                                  show_colorbar = TRUE,
                                  output_file = NULL,
                                  device = c("auto", "png", "eps", "pdf"),
                                  width = 7,
                                  height = 6,
                                  res = 220,
                                  ...) {
  Z <- as_lorentz_matrix(Z, "Z")
  check_hyperboloid(Z, stop_on_fail = TRUE)
  n <- nrow(Z)
  r <- acosh(pmax(Z[, 1], 1))
  if (is.null(values)) {
    values <- r
    if (is.null(value_label)) value_label <- "Hyperbolic radius"
  } else {
    if (length(values) != n) stop("`values` must have length nrow(Z).", call. = FALSE)
    values <- as.numeric(values)
    if (is.null(value_label)) value_label <- "Value"
  }
  if (!is.null(labels) && length(labels) != n) stop("`labels` must have length nrow(Z).", call. = FALSE)
  if (!is.null(sizes) && length(sizes) != n) stop("`sizes` must have length nrow(Z).", call. = FALSE)
  if (is.null(labels)) labels <- rep("", n)

  auto_label_cex <- is.character(label_cex) && length(label_cex) == 1L && identical(tolower(label_cex), "auto")
  if (is.character(label_cex) && !auto_label_cex) {
    stop('`label_cex` must be numeric, NULL, or "auto".', call. = FALSE)
  }
  ctl <- lhwr_resolve_plot_control(
    plot_control = plot_control,
    label_cex = if (auto_label_cex) NULL else label_cex,
    colorbar_cex = colorbar_cex,
    title_cex = title_cex,
    point_min = point_min,
    point_max = point_max
  )
  if (isTRUE(show_labels) && (auto_label_cex || (is.null(label_cex) && is.null(plot_control)))) {
    ctl$label_cex <- lhwr_hyperboloid_auto_label_cex(n)
  }

  if (is.null(sizes)) {
    point_cex <- rep(point_size, n)
    sizes <- rep(1, n)
  } else {
    sizes <- as.numeric(sizes)
    point_cex <- lhwr_map_scale(sizes, min_out = ctl$point_min, max_out = ctl$point_max)
  }

  if (!is.numeric(grid_n) || length(grid_n) != 1L || !is.finite(grid_n) || grid_n < 10) {
    stop("`grid_n` must be a finite number at least 10.", call. = FALSE)
  }
  grid_n <- as.integer(grid_n)
  if (!is.numeric(mar) || length(mar) != 4L || any(!is.finite(mar)) || any(mar < 0)) {
    stop("`mar` must be a numeric vector of length 4 with nonnegative finite values.", call. = FALSE)
  }
  auto_surface_quantile <- is.character(surface_quantile) && length(surface_quantile) == 1L &&
    identical(tolower(surface_quantile), "auto")
  if (!auto_surface_quantile) {
    if (!is.numeric(surface_quantile) || length(surface_quantile) != 1L ||
        !is.finite(surface_quantile) || surface_quantile <= 0 || surface_quantile > 1) {
      stop('`surface_quantile` must be a number in (0, 1] or "auto".', call. = FALSE)
    }
    surface_quantile <- as.numeric(surface_quantile)
  }

  close_device <- FALSE
  if (!is.null(output_file)) {
    close_device <- lhwr_open_output_device(output_file, device = device, width = width, height = height, res = res)
    on.exit(if (close_device) grDevices::dev.off(), add = TRUE)
  }

  old_par <- graphics::par(mar = mar, family = "sans", bg = "white")
  on.exit(graphics::par(old_par), add = TRUE)

  abs_coords <- abs(c(Z[, 2], Z[, 3]))
  abs_coords <- abs_coords[is.finite(abs_coords)]
  surface_quantile_used <- if (auto_surface_quantile) {
    lhwr_hyperboloid_auto_surface_quantile(abs_coords)
  } else {
    surface_quantile
  }
  lim <- if (length(abs_coords)) {
    as.numeric(stats::quantile(abs_coords, probs = surface_quantile_used, names = FALSE, type = 8)) * 1.12
  } else {
    NA_real_
  }
  if (!is.finite(lim) || lim <= 0) lim <- 1
  x1 <- seq(-lim, lim, length.out = grid_n)
  x2 <- seq(-lim, lim, length.out = grid_n)
  x0 <- outer(x1, x2, function(a, b) sqrt(1 + a^2 + b^2))
  cols <- plot_value_colors(values, palette = point_palette)
  zlim <- range(c(x0, Z[, 1]), finite = TRUE)

  pmat <- graphics::persp(
    x1,
    x2,
    x0,
    theta = theta,
    phi = phi,
    expand = expand,
    col = grDevices::adjustcolor(surface_col, alpha.f = surface_alpha),
    border = surface_border,
    shade = 0.22,
    xlab = "z1",
    ylab = "z2",
    zlab = "z0",
    zlim = zlim,
    ticktype = "detailed",
    box = FALSE,
    main = title,
    cex.main = ctl$title_cex,
    ...
  )
  pts <- grDevices::trans3d(Z[, 2], Z[, 3], Z[, 1], pmat = pmat)
  ord <- order(Z[, 1])
  graphics::points(pts$x[ord], pts$y[ord], pch = 21,
                   bg = grDevices::adjustcolor(cols[ord], alpha.f = 0.92),
                   col = "#263238", cex = point_cex[ord], lwd = 0.35)

  if (isTRUE(show_labels) && n_label > 0) {
    n_label <- min(as.integer(n_label), n)
    idx <- if (!is.null(sizes)) order(sizes, decreasing = TRUE)[seq_len(n_label)] else seq_len(n_label)
    graphics::text(pts$x[idx], pts$y[idx], labels = labels[idx], cex = ctl$label_cex,
                   col = "#263238", pos = 4, offset = 0.15, xpd = NA)
  }

  if (isTRUE(show_colorbar)) {
    lhwr_hyperboloid_colorbar(values, value_label, palette = point_palette, cex = ctl$colorbar_cex)
  }

  out <- data.frame(
    label = labels,
    size = as.numeric(sizes),
    value = as.numeric(values),
    hyperbolic_radius = r,
    z0 = Z[, 1],
    z1 = Z[, 2],
    z2 = Z[, 3]
  )
  attr(out, "metadata") <- list(
    n = n,
    title = title,
    value_label = value_label,
    surface = "upper Lorentz hyperboloid",
    distance = "Lorentz-hyperbolic",
    surface_quantile = surface_quantile_used
  )
  class(out) <- c("lhwr_hyperboloid_plot", class(out))
  invisible(out)
}

lhwr_open_output_device <- function(output_file,
                                    device = c("auto", "png", "eps", "pdf"),
                                    width = 7,
                                    height = 6,
                                    res = 220) {
  device <- match.arg(device)
  output_file <- path.expand(output_file)
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  if (device == "auto") {
    ext <- tolower(tools::file_ext(output_file))
    device <- if (ext %in% c("eps", "ps")) "eps" else if (ext == "pdf") "pdf" else "png"
  }
  if (device == "eps") {
    grDevices::cairo_ps(output_file, width = width, height = height, onefile = FALSE, fallback_resolution = res)
  } else if (device == "pdf") {
    grDevices::pdf(output_file, width = width, height = height, onefile = TRUE)
  } else {
    grDevices::png(output_file, width = width, height = height, units = "in", res = res)
  }
  TRUE
}


lhwr_hyperboloid_auto_surface_quantile <- function(abs_coords) {
  abs_coords <- abs_coords[is.finite(abs_coords) & abs_coords > 0]
  if (length(abs_coords) < 20L) return(1)
  q98 <- as.numeric(stats::quantile(abs_coords, probs = 0.98, names = FALSE, type = 8))
  q95 <- as.numeric(stats::quantile(abs_coords, probs = 0.95, names = FALSE, type = 8))
  mx <- max(abs_coords, na.rm = TRUE)
  if (is.finite(q98) && q98 > 0 && mx / q98 > 1.75) return(0.995)
  if (is.finite(q95) && q95 > 0 && mx / q95 > 2.25) return(0.99)
  1
}

lhwr_hyperboloid_auto_label_cex <- function(n) {
  n <- max(1L, as.integer(n))
  max(0.20, min(0.50, 0.55 / sqrt(max(1, n / 35))))
}

lhwr_hyperboloid_colorbar <- function(values, title, palette = "Blue-Red 3",
                                      cex = 0.55) {
  values <- as.numeric(values)
  rng <- range(values, finite = TRUE)
  if (!all(is.finite(rng)) || diff(rng) < .Machine$double.eps) return(invisible(NULL))
  usr <- graphics::par("usr")
  dx <- diff(usr[1:2])
  dy <- diff(usr[3:4])
  x0 <- usr[1] + 0.31 * dx
  x1 <- usr[1] + 0.69 * dx
  y <- usr[3] + 0.315 * dy
  h <- 0.019 * dy
  pal <- grDevices::hcl.colors(100, palette)
  xs <- seq(x0, x1, length.out = length(pal) + 1L)
  for (i in seq_along(pal)) {
    graphics::rect(xs[i], y, xs[i + 1L], y + h, col = pal[i], border = NA, xpd = NA)
  }
  graphics::rect(x0, y, x1, y + h, border = "#263238", lwd = 0.35, xpd = NA)
  graphics::text(x0, y - 0.018 * dy, labels = sprintf("%.2f", rng[1]), adj = c(0.5, 1), cex = cex * 0.86, col = "#263238", xpd = NA)
  graphics::text(x1, y - 0.018 * dy, labels = sprintf("%.2f", rng[2]), adj = c(0.5, 1), cex = cex * 0.86, col = "#263238", xpd = NA)
  graphics::text((x0 + x1) / 2, y + h + 0.018 * dy, labels = title, adj = c(0.5, 0), cex = cex, font = 2, col = "#263238", xpd = NA)
  invisible(NULL)
}

#' Plot Lorentz Coordinates in the Poincare Disk
#'
#' Visualize Lorentz hyperboloid coordinates in the Poincare disk. This is a
#' display tool: distances and model estimation in LHWR are still based on the
#' original Lorentz-hyperbolic geometry unless the user explicitly uses another
#' geometry. The function is useful for network, similarity, and trade-map
#' applications where nodes, edge backbones, labels, node sizes, and value-based
#' colors should be shown in a hyperbolic disk view.
#'
#' @param Z An `n x 3` Lorentz hyperboloid coordinate matrix with columns
#'   `(z0, z1, z2)`.
#' @param values Optional numeric vector used to color points. If omitted and
#'   `groups` is also omitted, points use a single color.
#' @param groups Optional grouping variable used for categorical point colors.
#'   Ignored when `values` or `significant` is supplied.
#' @param significant Optional logical vector. If supplied, points are colored as
#'   significant/non-significant using `significance_colors`; useful for local
#'   p-value maps.
#' @param significance_adjusted Optional logical vector used together with
#'   `significant`. If supplied, local significance is classified as not
#'   significant, raw-only significant, or significant after adjustment. This is
#'   useful for comparing raw local p-values with multiplicity-adjusted
#'   p-values without drawing a second map.
#' @param significance_direction Optional numeric vector used together with
#'   `significant`. Negative values are shown as negative significant, positive
#'   values as positive significant, and non-significant locations remain grey.
#'   This is useful for local coefficient maps where both significance and sign
#'   should be visible.
#' @param labels Optional labels for points, for example country ISO3 codes.
#' @param node_names Optional longer node names used in map summaries.
#' @param sizes Optional nonnegative numeric vector for point sizes. A common
#'   choice is original total trade or a node degree/strength measure. The map
#'   already applies `log1p` before rescaling sizes; do not log-transform twice
#'   unless that additional visual transformation is explicitly intended.
#' @param edges Optional edge list as a data frame or matrix with two columns.
#'   Entries may be integer row indices or labels matching `node_id`.
#' @param edge_weight Optional numeric vector of edge weights used for line
#'   widths and transparency.
#' @param edge_alpha_range Length-2 numeric vector giving the minimum and
#'   maximum edge transparency. Larger values make network links more visible.
#' @param edge_lwd_range Length-2 numeric vector giving the minimum and maximum
#'   edge line width.
#' @param node_id Optional node identifiers used to match `edges` when edge
#'   endpoints are not integer row indices.
#' @param title Plot title. Use `NULL` to suppress the title, which is useful when the title is supplied in a manuscript caption.
#' @param value_label Label for the colorbar when `values` is supplied.
#' @param group_label Label used in summaries and captions when `groups` is supplied.
#' @param size_label Label explaining point size.
#' @param edge_label Label explaining edges.
#' @param disk_radius Visual radius used by the largest point after optional
#'   radial scaling. This changes only the display, not `Z` or distances.
#' @param radial_rescale If `TRUE`, rescales Poincare coordinates for display.
#'   This changes only the graphic, not `Z` or the distances used by LHWR.
#' @param radial_quantile Quantile of the Poincare radii used as the visual
#'   reference radius when `radial_rescale = TRUE`. The default `1` preserves
#'   the previous behavior. Values below 1 are useful when a few peripheral
#'   points would otherwise compress the rest of the map.
#' @param radial_power Positive power applied to rescaled radii for display.
#'   Values below 1 spread central points outward while preserving angles.
#' @param show_all_labels If `TRUE`, all labels are drawn. If `FALSE`, only the
#'   largest `n_label` points according to `sizes` are labeled.
#' @param n_label Number of labels to draw when `show_all_labels = FALSE`.
#' @param label_cex Label size. If `NULL`, the value in `plot_control` is used.
#' @param degree_cex Character size for degree labels on the disk boundary.
#' @param legend_cex Character size for categorical, significance, and size legends.
#' @param colorbar_cex Character size for numeric colorbar labels and title.
#' @param title_cex Character size for the plot title.
#' @param info_cex Character size for the optional information box.
#' @param point_min,point_max Minimum and maximum point sizes after scaling.
#' @param plot_control Optional list created by [lhwr_plot_control()]. Direct
#'   arguments such as `label_cex` or `legend_cex` override values in this list.
#' @param point_palette HCL palette for numeric `values`.
#' @param group_palette Optional vector of colors for `groups`.
#' @param significance_colors Named or unnamed length-2 color vector for
#'   `significant = FALSE` and `significant = TRUE`, respectively.
#' @param signed_significance_colors Named or unnamed length-3 color vector for
#'   not significant, negative significant, and positive significant locations
#'   when both `significant` and `significance_direction` are supplied.
#' @param adjusted_significance_colors Named or unnamed length-3 color vector
#'   for not significant, raw-only significant, and raw-and-adjusted significant
#'   locations when both `significant` and `significance_adjusted` are supplied.
#' @param significance_label Legend title used when `significant` is supplied.
#' @param background_col Disk background color.
#' @param boundary_col Disk boundary color.
#' @param edge_col Edge color.
#' @param show_colorbar Logical; draw a colorbar for numeric `values`.
#' @param show_group_legend Logical; draw a legend for categorical `groups`.
#' @param group_legend_position Position for the categorical/significance
#'   legend. Use `"right"` for a side legend or `"bottom"` to keep the
#'   legend below the disk when side placement overlaps the plot.
#' @param show_size_legend Logical; draw point-size legend.
#' @param show_info Logical; draw an explanatory annotation box.
#' @param output_file Optional output path. PDF and EPS/PS extensions select
#'   the corresponding device; other extensions use PNG. Existing files may
#'   be overwritten: check the target first if preservation is required.
#' @param width,height,res Width and height in pixels and resolution in pixels
#'   per inch. PDF/EPS dimensions in inches are `width/res` and `height/res`.
#'   These differ from the inch-based hyperboloid output settings.
#' @param ... Additional graphical arguments currently ignored.
#' @return Invisibly returns a data frame containing the Poincare coordinates
#'   used for plotting.
#' @examples
#' set.seed(1)
#' H <- simulate_H2(40, r_max = 3)
#' Z <- H$Z
#' plot_lhwr_poincare_map(
#'   Z,
#'   values = H$r,
#'   labels = paste0("n", seq_len(nrow(Z))),
#'   title = "Simulated points in the Poincare disk",
#'   value_label = "radius"
#' )
#' @export
plot_lhwr_poincare_map <- function(Z,
                                   values = NULL,
                                   groups = NULL,
                                   significant = NULL,
                                   significance_adjusted = NULL,
                                   significance_direction = NULL,
                                   labels = NULL,
                                   node_names = NULL,
                                   sizes = NULL,
                                   edges = NULL,
                                   edge_weight = NULL,
                                   edge_alpha_range = c(0.08, 0.30),
                                   edge_lwd_range = c(0.25, 1.60),
                                   node_id = NULL,
                                   title = "LHWR Poincare disk map",
                                   value_label = "Value",
                                   group_label = "groups",
                                   size_label = "Node size",
                                   edge_label = "Edges",
                                   disk_radius = 0.94,
                                   radial_rescale = TRUE,
                                   radial_quantile = 1,
                                   radial_power = 1,
                                   show_all_labels = TRUE,
                                   n_label = 30L,
                                   label_cex = NULL,
                                   degree_cex = NULL,
                                   legend_cex = NULL,
                                   colorbar_cex = NULL,
                                   title_cex = NULL,
                                   info_cex = NULL,
                                   point_min = NULL,
                                   point_max = NULL,
                                   plot_control = NULL,
                                   point_palette = "Blue-Red 3",
                                   group_palette = NULL,
                                   significance_colors = c("Not significant" = "#BDBDBD", "Significant" = "#2E7D32"),
                                   signed_significance_colors = c(
                                     "Not significant" = "#C9C9C9",
                                     "Negative significant" = "#7B3294",
                                     "Positive significant" = "#2E7D32"
                                   ),
                                   adjusted_significance_colors = c(
                                     "Not significant" = "#C9C9C9",
                                     "Raw only" = "#E69F00",
                                     "Raw and BH" = "#1B7837"
                                   ),
                                   significance_label = "Local significance",
                                   background_col = "#F6F1E7",
                                   boundary_col = "#0A9396",
                                   edge_col = "#607D8B",
                                   show_colorbar = TRUE,
                                   show_group_legend = TRUE,
                                   group_legend_position = c("right", "bottom"),
                                   show_size_legend = TRUE,
                                   show_info = TRUE,
                                   output_file = NULL,
                                   width = 1500,
                                   height = 1500,
                                   res = 170,
                                   ...) {
  Z <- as_lorentz_matrix(Z, "Z")
  check_hyperboloid(Z, stop_on_fail = TRUE)
  n <- nrow(Z)
  if (!is.null(values) && length(values) != n) stop("`values` must have length nrow(Z).", call. = FALSE)
  if (!is.null(groups) && length(groups) != n) stop("`groups` must have length nrow(Z).", call. = FALSE)
  if (!is.null(significant) && length(significant) != n) stop("`significant` must have length nrow(Z).", call. = FALSE)
  if (!is.null(significance_adjusted) && length(significance_adjusted) != n) {
    stop("`significance_adjusted` must have length nrow(Z).", call. = FALSE)
  }
  if (!is.null(significance_direction) && length(significance_direction) != n) {
    stop("`significance_direction` must have length nrow(Z).", call. = FALSE)
  }
  if (!is.null(labels) && length(labels) != n) stop("`labels` must have length nrow(Z).", call. = FALSE)
  if (!is.null(node_names) && length(node_names) != n) stop("`node_names` must have length nrow(Z).", call. = FALSE)
  if (!is.null(sizes) && length(sizes) != n) stop("`sizes` must have length nrow(Z).", call. = FALSE)
  if (!is.null(node_id) && length(node_id) != n) stop("`node_id` must have length nrow(Z).", call. = FALSE)

  px <- Z[, 2] / (Z[, 1] + 1)
  py <- Z[, 3] / (Z[, 1] + 1)
  pr <- sqrt(px^2 + py^2)
  if (!is.numeric(radial_quantile) || length(radial_quantile) != 1L ||
      !is.finite(radial_quantile) || radial_quantile <= 0 || radial_quantile > 1) {
    stop("`radial_quantile` must be a number in (0, 1].", call. = FALSE)
  }
  if (!is.numeric(radial_power) || length(radial_power) != 1L ||
      !is.finite(radial_power) || radial_power <= 0) {
    stop("`radial_power` must be a positive finite number.", call. = FALSE)
  }

  scale_factor <- 1
  if (isTRUE(radial_rescale) && max(pr, na.rm = TRUE) > 0) {
    finite_pr <- pr[is.finite(pr) & pr > 0]
    reference_radius <- as.numeric(stats::quantile(finite_pr, probs = radial_quantile, names = FALSE, type = 8))
    if (!is.finite(reference_radius) || reference_radius <= 0) reference_radius <- max(finite_pr, na.rm = TRUE)
    scale_factor <- disk_radius / reference_radius
    theta <- atan2(py, px)
    pr_unit <- pmin(pr / reference_radius, 1)
    pr_plot <- disk_radius * pr_unit^radial_power
    px_plot <- pr_plot * cos(theta)
    py_plot <- pr_plot * sin(theta)
    px_plot[pr == 0 | !is.finite(pr)] <- 0
    py_plot[pr == 0 | !is.finite(pr)] <- 0
  } else {
    px_plot <- px
    py_plot <- py
    pr_plot <- pr
  }

  ctl <- lhwr_resolve_plot_control(
    plot_control = plot_control,
    label_cex = label_cex,
    degree_cex = degree_cex,
    legend_cex = legend_cex,
    colorbar_cex = colorbar_cex,
    title_cex = title_cex,
    info_cex = info_cex,
    point_min = point_min,
    point_max = point_max
  )

  if (is.null(labels)) labels <- rep("", n)
  if (is.null(node_names)) node_names <- labels
  if (is.null(sizes)) sizes <- rep(1, n)
  point_cex <- lhwr_map_scale(sizes, min_out = ctl$point_min, max_out = ctl$point_max)
  significance_factor <- NULL
  if (!is.null(significant)) {
    significance_factor <- if (!is.null(significance_adjusted)) {
      lhwr_adjusted_significance_factor(significant, significance_adjusted)
    } else if (is.null(significance_direction)) {
      lhwr_significance_factor(significant)
    } else {
      lhwr_signed_significance_factor(significant, significance_direction)
    }
  }
  point_cols <- lhwr_map_point_colors(
    values,
    groups,
    point_palette,
    group_palette,
    significant,
    significance_adjusted,
    significance_direction,
    significance_colors,
    signed_significance_colors,
    adjusted_significance_colors
  )

  edge_df <- lhwr_map_edges(
    edges, edge_weight, node_id, px_plot, py_plot,
    edge_alpha_range = edge_alpha_range,
    edge_lwd_range = edge_lwd_range
  )

  if (!is.null(output_file)) {
    output_file <- path.expand(output_file)
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    ext <- tolower(tools::file_ext(output_file))
    if (ext %in% c("eps", "ps")) {
      grDevices::cairo_ps(output_file, width = width / res, height = height / res,
                          onefile = FALSE, fallback_resolution = res)
    } else if (ext == "pdf") {
      grDevices::pdf(output_file, width = width / res, height = height / res,
                     onefile = TRUE)
    } else {
      grDevices::png(output_file, width = width, height = height, units = "px", res = res)
    }
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  old_par <- graphics::par(mar = c(1.5, 1.5, 3.5, 1.5), family = "sans", xaxs = "i", yaxs = "i")
  on.exit(graphics::par(old_par), add = TRUE)
  lhwr_map_disk_background(title, background_col = background_col, boundary_col = boundary_col,
                           degree_cex = ctl$degree_cex, title_cex = ctl$title_cex)
  lhwr_map_draw_edges(edge_df, edge_col = edge_col)
  graphics::points(px_plot, py_plot, pch = 21, bg = point_cols, col = "#263238", cex = point_cex, lwd = 0.35)

  if (show_all_labels || n_label > 0) {
    if (show_all_labels) {
      idx <- seq_len(n)
    } else {
      idx <- order(sizes, decreasing = TRUE)[seq_len(min(n_label, n))]
    }
    graphics::text(px_plot[idx], py_plot[idx], labels = labels[idx], cex = ctl$label_cex,
                   col = "#263238", pos = 4, offset = 0.12)
  }

  if (!is.null(values) && is.null(significant) && isTRUE(show_colorbar)) {
    lhwr_map_colorbar(values, value_label, palette = point_palette, cex = ctl$colorbar_cex)
  } else if (!is.null(significant) && (!is.null(significance_adjusted) || !is.null(significance_direction)) && isTRUE(show_group_legend)) {
    lhwr_map_group_legend(
      significance_factor,
      point_cols,
      position = group_legend_position,
      title = significance_label,
      cex = ctl$legend_cex
    )
  } else if (!is.null(significant) && isTRUE(show_group_legend)) {
    lhwr_map_significance_legend(significance_colors, title = significance_label, cex = ctl$legend_cex)
  } else if (is.null(values) && !is.null(groups) && isTRUE(show_group_legend)) {
    lhwr_map_group_legend(groups, point_cols, position = group_legend_position, title = group_label, cex = ctl$legend_cex)
  }
  if (isTRUE(show_size_legend)) {
    lhwr_map_size_legend(size_label, cex = ctl$legend_cex)
  }
  if (isTRUE(show_info)) {
    color_line <- if (!is.null(significant) && !is.null(significance_adjusted)) {
      "Color: raw and adjusted local significance"
    } else if (!is.null(significant) && !is.null(significance_direction)) {
      "Color: signed local significance"
    } else if (!is.null(significant)) {
      "Color: local significance"
    } else if (!is.null(values)) {
      paste("Color:", value_label)
    } else if (!is.null(groups)) {
      "Color: group"
    } else {
      "Color: fixed"
    }
    lhwr_map_info_box(c(
      color_line,
      paste("Size:", size_label),
      edge_label,
      "Disk view: Poincare projection of Lorentz Z"
    ), cex = ctl$info_cex)
  }

  out <- data.frame(
    label = labels,
    node_name = node_names,
    size = as.numeric(sizes),
    log_size = log1p(pmax(as.numeric(sizes), 0)),
    hyperbolic_radius = acosh(pmax(Z[, 1], 1)),
    poincare_radius = pr,
    p_x = px,
    p_y = py,
    p_x_plot = px_plot,
    p_y_plot = py_plot,
    p_r = pr,
    p_r_plot = pr_plot,
    visual_scale = scale_factor
  )
  edge_n <- if (is.null(edge_df)) 0L else nrow(edge_df)
  metadata <- list(
    n_nodes = n,
    n_edges = edge_n,
    title = title,
    value_label = if (!is.null(values)) value_label else NA_character_,
    group_label = if (!is.null(groups)) group_label else NA_character_,
    color_variable = if (!is.null(significant) && !is.null(significance_adjusted)) "raw and adjusted local significance" else if (!is.null(significant) && !is.null(significance_direction)) "signed local significance" else if (!is.null(significant)) "local significance" else if (!is.null(values)) value_label else if (!is.null(groups)) group_label else "fixed",
    size_label = size_label,
    edge_label = edge_label,
    disk_radius = disk_radius,
    radial_rescale = radial_rescale,
    visual_scale = scale_factor,
    radial_quantile = radial_quantile,
    radial_power = radial_power,
    cor_log_size_hyperbolic_radius = safe_cor(log1p(pmax(as.numeric(sizes), 0)), acosh(pmax(Z[, 1], 1))),
    cor_log_size_poincare_radius = safe_cor(log1p(pmax(as.numeric(sizes), 0)), pr),
    group_counts = if (!is.null(significance_factor)) sort(table(significance_factor), decreasing = TRUE) else if (!is.null(groups)) sort(table(groups), decreasing = TRUE) else NULL
  )
  attr(out, "metadata") <- metadata
  class(out) <- c("lhwr_poincare_map", class(out))
  invisible(out)
}

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) return(NA_real_)
  if (stats::sd(x[ok]) < .Machine$double.eps || stats::sd(y[ok]) < .Machine$double.eps) return(NA_real_)
  stats::cor(x[ok], y[ok])
}

#' Summarize a Poincare Disk Map
#'
#' Print a compact interpretation summary for an object returned by
#' [plot_lhwr_poincare_map()]. The summary keeps the graphic clean while making
#' the figure metadata reproducible, including the number of nodes and edges,
#' radial-size correlations, and the most central observations.
#'
#' @param object An object returned invisibly by [plot_lhwr_poincare_map()].
#' @param top_n Number of most central observations to print.
#' @param ... Unused.
#' @return Invisibly returns a list containing map metadata and the central-node table.
#' @export
summary.lhwr_poincare_map <- function(object, top_n = 10L, ...) {
  meta <- attr(object, "metadata")
  if (is.null(meta)) stop("`object` does not contain LHWR map metadata.", call. = FALSE)
  top_n <- max(1L, min(as.integer(top_n), nrow(object)))
  central <- object[order(object$hyperbolic_radius),
                    c("label", "node_name", "size", "log_size", "hyperbolic_radius", "poincare_radius")]
  central <- utils::head(central, top_n)
  rownames(central) <- NULL
  cat("------------------------------------------------------------------------
")
  cat(" LHWR Poincare Disk Map Summary
")
  cat("------------------------------------------------------------------------
")
  cat("Title        : ", meta$title, "
", sep = "")
  cat("Nodes        : ", meta$n_nodes, "
", sep = "")
  cat("Edges        : ", meta$n_edges, "
", sep = "")
  cat("Color        : ", meta$color_variable, "
", sep = "")
  cat("Node size    : ", meta$size_label, "
", sep = "")
  cat("Edges        : ", meta$edge_label, "
", sep = "")
  cat("Visual scale : ", signif(meta$visual_scale, 4), " (plot only)
", sep = "")
  cat("
Radial interpretation
")
  cat("------------------------------------------------------------------------
")
  cat("Smaller radius means closer to the disk center.
")
  cat("cor(log node size, hyperbolic radius): ", signif(meta$cor_log_size_hyperbolic_radius, 4), "
", sep = "")
  cat("cor(log node size, Poincare radius)  : ", signif(meta$cor_log_size_poincare_radius, 4), "
", sep = "")
  cat("
Most central observations
")
  cat("------------------------------------------------------------------------
")
  print(format_lhwr_map_central_table(central), row.names = FALSE)
  if (!is.null(meta$group_counts)) {
    cat("
Group counts
")
    cat("------------------------------------------------------------------------
")
    print(meta$group_counts)
  }
  invisible(list(metadata = meta, central = central))
}

format_lhwr_map_central_table <- function(x) {
  data.frame(
    label = x$label,
    name = x$node_name,
    log_size = round(x$log_size, 3),
    hyperbolic_radius = round(x$hyperbolic_radius, 3),
    poincare_radius = round(x$poincare_radius, 3),
    stringsAsFactors = FALSE
  )
}

#' Create a Caption for a Poincare Disk Map
#'
#' Generate a concise caption-style interpretation for a map returned by
#' [plot_lhwr_poincare_map()].
#'
#' @param object An object returned invisibly by [plot_lhwr_poincare_map()].
#' @param digits Number of digits for reported correlations.
#' @return A character string.
#' @export
caption_lhwr_map <- function(object, digits = 3L) {
  meta <- attr(object, "metadata")
  if (is.null(meta)) stop("`object` does not contain LHWR map metadata.", call. = FALSE)
  paste0(
    "Node sizes represent ", meta$size_label, ". ",
    "Smaller radial coordinates indicate observations closer to the center of the Poincare disk. ",
    "In this map, cor(log node size, hyperbolic radius) = ",
    formatC(meta$cor_log_size_hyperbolic_radius, digits = digits, format = "f"),
    " and cor(log node size, Poincare radius) = ",
    formatC(meta$cor_log_size_poincare_radius, digits = digits, format = "f"),
    ". Colors indicate ", meta$color_variable,
    "; ", meta$edge_label, "."
  )
}

lhwr_map_scale <- function(x, min_out = 0.45, max_out = 2.8) {
  x <- as.numeric(x)
  lx <- log1p(pmax(x, 0))
  rng <- range(lx, finite = TRUE)
  if (!all(is.finite(rng)) || diff(rng) < .Machine$double.eps) return(rep(1, length(x)))
  min_out + (lx - rng[1]) / diff(rng) * (max_out - min_out)
}

lhwr_map_point_colors <- function(values, groups, point_palette, group_palette,
                                  significant = NULL,
                                  significance_adjusted = NULL,
                                  significance_direction = NULL,
                                  significance_colors = NULL,
                                  signed_significance_colors = NULL,
                                  adjusted_significance_colors = NULL) {
  if (!is.null(significant)) {
    if (!is.null(significance_adjusted)) {
      f <- lhwr_adjusted_significance_factor(significant, significance_adjusted)
      cols <- rep(adjusted_significance_colors, length.out = 3L)
    } else if (is.null(significance_direction)) {
      f <- lhwr_significance_factor(significant)
      cols <- rep(significance_colors, length.out = 2L)
    } else {
      f <- lhwr_signed_significance_factor(significant, significance_direction)
      cols <- rep(signed_significance_colors, length.out = 3L)
    }
    names(cols) <- levels(f)
    return(unname(cols[as.character(f)]))
  }
  if (!is.null(values)) return(plot_value_colors(as.numeric(values), palette = point_palette))
  if (!is.null(groups)) {
    f <- as.factor(groups)
    lev <- levels(f)
    if (is.null(group_palette)) {
      group_palette <- grDevices::hcl.colors(length(lev), "Dark 3")
      names(group_palette) <- lev
    } else if (!is.null(names(group_palette)) && all(lev %in% names(group_palette))) {
      group_palette <- group_palette[lev]
    } else {
      group_palette <- rep(group_palette, length.out = length(lev))
      names(group_palette) <- lev
    }
    return(unname(group_palette[as.character(f)]))
  }
  rep("#2A9D8F", length.out = if (!is.null(values)) length(values) else if (!is.null(groups)) length(groups) else 1L)
}

lhwr_map_edges <- function(edges, edge_weight, node_id, px, py,
                           edge_alpha_range = c(0.08, 0.30),
                           edge_lwd_range = c(0.25, 1.60)) {
  if (is.null(edges)) return(NULL)
  edges <- as.data.frame(edges, stringsAsFactors = FALSE)
  if (ncol(edges) < 2L) stop("`edges` must have at least two columns.", call. = FALSE)
  from <- edges[[1L]]
  to <- edges[[2L]]
  if (is.null(node_id)) node_id <- seq_along(px)
  idx <- seq_along(px)
  names(idx) <- as.character(node_id)
  if (is.numeric(from) && is.numeric(to) && all(from %in% seq_along(px)) && all(to %in% seq_along(px))) {
    from_idx <- as.integer(from)
    to_idx <- as.integer(to)
  } else {
    from_idx <- unname(idx[as.character(from)])
    to_idx <- unname(idx[as.character(to)])
  }
  ok <- is.finite(from_idx) & is.finite(to_idx)
  from_idx <- from_idx[ok]
  to_idx <- to_idx[ok]
  if (is.null(edge_weight)) edge_weight <- rep(1, length(ok))
  edge_weight <- as.numeric(edge_weight)[ok]
  if (!is.numeric(edge_alpha_range) || length(edge_alpha_range) != 2L) {
    stop("`edge_alpha_range` must be a numeric vector of length 2.", call. = FALSE)
  }
  if (!is.numeric(edge_lwd_range) || length(edge_lwd_range) != 2L) {
    stop("`edge_lwd_range` must be a numeric vector of length 2.", call. = FALSE)
  }
  edge_alpha_range <- sort(edge_alpha_range)
  edge_lwd_range <- sort(edge_lwd_range)
  lwd <- lhwr_map_scale(edge_weight, min_out = edge_lwd_range[1L], max_out = edge_lwd_range[2L])
  alpha <- lhwr_map_scale(edge_weight, min_out = edge_alpha_range[1L], max_out = edge_alpha_range[2L])
  data.frame(
    x0 = px[from_idx], y0 = py[from_idx],
    x1 = px[to_idx], y1 = py[to_idx],
    lwd = lwd, alpha = alpha
  )
}

lhwr_map_disk_background <- function(title, background_col, boundary_col,
                                     degree_cex = 0.50, title_cex = 1.20) {
  graphics::plot.new()
  graphics::plot.window(xlim = c(-1.08, 1.08), ylim = c(-1.36, 1.08), asp = 1)
  theta <- seq(0, 2 * pi, length.out = 720)
  graphics::polygon(cos(theta), sin(theta), col = background_col, border = NA)
  graphics::symbols(0, 0, circles = 1, inches = FALSE, add = TRUE, bg = NA, fg = boundary_col, lwd = 11)
  graphics::symbols(0, 0, circles = 1, inches = FALSE, add = TRUE, bg = NA, fg = "#263238", lwd = 1.1)
  for (a in seq(0, 330, by = 30)) {
    aa <- a * pi / 180
    graphics::segments(0.965 * cos(aa), 0.965 * sin(aa), cos(aa), sin(aa), col = "#8D99AE", lwd = 0.5)
    graphics::text(1.035 * cos(aa), 1.035 * sin(aa), labels = paste0(a, " deg"), cex = degree_cex, col = "#455A64")
  }
  if (!is.null(title) && nzchar(title)) {
    graphics::title(main = title, cex.main = title_cex)
  }
}

lhwr_map_draw_edges <- function(edge_df, edge_col) {
  if (is.null(edge_df) || nrow(edge_df) == 0) return(invisible(NULL))
  for (i in seq_len(nrow(edge_df))) {
    graphics::segments(edge_df$x0[i], edge_df$y0[i], edge_df$x1[i], edge_df$y1[i],
                       col = grDevices::adjustcolor(edge_col, alpha.f = edge_df$alpha[i]),
                       lwd = edge_df$lwd[i])
  }
  invisible(NULL)
}

lhwr_map_colorbar <- function(values, title, palette = "Blue-Red 3",
                               x0 = -0.54, x1 = 0.54, y = -1.155, h = 0.035,
                               cex = 0.55) {
  values <- as.numeric(values)
  rng <- range(values, finite = TRUE)
  if (!all(is.finite(rng)) || diff(rng) < .Machine$double.eps) return(invisible(NULL))
  pal <- grDevices::hcl.colors(100, palette)
  xs <- seq(x0, x1, length.out = length(pal) + 1L)
  for (i in seq_along(pal)) {
    graphics::rect(xs[i], y, xs[i + 1L], y + h, col = pal[i], border = NA)
  }
  graphics::rect(x0, y, x1, y + h, border = "#263238", lwd = 0.35)
  graphics::text(x0, y - 0.018, labels = sprintf("%.2f", rng[1]), adj = c(0.5, 1), cex = cex * 0.86, col = "#263238")
  graphics::text(x1, y - 0.018, labels = sprintf("%.2f", rng[2]), adj = c(0.5, 1), cex = cex * 0.86, col = "#263238")
  graphics::text((x0 + x1) / 2, y + h + 0.018, labels = title, adj = c(0.5, 0), cex = cex, font = 2, col = "#263238")
  invisible(NULL)
}

lhwr_significance_factor <- function(significant) {
  significant <- as.logical(significant)
  out <- ifelse(significant, "Significant", "Not significant")
  out[is.na(significant)] <- "Not significant"
  factor(out, levels = c("Not significant", "Significant"))
}

lhwr_signed_significance_factor <- function(significant, direction) {
  significant <- as.logical(significant)
  significant[is.na(significant)] <- FALSE
  direction <- as.numeric(direction)
  out <- rep("Not significant", length(significant))
  out[significant & is.finite(direction) & direction < 0] <- "Negative significant"
  out[significant & is.finite(direction) & direction > 0] <- "Positive significant"
  factor(out, levels = c("Not significant", "Negative significant", "Positive significant"))
}

lhwr_adjusted_significance_factor <- function(significant, adjusted) {
  significant <- as.logical(significant)
  adjusted <- as.logical(adjusted)
  significant[is.na(significant)] <- FALSE
  adjusted[is.na(adjusted)] <- FALSE
  adjusted <- adjusted & significant
  out <- rep("Not significant", length(significant))
  out[significant & !adjusted] <- "Raw only"
  out[adjusted] <- "Raw and BH"
  factor(out, levels = c("Not significant", "Raw only", "Raw and BH"))
}

lhwr_map_significance_legend <- function(significance_colors,
                                         title = "Local significance",
                                         x0 = -0.42, x1 = 0.42,
                                         y = -1.155, h = 0.045,
                                         cex = 0.55) {
  cols <- rep(significance_colors, length.out = 2L)
  names(cols) <- c("Not significant", "Significant")
  graphics::text((x0 + x1) / 2, y + h + 0.022, labels = title,
                 adj = c(0.5, 0), cex = cex, font = 2, col = "#263238")
  mid <- (x0 + x1) / 2
  box_w <- 0.055
  yb <- y + 0.004
  graphics::rect(x0, yb, x0 + box_w, yb + h, col = cols[["Not significant"]], border = "#263238", lwd = 0.35)
  graphics::text(x0 + box_w + 0.025, yb + h / 2, labels = "Not significant",
                 adj = c(0, 0.5), cex = cex * 0.88, col = "#263238")
  graphics::rect(mid + 0.04, yb, mid + 0.04 + box_w, yb + h, col = cols[["Significant"]], border = "#263238", lwd = 0.35)
  graphics::text(mid + 0.04 + box_w + 0.025, yb + h / 2, labels = "Significant",
                 adj = c(0, 0.5), cex = cex * 0.88, col = "#263238")
  invisible(NULL)
}

lhwr_map_group_legend <- function(groups, cols, x = 1.08, y = 0.92,
                                  title = "Groups",
                                  position = c("right", "bottom"),
                                  cex = 0.55) {
  position <- match.arg(position)
  f <- as.factor(groups)
  lev <- levels(f)
  col_map <- tapply(cols, f, function(z) z[1L])
  missing_cols <- is.na(col_map)
  if (any(missing_cols)) {
    known_pal <- c(
      "Not significant" = "#C9C9C9",
      "Significant" = "#2E7D32",
      "Raw only" = "#E69F00",
      "Raw and BH" = "#1B7837",
      "Negative significant" = "#6A3D9A",
      "Positive significant" = "#1B7837",
      "High-High" = "#D7191C",
      "Low-Low" = "#2C7BB6",
      "Low-High" = "#ABD9E9",
      "High-Low" = "#FDAE61"
    )
    common <- intersect(names(col_map)[missing_cols], names(known_pal))
    col_map[common] <- known_pal[common]
  }
  missing_cols <- is.na(col_map)
  if (any(missing_cols) && exists("lisa_cluster_palette", mode = "function")) {
    lisa_pal <- lisa_cluster_palette()
    common <- intersect(names(col_map)[missing_cols], names(lisa_pal))
    col_map[common] <- lisa_pal[common]
  }
  col_map[is.na(col_map)] <- "#F7F7F7"

  if (position == "bottom") {
    n_col <- if (length(lev) <= 5L) length(lev) else min(4L, length(lev))
    rows <- ceiling(length(lev) / n_col)
    y0 <- if (rows > 1L) -1.315 else -1.225
    row_h <- 0.045
    text_cex <- cex * 0.88
    point_gap <- 0.030
    item_gap <- 0.135
    point_w <- 0.032
    label_w <- graphics::strwidth(lev, cex = text_cex, units = "user")
    if (rows == 1L) {
      item_w <- point_w + point_gap + label_w
      total_w <- sum(item_w) + item_gap * (length(lev) - 1L)
      x0 <- -total_w / 2
      box_w <- total_w
    } else {
      col_w <- if (n_col >= 5L) 0.36 else 0.39
      box_w <- n_col * col_w + 0.08
      x0 <- -box_w / 2
    }
    box_h <- 0.105 + rows * row_h
    graphics::text(x0 + box_w / 2, y0 + box_h - 0.012, labels = title,
                   adj = c(0.5, 0.5), cex = cex * 0.90, font = 2, col = "#263238")
    if (rows == 1L) {
      cursor <- x0
      for (i in seq_along(lev)) {
        px <- cursor + point_w / 2
        py <- y0 + box_h - 0.085
        graphics::points(px, py, pch = 21, bg = col_map[[lev[i]]], col = "#263238", cex = max(1.25, cex * 1.35), lwd = 0.46)
        graphics::text(px + point_gap, py, labels = lev[i], adj = c(0, 0.5), cex = text_cex, col = "#263238")
        cursor <- cursor + item_w[i] + item_gap
      }
      return(invisible(NULL))
    }
    col_w <- if (n_col >= 5L) 0.36 else 0.39
    for (i in seq_along(lev)) {
      row <- floor((i - 1L) / n_col)
      col <- (i - 1L) %% n_col
      px <- x0 + 0.038 + col * col_w
      py <- y0 + box_h - 0.085 - row * row_h
      graphics::points(px, py, pch = 21, bg = col_map[[lev[i]]], col = "#263238", cex = max(1.25, cex * 1.35), lwd = 0.46)
      graphics::text(px + 0.027, py, labels = lev[i], adj = c(0, 0.5), cex = text_cex, col = "#263238")
    }
    return(invisible(NULL))
  }

  n_col <- if (length(lev) >= 9L) 3L else if (length(lev) >= 5L) 2L else 1L
  rows <- ceiling(length(lev) / n_col)
  col_w <- 0.19
  box_w <- 0.08 + n_col * col_w
  box_h <- 0.12 + rows * 0.052
  graphics::text(x, y, labels = title, adj = c(0, 0.5),
                 cex = cex, font = 2, col = "#263238")
  for (i in seq_along(lev)) {
    row <- (i - 1L) %% rows
    col <- floor((i - 1L) / rows)
    px <- x + 0.03 + col * col_w
    py <- y - 0.075 - row * 0.052
    graphics::points(px, py, pch = 21, bg = col_map[[lev[i]]], col = "#263238", cex = 0.72, lwd = 0.35)
    graphics::text(px + 0.032, py, labels = lev[i], adj = c(0, 0.5), cex = cex * 0.90, col = "#263238")
  }
  invisible(NULL)
}

lhwr_map_size_legend <- function(size_label, x = 1.08, y = -0.76, cex = 0.55) {
  graphics::text(x, y + 0.04, labels = "Node size", adj = c(0, 0.5), cex = cex, font = 2, col = "#263238")
  graphics::points(c(x + 0.055, x + 0.19, x + 0.325), rep(y - 0.045, 3),
                   pch = 21, bg = "#BFD7EA", col = "#263238", cex = c(0.55, 1.05, 1.65), lwd = 0.35)
  graphics::text(x + 0.19, y - 0.105, labels = size_label, adj = c(0.5, 0.5), cex = cex * 0.90, col = "#263238")
  invisible(NULL)
}

lhwr_map_info_box <- function(lines, x = 1.08, y = -0.28, cex = 0.52) {
  line_h <- 0.065
  w <- 0.49
  h <- line_h * (length(lines) + 0.8)
  graphics::text(x, y - line_h * (seq_along(lines) - 1L), labels = lines,
                 adj = c(0, 0.5), cex = cex, col = "#263238")
  invisible(NULL)
}

#' @export
plot.hgwr <- plot.lhwr
#' @export
plot.bw.hgwr <- plot.bw.lhwr

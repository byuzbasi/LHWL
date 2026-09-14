#' Plot an LHWPR Coefficient or Focal-Selection Surface
#'
#' Display one fitted LHWPR slope over the Poincare-disk view of the Lorentz
#' representation. A coefficient surface colors nonzero focal coefficients
#' and leaves numerical zeros hollow. A selection surface shows the focal
#' active set for a positively penalized predictor. These displays summarize
#' penalized estimation and are not significance maps.
#'
#' @param x A fitted `"lhwpr"` object.
#' @param predictor Predictor name or one-based slope index. `"Intercept"` is
#'   available only for `type = "coefficient"`.
#' @param type Draw coefficient values or focal selection.
#' @param threshold Nonnegative numerical threshold used to distinguish stored
#'   zeros from active coefficients.
#' @param labels,node_names,sizes Optional point labels, longer names, and
#'   nonnegative point-size values; see [plot_lhwr_poincare_map()].
#' @param title Optional plot title. `NULL` creates a descriptive default.
#' @param point_palette HCL palette for nonzero coefficient values.
#' @param selection_colors Two colors for not-selected and selected locations.
#' @param show_legend Logical; show the selection legend or the
#'   filled-versus-hollow coefficient key.
#' @param radial_rescale,radial_quantile,radial_power Display-only Poincare
#'   radial controls; they do not alter fitting distances.
#' @param output_file Optional PNG, PDF, EPS, or PS output path.
#' @param width,height,res Output-device dimensions in inches and resolution
#'   for raster output.
#' @param ... Additional arguments passed to [plot_lhwr_poincare_map()].
#' @return Invisibly, a data frame containing plotted coordinates, the focal
#'   coefficient, the active indicator, predictor, threshold, and plot type.
#' @examples
#' set.seed(7)
#' H <- simulate_H2(35, r_max = 2)
#' dat <- data.frame(y = rnorm(35), x1 = rnorm(35), x2 = rnorm(35))
#' fit <- lhwpr(
#'   y ~ x1 + x2, dat, H$Z, bw = 18, lambda = c(0.2, 0.05),
#'   keep_path = "selected"
#' )
#' plot_lhwpr_surface(fit, "x1", type = "coefficient")
#' plot_lhwpr_surface(fit, "x1", type = "selection")
#' @export
plot_lhwpr_surface <- function(
    x,
    predictor,
    type = c("coefficient", "selection"),
    threshold = 1e-8,
    labels = NULL,
    node_names = NULL,
    sizes = NULL,
    title = NULL,
    point_palette = "Blue-Red 3",
    selection_colors = c(
      "Not selected" = "#D9D9D9", "Selected" = "#2166AC"
    ),
    show_legend = TRUE,
    radial_rescale = TRUE,
    radial_quantile = 1,
    radial_power = 1,
    output_file = NULL,
    width = 7,
    height = 7,
    res = 220,
    ...) {
  if (!inherits(x, "lhwpr")) {
    stop("`x` must inherit from class 'lhwpr'.", call. = FALSE)
  }
  type <- match.arg(type)
  if (!is.numeric(threshold) || length(threshold) != 1L ||
      !is.finite(threshold) || threshold < 0) {
    stop("`threshold` must be one non-negative finite value.", call. = FALSE)
  }
  show_legend <- lhwpr_logical_scalar(show_legend, "show_legend")
  coefficient_names <- colnames(x$beta_hat)
  if (is.null(coefficient_names)) {
    coefficient_names <- c(
      if (isTRUE(x$intercept)) "Intercept", x$predictor_names
    )
  }
  coefficient_index <- lhwpr_surface_coefficient_index(
    predictor, coefficient_names, x$intercept
  )
  coefficient_name <- coefficient_names[[coefficient_index]]
  slope_index <- coefficient_index - as.integer(isTRUE(x$intercept))
  if (type == "selection") {
    if (slope_index < 1L) {
      stop(
        "The unpenalized intercept does not have a focal-selection surface.",
        call. = FALSE
      )
    }
    if (!is.null(x$penalty_factor) && x$penalty_factor[[slope_index]] <= 0) {
      stop(
        "A focal-selection surface requires a positively penalized predictor; use `type = \"coefficient\"` for a forced predictor.",
        call. = FALSE
      )
    }
  }
  coefficient <- as.numeric(x$beta_hat[, coefficient_index])
  active <- is.finite(coefficient) & abs(coefficient) > threshold
  if (is.null(title)) {
    title <- if (type == "selection") {
      paste0("Focal selection: ", coefficient_name)
    } else {
      paste0("Local coefficient: ", coefficient_name)
    }
  }

  if (!is.null(output_file)) {
    device <- tolower(tools::file_ext(output_file))
    device <- if (device %in% c("eps", "ps")) {
      "eps"
    } else if (device == "pdf") {
      "pdf"
    } else {
      "png"
    }
    lhwr_open_output_device(
      output_file, device = device, width = width, height = height, res = res
    )
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  if (type == "selection") {
    status <- factor(
      ifelse(active, "Selected", "Not selected"),
      levels = c("Not selected", "Selected")
    )
    plotted <- plot_lhwr_poincare_map(
      x$Z,
      groups = status,
      labels = labels,
      node_names = node_names,
      sizes = sizes,
      title = title,
      group_label = "Focal selection",
      group_palette = selection_colors,
      show_group_legend = show_legend,
      show_size_legend = FALSE,
      show_info = FALSE,
      radial_rescale = radial_rescale,
      radial_quantile = radial_quantile,
      radial_power = radial_power,
      ...
    )
  } else {
    displayed <- coefficient
    displayed[!active] <- NA_real_
    plotted <- plot_lhwr_poincare_map(
      x$Z,
      values = displayed,
      labels = labels,
      node_names = node_names,
      sizes = sizes,
      title = title,
      value_label = paste0("Coefficient: ", coefficient_name),
      point_palette = point_palette,
      show_colorbar = show_legend,
      show_size_legend = FALSE,
      show_info = FALSE,
      radial_rescale = radial_rescale,
      radial_quantile = radial_quantile,
      radial_power = radial_power,
      ...
    )
    if (show_legend) {
      graphics::legend(
        "topleft",
        legend = c("Nonzero coefficient", "Numerical zero"),
        pch = 21,
        pt.bg = c("#4C78A8", NA),
        col = "#263238",
        bty = "n",
        cex = 0.72
      )
    }
  }
  plotted$coefficient <- coefficient
  plotted$active <- active
  plotted$predictor <- coefficient_name
  plotted$threshold <- threshold
  plotted$type <- type
  invisible(plotted)
}

lhwpr_surface_coefficient_index <- function(predictor, names, intercept) {
  if (is.character(predictor) && length(predictor) == 1L &&
      !is.na(predictor)) {
    index <- match(predictor, names)
    if (is.na(index)) {
      stop(
        "`predictor` is not a fitted coefficient name: ", predictor, ".",
        call. = FALSE
      )
    }
    return(index)
  }
  if (!is.numeric(predictor) || length(predictor) != 1L ||
      !is.finite(predictor) || predictor != as.integer(predictor)) {
    stop(
      "`predictor` must be one fitted name or one-based slope index.",
      call. = FALSE
    )
  }
  slope_index <- as.integer(predictor)
  n_slopes <- length(names) - as.integer(isTRUE(intercept))
  if (slope_index < 1L || slope_index > n_slopes) {
    stop("The slope index in `predictor` is outside the fitted model.", call. = FALSE)
  }
  slope_index + as.integer(isTRUE(intercept))
}

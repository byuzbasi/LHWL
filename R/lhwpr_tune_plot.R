#' Plot an LHWPR Tuning Object
#'
#' Display separate bandwidth and lambda profiles, the joint
#' bandwidth--lambda surface, the exact evaluated-candidate landscape, and the
#' alpha profile. The landscape distinguishes coarse, refined, ineligible, and
#' selected candidates without interpolating unevaluated locations. When
#' several alpha values were evaluated, `type = "surface"` or
#' `type = "landscape"` draws one panel per alpha.
#'
#' @param x A `"tune.lhwpr"` object.
#' @param type Plot type.
#' @param criterion Score to display. `NULL` uses the selection criterion.
#' @param alpha,gamma Optional surface slice. `alpha = NULL` facets the
#'   standalone surface plot over all evaluated alpha values.
#' @param include_lhwr Logical; add the same-fold unpenalized LHWR reference to
#'   the bandwidth profile when it is available.
#' @param ... Additional graphical arguments passed to the profile plots.
#' @return Invisibly, the data used by the requested plot.
#' @export
plot.tune.lhwpr <- function(
    x,
    type = c(
      "all", "bandwidth", "lambda", "landscape", "surface", "alpha"
    ),
    criterion = NULL,
    alpha = NULL,
    gamma = NULL,
    include_lhwr = TRUE,
    ...) {
  if (!inherits(x, "tune.lhwpr")) {
    stop("`x` must inherit from class 'tune.lhwpr'.", call. = FALSE)
  }
  type <- match.arg(type)
  include_lhwr <- lhwpr_logical_scalar(include_lhwr, "include_lhwr")
  criterion <- if (is.null(criterion)) {
    x$criterion
  } else {
    match.arg(
      tolower(criterion),
      c("cv", "aicc", "aic", "bic", "gcv", "ebic", "ebic_refit")
    )
  }
  score_column <- lhwpr_criterion_column(criterion)
  if (!any(is.finite(x$selection_table[[score_column]]))) {
    stop(
      "The requested criterion is not available for this tuning object.",
      call. = FALSE
    )
  }
  eligibility <- lhwpr_tuning_eligibility(
    x$selection_table,
    criterion = criterion,
    kkt_tol = x$eligibility_kkt_tol %||% 1e-5
  )
  x$selection_table$eligible <- eligibility$eligible
  x$selection_table$ineligibility_reason <- eligibility$reason

  if (type == "all") {
    old <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old), add = TRUE)
    graphics::par(mfrow = c(2, 2), mar = c(4.2, 4.4, 3.1, 1.2))
    bandwidth_data <- lhwpr_plot_tune_bandwidth(
      x, score_column, criterion, include_lhwr, ...
    )
    lambda_data <- lhwpr_plot_tune_lambda(
      x, score_column, criterion, ...
    )
    landscape_data <- lhwpr_plot_tune_landscape(
      x,
      score_column,
      criterion,
      alpha = x$selected$alpha,
      gamma = x$selected$gamma
    )
    alpha_data <- lhwpr_plot_tune_alpha(
      x, score_column, criterion, ...
    )
    return(invisible(list(
      bandwidth = bandwidth_data,
      lambda = lambda_data,
      landscape = landscape_data,
      alpha = alpha_data
    )))
  }
  if (type == "bandwidth") {
    return(invisible(lhwpr_plot_tune_bandwidth(
      x, score_column, criterion, include_lhwr, ...
    )))
  }
  if (type == "lambda") {
    return(invisible(lhwpr_plot_tune_lambda(
      x, score_column, criterion, ...
    )))
  }
  if (type == "alpha") {
    return(invisible(lhwpr_plot_tune_alpha(
      x, score_column, criterion, ...
    )))
  }

  alphas <- if (is.null(alpha)) {
    sort(unique(x$selection_table$alpha))
  } else {
    as.numeric(alpha)
  }
  old <- NULL
  if (length(alphas) > 1L) {
    old <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old), add = TRUE)
    columns <- min(3L, length(alphas))
    graphics::par(
      mfrow = c(ceiling(length(alphas) / columns), columns),
      mar = c(4.2, 4.4, 3.1, 1.2)
    )
  }
  result <- lapply(
    alphas,
    function(value) {
      if (type == "landscape") {
        lhwpr_plot_tune_landscape(
          x, score_column, criterion, alpha = value, gamma = gamma
        )
      } else {
        lhwpr_plot_tune_surface(
          x, score_column, criterion, alpha = value, gamma = gamma
        )
      }
    }
  )
  invisible(if (length(result) == 1L) result[[1L]] else result)
}

lhwpr_tune_score_label <- function(criterion) {
  switch(
    criterion,
    cv = "Cross-validation loss",
    aic = "AIC",
    aicc = "AICc",
    bic = "BIC",
    gcv = "GCV",
    ebic = "EBIC",
    ebic_refit = "Refit EBIC"
  )
}

lhwpr_tune_profile <- function(table, group, score_column) {
  if ("eligible" %in% names(table)) {
    table <- table[!is.na(table$eligible) & table$eligible, , drop = FALSE]
  }
  values <- split(table[[score_column]], table[[group]])
  data.frame(
    value = as.numeric(names(values)),
    score = vapply(
      values,
      function(score) {
        finite <- score[is.finite(score)]
        if (length(finite)) min(finite) else Inf
      },
      numeric(1)
    ),
    stringsAsFactors = FALSE
  )
}

lhwpr_plot_tune_boundary_note <- function(x) {
  if (isTRUE(x$bandwidth_search$selected_at_boundary)) {
    graphics::mtext(
      paste0(
        "Selected bandwidth is on the ",
        x$selected$bw_boundary,
        " boundary"
      ),
      side = 3,
      line = 0.15,
      col = "#B22222",
      cex = 0.75
    )
  }
}

lhwpr_plot_tune_bandwidth <- function(
    x, score_column, criterion, include_lhwr, ...) {
  profile <- lhwpr_tune_profile(
    x$selection_table, "bw", score_column
  )
  profile <- profile[order(profile$value), , drop = FALSE]
  graphics::plot(
    profile$value,
    profile$score,
    type = "b",
    pch = 19,
    col = "#0868AC",
    xlab = if (x$bandwidth_search$mode == "automatic" ||
                isTRUE(x$fit$adaptive)) {
      "Adaptive bandwidth (nearest neighbors)"
    } else {
      "Bandwidth"
    },
    ylab = lhwpr_tune_score_label(criterion),
    main = paste("Bandwidth selection by", toupper(criterion)),
    ...
  )
  graphics::abline(
    v = x$selected$bw, col = "#CB181D", lty = 2, lwd = 2
  )
  legend_labels <- "LHWPR profile"
  legend_colors <- "#0868AC"
  legend_types <- 1
  if (include_lhwr && !is.null(x$lhwr_comparison) &&
      any(is.finite(x$lhwr_comparison[[score_column]]))) {
    reference <- lhwpr_tune_profile(
      x$lhwr_comparison, "bw", score_column
    )
    reference <- reference[order(reference$value), , drop = FALSE]
    graphics::lines(
      reference$value,
      reference$score,
      type = "b",
      pch = 17,
      col = "#636363",
      lty = 3
    )
    legend_labels <- c(legend_labels, "LHWR reference")
    legend_colors <- c(legend_colors, "#636363")
    legend_types <- c(legend_types, 3)
  }
  graphics::legend(
    "topright",
    legend = legend_labels,
    col = legend_colors,
    lty = legend_types,
    pch = c(19, 17)[seq_along(legend_labels)],
    bty = "n",
    cex = 0.8
  )
  lhwpr_plot_tune_boundary_note(x)
  profile
}

lhwpr_same_gamma <- function(value, target) {
  (is.na(value) & is.na(target)) |
    (!is.na(value) & !is.na(target) &
       abs(value - target) <=
       sqrt(.Machine$double.eps) * (1 + abs(target)))
}

lhwpr_plot_tune_lambda <- function(x, score_column, criterion, ...) {
  table <- x$selection_table
  keep <- table$bw == x$selected$bw &
    table$alpha == x$selected$alpha &
    lhwpr_same_gamma(table$gamma, x$selected$gamma)
  path <- table[keep, , drop = FALSE]
  path <- path[order(path$lambda, decreasing = TRUE), , drop = FALSE]
  eligible <- !is.na(path$eligible) & path$eligible
  graphics::plot(
    log10(pmax(path$lambda, .Machine$double.xmin)),
    path[[score_column]],
    type = "n",
    pch = 19,
    col = "#0868AC",
    xlab = expression(log[10](lambda)),
    ylab = lhwpr_tune_score_label(criterion),
    main = paste("Lambda selection by", toupper(criterion)),
    ...
  )
  if (any(eligible)) {
    graphics::lines(
      log10(pmax(path$lambda[eligible], .Machine$double.xmin)),
      path[[score_column]][eligible],
      type = "b", pch = 19, col = "#0868AC"
    )
  }
  if (any(!eligible)) {
    graphics::points(
      log10(pmax(path$lambda[!eligible], .Machine$double.xmin)),
      path[[score_column]][!eligible],
      pch = 4, col = "#969696"
    )
  }
  graphics::abline(
    v = log10(max(x$selected$lambda, .Machine$double.xmin)),
    col = "#CB181D",
    lty = 2,
    lwd = 2
  )
  path
}

lhwpr_tune_surface_slice <- function(x, alpha, gamma = NULL) {
  table <- x$selection_table
  if (is.null(gamma)) {
    gamma <- if (alpha == x$selected$alpha) {
      x$selected$gamma
    } else {
      candidates <- table$gamma[table$alpha == alpha]
      candidates[1L]
    }
  }
  keep <- table$alpha == alpha &
    lhwpr_same_gamma(table$gamma, gamma) &
    table$lambda_ratio > 0
  list(
    table = table[keep, , drop = FALSE],
    gamma = gamma
  )
}

lhwpr_plot_tune_landscape <- function(
    x, score_column, criterion, alpha, gamma = NULL) {
  slice <- lhwpr_tune_surface_slice(x, alpha, gamma)
  landscape <- slice$table
  gamma <- slice$gamma
  finite <- is.finite(landscape[[score_column]])
  eligible <- finite & !is.na(landscape$eligible) & landscape$eligible
  if (!any(finite)) {
    graphics::plot.new()
    graphics::title(main = paste0("No finite landscape for alpha = ", alpha))
    return(landscape)
  }
  score <- landscape[[score_column]]
  delta <- rep(NA_real_, length(score))
  if (any(eligible)) {
    delta[eligible] <- score[eligible] - min(score[eligible])
  }
  cap <- if (sum(eligible) > 1L) {
    as.numeric(stats::quantile(delta[eligible], 0.95, names = FALSE, type = 8))
  } else {
    1
  }
  if (!is.finite(cap) || cap <= 0) cap <- 1
  color_index <- pmax(
    1L,
    pmin(100L, 1L + floor(99 * pmin(delta, cap) / cap))
  )
  palette <- grDevices::hcl.colors(100, "YlGnBu", rev = TRUE)
  stage <- as.character(landscape$stage)
  pch <- ifelse(stage == "refined", 17, 16)
  graphics::plot(
    landscape$bw,
    log10(landscape$lambda_ratio),
    type = "n",
    xlab = "Bandwidth",
    ylab = expression(log[10](lambda/lambda[max])),
    main = paste0(toupper(criterion), " candidate landscape, alpha = ", alpha)
  )
  if (any(eligible)) {
    graphics::points(
      landscape$bw[eligible],
      log10(landscape$lambda_ratio[eligible]),
      pch = pch[eligible],
      col = palette[color_index[eligible]],
      cex = 0.85
    )
  }
  ineligible <- finite & !eligible
  if (any(ineligible)) {
    graphics::points(
      landscape$bw[ineligible],
      log10(landscape$lambda_ratio[ineligible]),
      pch = 4, col = "#969696", cex = 0.7
    )
  }
  selected_slice <- alpha == x$selected$alpha &&
    lhwpr_same_gamma(gamma, x$selected$gamma)
  if (selected_slice) {
    graphics::points(
      x$selected$bw,
      log10(max(x$selected$lambda_ratio, .Machine$double.xmin)),
      pch = 23, bg = "#CB181D", col = "#7F0000", cex = 1.45, lwd = 1.1
    )
    lhwpr_plot_tune_boundary_note(x)
  }
  graphics::legend(
    "topright",
    legend = c("Eligible coarse", "Eligible refined", "Ineligible", "Selected"),
    pch = c(16, 17, 4, 23),
    col = c("#0868AC", "#0868AC", "#969696", "#7F0000"),
    pt.bg = c(NA, NA, NA, "#CB181D"),
    bty = "n", cex = 0.72
  )
  landscape$plot_eligible <- eligible
  landscape$score_delta <- delta
  landscape
}

lhwpr_plot_tune_surface <- function(
    x, score_column, criterion, alpha, gamma = NULL) {
  slice <- lhwpr_tune_surface_slice(x, alpha, gamma)
  gamma <- slice$gamma
  candidates <- slice$table
  keep <- is.finite(candidates[[score_column]]) &
    !is.na(candidates$eligible) & candidates$eligible
  surface <- candidates[keep, , drop = FALSE]
  if (!nrow(surface)) {
    graphics::plot.new()
    graphics::title(main = paste0("No finite surface for alpha = ", alpha))
    return(surface)
  }
  bandwidth <- sort(unique(surface$bw))
  ratio <- sort(unique(surface$lambda_ratio))
  complete <- nrow(surface) == length(bandwidth) * length(ratio)
  palette <- grDevices::hcl.colors(24, "YlGnBu", rev = TRUE)
  title <- paste0(
    toupper(criterion),
    " surface, alpha = ",
    format(alpha)
  )

  if (complete && length(bandwidth) > 1L && length(ratio) > 1L) {
    score <- matrix(
      NA_real_, nrow = length(bandwidth), ncol = length(ratio)
    )
    for (row in seq_len(nrow(surface))) {
      score[
        match(surface$bw[row], bandwidth),
        match(surface$lambda_ratio[row], ratio)
      ] <- surface[[score_column]][row]
    }
    graphics::image(
      bandwidth,
      log10(ratio),
      score,
      col = palette,
      xlab = "Bandwidth",
      ylab = expression(log[10](lambda/lambda[max])),
      main = title
    )
    if (length(unique(as.numeric(score))) > 1L) {
      graphics::contour(
        bandwidth,
        log10(ratio),
        score,
        add = TRUE,
        drawlabels = FALSE,
        col = grDevices::adjustcolor("white", 0.75)
      )
    }
  } else {
    score <- surface[[score_column]]
    breaks <- unique(stats::quantile(
      score,
      probs = seq(0, 1, length.out = length(palette) + 1),
      na.rm = TRUE
    ))
    scaled <- if (length(breaks) < 2L) {
      rep(1L, length(score))
    } else {
      cut(
        score,
        breaks = breaks,
        include.lowest = TRUE,
        labels = FALSE
      )
    }
    scaled[is.na(scaled)] <- 1L
    graphics::plot(
      surface$bw,
      log10(surface$lambda_ratio),
      pch = 15,
      cex = 1.4,
      col = palette[scaled],
      xlab = "Bandwidth",
      ylab = expression(log[10](lambda/lambda[max])),
      main = title
    )
  }
  if (alpha == x$selected$alpha &&
      lhwpr_same_gamma(gamma, x$selected$gamma)) {
    graphics::points(
      x$selected$bw,
      log10(max(x$selected$lambda_ratio, .Machine$double.xmin)),
      pch = 4,
      cex = 1.8,
      lwd = 2.5,
      col = "#CB181D"
    )
    lhwpr_plot_tune_boundary_note(x)
  }
  ineligible <- candidates[
    is.finite(candidates[[score_column]]) &
      (is.na(candidates$eligible) | !candidates$eligible),
    , drop = FALSE
  ]
  if (nrow(ineligible)) {
    graphics::points(
      ineligible$bw,
      log10(ineligible$lambda_ratio),
      pch = 4, col = "#969696", cex = 0.65
    )
  }
  surface
}

lhwpr_plot_tune_alpha <- function(x, score_column, criterion, ...) {
  profile <- lhwpr_tune_profile(
    x$selection_table, "alpha", score_column
  )
  profile <- profile[order(profile$value), , drop = FALSE]
  if (nrow(profile) == 1L) {
    table <- x$selection_table
    keep <- table$bw == x$selected$bw &
      table$alpha == x$selected$alpha &
      lhwpr_same_gamma(table$gamma, x$selected$gamma)
    path <- table[keep, , drop = FALSE]
    active <- if (any(is.finite(path$surface_size))) {
      path$surface_size
    } else {
      path$mean_active
    }
    graphics::plot(
      log10(pmax(path$lambda, .Machine$double.xmin)),
      active,
      type = "b",
      pch = 19,
      col = "#238B45",
      xlab = expression(log[10](lambda)),
      ylab = "Selected coefficient surfaces",
      main = "Sparsity path",
      ...
    )
    return(path)
  }
  graphics::plot(
    profile$value,
    profile$score,
    type = "b",
    pch = 19,
    col = "#6A51A3",
    xlab = expression(alpha),
    ylab = lhwpr_tune_score_label(criterion),
    main = paste("Alpha selection by", toupper(criterion)),
    ...
  )
  graphics::abline(
    v = x$selected$alpha, col = "#CB181D", lty = 2, lwd = 2
  )
  profile
}

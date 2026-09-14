#' Extract and Predict from LHWR Models
#'
#' Convenience S3 methods for fitted LHWR models. These methods make LHWR
#' objects behave like ordinary R model objects for common downstream tasks.
#' `predict.lhwr()` returns fitted values for the calibration data when
#' `newdata` is omitted. For new observations, supply both the Euclidean
#' covariates in `newdata` and their Lorentz-hyperboloid coordinates in `Znew`.
#' Prediction is then obtained by fitting a local weighted Gaussian or
#' binomial model at each new Lorentz location using only the original
#' calibration data.
#'
#' Tangent-E prediction reuses the stored training reference point for every
#' prediction batch. If an older object has no stored reference, it is
#' reconstructed from `object$Z`, never from `Znew`. Identical mapped points
#' have exactly zero distance; nearby distinct tangent coordinates are not
#' rounded to zero by a pairwise-distance tolerance. With unchanged model
#' settings, predicting the training rows reproduces the calibration fit up
#' to floating-point error.
#'
#' @param object An LHWR model object returned by [lhwr()].
#' @param newdata Optional data frame containing predictor variables for new
#'   observations. If omitted, `predict()` returns fitted values from the
#'   calibration data.
#' @param Znew Optional `m x 3` Lorentz-hyperboloid coordinate matrix for
#'   `newdata`. Required when `newdata` is supplied.
#' @param type Prediction type. `"response"` returns fitted responses
#'   (probabilities for binomial LHWR), `"link"` returns linear predictors, and
#'   `"coefficients"` returns the local coefficient matrix used for prediction.
#' @param dMat Optional `m x n` distance matrix from new observations to the
#'   original calibration observations. If supplied, it is used directly.
#' @param bw Optional prediction bandwidth. By default the bandwidth stored in
#'   the fitted model is used.
#' @param kernel Optional kernel name. By default the fitted model kernel is
#'   used.
#' @param adaptive Optional logical bandwidth flag. By default the fitted model
#'   setting is used.
#' @param bandwidth_type Optional bandwidth type. By default the fitted model
#'   setting is used.
#' @param geometry Optional prediction geometry. By default the fitted model
#'   geometry is used.
#' @param ... Unused.
#' @return `fitted()` returns the fitted response vector; `residuals()` returns
#'   residuals; `coef()` returns the `n x p` calibration local coefficient
#'   matrix. `predict()` returns either a numeric response vector or an `m x p`
#'   matrix of local coefficients for `newdata`.
#' @examples
#' set.seed(11)
#' H <- simulate_H2(50, r_max = 3)
#' x <- rnorm(50)
#' y <- 1 + 2 * x + rnorm(50, sd = 0.3)
#' fit <- lhwr(y ~ x, data.frame(y, x), Z = H$Z, bw = 20)
#' head(fitted(fit))
#' head(residuals(fit))
#' head(coef(fit))
#'
#' Znew <- simulate_H2(3, r_max = 3)$Z
#' newdata <- data.frame(x = c(-1, 0, 1))
#' predict(fit, newdata = newdata, Znew = Znew)
#' @export
fitted.lhwr <- function(object, ...) {
  if (!is_lhwr_object(object)) stop('`object` must be an LHWR model object.', call. = FALSE)
  object$y_hat
}

#' @rdname fitted.lhwr
#' @export
residuals.lhwr <- function(object, ...) {
  if (!is_lhwr_object(object)) stop('`object` must be an LHWR model object.', call. = FALSE)
  object$residuals
}

#' @rdname fitted.lhwr
#' @export
coef.lhwr <- function(object, ...) {
  if (!is_lhwr_object(object)) stop('`object` must be an LHWR model object.', call. = FALSE)
  object$beta_hat
}

#' @rdname fitted.lhwr
#' @export
predict.lhwr <- function(object, newdata = NULL, Znew = NULL,
                         type = c("response", "link", "coefficients"),
                         dMat = NULL, bw = NULL, kernel = NULL,
                         adaptive = NULL, bandwidth_type = NULL,
                         geometry = NULL, ...) {
  if (!is_lhwr_object(object)) stop('`object` must be an LHWR model object.', call. = FALSE)
  type <- match.arg(type)
  if (is.null(newdata)) {
    if (type == "response") return(object$y_hat)
    if (type == "link") {
      return(object$linear_predictor %||% object$y_hat)
    }
    return(object$beta_hat)
  }

  args <- object$LHWR.arguments %||% object$HGWR.arguments
  if (is.null(args)) {
    stop("The LHWR object does not contain model arguments required for prediction.", call. = FALSE)
  }
  kernel <- kernel %||% args$kernel
  adaptive <- adaptive %||% args$adaptive
  bandwidth_type <- bandwidth_type %||% args$bandwidth_type
  geometry <- geometry %||% args$geometry
  bw <- bw %||% args$bw
  family <- object$family %||% args$family %||% "gaussian"

  kernel <- match.arg(kernel, hgwr_kernel_names())
  bandwidth_type <- match.arg(bandwidth_type, c("knn", "volume"))
  geometry <- match.arg(geometry, c("lorentz", "tangent"))
  adaptive <- isTRUE(adaptive)

  X_train <- stats::model.matrix(object$lm)
  colnames(X_train) <- clean_coef_names(colnames(X_train))
  y_train <- if (family == "binomial") {
    as.numeric(object$SDF$y)
  } else {
    stats::model.response(stats::model.frame(object$lm))
  }
  X_new <- lhwr_predict_model_matrix(object, newdata)
  if (!identical(colnames(X_new), colnames(X_train))) {
    stop("`newdata` does not produce the same model matrix columns as the fitted LHWR model.", call. = FALSE)
  }

  if (is.null(dMat)) {
    if (is.null(Znew)) {
      stop("`Znew` is required when predicting at new observations unless `dMat` is supplied.", call. = FALSE)
    }
    Znew <- as_lorentz_matrix(Znew, "Znew")
    check_hyperboloid(Znew, stop_on_fail = TRUE)
    if (nrow(Znew) != nrow(X_new)) {
      stop("`Znew` must have the same number of rows as `newdata`.", call. = FALSE)
    }
    D_new <- lhwr_predict_distance(object, Znew, geometry = geometry)
  } else {
    D_new <- as.matrix(dMat)
    storage.mode(D_new) <- "double"
    if (nrow(D_new) != nrow(X_new) || ncol(D_new) != nrow(X_train)) {
      stop("`dMat` must be an m x n matrix of distances from new observations to calibration observations.", call. = FALSE)
    }
  }

  bw_new <- lhwr_predict_bandwidth(
    D_new,
    bw = bw,
    adaptive = adaptive,
    bandwidth_type = bandwidth_type,
    training_bandwidth = object$bandwidth
  )
  beta_new <- if (family == "binomial") {
    prediction_neighborhood <- list(
      engine = "dense",
      D = t(D_new),
      bandwidth = bw_new,
      index = NULL,
      distance = NULL
    )
    prediction_fit <- lhwr_binomial_fit_native(
      y = y_train,
      X = X_train,
      X_focus = X_new,
      neighborhood = prediction_neighborhood,
      intercept = isTRUE(object$intercept),
      kernel = kernel,
      n_threads = args$n_threads %||% 1L,
      control = args$control %||% stats::glm.control(maxit = 100L),
      return_inference = FALSE
    )
    if (any(!prediction_fit$converged)) {
      warning(
        sum(!prediction_fit$converged),
        " new-location binomial fits did not meet the convergence criterion.",
        call. = FALSE
      )
    }
    prediction_fit$beta_hat
  } else {
    lhwr_predict_local_coefficients(
      y_train, X_train, D_new, bw_new, kernel
    )
  }
  if (type == "coefficients") return(beta_new)
  eta <- as.numeric(rowSums(X_new * beta_new))
  if (type == "link" || family == "gaussian") return(eta)
  stats::plogis(eta)
}

#' @rdname fitted.lhwr
#' @export
fitted.hgwr <- fitted.lhwr

#' @rdname fitted.lhwr
#' @export
residuals.hgwr <- residuals.lhwr

#' @rdname fitted.lhwr
#' @export
coef.hgwr <- coef.lhwr

#' @rdname fitted.lhwr
#' @export
predict.hgwr <- predict.lhwr

lhwr_predict_model_matrix <- function(object, newdata) {
  tt <- stats::delete.response(stats::terms(object$lm))
  mf <- stats::model.frame(
    tt,
    data = newdata,
    na.action = stats::na.pass,
    xlev = object$lm$xlevels
  )
  X <- stats::model.matrix(tt, mf, contrasts.arg = object$lm$contrasts)
  colnames(X) <- clean_coef_names(colnames(X))
  X
}

lhwr_predict_distance <- function(object, Znew, geometry = c("lorentz", "tangent")) {
  geometry <- match.arg(geometry)
  if (geometry == "lorentz") {
    return(hyperbolic_distance_matrix(Znew, object$Z, check = FALSE))
  }

  mu <- object$tangent_ref
  if (is.null(mu)) {
    args <- object$LHWR.arguments %||% object$HGWR.arguments
    mu <- args$tangent_ref
  }
  if (is.null(mu)) {
    mu <- extrinsic_lorentz_mean(object$Z)
  }
  train_coords <- tangent_coordinates(object$Z, mu = mu)
  new_coords <- tangent_coordinates(Znew, mu = mu)
  # Coordinate differences avoid cancellation in ||a||^2 + ||b||^2 - 2<a,b>.
  # In particular, identical points must have exactly zero distance: a tiny
  # positive self-distance would count as an extra adaptive-kNN neighbour.
  dx <- abs(outer(new_coords[, 1L], train_coords[, 1L], "-"))
  dy <- abs(outer(new_coords[, 2L], train_coords[, 2L], "-"))
  scale <- pmax(dx, dy)
  D <- matrix(0, nrow(new_coords), nrow(train_coords))
  positive <- scale > 0
  D[positive] <- scale[positive] * sqrt(
    (dx[positive] / scale[positive])^2 +
      (dy[positive] / scale[positive])^2
  )
  D
}

lhwr_predict_bandwidth <- function(D_new, bw, adaptive, bandwidth_type,
                                   training_bandwidth = NULL) {
  if (length(bw) != 1L) {
    if (!is.null(training_bandwidth) && length(training_bandwidth) == 1L) {
      bw <- training_bandwidth
    } else {
      stop("Prediction requires a scalar bandwidth or adaptive neighbor count.", call. = FALSE)
    }
  }

  if (isTRUE(adaptive)) {
    k <- as.integer(round(bw))
    n_train <- ncol(D_new)
    if (k < 1L || k > n_train) {
      stop("Adaptive prediction bandwidth must be between 1 and the number of calibration observations.", call. = FALSE)
    }
    if (bandwidth_type == "volume") {
      warning("Prediction with `bandwidth_type = \"volume\"` uses kNN radii for new observations.", call. = FALSE)
    }
    apply(D_new, 1L, function(d) {
      positive <- sort(d[is.finite(d) & d > 0])
      if (length(positive) < k) {
        stop(
          "Fewer than `k` positive finite calibration distances are available for prediction.",
          call. = FALSE
        )
      }
      positive[k]
    })
  } else {
    if (bandwidth_type == "volume") {
      rep(hyperbolic_radius_from_volume(bw), nrow(D_new))
    } else {
      check_bandwidth(bw)
      rep(as.numeric(bw), nrow(D_new))
    }
  }
}

lhwr_predict_local_coefficients <- function(y, X, D_new, bw_new, kernel) {
  p <- ncol(X)
  beta <- matrix(NA_real_, nrow(D_new), p)
  colnames(beta) <- colnames(X)
  for (i in seq_len(nrow(D_new))) {
    w <- as.numeric(hgwr_kernel(D_new[i, ], bw_new[i], kernel = kernel))
    if (any(!is.finite(w))) {
      stop("Prediction produced non-finite kernel weights.", call. = FALSE)
    }
    XtWX <- crossprod(X, X * w)
    XtWy <- crossprod(X, w * y)
    beta_i <- tryCatch(
      solve(XtWX, XtWy),
      error = function(e) MASS::ginv(XtWX) %*% XtWy
    )
    if (any(!is.finite(beta_i))) {
      beta_i <- MASS::ginv(XtWX) %*% XtWy
    }
    beta[i, ] <- as.numeric(beta_i)
  }
  beta
}

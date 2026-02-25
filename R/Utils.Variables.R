#' @title
#' Construct determinant variables for [kpss_single]
#'
#' @details
#' Procedure to compute deterministic terms
#' for KPSS with 1 structural break.
#'
#' @param model A scalar equal to
#' * 1: Model with trend, break in const,
#' * 2: Model with const and trend, break in const,
#' * 3: Model with const and trend, break in trend,
#' * 4: Model with const and trend, break in const and trend.
#' @param n.obs Number of observations.
#' @param break.point Break point.
#'
#' @return Matrix of determinant variables.
#'
#' @keywords internal
trend_kpss_single <- function(model,
                              n_obs,
                              break_point) {
  if (!model %in% 1:4) {
    stop("ERROR: Try to specify the deterministic component again")
  }

  const <- matrix(data = 1, nrow = n_obs, ncol = 1)

  trend <- if (model != 1) {
    matrix(data = 1:n_obs, nrow = n_obs, ncol = 1)
  } else {
    NULL
  }

  du <- if (model != 3) .du(break_point, n_obs) else NULL

  dt <- if (model %in% c(3, 4)) .dt(break_point, n_obs) else NULL

  cbind(const, trend, du, dt)
}


#' @title
#' Construct determinant variables for [kpss_double]
#'
#' @details
#' Procedure to compute deterministic terms
#' for KPSS with 2 structural breaks.
#'
#' @param model A scalar equal to
#' * 1: for the AA (without trend) model,
#' * 2: for the AA (with trend) model,
#' * 3: for the BB model,
#' * 4: for the CC model,
#' * 5: for the AC-CA model,
#' * 6: for the AC-CA model,
#' * 7: for the AC-CA model.
#' @param n.obs Number of observations.
#' @param break.point Positions for the first and second structural breaks
#'            (respective to the origin which is 1).
#'
#' @return Matrix of deterministic terms.
#'
#' @keywords internal
trend_kpss_double <- function(model,
                              n_obs,
                              break_point) {
  if (any(!model %in% 1:7)) {
    stop("ERROR: Try to specify the deterministic component again")
  }

  const <- matrix(data = 1, nrow = n_obs, ncol = 1)

  trend <- if (model %in% c(2, 3, 4, 5)) {
    matrix(data = 1:n_obs, nrow = n_obs, ncol = 1)
  } else {
    NULL
  }

  du1 <- if (model %in% c(1, 2, 4, 5, 6, 7)) {
    .du(break_point[1], n_obs)
  } else {
    NULL
  }

  du2 <- if (model %in% c(1, 2, 4, 6)) {
    .du(break_point[2], n_obs)
  } else {
    NULL
  }

  dt1 <- if (model %in% c(3, 4, 6, 7)) {
    .dt(break_point[1], n_obs)
  } else {
    NULL
  }

  dt2 <- if (model %in% c(3, 4, 5, 7)) {
    .dt(break_point[2], n_obs)
  } else {
    NULL
  }

  cbind(const, trend, du1, dt1, du2, dt2)
}


#' @title
#' Deterministic terms for [kpss_multiple]
#'
#' @description
#' Procedure to compute deterministic terms for KPSS with \eqn{m}
#' structural breaks.
#'
#' @details
#' **model** should be either a scalar or a vector of the same size as the
#' **break.point**. If scalar **model** will be repeated till the length of
#' **break.point** is achieved.
#'
#' @param model A scalar or vector of
#' * 1: for the break in const,
#' * 2: for the break in trend,
#' * 3: for the break in const and trend.
#' @param n.obs Number of observations.
#' @param break.point Array of structural breaks.
#' @param const,trend Include constant and trend if `TRUE`.
#'
#' @return Matrix of deterministic terms.
#'
#' @keywords internal
trend_kpss_miltiple <- function(model,
                                n_obs,
                                break_point,
                                const = FALSE,
                                trend = FALSE) {
  n_breaks <- length(break_point)

  if (length(model) == 1) {
    model <- rep(model, n_breaks)
  } else if (length(model) != n_breaks) {
    stop("ERROR! Inconsistent sizes of model and break.point")
  }

  if (any(!model %in% 1:3)) {
    stop("ERROR: Try to specify the deterministic component again")
  }

  xt <- NULL
  if (const) {
    xt <- matrix(data = 1, nrow = n_obs, ncol = 1)
  }
  if (trend) {
    xt <- cbind(
      xt,
      matrix(data = 1:n_obs, nrow = n_obs, ncol = 1)
    )
  }

  for (i in 1:n_breaks) {
    xt <- switch(model[i],
      cbind(xt, .du(break_point[i], n_obs)),
      cbind(xt, .dt(break_point[i], n_obs)),
      cbind(xt, .du(break_point[i], n_obs), .dt(break_point[i], n_obs))
    )
  }

  xt
}


#' @title
#' Preparing variables for DOLS regression with multiple known break points
#'
#' @param y A time series of interest.
#' @param x A matrix of explanatory stochastic regressors.
#' @param model A scalar or vector of
#' * 1: for the break in const.
#' * 2: for the break in trend.
#' * 3: for the break in const and trend.
#' @param break.point An array of moments of structural breaks.
#' @param const,trend Whether a constant or trend are to be included.
#' @param k.lags,k.leads A number of lags and leads in DOLS regression.
#'
#' @return A list of LHS and RHS variables.
#'
#' @keywords internal
variables_dols_multiple <- function(
  y,
  x,
  model,
  break_point,
  const = FALSE,
  trend = FALSE,
  n_lags,
  n_leads
) {
  if (is.null(x)) {
    stop("ERROR! Explanatory variables needed for DOLS")
  }
  if (!is.matrix(y)) y <- as.matrix(y)
  if (!is.matrix(x)) x <- as.matrix(x)

  n_obs <- nrow(y)

  .dx_step <- x[2:n_obs, , drop = FALSE] - x[1:(n_obs - 1), , drop = FALSE]
  .dx_lags <- .dx_step
  .dx_leads <- .dx_step

  for (i in 1:n_lags) {
    .dx_lags <- cbind(
      .dx_lags,
      lagn(.dx_step, i)
    )
  }

  for (i in 1:n_leads) {
    .dx_leads <- cbind(
      .dx_leads,
      lagn(.dx_step, -i)
    )
  }

  if (n_lags != 0 && n_leads != 0) {
    lags <- .dx_lags
    leads <- .dx_leads[, (ncol(x) + 1):(ncol(.dx_leads)), drop = FALSE]
    .lags_leads <- cbind(lags, leads)
    .lags_leads <-
      .lags_leads[(n_lags + 1):(n_obs - 1 - n_leads), , drop = FALSE]
  } else if (n_lags != 0 && n_leads == 0) {
    lags <- .dx_lags
    .lags_leads <- lags[(n_lags + 1):(n_obs - 1), , drop = FALSE]
  } else if (n_lags == 0 && n_leads != 0) {
    lags <- .dx_lags
    leads <- .dx_leads[, (ncol(x) + 1):(ncol(.dx_leads)), drop = FALSE]
    .lags_leads <- cbind(lags, leads)
    .lags_leads <- .lags_leads[1:(n_obs - 1 - n_leads), , drop = FALSE]
  } else if (n_lags == 0 && n_leads == 0) {
    .lags_leads <- .dx_lags
  }
  deter <- trend_kpss_miltiple(model, n_obs, break_point, const, trend)

  list(
    yreg = y[(n_lags + 2):(n_obs - n_leads), 1, drop = FALSE],
    xreg = cbind(
      deter[(n_lags + 2):(n_obs - n_leads), , drop = FALSE],
      x[(n_lags + 2):(n_obs - n_leads), , drop = FALSE],
      .lags_leads
    )
  )
}


#' @title
#' Generating monthly seasonal dummy variables
#'
#' @param n.obs number of observations.
#'
#' @return The matrix of values od seasonal dummies.
#'
#' @keywords internal
seasonal_dummies <- function(n_obs) {
  s1 <- c(1 - 1 / 12, rep(-1 / 12, 11))

  result <- NULL
  for (i in 0:10) {
    result <- cbind(
      result,
      c(
        rep(-1 / 12, i),
        rep(s1, length.out = n_obs - i)
      )
    )
  }

  result
}


#' @title
#' Generating break in constant
#'
#' @param bp index of break point.
#' @param n_obs number of observations.
#'
#' @return The matrix of values od seasonal dummies.
#'
#' @keywords internal
.du <- function(bp, n_obs) {
  rbind(
    matrix(data = 0, nrow = bp, ncol = 1),
    matrix(data = 1, nrow = n_obs - bp, ncol = 1)
  )
}


#' @title
#' Generating break in trend
#'
#' @param bp index of break point.
#' @param n_obs number of observations.
#'
#' @return The matrix of values od seasonal dummies.
#'
#' @keywords internal
.dt <- function(bp, n_obs) {
  rbind(
    matrix(data = 0, nrow = bp, ncol = 1),
    matrix(data = 1:(n_obs - bp), nrow = n_obs - bp, ncol = 1)
  )
}

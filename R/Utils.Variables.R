#' @title
#' Construct determinant variables for [kpss.1br]
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
trend.kpss.single <- function(model, N, bp) {
  if (!model %in% 1:4) {
    stop("ERROR! kpss.single: Try to specify the deterministic component again")
  }

  cbind(
    1,
    1:N,
    if (model != 3) .du(bp, N) else NULL,
    if (model %in% c(3, 4)) .dt(bp, N) else NULL
  )
}


#' @title
#' Construct determinant variables for [kpss.2br]
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
trend.kpss.double <- function(model, N, bp) {
  if (any(!model %in% 1:7)) {
    stop("ERROR! kpss.double: Try to specify the deterministic component again")
  }

  cbind(
    1,
    1:N,
    if (model %in% c(1, 2, 4, 5, 6, 7)) .du(bp[1], N) else NULL,
    if (model %in% c(3, 4, 6, 7)) .dt(bp[1], N) else NULL,
    if (model %in% c(1, 2, 4, 6)) .du(bp[2], N) else NULL,
    if (model %in% c(3, 4, 5, 7)) .dt(bp[2], N) else NULL
  )
}


#' @title
#' Deterministic terms for [kpss.mlt]
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
trend.kpss.miltiple <- function(
  model,
  N,
  bps,
  const = FALSE,
  trend = FALSE
) {
  nb <- length(bps)

  if (length(model) == 1) {
    model <- rep(model, nb)
  } else if (length(model) != nb) {
    stop("ERROR! kpss.multiple: Inconsistent sizes of model and break.point")
  }

  if (any(!model %in% 1:3)) {
    stop("ERROR: kpss.multiple: Try to specify the deterministic component again")
  }

  xt <- cbind(
    if (const) 1 else NULL,
    if (trend) 1:N else NULL
  )

  for (i in 1:nb) {
    xt <- switch(model[i],
      cbind(xt, .du(bps[i], N)),
      cbind(xt, .dt(bps[i], N)),
      cbind(xt, .du(bps[i], N), .dt(bps[i], N))
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
variables.dols.multiple <- function(
  y,
  x,
  model,
  bp,
  const = FALSE,
  trend = FALSE,
  n.lags,
  n.leads
) {
  if (is.null(x)) {
    stop("ERROR! dols.multiple: Explanatory variables needed for DOLS")
  }
  if (!is.matrix(y)) y <- as.matrix(y)
  if (!is.matrix(x)) x <- as.matrix(x)

  N <- nrow(y)

  .dx_step <- x[2:N, , drop = FALSE] - x[1:(N - 1), , drop = FALSE]
  .dx_lags <- .dx_step
  .dx_leads <- .dx_step

  for (i in 1:n.lags) {
    .dx_lags <- cbind(.dx_lags, .lagn(.dx_step, i))
  }

  for (i in 1:n.leads) {
    .dx_leads <- cbind(
      .dx_leads,
      .lagn(.dx_step, -i)
    )
  }

  if (n.lags != 0 && n.leads != 0) {
    lags <- .dx_lags
    leads <- .dx_leads[, (ncol(x) + 1):(ncol(.dx_leads)), drop = FALSE]
    .lags_leads <- cbind(lags, leads)
    .lags_leads <-
      .lags_leads[(n.lags + 1):(N - 1 - n.leads), , drop = FALSE]
  } else if (n.lags != 0 && n.leads == 0) {
    lags <- .dx_lags
    .lags_leads <- lags[(n.lags + 1):(N - 1), , drop = FALSE]
  } else if (n.lags == 0 && n.leads != 0) {
    lags <- .dx_lags
    leads <- .dx_leads[, (ncol(x) + 1):(ncol(.dx_leads)), drop = FALSE]
    .lags_leads <- cbind(lags, leads)
    .lags_leads <- .lags_leads[1:(N - 1 - n.leads), , drop = FALSE]
  } else if (n.lags == 0 && n.leads == 0) {
    .lags_leads <- .dx_lags
  }
  deter <- trend.kpss.miltiple(model, N, bp, const, trend)

  list(
    yreg = y[(n.lags + 2):(N - n.leads), 1, drop = FALSE],
    xreg = cbind(
      deter[(n.lags + 2):(N - n.leads), , drop = FALSE],
      x[(n.lags + 2):(N - n.leads), , drop = FALSE],
      .lags_leads
    )
  )
}


#' @title
#' Generating monthly seasonal dummy variables
#'
#' @param N number of observations.
#'
#' @return The matrix of values od seasonal dummies.
#'
#' @keywords internal
seasonal.dummies <- function(N) {
  s1 <- c(1 - 1 / 12, rep(-1 / 12, 11))

  result <- NULL
  for (i in 0:10) {
    result <- cbind(
      result,
      c(
        rep(-1 / 12, i),
        rep(s1, length.out = N - i)
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
  matrix(
    as.numeric((1:n_obs) > bp),
    nrow = n_obs,
    ncol = 1
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
  matrix(
    as.numeric((1:n_obs) > bp) * ((1:n_obs) - bp),
    nrow = n_obs,
    ncol = 1
  )
}

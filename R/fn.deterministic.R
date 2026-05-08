#' @title
#' Deterministic terms for [coint.CSS]
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
trend.variables <- function(
  break.type,
  N,
  break.point,
  const = FALSE,
  trend = FALSE
) {
  Nb <- length(break.point)

  if (length(break.type) == 1) {
    break.type <- rep(break.type, Nb)
  } else if (length(break.type) != Nb) {
    stop("ERROR! kpss.multiple: Inconsistent sizes of model and break.point")
  }

  xt <- cbind(
    if (const) .const(N) else NULL,
    if (trend) .trend(N) else NULL
  )

  for (i in seq_len(Nb)) {
    xt <- switch(
      break.type[i],
      "c" = cbind(xt, .du(break.point[i], N)),
      "t" = cbind(xt, .dt(break.point[i], N)),
      "ct" = cbind(xt, .du(break.point[i], N), .dt(break.point[i], N)),
      stop("ERROR: kpss.multiple: unknown break value '", break.type[i], "'")
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
#' @importFrom Rfast rowAll
#'
#' @keywords internal
DOLS.mlt.regressors <- function(
  y,
  x,
  const = FALSE,
  trend = FALSE,
  break.type,
  break.point,
  break.coint = FALSE,
  n.lags,
  n.leads,
  max.ll = NULL
) {
  if (is.null(x)) {
    stop("ERROR! dols.multiple: Explanatory variables needed for DOLS")
  }
  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }

  if (!is.null(max.ll)) {
    max.lag <- max.ll[1]
    max.lead <- max.ll[2]
  } else {
    max.lag <- n.lags
    max.lead <- n.leads
  }

  N <- nrow(y)

  mDeter <- trend.variables(break.type, N, break.point, const, trend)

  mXdu <- NULL
  if (break.coint) {
    for (bp in break.point) {
      mXdu <- cbind(mXdu, sweep(x, 1, .du(bp, N), `*`))
    }
  }

  mDx <- .diffn(x)
  mLagLead <- mDx

  if (n.lags > 0) {
    for (i in 1:n.lags) {
      mLagLead <- cbind(mLagLead, .lagn(mDx, i))
    }
  }
  if (n.leads > 0) {
    for (i in 1:n.leads) {
      mLagLead <- cbind(mLagLead, .lagn(mDx, -i))
    }
  }

  rows <- rowAll(!is.na(y)) & rowAll(!is.na(x))
  rows[1:max.lag] <- FALSE
  rows[(N - max.lead + 1):N] <- FALSE

  mX <- cbind(
    mDeter,
    mXdu,
    x,
    mLagLead
  )

  list(
    yreg = y[rows, 1, drop = FALSE],
    xreg = mX[rows, , drop = FALSE]
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
#' Generating const and trend
#' @keywords internal
.const <- function(N) rep(1, N)
.trend <- function(N) 1:N


#' @title
#' Generating break in constant
#'
#' @param bp index of break point.
#' @param n_obs number of observations.
#'
#' @return The matrix of values od seasonal dummies.
#'
#' @keywords internal
.du <- function(bp, N) {
  matrix(
    as.numeric((1:N) > bp),
    nrow = N,
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
.dt <- function(bp, N) {
  matrix(
    as.numeric((1:N) > bp) * ((1:N) - bp),
    nrow = N,
    ncol = 1
  )
}

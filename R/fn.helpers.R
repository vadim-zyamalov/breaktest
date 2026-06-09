#' @title
#' Produce a vector lagged backward of forward
#'
#' @param x An initial vector.
#' @param i Size of lag (lead if negative).
#' @param na A value to fill missing observations, `NA` by default.
#'
#' @return Lagged or leaded vector.
#'
#' @keywords internal
.lagn <- function(x, i, na = NA) {
  .lag <- function(x, i, na) {
    N <- length(x)

    if (i > 0) {
      c(rep(na, i), x)[1:N]
    } else {
      i <- abs(i)
      c(x, rep(na, i))[(i + 1):(N + i)]
    }
  }

  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  apply(x, 2, (function(col) .lag(col, i, na)))
}


#' @title
#' Produce a vector or matrix of differences, keeping initial length
#'
#' @param x An initial vector.
#' @param lag Size of lag.
#' @param difference Order of differentiating.
#' @param na A value to fill missing observations, `NA` by default.
#'
#' @return Vector or matrix of differences.
#'
#' @keywords internal
.diffn <- function(x, lag = 1, differences = 1, na = NA) {
  .diff <- function(x, l, d, na) {
    N <- length(x)
    tmp <- diff(x, lag = l, differences = d)
    c(rep(na, N - length(tmp)), tmp)
  }

  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  apply(x, 2, (function(col) .diff(col, lag, differences, na)))
}


#' @title
#' Auxiliary function returning KPSS statistic value.
#'
#' @param resids A series of residuals.
#' @param variance A value of the long-run variance.
#'
#' @keywords internal
.kpss.statistic <- function(resids, variance) {
  N <- length(resids)
  sum(cumsum(resids)^2) / (N^2 * variance)
}


#' @title
#' Calculating M-statistics by Stock (1990) and Perron and Ng (1996).
#'
#' @param y A time series of interest.
#' @param l Number of lags for inner ADF test.
#' @param const,trend Whether a constant and trend are to be included.
#'
#' @return List of values of \eqn{MZ_\alpha}, \eqn{MZ_t} and \eqn{MSB}
#' statistics.
#'
#' @references
#' Perron, Pierre, and Serena Ng.
#' “Useful Modifications to Some Unit Root Tests with Dependent Errors
#' and Their Local Asymptotic Properties.”
#' The Review of Economic Studies 63, no. 3 (July 1, 1996): 435–63.
#' https://doi.org/10.2307/2297890.
#'
#' Stock, James H.
#' “A Class of Tests for Integration and Cointegration.”
#' Kennedy School of Government, Harvard University, 1990.
#'
#' @keywords internal
.mz.statistics <- function(y, l, const = FALSE, trend = FALSE) {
  N <- nrow(y)
  .adf <- uroot.ADF(y, const, trend, l, criterion = NULL)
  r_ADF <- .adf$model$residuals
  b_ADF <- .adf$model$coefficients
  a_ADF <- .adf$alpha

  denom <- if (l > 0) 1 - sum(b_ADF) + a_ADF else 1
  s_2 <- sum(r_ADF^2) / (length(r_ADF) - (1 + l)) / denom^2
  sum_y2 <- sum(y[1:(N - 1)]^2)

  mza <- (y[N]^2 / N - s_2) / (2 * sum_y2 / N^2)
  msb <- sqrt(sum_y2 / s_2 / N^2)
  mzt <- mza * msb

  list(
    mza = mza,
    msb = msb,
    mzt = mzt
  )
}


.nbcn.t.stat <- function(y, t1, t2) {
  lo <- min(t1, t2)
  hi <- max(t1, t2)

  dy <- y[lo:hi] - y[(lo - 1):(hi - 1)]
  yL <- y[(lo - 1):(hi - 1)]

  S_y1dy <- sum(yL * dy)
  S_yL2 <- sum(yL^2)
  S_dy2 <- sum(dy^2)

  list(
    S_y1dy = S_y1dy,
    S_yL2 = S_yL2,
    S_dy2 = S_dy2,
    t.stat = S_y1dy / sqrt(S_yL2)
  )
}

.nbcn.cval <- function(coef, lambda) {
  if (is.null(lambda)) {
    stop(".nbcn.cval: No lambda provided!")
  }

  d <- as.integer(lambda > 0.7)
  nc <- length(coef)

  x <- c(
    1,
    1 / lambda,
    lambda,
    lambda^2,
    lambda^3,
    d,
    d / lambda,
    d * lambda,
    d * lambda^2,
    d * lambda^3
  )[1:nc]

  sum(coef * x)
}

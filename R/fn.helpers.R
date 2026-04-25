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

  if (!is.matrix(x)) x <- as.matrix(x)
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
.diffn <- function(x,
                   lag = 1,
                   differences = 1,
                   na = NA) {
  .diff <- function(x, l, d, na) {
    N <- length(x)
    tmp <- diff(x, lag = l, differences = d)
    c(rep(na, N - length(tmp)), tmp)
  }

  if (!is.matrix(x)) x <- as.matrix(x)
  apply(x, 2, (function(col) .diff(col, lag, differences, na)))
}


#' @title
#' Auxiliary function returning KPSS statistic value.
#'
#' @param resids A series of residuals.
#' @param variance A value of the long-run variance.
#'
#' @keywords internal
.kpss.statistic <- function(resids,
                            variance) {
  if (!is.matrix(resids)) resids <- as.matrix(resids)
  N <- nrow(resids)
  s_t <- apply(resids, 2, cumsum)
  drop(t(s_t) %*% s_t) / (N^2 * variance)
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
.mz.statistics <- function(y,
                           l,
                           const = FALSE,
                           trend = FALSE) {
  n_obs <- nrow(y)
  .adf <- ADF.test(y, const, trend, l, criterion = NULL)

  denom <- 1 - sum(.adf$coefficients) + .adf$alpha
  s_2 <- drop(t(.adf$residuals) %*% .adf$residuals) /
    (nrow(.adf$residuals) - (1 + l)) / denom^2
  sum_y2 <- sum(.adf$yd[1:(n_obs - 1)]^2)

  mza <- (y[n_obs]^2 / n_obs - s_2) / (2 * sum_y2 / n_obs^2)
  msb <- sqrt(sum_y2 / s_2 / n_obs^2)
  mzt <- mza * msb

  list(
    mza = mza,
    msb = msb,
    mzt = mzt
  )
}

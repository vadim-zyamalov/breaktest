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
.lagn <- function(x,
                  i,
                  na = NA) {
  if (!is.matrix(x)) x <- as.matrix(x)
  n_obs <- nrow(x)
  n_var <- ncol(x)

  if (i > 0) {
    rbind(
      matrix(data = na, nrow = i, ncol = n_var),
      x[1:(n_obs - i), , drop = FALSE]
    )
  } else {
    rbind(
      x[(1 + abs(i)):n_obs, , drop = FALSE],
      matrix(data = na, nrow = abs(i), ncol = n_var)
    )
  }
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
  if (!is.matrix(x)) x <- as.matrix(x)
  n_obs <- nrow(x)
  n_var <- ncol(x)
  .diff <- diff(x, lag = lag, differences = differences)

  rbind(
    matrix(
      data = na,
      nrow = n_obs - nrow(.diff),
      ncol = n_var
    ),
    .diff
  )
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
  n_obs <- nrow(resids)
  s_t <- apply(resids, 2, cumsum)

  drop(t(s_t) %*% s_t) / (n_obs^2 * variance)
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
.mz_stats <- function(y,
                      l,
                      const = FALSE,
                      trend = FALSE) {
  n_obs <- nrow(y)
  .adf <- ADF.test(y, const, trend, l, criterion = NULL)

  denom <- 1 - sum(.adf$beta) + .adf$alpha
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

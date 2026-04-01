#' @title
#' A simple implementation of ADF test
#'
#' @description
#' A function for ADF test with the ability to select the number of lags.
#' Lags are selected by informational criterions which can be modified as in
#' Ng and Perron (2001) and Cavaliere et al. (2015).
#'
#' @details
#' Due to the Frisch-Waugh-Lovell theorem we first detrend `y` and then apply
#' the test to the detrended series.
#'
#' @param y A time series of interest.
#' @param const,trend Whether a constand and trend are to be included.
#' @param max.lag Maximum lag number.
#' @param criterion A criterion used to select number of lags.
#' If lag selection is not needed keep this NULL.
#' @param modified.criterion Whether the unit-root test modificaton is needed.
#' @param rescale.criterion Whether the rescaling informational criterion
#' is needed. Designed to cope with heteroscedasticity in residuals.
#'
#' @return A list containing:
#' * y,
#' * const,
#' * trend,
#' * residuals,
#' * coefficient estimates,
#' * t-statistic value,
#' * critical value,
#' * Number of lags,
#' * indicator of stationarity.
#'
#' @references
#' Cavaliere, Giuseppe, Peter C. B. Phillips, Stephan Smeekes,
#' and A. M. Robert Taylor. “Lag Length Selection for Unit Root Tests
#' in the Presence of Nonstationary Volatility.”
#' Econometric Reviews 34, no. 4 (April 21, 2015): 512–36.
#' https://doi.org/10.1080/07474938.2013.808065.
#'
#' Ng, Serena, and Pierre Perron. “Lag Length Selection and the Construction of
#' Unit Root Tests with Good Size and Power.”
#' Econometrica 69, no. 6 (2001): 1519–54.
#' https://doi.org/10.1111/1468-0262.00256.
#'
#' @export
ADF.test <- function(y,
                     const = TRUE,
                     trend = FALSE,
                     max.lag = 0,
                     criterion = NULL,
                     modified.criterion = FALSE,
                     rescale.criterion = FALSE) {
  if (!is.null(criterion)) {
    if (!criterion %in% c("bic", "aic", "lwz", "hq")) {
      stop("ERROR! Unknown criterion, none is used")
    }
  }

  if (!is.matrix(y)) y <- as.matrix(y)

  N <- nrow(y)
  rows <- (2 + max.lag):N

  deter <- cbind(
    if (const) 1 else NULL,
    if (trend) 1:N else NULL
  )

  ## Detrending
  if (!is.null(deter)) {
    yd <- .OLS(y, deter)$residuals
  } else {
    yd <- y
  }

  d.y <- .diffn(yd)
  d.y[1] <- yd[1]

  x <- .lagn(yd, 1, na = 0)
  if (max.lag > 0) {
    x <- cbind(
      x,
      apply(as.array(1:max.lag), 1, function(l) .lagn(d.y, l, na = 0))
    )
  }

  if (is.null(criterion)) {
    res.lag <- max.lag
  } else {
    if (rescale.criterion) {
      tmp.rescale <- rescale.CPST(d.y, x, deter, 0, max.lag)
      d.yr <- tmp.rescale$d.y
      xr <- tmp.rescale$x
      rm(tmp.rescale)
    } else {
      d.yr <- d.y
      xr <- x
    }

    tmp.ols <- .OLS(
      d.yr[rows, , drop = FALSE],
      xr[rows, 1, drop = FALSE]
    )
    b <- tmp.ols$beta
    e <- tmp.ols$residuals
    rm(tmp.ols)

    res.ic <- .ic.values(
      e, 0,
      modification = modified.criterion,
      alpha = b[1], y = xr[rows, 1, drop = FALSE]
    )[[criterion]]
    res.lag <- 0

    for (l in 1:max.lag) {
      if (max.lag == 0) break

      if (rescale.criterion) {
        tmp.rescale <- rescale.CPST(d.y, x, deter, l, max.lag)
        d.yr <- tmp.rescale$d.y
        xr <- tmp.rescale$x
        rm(tmp.rescale)
      } else {
        d.yr <- d.y
        xr <- x
      }

      tmp.ols <- .OLS(
        d.yr[rows, , drop = FALSE],
        xr[rows, 1:(1 + l), drop = FALSE]
      )
      b <- tmp.ols$beta
      e <- tmp.ols$residuals
      rm(tmp.ols)

      tmp.ic <- .ic.values(
        e, l,
        modification = modified.criterion,
        alpha = b[1], y = xr[rows, 1, drop = FALSE]
      )[[criterion]]

      if (tmp.ic < res.ic) {
        res.ic <- tmp.ic
        res.lag <- l
      }
    }
  }

  res.OLS <- .OLS(
    d.y[(2 + res.lag):N, , drop = FALSE],
    x[(2 + res.lag):N, 1:(1 + res.lag), drop = FALSE]
  )

  Z.stat <- (N - res.lag - 1) * drop(res.OLS$beta[1] - 1)

  list(
    y = drop(y),
    yd = drop(yd),
    const = const,
    trend = trend,
    beta = res.OLS$beta,
    t.beta = drop(res.OLS$t.beta),
    alpha = drop(res.OLS$beta[1]),
    t.alpha = drop(res.OLS$t.beta[1]),
    Z.stat = Z.stat,
    residuals = res.OLS$residuals,
    lag = res.lag
  )
}


#' @title
#' Generating rescaled series as in Cavaliere et al. (2015)
#'
#' @description
#' This rescaling procedure is needed to cope with possible heteroscedasticity
#' in the data. Simply it's achieved by taking a cumulative sum of the
#' first difference normalized by the non-parametric local estimate of the
#' variance.
#'
#' @param d.y A series of first differences.
#' @param x A matrix of ADF variables.
#' @param deter A matrix of deterministic variables for detrending.
#' @param adf.lag A lag of the corresponding ADF model.
#' @param max.lag The maximum possible lag.
#'
#' @return A rescaled series.
#'
#' @references
#' Cavaliere, Giuseppe, Peter C. B. Phillips, Stephan Smeekes,
#' and A. M. Robert Taylor. “Lag Length Selection for Unit Root Tests
#' in the Presence of Nonstationary Volatility.”
#' Econometric Reviews 34, no. 4 (April 21, 2015): 512–36.
#' https://doi.org/10.1080/07474938.2013.808065.
#'
#' @keywords internal
rescale.CPST <- function(d.y,
                         x,
                         deter,
                         adf.lag,
                         max.lag) {
  e <- .OLS(d.y, x[, 1:(1 + adf.lag), drop = FALSE])$residuals

  NW.se <- .NW.variance(
    e,
    .NW.bandwidth(e^2, rep(1, nrow(e)))$h
  )$se

  yr <- cumsum(d.y / NW.se)

  if (!is.null(deter)) {
    yr <- .OLS(yr, deter)$residuals
  }

  d.yr <- as.matrix(c(0, diff(yr)))

  xr <- .lagn(yr, 1, na = 0)
  if (max.lag > 0) {
    for (l in 1:max.lag) {
      xr <- cbind(xr, .lagn(d.yr, l, na = 0))
    }
  }

  list(
    d.y = d.yr,
    x = xr
  )
}

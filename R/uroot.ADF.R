#' @title
#' A simple implementation of ADF test
#'
#' @description
#' A function for ADF test with the ability to select the number of lags.
#' Lags are selected by informational criterions which can be modified as in
#' Ng and Perron (2001) and Cavaliere et al. (2015).
#'
#' @details
#' A function for ADF test with the ability to select the number of lags.
#' Lags are selected by informational criterions which can be modified as in
#' Ng and Perron (2001) and Cavaliere et al. (2015).
#'
#' It's an ordinary Augmented Dickey-Fuller test using OLS estimation under the hood.
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
#' @param recursive Whether a recursive detrending should be applied.
#' See Smeekes (2013) for details.
#' @param cc A filtration parameter used to construct an autocorrelation
#' coefficient.
#' @param gamma Detrending type selection parameter. If 0 the OLS detrending
#' is applied, if 1 the GLS detrending is applied, otherwise the autocorrelation
#' coefficient is calculated as \eqn{1 + c^{\gamma} T^{-\gamma}}.
#' @param trim A trimming parameter.
#' @param boot.p Whether bootstrapped p-values should be returned.
#' @param boot.iter The number of bootstrap iterations.
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
#' Taylor, A. M. Robert.
#' “Regression-Based Unit Root Tests With Recursive Mean Adjustment for
#' Seasonal and Nonseasonal Time Series.”
#' Journal of Business & Economic Statistics 20, no. 2 (April 2002): 269–81.
#' https://doi.org/10.1198/073500102317352001.
#'
#' MacKinnon, James G.
#' “Critical Values for Cointegration Tests.”
#' Working Paper. Economics Department, Queen’s University, January 2010.
#' https://ideas.repec.org/p/qed/wpaper/1227.html.
#'
#' Smeekes, Stephan.
#' “Detrending Bootstrap Unit Root Tests.”
#' Econometric Reviews 32, no. 8 (July 2013): 869–91.
#' https://doi.org/10.1080/07474938.2012.690693.
#'
#' Elliott, Graham, Thomas J. Rothenberg, and James H. Stock.
#' “Efficient Tests for an Autoregressive Unit Root.”
#' Econometrica 64, no. 4 (1996): 813–36.
#' https://doi.org/10.2307/2171846.
#'
#' @import doSNOW
#' @import foreach
#' @import parallel
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#'
#' @export
uroot.ADF <- function(
  y,
  const = TRUE,
  trend = FALSE,
  max.lag = 0,
  criterion = NULL,
  modified.criterion = FALSE,
  rescale.criterion = FALSE,
  recursive = FALSE,
  cc = 0,
  gamma = 0,
  trim = 0.15,
  boot.p = FALSE,
  boot.iter = 999
) {
  if (!is.null(criterion)) {
    if (!criterion %in% c("bic", "aic", "lwz", "hq")) {
      stop("ERROR! Unknown criterion, none is used")
    }
  }

  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }

  N <- nrow(y)
  rows <- (2 + max.lag):N

  mDeter <- cbind(
    if (const) .const(N) else NULL,
    if (trend) .trend(N) else NULL
  )

  ## Detrending
  if (!is.null(mDeter)) {
    y <- if (recursive) {
      detrend.recursively(y, mDeter, cc, gamma, trim)
    } else {
      OLS.reg(y, mDeter)$residuals
    }
    y <- as.matrix(y)
  }

  diffY <- .diffn(y)
  mX <- .lagn(y, 1)
  for (l in seq_len(max.lag)) {
    mX <- cbind(mX, .lagn(diffY, l))
  }

  if (is.null(criterion)) {
    rLag <- max.lag
  } else {
    if (rescale.criterion) {
      tmp.rescale <- rescale.CPST(diffY, mX, mDeter, 0, max.lag)
      diffYr <- tmp.rescale$d.y
      mXr <- tmp.rescale$x
    } else {
      diffYr <- diffY
      mXr <- mX
    }

    tmp.ols <- OLS.reg(
      diffYr[rows, , drop = FALSE],
      mXr[rows, 1, drop = FALSE]
    )
    b <- tmp.ols$coefficients

    rIC <- info.criterions(
      tmp.ols,
      criterion,
      alpha = b[1],
      y = mXr[rows, 1, drop = FALSE]
    )
    rLag <- 0

    for (l in seq_len(max.lag)) {
      if (rescale.criterion) {
        tmp.rescale <- rescale.CPST(diffY, mX, mDeter, l, max.lag)
        diffYr <- tmp.rescale$d.y
        mXr <- tmp.rescale$x
      } else {
        diffYr <- diffY
        mXr <- mX
      }

      tmp.ols <- OLS.reg(
        diffYr[rows, , drop = FALSE],
        mXr[rows, 1:(1 + l), drop = FALSE]
      )
      b <- tmp.ols$coefficients

      tmp.ic <- info.criterions(
        tmp.ols,
        criterion,
        alpha = b[1],
        y = mXr[rows, 1, drop = FALSE]
      )

      if (tmp.ic < rIC) {
        rIC <- tmp.ic
        rLag <- l
      }
    }
  }

  res.OLS <- OLS.reg(
    diffY[rows, , drop = FALSE],
    mX[rows, 1:(1 + rLag), drop = FALSE]
  )

  dZstat <- (N - rLag - 1) * drop(res.OLS$coefficients[1] - 1)

  result <- list(
    const = const,
    trend = trend,
    model = res.OLS,
    alpha = as.numeric(res.OLS$coefficients[1]),
    t.alpha = as.numeric(res.OLS$t.stats[1]),
    Z.stat = dZstat,
    lag = rLag,
    recursive = recursive
  )
  class(result) <- "bt_adf"

  if (recursive) {
    result$recursive.params <- list(cc = cc, gamma = gamma, trim = trim)
  }

  if (boot.p) {
    result$p.value <- bootstrap(result, boot.iter)
  }

  result
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
rescale.CPST <- function(d.y, x, deter, adf.lag, max.lag) {
  e <- OLS.reg(d.y, x[, 1:(1 + adf.lag), drop = FALSE])$residuals

  dNWse <- NW.variance(
    e,
    NW.bandwidth(e^2, rep(1, nrow(e)))$h
  )$se

  vYresc <- cumsum(d.y[-1] / dNWse)

  if (!is.null(deter)) {
    vYresc <- OLS.reg(vYresc, deter)$residuals
  }

  diffYresc <- .diffn(vYresc, na = 0)

  xr <- .lagn(vYresc, 1, na = 0)
  for (l in seq_len(max.lag)) {
    xr <- cbind(xr, .lagn(diffYresc, l, na = 0))
  }

  list(
    d.y = diffYresc,
    x = xr
  )
}


#' @title
#' Detrending the data recursively
#'
#' @description
#' This procedure is aimed to provide a recursively detrended series. More or
#' less classical approach of full-sample detrending may lead to the regressors
#' correlated with the error term.
#'
#' @details
#' Elliott et al (1996) recommend using \eqn{c = -7} for the model with only
#' an intercept, and \eqn{c = -13.5} for the model with a linear trend.
#'
#' @param y A time series of interest.
#' @param x A matrix of explanatory variables.
#' @param c A filtration parameter used to construct an autocorrelation
#' coefficient.
#' @param gamma A detrending type selection parameter. If 0 the OLS detrending
#' is applied, if 1 the GLS detrending is applied, otherwise the autocorrelation
#' coefficient is calculated as \eqn{1 + c^{\gamma} T^{-\gamma}}.
#' @param trim A trimming parameter. It's used to find the minimum size of
#' subsamples while calculating recursive estimates. The ending point of the
#' subsample for the \eqn{t} is \eqn{max(t, trim \times T)}.
#'
#' @return A detrended series.
#'
#' @references
#' Elliott, Graham, Thomas J. Rothenberg, and James H. Stock.
#' “Efficient Tests for an Autoregressive Unit Root.”
#' Econometrica 64, no. 4 (1996): 813–36.
#' https://doi.org/10.2307/2171846.
#'
#' Taylor, A. M. Robert.
#' “Regression-Based Unit Root Tests With Recursive Mean Adjustment for
#' Seasonal and Nonseasonal Time Series.”
#' Journal of Business & Economic Statistics 20, no. 2 (April 2002): 269–81.
#' https://doi.org/10.1198/073500102317352001.
#'
#' @keywords internal
detrend.recursively <- function(y, x, cc, gamma, trim) {
  if (is.null(x)) {
    return(y)
  }

  n.obs <- nrow(y)
  beg <- trunc(trim * n.obs)
  ct <- (cc / n.obs)^gamma

  yt <- y - (1 + ct) * .lagn(y, 1, na = 0)
  xt <- x - (1 + ct) * .lagn(x, 1, na = 0)

  yd <- OLS.reg(
    yt[1:beg, , drop = FALSE],
    xt[1:beg, , drop = FALSE]
  )$residuals

  for (lstar in (beg + 1):n.obs) {
    ystar <- OLS.reg(
      yt[1:lstar, , drop = FALSE],
      xt[1:lstar, , drop = FALSE]
    )$residuals
    yd <- c(yd, ystar[lstar])
  }
  as.matrix(yd)
}

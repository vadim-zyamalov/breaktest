#' @title
#' Detrending bootstrap test by Smeekes (2013)
#'
#' @description
#' This bootstrap test is based on the recursive detrending procedure of Taylor
#' (2002). The main idea is to apply the standard ADF test to the series
#' with nuissanse parameters eliminated.
#'
#' @details
#' Critical values are calculated via a bootstrapping using MacKinnon-like
#' regressions. For each number of observations and each number of variables
#' obtained were 1999 values of test statistics. After that 1st, 2.5-th, 5-th,
#' 10-th, and 97.5-th percentiles were calculated and saved along with the
#' corresponding number of observations. This step was repeated 5 times to cope
#' with possible biases. After that MacKinnon-like regressions were estimated.
#'
#' @param y A time series of interest.
#' @param const,trend Whether the constant and trend are to be included.
#' @param c A filtration parameter used to construct an autocorrelation
#' coefficient.
#' @param gamma Detrending type selection parameter. If 0 the OLS detrending
#' is applied, if 1 the GLS detrending is applied, otherwise the autocorrelation
#' coefficient is calculated as \eqn{1 + c^{\gamma} T^{-\gamma}}.
#' @param trim A trimming parameter.
#' @param max.lag The maximum lag for inner ADF testing.
#' @param criterion A criterion used to select number of lags.
#' If lag selection is not needed keep this NULL.
#' @param modified.criterion Whether the unit-root test modificaton is needed.
#' @param iter The number of bootstrap iterations.
#'
#' @references
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
ADF.test.S <- function(y,
                       const = TRUE,
                       trend = FALSE,
                       c = 0,
                       gamma = 0,
                       trim = 0.15,
                       max.lag = 0,
                       criterion = NULL,
                       modified.criterion = FALSE,
                       iter = 999) {
  if (!is.matrix(y)) y <- as.matrix(y)

  N <- nrow(y)

  x <- cbind(
    if (const) .const(N) else NULL,
    if (trend) .trend(N) else NULL,
  )

  yd <- detrend.recursively(y, x, c, gamma, trim)

  res.ADF <- ADF.test(
    yd, const, trend, max.lag,
    criterion, modified.criterion
  )

  res.lag <- res.ADF$lag
  res.stat <- res.ADF$t.alpha
  res.beta <- res.ADF$beta[-1]

  e <- res.ADF$residuals

  progress.bar <- txtProgressBar(max = iter, style = 3)
  progress <- function(n) setTxtProgressBar(progress.bar, n)

  cores <- detectCores()
  cluster <- makeCluster(max(cores - 1, 1), type = "SOCK")
  registerDoSNOW(cluster)

  tmp.stats <- foreach(
    i = 1:iter,
    .combine = c,
    .inorder = FALSE,
    .errorhandling = "remove",
    .packages = c("breaktest"),
    .options.snow = list(progress = progress)
  ) %dopar% {
    u <- rep(0, res.lag + N)
    eps <- sample(e, N, replace = TRUE)

    if (res.lag > 0) {
      for (s in 1:N) {
        u[res.lag + s] <- u[(res.lag + s - 1):s] %*% res.beta + eps[s]
      }
      u <- u[-(1:res.lag)]
    } else {
      for (s in 1:N) {
        u[s] <- eps[s]
      }
    }

    tmp.y <- as.matrix(rep(0, N))

    tmp.y[1] <- u[1]
    for (s in 2:N) {
      tmp.y[s] <- tmp.y[s - 1] + u[s]
    }

    tmp.yd <- detrend.recursively(tmp.y, x, c, gamma, trim)

    tmp.res <- ADF.test(
      tmp.yd, const, trend, res.lag,
      criterion = NULL
    )

    tmp.res$t.alpha
  }

  stopCluster(cluster)

  p.value <- sum(tmp.stats < res.stat) / iter

  list(
    stat = res.stat,
    p.value = p.value
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
detrend.recursively <- function(y,
                                x,
                                c,
                                gamma,
                                trim) {
  if (is.null(x)) {
    return(y)
  }

  n.obs <- nrow(y)
  beg <- trunc(trim * n.obs)
  ct <- (c / n.obs)^gamma

  yt <- y - (1 + ct) * .lagn(y, 1, na = 0)
  xt <- x - (1 + ct) * .lagn(x, 1, na = 0)

  yd <- .OLS(
    yt[1:beg, , drop = FALSE],
    xt[1:beg, , drop = FALSE]
  )$residuals

  yd <- rbind(
    yd,
    as.matrix(rep(0, n.obs - beg))
  )

  for (lstar in (beg + 1):n.obs) {
    ystar <- .OLS(
      yt[1:lstar, , drop = FALSE],
      xt[1:lstar, , drop = FALSE]
    )$residuals
    yd[lstar, ] <- ystar[lstar, ]
  }

  yd
}

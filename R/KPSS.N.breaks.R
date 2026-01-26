#' @title
#' KPSS-test with multiple known structural breaks
#'
#' @description
#' Procedure to compute the KPSS test with multiple known structural breaks
#'
#' @details
#' The code provided is based on the original code
#' by Carrion-i-Silvestre and Sansó.
#'
#' @param y A time series of interest.
#' @param x A matrix of explanatory stochastic regressors.
#' @param model A scalar or vector of
#' * 1: for the break in const,
#' * 2: for the break in trend,
#' * 3: for the break in const and trend.
#' @param break.point Array of structural breaks.
#' @param const,trend Whether a constant or trend should be included.
#' @param weakly.exog Boolean where we specify
#' whether the stochastic regressors are exogenous or not
#' * `TRUE`: if the regressors are weakly exogenous,
#' * `FALSE`: if the regressors are not weakly exogenous
#' (DOLS is used in this case).
#' @param lags.init,leads.init Scalars defininig the initial number of lags and
#' leads for DOLS.
#' @param max.lag scalar, with the maximum order of the parametric correction.
#' The final order of the parametric correction is selected
#' using the BIC information criterion.
#' @param kernel Kernel for calculating long-run variance
#' * `bartlett`: for Bartlett kernel,
#' * `quadratic`: for Quadratic Spectral kernel,
#' * `NULL` for the Kurozumi's proposal, using Bartlett kernel.
#' @param criterion Information criterion for DOLS lags and leads selection:
#' aic, bic, hq, or lwz.
#'
#' @return A list of
#' * `beta`: DOLS estimates of the coefficients,
#' * `tests`: SC test (coinKPSS-test),
#' * `resid`: Residuals of the model,
#' * `t.beta`: \eqn{t}-statistics for `beta`,
#' * `DOLS.lags`: The estimated number of lags and leads in DOLS,
#' * `break_point`: Break points.
#'
#' @references
#' Carrion-i-Silvestre, Josep Lluís, and Andreu Sansó.
#' “Testing the Null of Cointegration with Structural Breaks.”
#' Oxford Bulletin of Economics and Statistics 68, no. 5 (October 2006): 623–46.
#' https://doi.org/10.1111/j.1468-0084.2006.00180.x.
#'
#' Carrion-i-Silvestre, Josep Lluís, and Andreu Sansó.
#' “The KPSS Test with Two Structural Breaks.”
#' Spanish Economic Review 9, no. 2 (May 16, 2007): 105–27.
#' https://doi.org/10.1007/s10108-006-9017-8.
#'
#' @export
KPSS.N.breaks <- function(y,
                          x,
                          model,
                          break.point,
                          const = FALSE,
                          trend = FALSE,
                          weakly.exog = TRUE,
                          lags.init,
                          leads.init,
                          max.lag,
                          kernel,
                          criterion = "bic") {
  if (!is.matrix(y)) y <- as.matrix(y)
  if (!is.null(x)) {
    if (!is.matrix(x)) x <- as.matrix(x)
  }

  n.obs <- nrow(y)

  if (weakly.exog) {
    xt <- cbind(
      x,
      trend_kpss_miltiple(model, n.obs, break.point, const, trend)
    )

    tmp.OLS <- .estimate_ols(y, xt)
    beta <- tmp.OLS$beta
    resids <- tmp.OLS$residuals
    t.beta <- tmp.OLS$t.beta
    DOLS.lags <- 0
    DOLS.leads <- 0
    rm(tmp.OLS)
  } else {
    info.crit.min <- Inf
    for (i in lags.init:1) {
      for (j in leads.init:1) {
        tmp.DOLS <- DOLS.N.breaks(
          y, x, model, break.point, const, trend, i, j
        )
        beta <- tmp.DOLS$beta
        resids <- tmp.DOLS$residuals
        t.beta <- tmp.DOLS$t.beta
        info.crit <- tmp.DOLS$criterions
        if (info.crit[[criterion]] < info.crit.min) {
          info.crit.min <- info.crit[[criterion]]
          beta.min <- beta
          t.beta.min <- t.beta
          resid.min <- resids
          DOLS.lags <- i
          DOLS.leads <- j
        }
        rm(tmp.DOLS)
      }
    }
    resids <- resid.min
    beta <- beta.min
    t.beta <- t.beta.min
  }

  if (!is.null(kernel)) {
    test <- .kpss_stat(resids, lr.var.SPC(resids, max.lag, kernel))
  } else {
    test <- .kpss_stat(resids, lr.var.bartlett.AK(resids))
  }

  return(
    list(
      beta = beta,
      test = test,
      residuals = resids,
      t.beta = t.beta,
      DOLS.lags = DOLS.lags,
      DOLS.leads = DOLS.leads,
      break.point = break.point
    )
  )
}


#' @title
#' KPSS-test with multiple unknown structural breaks
#'
#' @description
#' Procedure to compute the KPSS test with multiple unknown structural breaks
#'
#' @param y A time series of interest.
#' @param x A matrix of explanatory stochastic regressors.
#' @param model A scalar or vector of
#' * 1: for the break in const,
#' * 2: for the break in trend,
#' * 3: for the break in const and trend.
#' @param break.point Array of structural breaks.
#' @param const Include constant if **TRUE**.
#' @param trend Include trend if **TRUE**.
#' @param weakly.exog Boolean where we specify
#' whether the stochastic regressors are exogenous or not
#' * `TRUE`: if the regressors are weakly exogenous,
#' * `FALSE`: if the regressors are not weakly exogenous
#' (DOLS is used in this case).
#' @param lags.init,leads.init Scalars defininig the initial number of lags and
#' leads for DOLS.
#' @param max.lag scalar, with the maximum order of the parametric correction.
#' The final order of the parametric correction is selected using the BIC
#' information criterion.
#' @param kernel Kernel for calculating long-run variance
#' * `bartlett`: for Bartlett kernel,
#' * `quadratic`: for Quadratic Spectral kernel,
#' * `NULL` for the Kurozumi's proposal, using Bartlett kernel.
#' @param iter Number of bootstrap iterations.
#' @param bootstrap Type of bootstrapping:
#' * `"sample"`: sampling from residuals with replacement,
#' * `"Cavaliere-Taylor"`: multiplying residuals by \eqn{N(0, 1)}-distributed
#' variable,
#' * `"Rademacher"`: multiplying residuals by Rademacher-distributed variable.
#' @param criterion Information criterion for DOLS lags and leads selection:
#' aic, bic or lwz.
#'
#' @return A list of:
#' * `test`: The value of KPSS test statistic,
#' * `p.value`: The estimates p-value,
#' * `bootstrapped`: Bootstrapped auxiliary statistics.
#'
#' @import doSNOW
#' @import foreach
#' @import parallel
#' @importFrom stats rnorm
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#'
#' @export
KPSS.N.breaks.bootstrap <- function(y,
                                    x,
                                    model,
                                    break.point,
                                    const = FALSE,
                                    trend = FALSE,
                                    weakly.exog = TRUE,
                                    lags.init,
                                    leads.init,
                                    max.lag,
                                    kernel,
                                    iter = 9999,
                                    bootstrap = "sample",
                                    criterion = "bic") {
  if (!is.matrix(y)) y <- as.matrix(y)
  if (!is.null(x)) {
    if (!is.matrix(x)) x <- as.matrix(x)
  }

  n.obs <- nrow(y)

  tmp.kpss <- KPSS.N.breaks(
    y, x,
    model, break.point,
    const, trend,
    weakly.exog,
    lags.init, leads.init,
    max.lag, kernel,
    criterion
  )
  test <- tmp.kpss$test
  u <- tmp.kpss$residuals
  DOLS.lags <- tmp.kpss$DOLS.lags
  DOLS.leads <- tmp.kpss$DOLS.leads
  rm(tmp.kpss)

  if (weakly.exog) {
    xreg <- cbind(
      x,
      trend_kpss_miltiple(model, n.obs, break.point, const, trend)
    )
  } else {
    xreg <- DOLS.vars.N.breaks(
      y, x,
      model, break.point,
      const, trend,
      DOLS.lags, DOLS.leads
    )$xreg
  }

  cores <- detectCores()

  progress.bar <- txtProgressBar(max = iter, style = 3)
  progress <- function(n) setTxtProgressBar(progress.bar, n)

  cluster <- makeCluster(max(cores - 1, 1))
  registerDoSNOW(cluster)

  result <- foreach(
    i = 1:iter,
    .combine = rbind,
    .options.snow = list(progress = progress)
  ) %dopar% {
    if (bootstrap == "sample") {
      temp.y <- sample(u, length(u), replace = TRUE)
    } else if (bootstrap == "Cavaliere-Taylor") {
      z <- rnorm(length(u))
      temp.y <- z * u
    } else if (bootstrap == "Rademacher") {
      z <- sample(c(-1, 1), length(u), replace = TRUE)
      temp.y <- z * u
    }

    resids <- .estimate_ols(temp.y, xreg)$residuals

    if (!is.null(kernel)) {
      temp.test <- .kpss_stat(resids, lr.var.SPC(resids, max.lag, kernel))
    } else {
      temp.test <- .kpss_stat(resids, lr.var.bartlett.AK(resids))
    }

    temp.test
  }

  stopCluster(cluster)

  p.value <- (1 / iter) * sum(I(test <= result))

  return(
    list(
      test = test,
      p.value = p.value,
      bootstrapped = result
    )
  )
}


#' @title
#' Estimating DOLS regression for multiple known break points
#'
#' @param y A time series of interest.
#' @param x A matrix of explanatory stochastic regressors.
#' @param model A scalar or vector of break types:
#' * 1: for the break in const.
#' * 2: for the break in trend.
#' * 3: for the break in const and trend.
#' @param break.point An array of moments of structural breaks.
#' @param const,trend Whether a constant or trend are to be included.
#' @param k.lags,k.leads A number of lags and leads in DOLS regression.
#'
#' @return A list of:
#' * Estimates of coefficients,
#' * Estimates of residuals,
#' * A set of informational criterions values,
#' * \eqn{t}-statistics for the estimates of coefficients.
#'
#' @keywords internal
DOLS.N.breaks <- function(y,
                          x,
                          model,
                          break.point,
                          const = FALSE,
                          trend = FALSE,
                          k.lags,
                          k.leads) {
  if (!is.matrix(y)) y <- as.matrix(y)
  if (is.null(x)) {
    stop("ERROR! Explanatory variables needed for DOLS")
  }
  if (!is.matrix(x)) x <- as.matrix(x)

  dols.vars <- DOLS.vars.N.breaks(
    y, x,
    model, break.point,
    const, trend,
    k.lags, k.leads
  )

  res.OLS <- .estimate_ols(dols.vars$yreg, dols.vars$xreg)

  criterions <- info.criterion(res.OLS$residuals, ncol(dols.vars$xreg))

  return(
    list(
      beta = res.OLS$beta,
      residuals = res.OLS$residuals,
      criterions = criterions,
      t.beta = res.OLS$t.beta
    )
  )
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
DOLS.vars.N.breaks <- function(y,
                               x,
                               model,
                               break.point,
                               const = FALSE,
                               trend = FALSE,
                               k.lags,
                               k.leads) {
  if (is.null(x)) {
    stop("ERROR! Explanatory variables needed for DOLS")
  }
  if (!is.matrix(y)) y <- as.matrix(y)
  if (!is.matrix(x)) x <- as.matrix(x)

  n.obs <- nrow(y)

  d.x.step <- x[2:n.obs, , drop = FALSE] - x[1:(n.obs - 1), , drop = FALSE]
  d.x.lag <- d.x.step
  d.x.lead <- d.x.step

  for (i in 1:k.lags) {
    d.x.lag <- cbind(
      d.x.lag,
      lagn(d.x.step, i)
    )
  }

  for (i in 1:k.leads) {
    d.x.lead <- cbind(
      d.x.lead,
      lagn(d.x.step, -i)
    )
  }

  if (k.lags != 0 && k.leads != 0) {
    lags <- d.x.lag
    leads <- d.x.lead[, (ncol(x) + 1):(ncol(d.x.lead)), drop = FALSE]
    lags.leads <- cbind(lags, leads)
    lags.leads <-
      lags.leads[(k.lags + 1):(n.obs - 1 - k.leads), , drop = FALSE]
  } else if (k.lags != 0 && k.leads == 0) {
    lags <- d.x.lag
    lags.leads <- lags[(k.lags + 1):(n.obs - 1), , drop = FALSE]
  } else if (k.lags == 0 && k.leads != 0) {
    lags <- d.x.lag
    leads <- d.x.lead[, (ncol(x) + 1):(ncol(d.x.lead)), drop = FALSE]
    lags.leads <- cbind(lags, leads)
    lags.leads <- lags.leads[1:(n.obs - 1 - k.leads), , drop = FALSE]
  } else if (k.lags == 0 && k.leads == 0) {
    lags.leads <- d.x.lag
  }
  deter <- trend_kpss_miltiple(model, n.obs, break.point, const, trend)

  xreg <- cbind(
    deter[(k.lags + 2):(n.obs - k.leads), , drop = FALSE],
    x[(k.lags + 2):(n.obs - k.leads), , drop = FALSE],
    lags.leads
  )

  yreg <- y[(k.lags + 2):(n.obs - k.leads), 1, drop = FALSE]

  return(
    list(
      yreg = yreg,
      xreg = xreg
    )
  )
}

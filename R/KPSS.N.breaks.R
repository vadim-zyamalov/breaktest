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
#' @param bp Array of structural breaks.
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
KPSS.mlt <- function(
  y,
  x,
  model,
  bp,
  const = FALSE,
  trend = FALSE,
  weakly.exog = TRUE,
  lags.init,
  leads.init,
  max.lag,
  kernel,
  criterion = "bic"
) {
  if (!is.matrix(y)) y <- as.matrix(y)
  if (!is.null(x) && !is.matrix(x)) x <- as.matrix(x)

  N <- nrow(y)

  if (weakly.exog) {
    xt <- cbind(
      x,
      trend.kpss.miltiple(model, N, bp, const, trend)
    )

    model.est <- .OLS(y, xt)

    beta <- model.est$beta
    resids <- model.est$residuals
    t.beta <- model.est$t.beta
    dols.lags <- 0
    dols.leads <- 0
  } else {
    min.ic <- Inf
    for (nL in lags.init:1) {
      for (nF in leads.init:1) {
        model.est <- .DOLS.multiple(
          y, x, model, bp, const, trend, nL, nF
        )
        .ic <- model.est$criterions
        if (.ic[[criterion]] < min.ic) {
          min.ic <- .ic[[criterion]]
          .beta <- model.est$beta
          .t_beta <- model.est$t.beta
          .resids <- model.est$residuals
          dols.lags <- nL
          dols.leads <- nF
        }
        rm(model.est)
      }
    }
    resids <- .resids
    beta <- .beta
    t.beta <- .t_beta
  }

  test <- ifelse(is.null(kernel),
    .kpss.statistic(resids, .lr.var.kurozumi(resids)),
    .kpss.statistic(resids, .lr.var.spc(resids, max.lag, kernel))
  )

  list(
    beta        = beta,
    statistic   = test,
    residuals   = resids,
    t.beta      = t.beta,
    DOLS.lags   = dols.lags,
    DOLS.leads  = dols.leads,
    break.point = bp
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
kpss.mlt.bootstrap <- function(y,
                               x,
                               model,
                               bp,
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
  if (!is.null(x) && !is.matrix(x)) x <- as.matrix(x)

  N <- nrow(y)

  kpss.result <- KPSS.mlt(
    y, x,
    model, bp,
    const, trend,
    weakly.exog,
    lags.init, leads.init,
    max.lag, kernel,
    criterion
  )
  test <- kpss.result$statistic
  u <- kpss.result$residuals
  .lags.dols <- kpss.result$DOLS.lags
  .leads.dols <- kpss.result$DOLS.leads

  if (weakly.exog) {
    xreg <- cbind(
      x,
      trend.kpss.miltiple(model, N, bp, const, trend)
    )
  } else {
    xreg <- variables.dols.multiple(
      y, x,
      model, bp,
      const, trend,
      .lags.dols, .leads.dols
    )$xreg
  }

  cores <- detectCores()

  .progress <- txtProgressBar(max = iter, style = 3)
  progress <- function(n) setTxtProgressBar(.progress, n)

  cluster <- makeCluster(max(cores - 1, 1))
  registerDoSNOW(cluster)

  result <- foreach(
    i = 1:iter,
    .combine = rbind,
    .options.snow = list(progress = progress)
  ) %dopar% {
    y.loop <- switch(bootstrap,
      sample = sample(u, length(u), replace = TRUE),
      "Cavaliere-Taylor" = rnorm(length(u)) * u,
      "Rademacher" = sample(c(-1, 1), length(u), replace = TRUE) * u
    )

    resids <- .OLS(y.loop, xreg)$residuals

    ifelse(is.null(kernel),
      kpss.result(resids, .lr.var.kurozumi(resids)),
      kpss.result(resids, .lr.var.spc(resids, max.lag, kernel))
    )
  }

  stopCluster(cluster)

  p.value <- (1 / iter) * sum(I(test <= result))

  list(
    statistic    = test,
    p.value      = p.value,
    bootstrapped = result
  )
}

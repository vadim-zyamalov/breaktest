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
kpss.multiple <- function(
  y,
  x,
  model,
  break_point,
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
  if (!is.null(x)) {
    if (!is.matrix(x)) x <- as.matrix(x)
  }

  n_obs <- nrow(y)

  if (weakly.exog) {
    xt <- cbind(
      x,
      trend.kpss.miltiple(model, n_obs, break_point, const, trend)
    )

    .model <- .OLS(y, xt)
    beta <- .model$beta
    resids <- .model$residuals
    t_beta <- .model$t.beta
    .lags_dols <- 0
    .leads_dols <- 0
  } else {
    .ic.min <- Inf
    for (.lags in lags.init:1) {
      for (.leads in leads.init:1) {
        .model <- .DOLS.multiple(
          y, x, model, break_point, const, trend, .lags, .leads
        )
        beta <- .model$beta
        resids <- .model$residuals
        t.beta <- .model$t.beta
        .ic <- .model$criterions
        if (.ic[[criterion]] < .ic.min) {
          .ic.min <- .ic[[criterion]]
          .beta <- beta
          .t_beta <- t.beta
          .resids <- resids
          .lags_dols <- .lags
          .leads_dols <- .leads
        }
        rm(.model)
      }
    }
    resids <- .resids
    beta <- .beta
    t_beta <- .t_beta
  }

  test <- ifelse(is.null(kernel),
    .kpss.statistic(resids, .lr.var.kurozumi(resids)),
    .kpss.statistic(resids, .lr.var.spc(resids, max.lag, kernel))
  )

  list(
    beta = beta,
    test = test,
    residuals = resids,
    t.beta = t_beta,
    DOLS.lags = .lags_dols,
    DOLS.leads = .leads_dols,
    break.point = break_point
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
kpss.multiple.bootstrap <- function(y,
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

  .kpss.stat <- kpss.multiple(
    y, x,
    model, break.point,
    const, trend,
    weakly.exog,
    lags.init, leads.init,
    max.lag, kernel,
    criterion
  )
  test <- .kpss.stat$test
  u <- .kpss.stat$residuals
  .lags.dols <- .kpss.stat$DOLS.lags
  .leads.dols <- .kpss.stat$DOLS.leads

  if (weakly.exog) {
    xreg <- cbind(
      x,
      trend.kpss.miltiple(model, n.obs, break.point, const, trend)
    )
  } else {
    xreg <- variables.dols.multiple(
      y, x,
      model, break.point,
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
    .y.loop <- switch(bootstrap,
      sample = sample(u, length(u), replace = TRUE),
      "Cavaliere-Taylor" = rnorm(length(u)) * u,
      "Rademacher" = sample(c(-1, 1), length(u), replace = TRUE) * u
    )

    resids <- .OLS(.y.loop, xreg)$residuals

    ifelse(is.null(kernel),
      .kpss.stat(resids, .lr.var.kurozumi(resids)),
      .kpss.stat(resids, .lr.var.spc(resids, max.lag, kernel))
    )
  }

  stopCluster(cluster)

  .p.value <- (1 / iter) * sum(I(test <= result))

  list(
    test = test,
    p.value = .p.value,
    bootstrapped = result
  )
}

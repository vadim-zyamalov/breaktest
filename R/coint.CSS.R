#' @title
#' KPSS-based cointegration test with multiple known structural breaks
#'
#' @description
#' Procedure to test the presence of cointegration in the case of  multiple known structural breaks using KPSS test
#'
#' @param y A time series of interest.
#' @param x A matrix of explanatory stochastic regressors.
#' If it's NULL then the result of this function is just a KPSS test for y.
#' @param const,trend Whether a constant or trend should be included.
#' @param break.type A single value or vector of
#' * "c": for the break in const,
#' * "t": for the break in trend,
#' * "ct": for the break in const and trend.
#' @param break.point Array of structural break moments.
#' @param break.coint Whether breaks in cointegrating relation are needed to be included.
#' See Carrion-i-Silvestre & Sansó (2006) for details.
#' @param weakly.exog Boolean where we specify
#' whether the stochastic regressors are exogenous or not
#' * `TRUE`: if the regressors are weakly exogenous,
#' * `FALSE`: if the regressors are not weakly exogenous
#' (DOLS is used in this case).
#' @param lags.init,leads.init Scalars defininig the initial number of lags and leads for DOLS.
#' @param max.lag scalar, with the maximum order of the parametric correction.
#' The final order of the parametric correction is selected using the BIC information criterion.
#' @param kernel Kernel for calculating long-run variance
#' * `bartlett`: for Bartlett kernel,
#' * `quadratic`: for Quadratic Spectral kernel,
#' * `NULL` for the Kurozumi's proposal, using Bartlett kernel.
#' @param criterion Information criterion for DOLS lags and leads selection: aic, bic, hq, or lwz,
#' @param boot.p Whether bootstrapped p-values should be returned.
#' @param boot.iter The number of bootstrap iterations,
#' @param boot.type Method of building synthetic \eqn{y}:
#' * `sample` --- using a sample from residuals,
#' * `Cavaliere-Taylor` --- multiplying residuals by \eqn{N(0, 1)},
#' * `Rademacher` --- multiplying residuals by \eqn{\pm 1}.
#' @param ... A dummy parameter for technical purposes. Just ignore it.
#'
#' @details
#' The code provided is based on the original code
#' by Carrion-i-Silvestre and Sansó.
#'
#' If you want to get results like in Carrion-i-Silvestre & Sansó (2006) then you should
#' * provide a `break.point` of length 1;
#' * for the specification of `An` set `const=TRUE`, `trend=FALSE`, `break.type=c("c")`;
#' * for the specification of `A` set `const=TRUE`, `trend=TRUE`, `break.type=c("c")`;
#' * for the specification of `B` set `const=TRUE`, `trend=TRUE`, `break.type=c("t")`;
#' * for the specification of `C` set `const=TRUE`, `trend=TRUE`, `break.type=c("ct")`;
#' * for the specification of `D` set `const=TRUE`, `trend=FALSE`, `break.type=c("c")`;
#' * for the specification of `E` set `const=TRUE`, `trend=TRUE`, `break.type=c("ct")`, and `break.coint=TRUE`.
#'
#' If you want to get results like in Carrion-i-Silvestre & Sansó (2007) then you should
#' * provide a `break.point` of length 2;
#' * for the specification of `AAn` set `const=TRUE`, `trend=FALSE`, `break.type=c("c", "c")`;
#' * for the specification of `AA` set `const=TRUE`, `trend=TRUE`, `break.type=c("c", "c")`;
#' * for the specification of `BB` set `const=TRUE`, `trend=TRUE`, `break.type=c("t", "t")`;
#' * for the specification of `CC` set `const=TRUE`, `trend=TRUE`, `break.type=c("ct", "ct")`;
#' * for the specification of `AB-BA` set `const=TRUE`, `trend=FALSE`, `break.type=c("c", "t")`;
#' * for the specification of `AC-CA` set `const=TRUE`, `trend=TRUE`, `break.type=c("c", "ct")`;
#' * for the specification of `BC-CB` set `const=TRUE`, `trend=TRUE`, `break.type=c("t", "ct")`.
#'
#' @return An object of class `bt_kpss` containing
#' * `statistic`: the value of test statistic,
#' * `coefficients`: OLS/DOLS estimates of the coefficients,
#' * `se.coefs`: standard errors of the coefficients above,
#' * `t.stats`: \eqn{t}-statistics for the coefficients above,
#' * `residuals`: Residuals of the model,
#' * `t.beta`: \eqn{t}-statistics for `beta`,
#' * `fitted.values`: fitted values of the estimated model,
#' * `endog`: final \eqn{y} vector used in the estimated model,
#' * `exog`: final \eqn{x} matrix used in the estimated model,
#' * `DOLS.lags`: The estimated number of lags and leads in DOLS,
#' * `break.type`, `break.point`, `break.coint`: breaks specification,
#' * `criterions`: values of the information criterions for the estimated model,
#' * `lags`, `leads`: number of lags and leads in the estimated model.
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
coint.CSS <- function(
  y,
  x = NULL,
  const = FALSE,
  trend = FALSE,
  break.type,
  break.point,
  break.coint = FALSE,
  weakly.exog = TRUE,
  lags.init,
  leads.init,
  max.lag,
  kernel,
  criterion = "bic",
  boot.p = FALSE,
  boot.type = "sample",
  boot.iter = 999
) {
  if (!is.matrix(y)) y <- as.matrix(y)
  if (!is.null(x) && !is.matrix(x)) x <- as.matrix(x)

  N <- nrow(y)

  if (weakly.exog || is.null(x)) {
    mXdu <- NULL
    if (!is.null(x) && break.coint) {
      for (bp in break.point) {
        mXdu <- cbind(mXdu, sweep(x, 1, .du(bp, N), `*`))
      }
    }

    mX <- cbind(
      x,
      trend.variables(break.type, N, break.point, const, trend),
      mXdu
    )

    result <- OLS.reg(y, mX)
    result <- c(
      result,
      list(
        break.type = break.type,
        break.point = break.point,
        break.coint = break.coint,
        criterions = info.criterions(result$residuals, ncol(mX)),
        lags = 0,
        leads = 0
      )
    )
  } else {
    minIC <- Inf
    for (nL in rev(seq_len(lags.init))) {
      for (nF in rev(seq_len(leads.init))) {
        model.est <- DOLS.many(
          y,
          x,
          const,
          trend,
          break.type,
          break.point,
          break.coint,
          nL,
          nF
        )
        .ic <- model.est$criterions
        if (.ic[[criterion]] < minIC) {
          minIC <- .ic[[criterion]]
          result <- model.est
        }
      }
    }
  }

  test <- ifelse(is.null(kernel),
    .kpss.statistic(result$residuals, .lr.var.kurozumi(result$residuals)),
    .kpss.statistic(result$residuals, .lr.var.spc(result$residuals, max.lag, kernel))
  )

  result$statistic <- test
  result$lr.var.max.lag <- max.lag
  result$lr.var.kernel <- kernel
  class(result) <- "bt_kpss"

  if (boot.p) {
    result$p.value <- bootstrap(result, boot.iter, boot.type)
  }

  result
}

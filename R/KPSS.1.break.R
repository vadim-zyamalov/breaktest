#' @title
#' KPSS-test with known structural break
#'
#' @description
#' Computes the cointegration test with one known structural break.
#'
#' @details
#' The code provided is the original GAUSS code by Carrion-i-Silvestre and Sansó
#' ported to R.
#'
#' See Carrion-i-Silvestre and Sansó (2006) for further details.
#'
#' @param y A time series of interest.
#' @param x A matrix of explanatory stochastic regressors.
#' @param model A scalar equal to
#' * 1: for model An,
#' * 2: for model A,
#' * 3: for model B,
#' * 4: for model C,
#' * 5: for model D,
#' * 6: for model E.
#' @param bp A position of the break point.
#' @param weakly.exog Exogeneity of the stochastic regressors
#' * `TRUE`: if the regressors are weakly exogenous,
#' * `FALSE`: if the regressors are not weakly exogenous
#' (DOLS is used in this case).
#' @param ll.init A scalar, defines the initial number of leads and lags
#' for DOLS.
#'
#' @return A list of:
#' * `beta`: DOLS estimates of the coefficients,
#' * `tests`: SC test (coinKPSS-test),
#' * `resid`: Residuals of the model,
#' * `t.beta`: Individual significance t-statistics,
#' * `break_point`: Break point.
#'
#' @references
#' Carrion-i-Silvestre, Josep Lluís, and Andreu Sansó.
#' “Testing the Null of Cointegration with Structural Breaks.”
#' Oxford Bulletin of Economics and Statistics 68, no. 5 (October 2006): 623–46.
#' https://doi.org/10.1111/j.1468-0084.2006.00180.x.
#'
#' @export
KPSS.1br <- function(y,
                     x,
                     model,
                     bp,
                     weakly.exog = TRUE,
                     ll.init) {
  if (!is.matrix(y)) y <- as.matrix(y)
  if (!is.null(x) && !is.matrix(x)) x <- as.matrix(x)

  N <- nrow(y)
  dmodel <- c(1, 2, 3, 4, 1, 4)

  if (model < 0 && model > 6) {
    stop("ERROR! kpss.single: try to specify the deterministic component again")
  }

  if (weakly.exog) {
    deter <- if (model != 0) trend.kpss.single(dmodel[model], N, bp) else NULL
    xdu <- if (model %in% c(5, 6)) sweep(x, 1, deter[, 2, drop = FALSE], `*`) else NULL
    xt <- cbind(deter, x, xdu)

    model.est <- OLS.reg(y, xt)
    beta <- model.est$beta
    resids <- model.est$residuals
    t.beta <- model.est$t.beta
  } else {
    bic <- Inf
    for (i in ll.init:1) {
      model.est <- DOLS.1br(y, x, model, bp, i, i)
      if (model.est$bic < bic) {
        bic <- model.est$bic
        beta <- model.est$beta
        t.beta <- model.est$t.beta
        resids <- model.est$residuals
      }
    }
  }

  list(
    beta = beta,
    statistic = .kpss.statistic(resids, .lr.var.kurozumi(resids)),
    residuals = resids,
    t.beta = t.beta,
    break.point = bp
  )
}


#' @title
#' KPSS-test of cointegration
#'
#' @description
#' Procedure for testing the null of cointegration in the possible presence of
#' structural breaks.
#'
#' @details
#' Computes the cointegration test with one unknown structural break
#' where the break point is estimated either minimizing the value of
#' the statistic or the sum of the squared residuals.
#' The estimation of the cointegrating relationship bases on DOLS.
#'
#' The code provided is the original GAUSS code ported to R.
#'
#' See Carrion-i-Silvestre and Sansó (2006) for further details.
#'
#' @param y A time series of interest.
#' @param x A matrix of explanatory stochastic regressors.
#' @param model A scalar equal to
#' * 1: for model An,
#' * 2: for model A,
#' * 3: for model B,
#' * 4: for model C,
#' * 5: for model D,
#' * 6: for model E.
#' @param weakly.exog Exogeneity of the stochastic regressors
#' * `TRUE`: if the regressors are weakly exogenous,
#' * `FALSE`: if the regressors are not weakly exogenous
#' (DOLS is used in this case).
#' @param ll.init Scalar, defines the initial number of leads and lags for DOLS.
#'
#' @return (2x2)-matrix, where the first rows gives the value of
#' the min(SC) test and the estimated break point;
#' the second row gives the value of the SC statistic,
#' where the break point is estimated as min(SSR).
#'
#' @references
#' Carrion-i-Silvestre, Josep Lluís, and Andreu Sansó.
#' “Testing the Null of Cointegration with Structural Breaks.”
#' Oxford Bulletin of Economics and Statistics 68, no. 5 (October 2006): 623–46.
#' https://doi.org/10.1111/j.1468-0084.2006.00180.x.
#'
#' @export
KPSS.1br.unknown <- function(y,
                             x,
                             model,
                             weakly.exog,
                             ll.init) {
  if (!is.matrix(y)) y <- as.matrix(y)
  if (!is.null(x) && !is.matrix(x)) x <- as.matrix(x)

  N <- nrow(y)

  min.test <- Inf
  idx.test <- NULL
  min.rss <- Inf
  idx.rss <- NULL

  for (i in 3:(N - 3)) {
    if (ll.init + 2 < i && i < N - 5 - ll.init) {
      .result <- KPSS.1br(y, x, model, i, weakly.exog, ll.init)
      .rss <- drop(t(.result$residuals) %*% .result$residuals)

      if (.result$test < min.test) {
        min.test <- .result$test
        idx.test <- i
      }

      if (.rss < min.rss) {
        min.rss <- .rss
        idx.rss <- i
      }
    }
  }

  result <- matrix(c(min.test, min.rss, idx.test, idx.rss), ncol = 2)
  colnames(result) <- c("stat", "tb")
  rownames(result) <- c("min(stat)", "min(RSS)")

  result
}

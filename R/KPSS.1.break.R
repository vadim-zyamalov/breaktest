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
#' @param break.point A position of the break point.
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
kpss_single <- function(y,
                        x,
                        model,
                        tb,
                        weakly_exog = TRUE,
                        n_lag_lead) {
  if (!is.matrix(y)) y <- as.matrix(y)
  if (!is.null(x)) {
    if (!is.matrix(x)) x <- as.matrix(x)
  }

  n_obs <- nrow(y)

  if (model < 0 && model > 6) {
    stop("ERROR: Try to specify the deterministic component again")
  }

  if (weakly_exog) {
    if (model == 0) {
      xt <- x
    } else if (1 <= model && model <= 4) {
      deter <- trend_kpss_single(model, n_obs, tb)
      xt <- cbind(deter, x)
    } else if (model == 5) {
      deter <- trend_kpss_single(1, n_obs, tb)
      xdu <- sweep(x, 1, deter[, 2, drop = FALSE], `*`)
      xt <- cbind(deter, x, xdu)
    } else if (model == 6) {
      deter <- trend_kpss_single(4, n_obs, tb)
      xdu <- sweep(x, 1, deter[, 2, drop = FALSE], `*`)
      xt <- cbind(deter, x, xdu)
    }

    .res_ols <- .estimate_ols(y, xt)
    beta <- .res_ols$beta
    resids <- .res_ols$residuals
    t_beta <- .res_ols$t.beta
  } else {
    bic <- Inf
    for (i in n_lag_lead:1) {
      .res_dols <- .estimate_dols_single(y, x, model, tb, i, i)
      if (.res_dols$bic < bic) {
        bic <- .res_dols$bic
        beta <- .res_dols$beta
        t_beta <- .res_dols$t.beta
        resids <- .res_dols$residuals
      }
    }
  }

  test <- .kpss_stat(resids, .variance_lr_kurozumi(resids))

  list(
    beta = beta,
    test = test,
    residuals = resids,
    t.beta = t_beta,
    break.point = tb
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
kpss_single_unknown <- function(y,
                                x,
                                model,
                                weakly_exog,
                                lag_lead) {
  if (!is.matrix(y)) y <- as.matrix(y)
  if (!is.matrix(x)) x <- as.matrix(x)

  n_obs <- nrow(y)

  min_test <- Inf
  idx_test <- NULL
  min_rss <- Inf
  idx_rss <- NULL

  for (i in 3:(n_obs - 3)) {
    if (lag_lead + 2 < i && i < n_obs - 5 - lag_lead) {
      .result <- kpss_single(y, x, model, i, weakly_exog, lag_lead)
      .rss <- drop(t(.result$residuals) %*% .result$residuals)

      if (.result$test < min_test) {
        min_test <- .result$test
        idx_test <- i
      }

      if (.rss < min_rss) {
        min_rss <- .rss
        idx_rss <- i
      }
    }
  }

  result <- matrix(c(min_test, min_rss, idx_test, idx_rss), ncol = 2)
  colnames(result) <- c("stat", "tb")
  rownames(result) <- c("min(stat)", "min(RSS)")
  result
}

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
    rm(.res_ols)
  } else {
    bic <- Inf
    for (i in n_lag_lead:1) {
      .res_dols <- .dols_single(y, x, model, tb, i, i)
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
  result
}


#' @title
#' Estimating DOLS regression for multiple known break points
#'
#' @param y A time series of interest.
#' @param x A matrix of explanatory stochastic regressors.
#' @param model See Carrion-i-Silvestre and Sansó (2006)
#' * 1: for model An,
#' * 2: for model A,
#' * 3: for model B,
#' * 4: for model C,
#' * 5: for model D,
#' * 6: for model E.
#' @param break.point A position of the break point.
#' @param k.lags,k.leads A number of lags and leads in DOLS regression.
#'
#' @return A list of:
#' * Estimates of coefficients,
#' * Estimates of residuals,
#' * A value of BIC,
#' * \eqn{t}-statistics for the estimates of coefficients.
#'
#' @references
#' Carrion-i-Silvestre, Josep Lluís, and Andreu Sansó.
#' “Testing the Null of Cointegration with Structural Breaks.”
#' Oxford Bulletin of Economics and Statistics 68, no. 5 (October 2006): 623–46.
#' https://doi.org/10.1111/j.1468-0084.2006.00180.x.
#'
#' @keywords internal
.dols_single <- function(y,
                         x,
                         model,
                         tb,
                         k_lags,
                         k_leads) {
  if (is.null(x)) {
    stop("ERROR! Explanatory variables needed for DOLS")
  }
  if (!is.matrix(y)) y <- as.matrix(y)
  if (!is.matrix(x)) x <- as.matrix(x)

  n_obs <- nrow(y)

  .diff_x <- diff(x)
  .diff_x_lag <- .diff_x
  .diff_x_lead <- .diff_x

  for (i in 1:k_lags) {
    .diff_x_lag <- cbind(
      .diff_x_lag,
      lagn(.diff_x, i)
    )
  }

  for (i in 1:k_leads) {
    .diff_x_lead <- cbind(
      .diff_x_lead,
      lagn(.diff_x, -i)
    )
  }

  if (k_lags != 0 && k_leads != 0) {
    lags <- .diff_x_lag
    leads <- .diff_x_lead[, (ncol(x) + 1):(ncol(.diff_x_lead)), drop = FALSE]
    lags_leads <- cbind(lags, leads)
    lags_leads <-
      lags_leads[(k_lags + 1):(n_obs - 1 - k_leads), , drop = FALSE]
  } else if (k_lags != 0 && k_leads == 0) {
    lags <- .diff_x_lag
    lags_leads <- lags[(k_lags + 1):(n_obs - 1), , drop = FALSE]
  } else if (k_lags == 0 && k_leads != 0) {
    lags <- .diff_x_lag
    leads <- .diff_x_lead[, (ncol(x) + 1):(ncol(.diff_x_lead)), drop = FALSE]
    lags_leads <- cbind(lags, leads)
    lags_leads <- lags_leads[1:(n_obs - 1 - k_leads), , drop = FALSE]
  } else if (k_lags == 0 && k_leads == 0) {
    lags_leads <- .diff_x_lag
  }

  if (model == 0) {
    xreg <- cbind(
      x[(k_lags + 2):(n_obs - k_leads), , drop = FALSE],
      lags_leads
    )
  } else if (model >= 1 && model <= 4) {
    deter <- trend_kpss_single(model, n_obs, tb)
    xreg <- cbind(
      deter[(k_lags + 2):(n_obs - k_leads), , drop = FALSE],
      x[(k_lags + 2):(n_obs - k_leads), , drop = FALSE],
      lags_leads
    )
  } else if (model == 5) {
    deter <- trend_kpss_single(1, n_obs, tb)
    xdu <- sweep(x, 1, deter[, 2, drop = FALSE], `*`)
    xreg <- cbind(
      deter[(k_lags + 2):(n_obs - k_leads), , drop = FALSE],
      x[(k_lags + 2):(n_obs - k_leads), , drop = FALSE],
      xdu[(k_lags + 2):(n_obs - k_leads), , drop = FALSE],
      lags_leads
    )
  } else if (model == 6) {
    deter <- trend_kpss_single(4, n_obs, tb)
    xdu <- sweep(x, 1, deter[, 2, drop = FALSE], `*`)
    xreg <- cbind(
      deter[(k_lags + 2):(n_obs - k_leads), , drop = FALSE],
      x[(k_lags + 2):(n_obs - k_leads), , drop = FALSE],
      xdu[(k_lags + 2):(n_obs - k_leads), , drop = FALSE],
      lags_leads
    )
  }

  .res_ols <- .estimate_ols(
    y[(k_lags + 2):(n_obs - k_leads), 1, drop = FALSE],
    xreg
  )

  bic <- log(drop(t(.res_ols$residuals) %*% .res_ols$residuals) / nrow(xreg)) +
    ncol(xreg) * log(nrow(xreg)) / nrow(xreg)

  list(
    beta   = .res_ols$beta,
    resid  = .res_ols$residuals,
    bic    = bic,
    t.beta = .res_ols$t.beta
  )
}

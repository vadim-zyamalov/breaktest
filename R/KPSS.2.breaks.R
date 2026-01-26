#' @title
#' KPSS-test with 2 known structural breaks
#'
#' @description
#' Procedure to compute the KPSS test with two structural breaks
#'
#' @details
#' The break points are known.
#'
#' The code provided is the original GAUSS code by Carrion-i-Silvestre and Sansó
#' ported to R.
#'
#' See Carrion-i-Silvestre and Sansó (2007) for further details.
#'
#' @param y A time series of interest.
#' @param model A scalar equal to
#' * 1: for the AA (without trend) model,
#' * 2: for the AA (with trend) model,
#' * 3: for the BB model,
#' * 4: for the CC model,
#' * 5: for the AC-CA model.
#' @param break.point Positions for the first and second structural breaks
#' (respective to the origin which is 1).
#' @param max.lag A scalar, with the maximum order of the parametric correction.
#' The final order of the parametric correction is selected using the
#' BIC information criterion.
#' @param kernel Kernel for calculating long-run variance
#' * `bartlett`: for Bartlett kernel,
#' * `quadratic`: for Quadratic Spectral kernel,
#' * `NULL` for the Kurozumi's proposal, using Bartlett kernel.
#'
#' @return A list of:
#' * `beta`: DOLS estimates of the coefficients,
#' * `tests`: SC test (coinKPSS-test),
#' * `resid`: Residuals of the model,
#' * `t.beta`: \eqn{t}-statistics for `beta`,
#' * `break_point`: Break points.
#'
#' @references
#' Carrion-i-Silvestre, Josep Lluís, and Andreu Sansó.
#' “The KPSS Test with Two Structural Breaks.”
#' Spanish Economic Review 9, no. 2 (May 16, 2007): 105–27.
#' https://doi.org/10.1007/s10108-006-9017-8.
#'
#' @export
kpss_double <- function(y,
                        model,
                        break_point,
                        max_lag,
                        kernel) {
  if (!is.matrix(y)) y <- as.matrix(y)

  n_obs <- nrow(y)

  z <- trend_kpss_double(model, n_obs, break_point)

  .model <- .estimate_ols(y, z)

  test <- ifelse(!is.null(kernel),
    .kpss_stat(
      .model$residuals,
      .variance_lr_spc(.model$residuals, max_lag, kernel)
    ),
    .kpss_stat(
      .model$residuals,
      .variance_lr_kurozumi(.model$residuals)
    )
  )

  list(
    beta        = .model$beta,
    test        = test,
    residuals   = .model$residuals,
    t.beta      = .model$t.beta,
    break_point = break_point
  )
}


#' @title
#' KPSS-test with 2 unknown structural breaks
#'
#' @description
#' Procedure to compute the KPSS test with two structural breaks
#'
#' @details
#' The break points are known
#'
#' The code provided is the original GAUSS code ported to R.
#'
#' See Carrion-i-Silvestre and Sansó (2007) for further details.
#'
#' @param y A time series of interest.
#' @param model A scalar equal to
#' * 1: for the AA (without trend) model,
#' * 2: for the AA (with trend) model,
#' * 3: for the BB model,
#' * 4: for the CC model,
#' * 5: for the AC-CA model.
#' @param max.lag A scalar, with the maximum order of the parametric correction.
#' The final order of the parametric correction is selected using
#' the BIC information criterion.
#' @param kernel Kernel for calculating long-run variance
#' * `bartlett`: for Bartlett kernel,
#' * `quadratic`: for Quadratic Spectral kernel,
#' * `NULL` for the Kurozumi's proposal, using Bartlett kernel.
#'
#' @return Value of test statistic.
#'
#' @references
#' Carrion-i-Silvestre, Josep Lluís, and Andreu Sansó.
#' “Testing the Null of Cointegration with Structural Breaks.”
#' Oxford Bulletin of Economics and Statistics 68, no. 5 (October 2006): 623–46.
#' https://doi.org/10.1111/j.1468-0084.2006.00180.x.
#'
#' @export
kpww_double_unknown <- function(y,
                                model,
                                max_lag = 0,
                                kernel = "bartlett") {
  if (!is.matrix(y)) y <- as.matrix(y)

  .segments <- .segments_ols_double(y, model)

  test <- ifelse(!is.null(kernel),
    .kpss_stat(
      .segments$residuals,
      .variance_lr_spc(.segments$residuals, max_lag, kernel)
    ),
    .kpss_stat(
      .segments$residuals,
      .variance_lr_kurozumi(.segments$residuals)
    )
  )

  list(
    test = test,
    tb1  = .segments$tb1,
    tb2  = .segments$tb2
  )
}

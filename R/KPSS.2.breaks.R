#' KPSS-test with 2 known structural breaks
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
#' @param y (Tx1)-vector of time series
#' @param model \describe{
#' \item{1}{for the AA (without trend) model.}
#' \item{2}{for the AA (with trend) model.}
#' \item{3}{for the BB model.}
#' \item{4}{for the CC model.}
#' \item{5}{for the AC-CA model.}
#' }
#' @param tb1 The first break point.
#' @param tb2 The second break point.
#' @param kmax scalar, with the maximum order of the parametric correction.
#' The final order of the parametric correction is selected using the
#' BIC information criterion.
#' @param kernel Kernel for calculating long-run variance
#' \describe{
#' \item{bartlett}{for Bartlett kernel.}
#' \item{quadratic}{for Quadratic Spectral kernel.}
#' \item{NULL}{for the Kurozumi's proposal, using Bartlett kernel.}
#' }
#'
#' @return \describe{
#' \item{beta}{DOLS estimates of the coefficients regressors.}
#' \item{tests}{SC test (coinKPSS-test).}
#' \item{resid}{Residuals of the model.}
#' \item{break_point}{Break points.}
#' }
#'
#' @importFrom zeallot %<-%
#' @export
KPSS.2.breaks <- function(y, model, break.point, kmax, kernel) {
    if (!is.matrix(y)) y <- as.matrix(y)

    N <- nrow(y)

    z <- determinants.KPSS.2.breaks(model, N, break.point)
    c(beta, resid, ., t.beta) %<-% OLS(y, z)

    S.t <- apply(resid, 2, cumsum)

    if (!is.null(kernel))
        test <- N^(-2) * drop(t(S.t) %*% S.t) /
            alrvr.kernel(resid, corr.max, kernel)
    else
        test <- N^(-2) * drop(t(S.t) %*% S.t) / alrvr(resid)

    return(
        list(
            beta = beta,
            test = test,
            resid = resid,
            t.beta = t.beta,
            break_point = break.point
        )
    )
}
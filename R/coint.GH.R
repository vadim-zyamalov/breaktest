#' Gregory-Hansen test for the absense of cointegration
#'
#' @description
#' Gregory and Hansen (1996) test for the null hypothesis of no cointegration
#' under a possible structural break at the unknown moment of time.
#'
#' The authors proposed ADF- and Z-type tests, slightly modified to allow
#' the presence of a possible regime shift. Three type of shifts are allowed:
#' * a shift in the constant,
#' * a shift in the constand with the trend included,
#' * and a shift in the constant and the cointegrating vector.
#'
#' Critical values are calculated via the adopted MacKinnon procedure of
#' estimating the model for the response surface.
#'
#' @param ... Variables of interest.
#' @param shift Expected break type.
#' @param trim The trimming parameter to calculate break moment bounds.
#' @param max.lag The maximum number of lags for the internal ADF testing.
#' @param criterion The criterion for lag selection.
#' @param add.cvals Whether critical values are to be returned.
#' This argument is needed to suppress the calculation of critical values
#' during the precalculation of tables needed for the p-values estimating.
#'
#' @return An object of type `cointGH`. It's a list of
#' * `shift`: shift type,
#' * `Za`: \eqn{MZ_\alpha} statistic and c.v.,
#' * `Zt`: \eqn{MZ_t} statistic and c.v.,
#' * `ADF`: \eqn{ADF} statistic and c.v..
#'
#' @references
#' MacKinnon, James G.
#' “Critical Values for Cointegration Tests.”
#' Working Paper. Working Paper.
#' Economics Department, Queen’s University, January 2010.
#' https://ideas.repec.org/p/qed/wpaper/1227.html.
#'
#' Gregory, Allan W., and Bruce E. Hansen.
#' “Residual-Based Tests for Cointegration in Models with Regime Shifts.”
#' Journal of Econometrics 70, no. 1 (January 1, 1996): 99–126.
#' https://doi.org/10.1016/0304-4076(69)41685-7.
#'
#' @export
coint.GH <- function(...,
                     shift = "level",
                     trim = 0.15,
                     max.lag = 10,
                     criterion = "aic",
                     add.cvals = TRUE) {
  if (...length() < 2) {
    stop("ERROR! coint.test.GH: Two or more variables are needed")
  }

  if (!shift %in% c("level", "level-trend", "regime")) {
    stop("ERROR! coint.test.GH: Unknown model specification")
  }

  if (length(unique(sapply(list(...), length))) != 1) {
    stop("ERROR! coint.test.GH: Series of different length")
  }

  y1 <- as.matrix(...elt(1))
  y2 <- NULL
  for (i in 2:...length()) {
    y2 <- cbind(y2, ...elt(i))
  }

  N <- nrow(y1)
  first.break <- trunc(trim * N)
  last.break <- trunc((1 - trim) * N)

  res.Za <- Inf
  res.Zt <- Inf
  res.ADF <- Inf

  for (tb in first.break:last.break) {
    phi <- .du(tb, N)

    x <- switch(shift,
      "level" = cbind(.const(N), phi, y2),
      "level-trend" = cbind(.const(N), phi, .trend(N), y2),
      "regime" = cbind(.const(N), phi, y2, phi * y2)
    )

    e <- OLS.reg(y1, x)$residuals
    Le <- .lagn(e, 1)

    rho <- sum(e * Le, na.rm = TRUE) /
      sum(Le^2, na.rm = TRUE)

    nu <- e - rho * .lagn(e, 1, na = 0)

    lrv <- .lr.var.bartlett(nu)
    lambda <- (lrv - drop(t(nu) %*% nu) / N) / 2

    rho.star <- sum(e * Le - lambda, na.rm = TRUE) /
      sum(Le^2, na.rm = TRUE)

    res.Za <- min(
      N * (rho.star - 1),
      res.Za
    )
    res.Zt <- min(
      (rho.star - 1) * sqrt(sum(e[1:(N - 1)]^2) / lrv),
      res.Zt
    )

    res.ADF <- min(
      ADF.test(
        e,
        const = FALSE, trend = FALSE,
        max.lag = max.lag,
        criterion = criterion,
        modified.criterion = TRUE
      )$t.alpha,
      res.ADF
    )
  }

  result <- list(
    shift = shift,
    Za = list(
      statistic = res.Za
    ),
    Zt = list(
      statistic = res.Zt
    ),
    ADF = list(
      statistic = res.ADF
    )
  )

  ## Critical values
  if (add.cvals) {
    m <- min(6, ncol(y2))

    result$Za$asy.cv <- .cval_coint_gh[[shift]]$Za$b0[m, 3]
    result$Za$cv <- .cval_coint_gh[[shift]]$Za$b0[m, 3] +
      .cval_coint_gh[[shift]]$Za$b1[m, 3] / N

    result$Zt$asy.cv <- .cval_coint_gh[[shift]]$Zt$b0[m, 3]
    result$Zt$cv <- .cval_coint_gh[[shift]]$Zt$b0[m, 3] +
      .cval_coint_gh[[shift]]$Zt$b1[m, 3] / N

    result$ADF$asy.cv <- .cval_coint_gh[[shift]]$ADF$b0[m, 3]
    result$ADF$cv <- .cval_coint_gh[[shift]]$ADF$b0[m, 3] +
      .cval_coint_gh[[shift]]$ADF$b1[m, 3] / N
  }

  class(result) <- "cointGH"
  result
}

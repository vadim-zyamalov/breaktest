#' @title
#' MDF procedure for multiple unknown breaks.
#'
#' @details
#' The code provided is the original Ox code by Skrobotov (2018)
#' ported to R.
#'
#' @param y A time series of interest.
#' @param const Whether the constant term should be included.
#' @param breaks Number of breaks.
#' @param breaks.star Number of breaks got from the Kejrival-Perron procedure.
#' @param trim Trimming value for a possible break date bounds.
#' @param ZA Whether ZA variant should be used.
#'
#' @return A list of sublists each containing
#' * The value of statistic: \eqn{MDF-GLS}, \eqn{MDF-OLS},
#' * The asymptotic critical values.
#' \eqn{UR} values are included as well.
#'
#' @export
MDF.mlt <- function(y,
                    const = FALSE,
                    breaks = 1,
                    breaks.star = 1,
                    trim = 0.15,
                    ZA = FALSE) {
  if (!is.matrix(y)) y <- as.matrix(y)

  N <- nrow(y)

  ## Critical values ##
  model <- if (const && ZA) {
    "cz"
  } else if (const && !ZA) {
    "c"
  } else if (!const %% ZA) {
    "nz"
  } else {
    "n"
  }

  cv.DF.OLS.t <- .cval_MDF_multiple[[model]]$cv.DF.OLS.t
  cv.MDF.OLS1 <- .cval_MDF_multiple[[model]]$cv.MDF.OLS1
  cv.MDF.OLS2 <- .cval_MDF_multiple[[model]]$cv.MDF.OLS2
  cv.MDF.OLS3 <- .cval_MDF_multiple[[model]]$cv.MDF.OLS3
  cv.DF.GLS.t <- .cval_MDF_multiple[[model]]$cv.DF.GLS.t
  cv.MDF.GLS1 <- .cval_MDF_multiple[[model]]$cv.MDF.GLS1
  cv.MDF.GLS2 <- .cval_MDF_multiple[[model]]$cv.MDF.GLS2
  cv.MDF.GLS3 <- .cval_MDF_multiple[[model]]$cv.MDF.GLS3
  sap.ur3.ols <- .cval_MDF_multiple[[model]]$sap.ur3.ols
  sap.ur3.olsgls <- .cval_MDF_multiple[[model]]$sap.ur3.olsgls
  sap.ur3.1ols <- .cval_MDF_multiple[[model]]$sap.ur3.1ols
  sap.ur3.1olsgls <- .cval_MDF_multiple[[model]]$sap.ur3.1olsgls
  sap.ur3.2ols <- .cval_MDF_multiple[[model]]$sap.ur3.2ols
  sap.ur3.2olsgls <- .cval_MDF_multiple[[model]]$sap.ur3.2olsgls
  sap.ur3.3ols <- .cval_MDF_multiple[[model]]$sap.ur3.3ols
  sap.ur3.3olsgls <- .cval_MDF_multiple[[model]]$sap.ur3.3olsgls
  sap.ur2.ols <- .cval_MDF_multiple[[model]]$sap.ur2.ols
  sap.ur2.olsgls <- .cval_MDF_multiple[[model]]$sap.ur2.olsgls
  sap.ur2.1ols <- .cval_MDF_multiple[[model]]$sap.ur2.1ols
  sap.ur2.1olsgls <- .cval_MDF_multiple[[model]]$sap.ur2.1olsgls
  sap.ur2.2ols <- .cval_MDF_multiple[[model]]$sap.ur2.2ols
  sap.ur2.2olsgls <- .cval_MDF_multiple[[model]]$sap.ur2.2olsgls
  sap.cv.ur.2 <- .cval_MDF_multiple[[model]]$sap.cv.ur.2
  sap.cv.ur.3 <- .cval_MDF_multiple[[model]]$sap.cv.ur.3


  ## Start ##
  max.lag <- trunc(12 * (N / 100)^(1 / 4))

  first.break <- trunc(trim * N) + 1
  width <- first.break - 1
  last.break <- trunc((1 - trim) * N) + 1

  x <- cbind(1, 1:N)

  ## GLS case
  resid.GLS.t <- .GLS(y, x, -13.5)$residuals

  resid.OLS.t <- .OLS(y, x)$residuals
  DF.OLS.t <- ADF.test(resid.OLS.t,
    const = FALSE, trend = FALSE,
    max.lag = max.lag,
    criterion = "aic",
    modified.criterion = TRUE
  )
  k.t <- max(1, DF.OLS.t$lag)

  DF.GLS.t <- ADF.test(resid.GLS.t,
    const = FALSE, trend = FALSE,
    max.lag = k.t,
    criterion = NULL
  )
  DF.GLS.t <- DF.GLS.t$t.alpha

  ## OLS case
  DF.OLS.t <- ADF.test(resid.OLS.t,
    const = FALSE, trend = FALSE,
    max.lag = k.t,
    criterion = NULL
  )
  DF.OLS.t <- DF.OLS.t$t.alpha

  ## OLS-GLS ##
  ## One break
  MDF.OLS1 <- Inf
  MDF.GLS1 <- Inf

  for (tb1 in first.break:last.break) {
    x <- cbind(
      1,
      1:N,
      if (const) .du(tb1, N) else NULL,
      .dt(tb1, N)
    )

    resid.OLS <- .OLS(y, x)$residuals
    DF1.tb <- ADF.test(resid.OLS,
      const = FALSE, trend = FALSE,
      max.lag = max.lag,
      criterion = "aic",
      modified.criterion = TRUE
    )
    k.t <- max(1, DF1.tb$lag)

    DF1.tb <- ADF.test(resid.OLS,
      const = FALSE, trend = FALSE,
      max.lag = k.t,
      criterion = NULL
    )
    if (!ZA) {
      denom <- 1 - sum(DF1.tb$beta) + DF1.tb$alpha
      stat.OLS <- N * DF1.tb$alpha / denom
    } else {
      stat.OLS <- DF1.tb$t.alpha
    }

    resid.GLS <- .GLS(y, x, -17.6)$residuals
    DF1.tb <- ADF.test(resid.GLS,
      const = FALSE, trend = FALSE,
      max.lag = k.t,
      criterion = NULL
    )

    if (stat.OLS < MDF.OLS1) MDF.OLS1 <- stat.OLS
    if (DF1.tb$t.alpha < MDF.GLS1) MDF.GLS1 <- DF1.tb$t.alpha
  }

  ## Two breaks
  MDF.OLS2 <- Inf
  MDF.GLS2 <- Inf

  for (tb1 in first.break:(last.break - width)) {
    for (tb2 in (tb1 + width):last.break) {
      DU1 <- as.numeric(1:N > tb1)
      DT1 <- DU1 * (1:N - tb1)
      DU2 <- as.numeric(1:N > tb2)
      DT2 <- DU2 * (1:N - tb2)

      x <- cbind(
        1,
        1:N,
        if (const) .du(tb1, N) else NULL,
        .dt(tb1, N),
        if (const) .du(tb2, N) else NULL,
        .dt(tb2, N)
      )

      resid.OLS <- .OLS(y, x)$residuals
      DF2.tb <- ADF.test(resid.OLS,
        const = FALSE, trend = FALSE,
        max.lag = max.lag,
        criterion = "aic",
        modified.criterion = TRUE
      )
      k.t <- max(1, DF2.tb$lag)

      DF2.tb <- ADF.test(resid.OLS,
        const = FALSE, trend = FALSE,
        max.lag = k.t,
        criterion = NULL
      )
      if (!ZA) {
        denom <- 1 - sum(DF2.tb$beta) + DF2.tb$alpha
        stat.OLS <- N * DF2.tb$alpha / denom
      } else {
        stat.OLS <- DF2.tb$t.alpha
      }

      resid.GLS <- .GLS(y, x, -21.5)$residuals
      DF2.tb <- ADF.test(resid.GLS,
        const = FALSE, trend = FALSE,
        max.lag = k.t,
        criterion = NULL
      )

      if (stat.OLS < MDF.OLS2) MDF.OLS2 <- stat.OLS
      if (DF2.tb$t.alpha < MDF.GLS2) MDF.GLS2 <- DF2.tb$t.alpha
    }
  }

  ## Three breaks
  MDF.OLS3 <- Inf
  MDF.GLS3 <- Inf

  for (tb1 in first.break:(last.break - 2 * width)) {
    for (tb2 in (tb1 + width):(last.break - width)) {
      for (tb3 in (tb2 + width):last.break) {
        x <- cbind(
          1,
          1:N,
          if (const) .du(tb1, N) else NULL,
          .dt(tb1, N),
          if (const) .du(tb2, N) else NULL,
          .dt(tb2, N),,
          if (const) .du(tb3, N) else NULL,
          .dt(tb3, N),
        )

        resid.OLS <- .OLS(y, x)$residuals
        DF3.tb <- ADF.test(resid.OLS,
          const = FALSE, trend = FALSE,
          max.lag = max.lag,
          criterion = "aic",
          modified.criterion = TRUE
        )
        k.t <- max(1, DF3.tb$lag)

        DF3.tb <- ADF.test(resid.OLS,
          const = FALSE, trend = FALSE,
          max.lag = k.t,
          criterion = NULL
        )
        if (!ZA) {
          denom <- 1 - sum(DF3.tb$beta) + DF3.tb$alpha
          stat.OLS <- N * DF3.tb$alpha / denom
        } else {
          stat.OLS <- DF3.tb$t.alpha
        }

        resid.GLS <- .GLS(y, x, -25.5)$residuals
        DF3.tb <- ADF.test(resid.GLS,
          const = FALSE, trend = FALSE,
          max.lag = k.t,
          criterion = NULL
        )

        if (stat.OLS < MDF.OLS3) MDF.OLS3 <- stat.OLS
        if (DF3.tb$t.alpha < MDF.GLS3) MDF.GLS3 <- DF3.tb$t.alpha
      }
    }
  }

  ## Alternative break selection
  if (breaks == 2) {
    tbs <- segments.GLS(y, const, TRUE, 2)

    x <- cbind(
      1,
      1:N,
      if (const) .du(tbs[1], N) else NULL,
      .dt(tbs[1], N),
      if (const) .du(tbs[2], N) else NULL,
      .dt(tbs[2], N)
    )

    tmp.OLS <- .OLS(y, x)
    bb <- tmp.OLS$beta
    rr <- tmp.OLS$residuals
    rm(tmp.OLS)
    t.alpha <- bb[1] / sqrt(drop(t(rr) %*% rr) / N)
    t.alpha.2.id <- as.numeric(t.alpha > 1)
  }
  if (breaks == 3) {
    tbs <- segments.GLS(y, const, TRUE, 3)

    x <- cbind(
      1,
      1:N,
      if (const) .du(tbs[1], N) else NULL,
      .dt(tbs[1], N),
      if (const) .du(tbs[2], N) else NULL,
      .dt(tbs[2], N),
      if (const) .du(tbs[3], N) else NULL,
      .dt(tbs[3], N)
    )

    tmp.OLS <- .OLS(y, x)
    bb <- tmp.OLS$beta
    rr <- tmp.OLS$residuals
    rm(tmp.OLS)
    t.alpha <- bb[1] / sqrt(drop(t(rr) %*% rr) / N)
    t.alpha.3.id <- as.numeric(t.alpha > 1)
  }

  ## breaks.star
  if (breaks.star == 0) {
    tbb <- 0
  } else {
    tbb <- segments.GLS(y, const, TRUE, breaks.star)
  }

  if (breaks == 2) {
    ur2.ols.sa <- as.numeric(
      (DF.OLS.t < (sap.cv.ur.2 * sap.ur2.ols * cv.DF.OLS.t)) ||
        (MDF.OLS1 < (sap.cv.ur.2 * sap.ur2.ols * cv.MDF.OLS1)) ||
        (MDF.OLS2 < (sap.cv.ur.2 * sap.ur2.ols * cv.MDF.OLS2))
    )

    ur2.olsgls.sa <- as.numeric(
      (DF.GLS.t < (sap.cv.ur.2 * sap.ur2.olsgls * cv.DF.GLS.t)) ||
        (MDF.GLS1 < (sap.cv.ur.2 * sap.ur2.olsgls * cv.MDF.GLS1)) ||
        (DF.OLS.t < (sap.cv.ur.2 * sap.ur2.olsgls * cv.DF.OLS.t)) ||
        (MDF.OLS1 < (sap.cv.ur.2 * sap.ur2.olsgls * cv.MDF.OLS1)) ||
        (MDF.GLS2 < (sap.cv.ur.2 * sap.ur2.olsgls * cv.MDF.GLS2)) ||
        (MDF.OLS2 < (sap.cv.ur.2 * sap.ur2.olsgls * cv.MDF.OLS2))
    )

    ur2.1ols.sa <- as.numeric(
      (MDF.OLS1 < (sap.cv.ur.2 * sap.ur2.1ols * cv.MDF.OLS1)) ||
        (MDF.OLS2 < (sap.cv.ur.2 * sap.ur2.1ols * cv.MDF.OLS2))
    )

    ur2.1olsgls.sa <- as.numeric(
      (MDF.GLS1 < (sap.cv.ur.2 * sap.ur2.1olsgls * cv.MDF.GLS1)) ||
        (MDF.OLS1 < (sap.cv.ur.2 * sap.ur2.1olsgls * cv.MDF.OLS1)) ||
        (MDF.GLS2 < (sap.cv.ur.2 * sap.ur2.1olsgls * cv.MDF.GLS2)) ||
        (MDF.OLS2 < (sap.cv.ur.2 * sap.ur2.1olsgls * cv.MDF.OLS2))
    )

    ur2.2ols.sa <- as.numeric(
      MDF.OLS2 < (sap.cv.ur.2 * sap.ur2.2ols * cv.MDF.OLS2)
    )

    ur2.2olsgls.sa <- as.numeric(
      (MDF.GLS2 < (sap.cv.ur.2 * sap.ur2.2olsgls * cv.MDF.GLS2)) ||
        (MDF.OLS2 < (sap.cv.ur.2 * sap.ur2.2olsgls * cv.MDF.OLS2))
    )

    UR <- as.numeric(breaks.star == 2) *
      ((1 - t.alpha.2.id) * ur2.olsgls.sa + t.alpha.2.id * ur2.ols.sa) +
      as.numeric(breaks.star == 1) *
        ((1 - t.alpha.2.id) * ur2.1olsgls.sa + t.alpha.2.id * ur2.1ols.sa) +
      as.numeric(breaks.star == 0) *
        ((1 - t.alpha.2.id) * ur2.2olsgls.sa + t.alpha.2.id * ur2.2ols.sa)

    ## without pre-test    for breaks
    ur2.ols.sa <- as.numeric(
      (DF.OLS.t < (sap.cv.ur.2 * sap.ur2.ols * cv.DF.OLS.t)) ||
        (MDF.OLS1 < (sap.cv.ur.2 * sap.ur2.ols * cv.MDF.OLS1)) ||
        (MDF.OLS2 < (sap.cv.ur.2 * sap.ur2.ols * cv.MDF.OLS2))
    )

    ur2.olsgls.sa <- as.numeric(
      (DF.GLS.t < (sap.cv.ur.2 * sap.ur2.olsgls * cv.DF.GLS.t)) ||
        (MDF.GLS1 < (sap.cv.ur.2 * sap.ur2.olsgls * cv.MDF.GLS1)) ||
        (DF.OLS.t < (sap.cv.ur.2 * sap.ur2.olsgls * cv.DF.OLS.t)) ||
        (MDF.OLS1 < (sap.cv.ur.2 * sap.ur2.olsgls * cv.MDF.OLS1)) ||
        (MDF.GLS2 < (sap.cv.ur.2 * sap.ur2.olsgls * cv.MDF.GLS2)) ||
        (MDF.OLS2 < (sap.cv.ur.2 * sap.ur2.olsgls * cv.MDF.OLS2))
    )

    UR1 <- (1 - t.alpha.2.id) * ur2.olsgls.sa + t.alpha.2.id * ur2.ols.sa
  } else if (breaks == 3) {
    ur3.ols.sa <- as.numeric(
      (DF.OLS.t < (sap.cv.ur.3 * sap.ur3.ols * cv.DF.OLS.t)) ||
        (MDF.OLS1 < (sap.cv.ur.3 * sap.ur3.ols * cv.MDF.OLS1)) ||
        (MDF.OLS2 < (sap.cv.ur.3 * sap.ur3.ols * cv.MDF.OLS2)) ||
        (MDF.OLS3 < (sap.cv.ur.3 * sap.ur3.ols * cv.MDF.OLS3))
    )

    ur3.olsgls.sa <- as.numeric(
      (DF.GLS.t < (sap.cv.ur.3 * sap.ur3.olsgls * cv.DF.GLS.t)) ||
        (MDF.GLS1 < (sap.cv.ur.3 * sap.ur3.olsgls * cv.MDF.GLS1)) ||
        (DF.OLS.t < (sap.cv.ur.3 * sap.ur3.olsgls * cv.DF.OLS.t)) ||
        (MDF.OLS1 < (sap.cv.ur.3 * sap.ur3.olsgls * cv.MDF.OLS1)) ||
        (MDF.GLS2 < (sap.cv.ur.3 * sap.ur3.olsgls * cv.MDF.GLS2)) ||
        (MDF.OLS2 < (sap.cv.ur.3 * sap.ur3.olsgls * cv.MDF.OLS2)) ||
        (MDF.GLS3 < (sap.cv.ur.3 * sap.ur3.olsgls * cv.MDF.GLS3)) ||
        (MDF.OLS3 < (sap.cv.ur.3 * sap.ur3.olsgls * cv.MDF.OLS3))
    )

    ur3.1ols.sa <- as.numeric(
      (MDF.OLS1 < (sap.cv.ur.3 * sap.ur3.1ols * cv.MDF.OLS1)) ||
        (MDF.OLS2 < (sap.cv.ur.3 * sap.ur3.1ols * cv.MDF.OLS2)) ||
        (MDF.OLS3 < (sap.cv.ur.3 * sap.ur3.1ols * cv.MDF.OLS3))
    )

    ur3.1olsgls.sa <- as.numeric(
      (MDF.GLS1 < (sap.cv.ur.3 * sap.ur3.1olsgls * cv.MDF.GLS1)) ||
        (MDF.OLS1 < (sap.cv.ur.3 * sap.ur3.1olsgls * cv.MDF.OLS1)) ||
        (MDF.GLS2 < (sap.cv.ur.3 * sap.ur3.1olsgls * cv.MDF.GLS2)) ||
        (MDF.OLS2 < (sap.cv.ur.3 * sap.ur3.1olsgls * cv.MDF.OLS2)) ||
        (MDF.GLS3 < (sap.cv.ur.3 * sap.ur3.1olsgls * cv.MDF.GLS3)) ||
        (MDF.OLS3 < (sap.cv.ur.3 * sap.ur3.1olsgls * cv.MDF.OLS3))
    )

    ur3.2ols.sa <- as.numeric(
      (MDF.OLS2 < (sap.cv.ur.3 * sap.ur3.2ols * cv.MDF.OLS2)) ||
        (MDF.OLS3 < (sap.cv.ur.3 * sap.ur3.2ols * cv.MDF.OLS3))
    )

    ur3.2olsgls.sa <- as.numeric(
      (MDF.GLS2 < (sap.cv.ur.3 * sap.ur3.2olsgls * cv.MDF.GLS2)) ||
        (MDF.OLS2 < (sap.cv.ur.3 * sap.ur3.2olsgls * cv.MDF.OLS2)) ||
        (MDF.GLS3 < (sap.cv.ur.3 * sap.ur3.2olsgls * cv.MDF.GLS3)) ||
        (MDF.OLS3 < (sap.cv.ur.3 * sap.ur3.2olsgls * cv.MDF.OLS3))
    )

    ur3.3ols.sa <- as.numeric(
      MDF.OLS3 < (sap.cv.ur.3 * sap.ur3.3ols * cv.MDF.OLS3)
    )

    ur3.3olsgls.sa <- as.numeric(
      (MDF.GLS3 < (sap.cv.ur.3 * sap.ur3.3olsgls * cv.MDF.GLS3)) ||
        (MDF.OLS3 < (sap.cv.ur.3 * sap.ur3.3olsgls * cv.MDF.OLS3))
    )

    UR <- as.numeric(breaks.star == 3) *
      ((1 - t.alpha.3.id) * ur3.olsgls.sa + t.alpha.3.id * ur3.ols.sa) +
      as.numeric(breaks.star == 2) *
        ((1 - t.alpha.3.id) * ur3.1olsgls.sa + t.alpha.3.id * ur3.1ols.sa) +
      as.numeric(breaks.star == 1) *
        ((1 - t.alpha.3.id) * ur3.2olsgls.sa + t.alpha.3.id * ur3.2ols.sa) +
      as.numeric(breaks.star == 0) *
        ((1 - t.alpha.3.id) * ur3.3olsgls.sa + t.alpha.3.id * ur3.3ols.sa)

    ## without pre-test for breaks

    ur3.ols.sa <- as.numeric(
      (DF.OLS.t < (sap.cv.ur.3 * sap.ur3.ols * cv.DF.OLS.t)) ||
        (MDF.OLS1 < (sap.cv.ur.3 * sap.ur3.ols * cv.MDF.OLS1)) ||
        (MDF.OLS2 < (sap.cv.ur.3 * sap.ur3.ols * cv.MDF.OLS2)) ||
        (MDF.OLS3 < (sap.cv.ur.3 * sap.ur3.ols * cv.MDF.OLS3))
    )

    ur3.olsgls.sa <- as.numeric(
      (DF.GLS.t < (sap.cv.ur.3 * sap.ur3.olsgls * cv.DF.GLS.t)) ||
        (MDF.GLS1 < (sap.cv.ur.3 * sap.ur3.olsgls * cv.MDF.GLS1)) ||
        (DF.OLS.t < (sap.cv.ur.3 * sap.ur3.olsgls * cv.DF.OLS.t)) ||
        (MDF.OLS1 < (sap.cv.ur.3 * sap.ur3.olsgls * cv.MDF.OLS1)) ||
        (MDF.GLS2 < (sap.cv.ur.3 * sap.ur3.olsgls * cv.MDF.GLS2)) ||
        (MDF.OLS2 < (sap.cv.ur.3 * sap.ur3.olsgls * cv.MDF.OLS2)) ||
        (MDF.GLS3 < (sap.cv.ur.3 * sap.ur3.olsgls * cv.MDF.GLS3)) ||
        (MDF.OLS3 < (sap.cv.ur.3 * sap.ur3.olsgls * cv.MDF.OLS3))
    )

    UR1 <- (1 - t.alpha.3.id) * ur3.olsgls.sa + t.alpha.3.id * ur3.ols.sa
  }

  result <- list()
  result$const <- const
  result$breaks.star <- breaks.star
  result$breaks.tbb <- tbb
  result$breaks <- breaks

  result$MDF.GLS.1 <- list(
    statistic = MDF.GLS1,
    cv = cv.MDF.GLS1
  )
  result$MDF.GLS.2 <- list(
    statistic = MDF.GLS2,
    cv = cv.MDF.GLS2
  )
  if (breaks == 3) {
    result$MDF.GLS.3 <- list(
      statistic = MDF.GLS3,
      cv = cv.MDF.GLS3
    )
  }

  result$MDF.OLS.1 <- list(
    statistic = MDF.OLS1,
    cv = cv.MDF.OLS1
  )
  result$MDF.OLS.2 <- list(
    statistic = MDF.OLS2,
    cv = cv.MDF.OLS2
  )
  if (breaks == 3) {
    result$MDF.OLS.3 <- list(
      statistic = MDF.OLS3,
      cv = cv.MDF.OLS3
    )
  }

  result$UR1 <- UR1
  result$UR <- UR

  class(result) <- "mdfHLTN"

  result
}

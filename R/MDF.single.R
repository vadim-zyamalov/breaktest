#' @title
#' MDF procedure for a single unknown break.
#'
#' @details
#' The code provided is the original Ox code by Skrobotov (2018)
#' ported to R.
#'
#' @param y A time series of interest.
#' @param const Whether the constant term should be included.
#' @param trend Whether the trend term should be included.
#' @param trim Trimming value for a possible break date bounds.
#'
#' @return A list of sublists each containing
#' * The value of statistic: \eqn{MDF-GLS}, \eqn{MDF-OLS},
#' * The asymptotic critical values.
#' \eqn{UR} values are included as well.
#'
#' @export
MDF.1br <- function(y,
                    const = FALSE, trend = FALSE,
                    trim = 0.15) {
  if (!is.matrix(y)) y <- as.matrix(y)
  model <- ifelse(const && trend, "ct", "t")

  result <- list(
    const = const,
    trend = trend
  )

  ## Magic numbers ##
  cv.DF.GLS.m <- .cval_MDF_single[[model]]$cv.DF.GLS.m
  cv.DF.OLS.m <- .cval_MDF_single[[model]]$cv.DF.OLS.m
  cv.DF.GLS.t <- .cval_MDF_single[[model]]$cv.DF.GLS.t
  cv.MDF.GLS <- .cval_MDF_single[[model]]$cv.MDF.GLS
  cv.DF.OLS.t <- .cval_MDF_single[[model]]$cv.DF.OLS.t
  cv.MDF.OLS <- .cval_MDF_single[[model]]$cv.MDF.OLS
  cv.MDF.t <- .cval_MDF_single[[model]]$cv.MDF.t
  cv.HLT <- .cval_MDF_single[[model]]$cv.HLT
  sap.ur <- .cval_MDF_single[[model]]$sap.ur
  sap.ur1 <- .cval_MDF_single[[model]]$sap.ur1
  sap.ur2 <- .cval_MDF_single[[model]]$sap.ur2
  sap.ur3 <- .cval_MDF_single[[model]]$sap.ur3
  sap.cv.ur.HLT.k0 <- .cval_MDF_single[[model]]$sap.cv.ur.HLT.k0
  sap.cv.A.HLT.k0 <- .cval_MDF_single[[model]]$sap.cv.A.HLT.k0
  sap.cv.ur.PY.k0 <- .cval_MDF_single[[model]]$sap.cv.ur.PY.k0
  sap.cv.A.PY.k0 <- .cval_MDF_single[[model]]$sap.cv.A.PY.k0

  ## Start ##
  N <- nrow(y)

  max.lag <- trunc(12 * (N / 100)^(1 / 4))

  first.break <- trunc(trim * N) + 1
  last.break <- trunc((1 - trim) * N) + 1

  tb <- segments.GLS(
    y, const, trend, 1,
    first.break, last.break,
    trim
  )
  tb <- drop(tb)
  result$break.time <- tb

  x <- cbind(
    .const(N),
    .trend(N),
    if (const) .du(tb, N) else NULL,
    if (trend) .dt(tb, N) else NULL
  )

  ## OLS/GLS Part ##
  ## Mean case
  resid.OLS.m <- .OLS(y, x[, 1, drop = FALSE])$residuals
  DF.OLS.m <- ADF.test(resid.OLS.m,
    const = FALSE, trend = FALSE,
    max.lag = max.lag,
    criterion = "aic",
    modified.criterion = TRUE
  )
  k.m <- max(1, DF.OLS.m$lag)

  DF.OLS.m <- ADF.test(resid.OLS.m,
    const = FALSE, trend = FALSE,
    max.lag = k.m,
    criterion = NULL
  )

  resid.GLS.m <- .GLS(y, x[, 1, drop = FALSE], -7)$residuals
  DF.GLS.m <- ADF.test(resid.GLS.m,
    const = FALSE, trend = FALSE,
    max.lag = k.m,
    criterion = NULL
  )

  ## Trend case
  resid.OLS.t <- .OLS(y, x[, 1:2])$residuals
  DF.OLS.t <- ADF.test(resid.OLS.t,
    const = FALSE, trend = FALSE,
    max.lag = max.lag,
    criterion = "aic",
    modified.criterion = TRUE
  )
  k.t <- max(1, DF.OLS.t$lag)

  DF.OLS.t <- ADF.test(resid.OLS.t,
    const = FALSE, trend = FALSE,
    max.lag = k.t,
    criterion = NULL
  )

  resid.GLS.t <- .GLS(y, x[, 1:2], -13.5)$residuals
  DF.GLS.t <- ADF.test(resid.GLS.t,
    const = FALSE, trend = FALSE,
    max.lag = k.t,
    criterion = NULL
  )

  ## ADF-OLS (lambda) ##
  resid.OLS <- .OLS(y, x)$residuals
  DF1 <- ADF.test(resid.OLS,
    const = FALSE, trend = FALSE,
    max.lag = max.lag,
    criterion = "aic",
    modified.criterion = TRUE
  )
  k.tb <- max(1, DF1$lag)
  DF1 <- ADF.test(resid.OLS,
    const = FALSE, trend = FALSE,
    max.lag = k.tb,
    criterion = NULL
  )
  MDF.OLS <- DF1$t.alpha

  ## MDF ##
  ## One break ##
  MDF.GLS <- Inf
  for (tb1 in first.break:last.break) {
    z <- cbind(
      .const(N),
      .trend(N),
      if (const) .du(tb1, N) else NULL,
      if (trend) .dt(tb1, N) else NULL
    )

    resid.OLS <- .OLS(y, z)$residuals
    DF1.tb <- ADF.test(resid.OLS,
      const = FALSE, trend = FALSE,
      max.lag = max.lag,
      criterion = "aic",
      modified.criterion = TRUE
    )
    k.tb <- max(1, DF1.tb$lag)

    resid.GLS <- .GLS(y, z, -17.6)$residuals
    DF1.tb <- ADF.test(resid.GLS,
      const = FALSE, trend = FALSE,
      max.lag = k.tb,
      criterion = NULL
    )

    if (DF1.tb$t.alpha < MDF.GLS) MDF.GLS <- DF1.tb$t.alpha
  }

  MDF.t <- Inf
  for (tb1 in first.break:last.break) {
    z <- cbind(
      .const(N),
      .trend(N),
      .dt(tb1, N)
    )

    resid.OLS <- .OLS(y, z)$residuals
    DF2 <- ADF.test(resid.OLS,
      const = FALSE, trend = FALSE,
      max.lag = max.lag,
      criterion = "aic",
      modified.criterion = TRUE
    )
    k.tb <- max(1, DF2$lag)
    DF2 <- ADF.test(resid.OLS,
      const = FALSE, trend = FALSE,
      max.lag = k.tb,
      criterion = NULL
    )
    DF2.tb <- N * DF2$alpha / (1 - sum(DF2$beta) + DF2$alpha)

    ######### TODO: Проверить!!!
    if (DF2$t.alpha < MDF.t) MDF.t <- DF2$t.alpha
  }

  t.HLT <- KPSS.HLT(y, const, trim)
  tmp.PY <- PY.single(y, const, trend, "aic", trim, max.lag)
  t.PY <- tmp.PY$wald
  cv.PY <- tmp.PY$critical.value
  rm(tmp.PY)

  tmp.OLS <- .OLS(y, x)
  t.alpha <- tmp.OLS$beta[1] /
    sqrt(drop(t(tmp.OLS$residuals) %*% tmp.OLS$residuals) / N)
  t.alpha.id <- as.numeric(abs(t.alpha) > 1)

  ## UR-HLT
  t.lambda.id <- as.numeric(t.HLT > cv.HLT)
  ur.id.sa <- as.numeric(
    (DF.GLS.t$t.alpha < (sap.ur * sap.cv.ur.HLT.k0 * cv.DF.OLS.t)) ||
      (MDF.GLS < (sap.ur * sap.cv.ur.HLT.k0 * cv.MDF.GLS)) ||
      (DF.OLS.t$t.alpha < (sap.ur * sap.cv.ur.HLT.k0 * cv.DF.OLS.t)) ||
      (MDF.t < (sap.ur * sap.cv.ur.HLT.k0 * cv.MDF.t))
  )
  ur1.id.sa <- as.numeric(
    (DF.OLS.t$t.alpha < (sap.ur1 * sap.cv.ur.HLT.k0 * cv.DF.OLS.t)) ||
      (MDF.t < (sap.ur1 * sap.cv.ur.HLT.k0 * cv.MDF.t))
  )
  ur2.id.sa <- as.numeric(
    (MDF.GLS < (sap.ur2 * sap.cv.ur.HLT.k0 * cv.MDF.GLS)) ||
      (MDF.OLS < (sap.ur2 * sap.cv.ur.HLT.k0 * cv.MDF.OLS))
  )
  ur3.id.sa <- as.numeric(MDF.OLS < sap.cv.ur.HLT.k0 * cv.MDF.OLS)
  UR.HLT <-
    (1 - t.lambda.id) * (1 - t.alpha.id) * ur.id.sa +
    (1 - t.lambda.id) * t.alpha.id * ur1.id.sa +
    t.lambda.id * (1 - t.alpha.id) * ur2.id.sa +
    t.lambda.id * t.alpha.id * ur3.id.sa

  ## UR-PY
  t.lambda.id <- as.numeric(t.PY > cv.PY[2])
  ur.id.sa <- as.numeric(
    (DF.GLS.t$t.alpha < (sap.ur * sap.cv.ur.PY.k0 * cv.DF.GLS.t)) ||
      (MDF.GLS < (sap.ur * sap.cv.ur.PY.k0 * cv.MDF.GLS)) ||
      (DF.OLS.t$t.alpha < (sap.ur * sap.cv.ur.PY.k0 * cv.DF.OLS.t)) ||
      (MDF.t < (sap.ur * sap.cv.ur.PY.k0 * cv.MDF.t))
  )
  ur1.id.sa <- as.numeric(
    (DF.OLS.t$t.alpha < (sap.ur1 * sap.cv.ur.PY.k0 * cv.DF.OLS.t)) ||
      (MDF.t < (sap.ur1 * sap.cv.ur.PY.k0 * cv.MDF.t))
  )
  ur2.id.sa <- as.numeric(
    (MDF.GLS < (sap.ur2 * sap.cv.ur.PY.k0 * cv.MDF.GLS)) ||
      (MDF.OLS < (sap.ur2 * sap.cv.ur.PY.k0 * cv.MDF.OLS))
  )
  ur3.id.sa <- as.numeric(MDF.OLS < sap.cv.ur.PY.k0 * cv.MDF.OLS)
  UR.PY <-
    (1 - t.lambda.id) * (1 - t.alpha.id) * ur.id.sa +
    (1 - t.lambda.id) * t.alpha.id * ur1.id.sa +
    t.lambda.id * (1 - t.alpha.id) * ur2.id.sa +
    t.lambda.id * t.alpha.id * ur3.id.sa

  ## A-HLT
  t.lambda.id <- as.numeric(t.HLT > cv.HLT)
  ur2.id.sa <- as.numeric(
    (MDF.GLS < (sap.ur2 * sap.cv.A.HLT.k0 * cv.MDF.GLS)) ||
      (MDF.OLS < (sap.ur2 * sap.cv.A.HLT.k0 * cv.MDF.OLS))
  )
  ur2a.id.sa <- as.numeric(MDF.OLS < sap.cv.A.HLT.k0 * cv.MDF.OLS)
  ur3.id.sa <- as.numeric(
    (MDF.GLS < (sap.ur3 * sap.cv.A.HLT.k0 * cv.MDF.GLS)) ||
      (MDF.t < (sap.ur3 * sap.cv.A.HLT.k0 * cv.MDF.t))
  )
  ur3a.id.sa <- as.numeric(MDF.t < sap.cv.A.HLT.k0 * cv.MDF.t)

  A.HLT <-
    (1 - t.lambda.id) * (1 - t.alpha.id) * ur3.id.sa +
    (1 - t.lambda.id) * t.alpha.id * ur3a.id.sa +
    t.lambda.id * (1 - t.alpha.id) * ur2.id.sa +
    t.lambda.id * t.alpha.id * ur2a.id.sa

  ## A-PY
  t.lambda.id <- as.numeric(t.PY > cv.PY[2])
  ur2.id.sa <- as.numeric(
    (MDF.GLS < (sap.ur2 * sap.cv.A.PY.k0 * cv.MDF.GLS)) ||
      (MDF.OLS < (sap.ur2 * sap.cv.A.PY.k0 * cv.MDF.OLS))
  )
  ur2a.id.sa <- as.numeric(MDF.OLS < sap.cv.A.PY.k0 * cv.MDF.OLS)
  ur3.id.sa <- as.numeric(
    (MDF.GLS < (sap.ur3 * sap.cv.A.PY.k0 * cv.MDF.GLS)) ||
      (MDF.t < (sap.ur3 * sap.cv.A.PY.k0 * cv.MDF.t))
  )
  ur3.id.sa <- as.numeric(MDF.t < sap.cv.A.PY.k0 * cv.MDF.t)

  A.PY <- (1 - t.lambda.id) * (1 - t.alpha.id) * ur3.id.sa +
    (1 - t.lambda.id) * t.alpha.id * ur3.id.sa +
    t.lambda.id * (1 - t.alpha.id) * ur2.id.sa +
    t.lambda.id * t.alpha.id * ur2a.id.sa

  result$HLT <- list(
    statistic = t.HLT,
    cv = cv.HLT
  )
  result$PY <- list(
    statistic = t.PY,
    cv = cv.PY[2]
  )
  result$DF.GLS <- list(
    statistic = DF.GLS.t$t.alpha,
    cv = cv.DF.GLS.t
  )
  result$DF.OLS <- list(
    statistic = DF.OLS.t$t.alpha,
    cv = cv.DF.OLS.t
  )
  result$MDF.GLS <- list(
    statistic = MDF.GLS,
    cv = cv.MDF.GLS
  )
  result$MDF.OLS <- list(
    statistic = MDF.OLS,
    cv = cv.MDF.OLS
  )
  result$MDF.t <- list(
    statistic = MDF.t,
    cv = cv.MDF.t
  )
  result$A.HLT <- A.HLT
  result$A.PY <- A.PY
  result$UR.HLT <- UR.HLT
  result$UR.PY <- UR.PY

  class(result) <- "mdfHLT"

  result
}

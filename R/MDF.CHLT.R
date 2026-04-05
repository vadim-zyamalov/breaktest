#' @title
#' MDF test for a single break and possible heteroscedasticity
#'
#' @details
#' The code provided is the original GAUSS code by Cavaliere et al.
#' ported to R.
#'
#' @param y A time series of interest.
#' @param max.lag The maximum possible lag.
#' @param trim Trimming parameter for lag selection
#' @param iter Number of bootstrap iterations.
#'
#' @return An object of type `mdfCHLT`. It's a list of four sublists
#' each containing:
#' * The value of \eqn{MZ_\alpha}, \eqn{MSB}, \eqn{MZ_t}, or \eqn{ADF},
#' * The asymptotic c.v.,
#' * The bootstrapped c.v.
#'
#' @references
#' Cavaliere, Giuseppe, David I. Harvey, Stephen J. Leybourne,
#' and A.M. Robert Taylor.
#' “Testing for Unit Roots in the Presence of a Possible Break in Trend and
#' Nonstationary Volatility.”
#' Econometric Theory 27, no. 5 (October 2011): 957–91.
#' https://doi.org/10.1017/S0266466610000605.
#'
#' @import doSNOW
#' @import foreach
#' @import parallel
#' @importFrom stats rnorm
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#'
#' @export
MDF.CHLT <- function(y,
                     max.lag = 10,
                     trim = 0.15,
                     iter = 499) {
  if (!is.matrix(y)) y <- as.matrix(y)

  ## CV ##
  g.brk.MZ <- .cval_MDF_CHLT$g.brk.MZ
  g.brk.ADF <- .cval_MDF_CHLT$g.brk.ADF
  cv.MZa.lim <- .cval_MDF_CHLT$cv.MZa.lim
  cv.MSB.lim <- .cval_MDF_CHLT$cv.MSB.lim
  cv.MZt.lim <- .cval_MDF_CHLT$cv.MZt.lim
  cv.ers.lim <- .cval_MDF_CHLT$cv.ers.lim
  MZa.cv.vals.lim <- .cval_MDF_CHLT$MZa.cv.vals.lim
  MSB.cv.vals.lim <- .cval_MDF_CHLT$MSB.cv.vals.lim
  MZt.cv.vals.lim <- .cval_MDF_CHLT$MZt.cv.vals.lim
  ADF.cv.vals.lim <- .cval_MDF_CHLT$ADF.cv.vals.lim
  cbar.vals <- .cval_MDF_CHLT$cbar.vals
  tau.cbar.MZa.cv.lim <- .cval_MDF_CHLT$tau.cbar.MZa.cv.lim
  tau.cbar.MSB.cv.lim <- .cval_MDF_CHLT$tau.cbar.MSB.cv.lim
  tau.cbar.MZt.cv.lim <- .cval_MDF_CHLT$tau.cbar.MZt.cv.lim
  tau.cbar.ADF.cv.lim <- .cval_MDF_CHLT$tau.cbar.ADF.cv.lim

  ## Start
  N <- nrow(y)

  first.break <- trunc(trim * N)
  last.break <- trunc((1 - trim) * N)

  d.y <- .diffn(y)

  ## Break date estimation
  res.ssr <- Inf
  tb.dy <- 0

  for (tb in first.break:last.break) {
    b1 <- (y[tb] - y[1]) / (tb - 1)
    b2 <- (y[N] - y[tb]) / (N - tb) - b1
    tmp.ssr <- (N - 1) * b1^2 + (N - tb) * b2^2 -
      2 * b1 * (y[N] - y[1]) - 2 * b2 * (y[N] - y[tb]) +
      2 * b1 * b2 * (N - tb)
    if (tmp.ssr < res.ssr) {
      res.ssr <- tmp.ssr
    }
    tb.dy <- tb
  }

  tau.dy <- tb.dy / N
  DT.tb.dy <- .dt(tb.dy, N)

  z <- cumsum(y)
  x <- cbind(
    .trend(N),
    cumsum(.trend(N)),
    cumsum(DT.tb.dy)
  )
  u.resid <- .OLS(z, x)$residuals
  x <- cbind(.trend(N), cumsum(.trend(N)))
  r.resid <- .OLS(z, x)$residuals
  W.stat.dy <- drop(t(r.resid) %*% r.resid) /
    drop(t(u.resid) %*% u.resid) - 1

  lam.MZ.brk.tau.dy <- exp(-g.brk.MZ * W.stat.dy / sqrt(N))
  tau.lam.MZ <- (1 - lam.MZ.brk.tau.dy) * tau.dy

  lam.ADF.brk.tau.dy <- exp(-g.brk.ADF * W.stat.dy / sqrt(N))
  tau.lam.ADF <- (1 - lam.ADF.brk.tau.dy) * tau.dy

  ## Unit root test
  x <- cbind(.const(N), .trend(N))

  resid.GLS <- .GLS(y, x, -13.5)$residuals
  resid.OLS <- .OLS(y, x)$residuals
  k.t <- ADF.test(
    resid.OLS,
    const = FALSE, trend = FALSE,
    max.lag = max.lag,
    criterion = "aic", modified.criterion = TRUE
  )$lag
  ers.DF <- ADF.test(
    resid.GLS,
    const = FALSE, trend = FALSE,
    max.lag = k.t,
    criterion = NULL
  )$t.alpha

  tmp.MZ <- .mz.statistics(resid.GLS, k.t)
  MZa <- tmp.MZ$mza
  MSB <- tmp.MZ$msb
  MZt <- tmp.MZ$mzt
  rm(tmp.MZ)

  ## MZ
  if (tau.lam.MZ < trim) {
    HHLT.MZa.kMAIC <- MZa
    HHLT.MZa.kMAIC.cv <- cv.MZa.lim
    HHLT.MSB.kMAIC <- MSB
    HHLT.MSB.kMAIC.cv <- cv.MSB.lim
    HHLT.MZt.kMAIC <- MZt
    HHLT.MZt.kMAIC.cv <- cv.MZt.lim
  } else {
    if (tau.lam.MZ == trim) {
      cbar.tau.lam.MZ <- tau.cbar.MZa.cv.lim[1, 2]
      cv.MZa.tau.lam.MZ.lim <- tau.cbar.MZa.cv.lim[1, 3]
      cv.MSB.tau.lam.MZ.lim <- tau.cbar.MSB.cv.lim[1, 3]
      cv.MZt.tau.lam.MZ.lim <- tau.cbar.MZt.cv.lim[1, 3]
    } else if (tau.lam.MZ == 1 - trim) {
      cbar.tau.lam.MZ <- tau.cbar.MZa.cv.lim[15, 2]
      cv.MZa.tau.lam.MZ.lim <- tau.cbar.MZa.cv.lim[15, 3]
      cv.MSB.tau.lam.MZ.lim <- tau.cbar.MSB.cv.lim[15, 3]
      cv.MZt.tau.lam.MZ.lim <- tau.cbar.MZt.cv.lim[15, 3]
    } else {
      tau.near.index <- which.min(
        abs(tau.lam.MZ - tau.cbar.MZa.cv.lim[, 1])
      )
      tau.near <- tau.cbar.MZa.cv.lim[tau.near.index, 1]

      if (tau.lam.MZ > tau.near) {
        tau.l <- tau.near
        tau.l.index <- tau.near.index
        tau.u <- tau.near + 0.05
        tau.u.index <- tau.near.index + 1
      } else {
        tau.l <- tau.near - 0.05
        tau.l.index <- tau.near.index - 1
        tau.u <- tau.near
        tau.u.index <- tau.near.index
      }

      weight.l <- 1 - (tau.lam.MZ - tau.l) / 0.05
      weight.u <- 1 - (tau.u - tau.lam.MZ) / 0.05

      cbar.tau.lam.MZ <-
        weight.l * tau.cbar.MZa.cv.lim[tau.l.index, 2] +
        weight.u * tau.cbar.MZa.cv.lim[tau.u.index, 2]
      cv.MZa.tau.lam.MZ.lim <-
        weight.l * tau.cbar.MZa.cv.lim[tau.l.index, 3] +
        weight.u * tau.cbar.MZa.cv.lim[tau.u.index, 3]
      cv.MSB.tau.lam.MZ.lim <-
        weight.l * tau.cbar.MSB.cv.lim[tau.l.index, 3] +
        weight.u * tau.cbar.MSB.cv.lim[tau.u.index, 3]
      cv.MZt.tau.lam.MZ.lim <-
        weight.l * tau.cbar.MZt.cv.lim[tau.l.index, 3] +
        weight.u * tau.cbar.MZt.cv.lim[tau.u.index, 3]
    }

    resid.GLS.bt <- GLS.bt(y, tau.lam.MZ, cbar.tau.lam.MZ)$residuals

    tb.lam.MZ <- trunc(tau.lam.MZ * N)
    resid.OLS.bt <- .OLS(
      y,
      cbind(.const(N), .trend(N), .dt(tb.lam.MZ, N))
    )$residuals

    k.bt <- ADF.test(
      resid.OLS.bt,
      const = FALSE, trend = FALSE,
      max.lag = max.lag,
      criterion = NULL
    )$lag

    tmp.MZ <- .mz.statistics(resid.GLS.bt, k.bt)
    MZa.tau.lam.MZ <- tmp.MZ$mza
    MSB.tau.lam.MZ <- tmp.MZ$msb
    MZt.tau.lam.MZ <- tmp.MZ$mzt
    rm(tmp.MZ)

    HHLT.MZa.kMAIC <- MZa.tau.lam.MZ
    HHLT.MZa.kMAIC.cv <- cv.MZa.tau.lam.MZ.lim
    HHLT.MSB.kMAIC <- MSB.tau.lam.MZ
    HHLT.MSB.kMAIC.cv <- cv.MSB.tau.lam.MZ.lim
    HHLT.MZt.kMAIC <- MZt.tau.lam.MZ
    HHLT.MZt.kMAIC.cv <- cv.MZt.tau.lam.MZ.lim
  }

  ## ADF
  if (tau.lam.ADF < trim) {
    HHLT.ADF.kMAIC <- ers.DF
    HHLT.ADF.kMAIC.cv <- cv.ers.lim
  } else {
    if (tau.lam.ADF == trim) {
      cbar.tau.lam.ADF <- tau.cbar.ADF.cv.lim[1, 2]
      cv.ADF.tau.lam.ADF.lim <- tau.cbar.ADF.cv.lim[1, 3]
    } else if (tau.lam.ADF == 1 - trim) {
      cbar.tau.lam.ADF <- tau.cbar.ADF.cv.lim[15, 2]
      cv.ADF.tau.lam.ADF.lim <- tau.cbar.ADF.cv.lim[15, 3]
    } else {
      tau.near.index <- which.min(
        abs(tau.lam.ADF - tau.cbar.ADF.cv.lim[, 1])
      )
      tau.near <- tau.cbar.ADF.cv.lim[tau.near.index, 1]

      if (tau.lam.ADF > tau.near) {
        tau.l <- tau.near
        tau.l.index <- tau.near.index
        tau.u <- tau.near + 0.05
        tau.u.index <- tau.near.index + 1
      } else {
        tau.l <- tau.near - 0.05
        tau.l.index <- tau.near.index - 1
        tau.u <- tau.near
        tau.u.index <- tau.near.index
      }

      weight.l <- 1 - (tau.lam.ADF - tau.l) / 0.05
      weight.u <- 1 - (tau.u - tau.lam.ADF) / 0.05

      cbar.tau.lam.ADF <-
        weight.l * tau.cbar.ADF.cv.lim[tau.l.index, 2] +
        weight.u * tau.cbar.ADF.cv.lim[tau.u.index, 2]
      cv.ADF.tau.lam.ADF.lim <-
        weight.l * tau.cbar.ADF.cv.lim[tau.l.index, 3] +
        weight.u * tau.cbar.ADF.cv.lim[tau.u.index, 3]
    }

    resid.GLS.bt <- GLS.bt(y, tau.lam.ADF, cbar.tau.lam.ADF)$residuals

    tb.lam.ADF <- trunc(tau.lam.ADF * N)
    resid.OLS.bt <- .OLS(
      y,
      cbind(.const(N), .trend(N), .dt(tb.lam.ADF, N))
    )$residuals

    k.bt <- ADF.test(
      resid.OLS.bt,
      const = FALSE, trend = FALSE,
      max.lag = max.lag,
      criterion = NULL
    )$lag

    ers.tau.lam.ADF.0 <- ADF.test(
      resid.GLS.bt,
      const = FALSE, trend = FALSE,
      max.lag = 0,
      criterion = NULL
    )$t.alpha

    ers.tau.lam.ADF.k <- ADF.test(
      resid.GLS.bt,
      const = FALSE, trend = FALSE,
      max.lag = k.bt,
      criterion = NULL
    )$t.alpha

    HHLT.ADF.kMAIC <- ers.tau.lam.ADF.k
    HHLT.ADF.kMAIC.cv <- cv.ADF.tau.lam.ADF.lim
  }

  ## Bootstrap
  eps <- .OLS(
    d.y,
    cbind(.const(N), .du(tb.dy, N))
  )$residuals
  eps <- as.matrix(eps)
  eps[1] <- 0

  cores <- detectCores()

  progress.bar <- txtProgressBar(max = iter, style = 3)
  progress <- function(n) setTxtProgressBar(progress.bar, n)

  cluster <- makeCluster(max(cores - 1, 1))
  registerDoSNOW(cluster)

  tmp.result <- foreach(
    i = 1:iter,
    .combine = rbind,
    .options.snow = list(progress = progress)
  ) %dopar% {
    z <- rnorm(N)
    y.wb <- cumsum(eps * z)


    if (tau.lam.MZ < trim) {
      resid.wb <- .GLS(
        y.wb,
        cbind(.const(N), .trend(N)),
        -13.5
      )$residuals
      MZ.wb <- .mz.statistics(resid.wb, 0)
      MZa.wb <- MZ.wb$mza
      MSB.wb <- MZ.wb$msb
      MZt.wb <- MZ.wb$mzt
      rm(MZ.wb)
    } else {
      resid.wb <- GLS.bt(y, tau.lam.MZ, cbar.tau.lam.MZ)$residuals
      MZ.wb <- .mz.statistics(resid.wb, 0)
      MZa.wb <- MZ.wb$mza
      MSB.wb <- MZ.wb$msb
      MZt.wb <- MZ.wb$mzt
      rm(MZ.wb)
    }

    if (tau.lam.ADF < trim) {
      resid.wb <- .GLS(
        y.wb,
        cbind(.const(N), .trend(N)),
        -13.5
      )$residuals
      ers.ADF.wb <- ADF.test(
        resid.wb,
        const = FALSE, trend = FALSE,
        max.lag = 0,
        criterion = NULL
      )$t.alpha
    } else {
      resid.wb <- GLS.bt(y, tau.lam.ADF, cbar.tau.lam.ADF)$residuals
      ers.ADF.wb <- ADF.test(
        resid.wb,
        const = FALSE, trend = FALSE,
        max.lag = 0,
        criterion = NULL
      )$t.alpha
    }

    c(MZa.wb, MSB.wb, MZt.wb, ers.ADF.wb)
  }

  stopCluster(cluster)

  s.stat <- sort(drop(tmp.result[, 1]))
  cv.MZa.0.wb <- s.stat[trunc(0.05 * iter)]
  s.stat <- sort(drop(tmp.result[, 2]))
  cv.MSB.0.wb <- s.stat[trunc(0.05 * iter)]
  s.stat <- sort(drop(tmp.result[, 3]))
  cv.MZt.0.wb <- s.stat[trunc(0.05 * iter)]
  s.stat <- sort(drop(tmp.result[, 4]))
  cv.ADF.0.wb <- s.stat[trunc(0.05 * iter)]

  result <- list(
    MZa = list(
      statistic = HHLT.MZa.kMAIC,
      cv = HHLT.MZa.kMAIC.cv,
      cv.bootstrap = cv.MZa.0.wb
    ),
    MSB = list(
      statistic = HHLT.MSB.kMAIC,
      cv = HHLT.MSB.kMAIC.cv,
      cv.bootstrap = cv.MSB.0.wb
    ),
    MZt = list(
      statistic = HHLT.MZt.kMAIC,
      cv = HHLT.MZt.kMAIC.cv,
      cv.bootstrap = cv.MZt.0.wb
    ),
    ADF = list(
      statistic = HHLT.ADF.kMAIC,
      cv = HHLT.ADF.kMAIC.cv,
      cv.bootstrap = cv.ADF.0.wb
    )
  )

  class(result) <- "mdfCHLT"

  result
}

#' @title
#' GLS fitering with a break
#'
#' @param y A time series of interest.
#' @param lambda A break relative position.
#' @param c A coefficient for \eqn{\rho} calculation.
#'
#' @keywords internal
GLS.bt <- function(y,
                   lambda,
                   c) {
  N <- nrow(y)
  tb <- trunc(lambda * N)
  x <- cbind(.const(N), .dt(tb, N))
  .GLS(y, x, c)
}

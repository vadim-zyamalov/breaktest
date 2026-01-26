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

  g.brk.MZ <- 6
  g.brk.ADF <- 3

  cv.MZa.lim <- -16.62
  cv.MSB.lim <- 0.171
  cv.MZt.lim <- -2.85
  cv.ers.lim <- -2.85

  MZa.cv.vals.lim <- c(
    -23.06, -23.49, -23.83, -23.95, -23.98,
    -23.95, -23.90, -23.77, -23.52, -23.25,
    -22.90, -22.46, -21.83, -21.16, -20.24
  )
  MSB.cv.vals.lim <- c(
    0.146, 0.144, 0.143, 0.143, 0.143,
    0.143, 0.143, 0.144, 0.145, 0.146,
    0.147, 0.148, 0.151, 0.153, 0.157
  )
  MZt.cv.vals.lim <- c(
    -3.37, -3.40, -3.42, -3.43, -3.43,
    -3.44, -3.44, -3.42, -3.41, -3.39,
    -3.37, -3.34, -3.29, -3.24, -3.17
  )
  ADF.cv.vals.lim <- MZt.cv.vals.lim

  cbar.vals <- c(
    -17.6, -17.8, -18.2, -18.4, -18.6,
    -18.4, -18.4, -18.2, -18.0, -17.6,
    -17.4, -17.0, -16.6, -16.0, -15.2
  )

  tau.cbar.MZa.cv.lim <- cbind(
    seq(from = 0.15, by = 0.05, length.out = 15),
    cbar.vals,
    MZa.cv.vals.lim
  )
  tau.cbar.MSB.cv.lim <- cbind(
    seq(from = 0.15, by = 0.05, length.out = 15),
    cbar.vals,
    MSB.cv.vals.lim
  )
  tau.cbar.MZt.cv.lim <- cbind(
    seq(from = 0.15, by = 0.05, length.out = 15),
    cbar.vals,
    MZt.cv.vals.lim
  )
  tau.cbar.ADF.cv.lim <- cbind(
    seq(from = 0.15, by = 0.05, length.out = 15),
    cbar.vals,
    ADF.cv.vals.lim
  )

  ## Start
  n.obs <- nrow(y)

  x.const <- rep(1, n.obs)
  x.trend <- 1:n.obs

  first.break <- trunc(trim * n.obs)
  last.break <- trunc((1 - trim) * n.obs)

  d.y <- as.matrix(diff(y))

  ## Break date estimation
  res.ssr <- Inf
  tb.dy <- 0

  for (tb in first.break:last.break) {
    b1 <- (y[tb] - y[1]) / (tb - 1)
    b2 <- (y[n.obs] - y[tb]) / (n.obs - tb) - b1
    tmp.ssr <- (n.obs - 1) * b1^2 + (n.obs - tb) * b2^2 -
      2 * b1 * (y[n.obs] - y[1]) - 2 * b2 * (y[n.obs] - y[tb]) +
      2 * b1 * b2 * (n.obs - tb)
    if (tmp.ssr < res.ssr) {
      res.ssr <- tmp.ssr
    }
    tb.dy <- tb
  }

  tau.dy <- tb.dy / n.obs

  DU.tb.dy <- as.numeric(x.trend > tb.dy)
  DT.tb.dy <- DU.tb.dy * (x.trend - tb.dy)

  z <- cumsum(y)
  x <- cbind(
    x.trend,
    cumsum(x.trend),
    cumsum(DT.tb.dy)
  )
  u.resid <- .estimate_ols(z, x)$residuals
  x <- cbind(
    x.trend,
    cumsum(x.trend)
  )
  r.resid <- .estimate_ols(z, x)$residuals
  W.stat.dy <- drop(t(r.resid) %*% r.resid) /
    drop(t(u.resid) %*% u.resid) - 1

  lam.MZ.brk.tau.dy <- exp(-g.brk.MZ * W.stat.dy / sqrt(n.obs))
  tau.lam.MZ <- (1 - lam.MZ.brk.tau.dy) * tau.dy

  lam.ADF.brk.tau.dy <- exp(-g.brk.ADF * W.stat.dy / sqrt(n.obs))
  tau.lam.ADF <- (1 - lam.ADF.brk.tau.dy) * tau.dy

  ## Unit root test
  x <- cbind(
    x.const,
    x.trend
  )

  resid.GLS <- .estimate_gls(y, x, -13.5)$residuals
  resid.OLS <- .estimate_ols(y, x)$residuals
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

  tmp.MZ <- .mz_stats(resid.GLS, k.t)
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

    tb.lam.MZ <- trunc(tau.lam.MZ * n.obs)
    DTr <- as.numeric(x.trend > tb.lam.MZ) *
      (x.trend - tb.lam.MZ)

    resid.OLS.bt <- .estimate_ols(
      y,
      cbind(
        x.const,
        x.trend,
        DTr
      )
    )$residuals

    k.bt <- ADF.test(
      resid.OLS.bt,
      const = FALSE, trend = FALSE,
      max.lag = max.lag,
      criterion = NULL
    )$lag

    tmp.MZ <- .mz_stats(resid.GLS.bt, k.bt)
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

    tb.lam.ADF <- trunc(tau.lam.ADF * n.obs)
    DTr <- as.numeric(x.trend > tb.lam.ADF) *
      (x.trend - tb.lam.ADF)

    resid.OLS.bt <- .estimate_ols(
      y,
      cbind(
        x.const,
        x.trend,
        DTr
      )
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
  eps <- .estimate_ols(
    d.y,
    cbind(
      x.const,
      DU.tb.dy
    )[2:n.obs, ]
  )$residuals
  eps <- as.matrix(c(0, eps))

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
    z <- rnorm(n.obs)
    y.wb <- cumsum(eps * z)


    if (tau.lam.MZ < trim) {
      resid.wb <- .estimate_gls(
        y.wb,
        cbind(
          x.const,
          x.trend
        ),
        -13.5
      )$residuals
      MZ.wb <- .mz_stats(resid.wb, 0)
      MZa.wb <- MZ.wb$mza
      MSB.wb <- MZ.wb$msb
      MZt.wb <- MZ.wb$mzt
      rm(MZ.wb)
    } else {
      resid.wb <- GLS.bt(y, tau.lam.MZ, cbar.tau.lam.MZ)$residuals
      MZ.wb <- .mz_stats(resid.wb, 0)
      MZa.wb <- MZ.wb$mza
      MSB.wb <- MZ.wb$msb
      MZt.wb <- MZ.wb$mzt
      rm(MZ.wb)
    }

    if (tau.lam.ADF < trim) {
      resid.wb <- .estimate_gls(
        y.wb,
        cbind(
          x.const,
          x.trend
        ),
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
      stat = HHLT.MZa.kMAIC,
      cv = HHLT.MZa.kMAIC.cv,
      cv.bootstrap = cv.MZa.0.wb
    ),
    MSB = list(
      stat = HHLT.MSB.kMAIC,
      cv = HHLT.MSB.kMAIC.cv,
      cv.bootstrap = cv.MSB.0.wb
    ),
    MZt = list(
      stat = HHLT.MZt.kMAIC,
      cv = HHLT.MZt.kMAIC.cv,
      cv.bootstrap = cv.MZt.0.wb
    ),
    ADF = list(
      stat = HHLT.ADF.kMAIC,
      cv = HHLT.ADF.kMAIC.cv,
      cv.bootstrap = cv.ADF.0.wb
    )
  )

  class(result) <- "mdfCHLT"

  return(result)
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
  n.obs <- nrow(y)

  tb <- trunc(lambda * n.obs)

  x <- cbind(
    rep(1, n.obs),
    c(
      rep(0, tb),
      1:(n.obs - tb)
    )
  )

  return(.estimate_gls(y, x, c))
}

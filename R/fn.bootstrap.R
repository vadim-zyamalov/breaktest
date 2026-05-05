#' @title
#' Functions to receive bootstrapped p-values.
#'
#' @param obj An object containing test results.
#' @param ... Any additional arguments for [bootstrap] function.
#'
#' @return A bootstrapped \eqn{p}-value.
#'
#' @keywords internal
#' @export
bootstrap <- function(obj, ...) UseMethod("bootstrap")


#' @rdname bootstrap
#'
#' @import doSNOW
#' @import foreach
#' @import parallel
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#'
#' @keywords internal
#' @export
bootstrap.bt_adf <- function(obj, iter = 999, ...) {
  vCoefs <- obj$model$coefficients[-1]
  vEps <- obj$model$residuals

  cLag <- obj$lag

  cN <- length(vEps)
  mDeter <- cbind(
    if (obj$const) .const(cN) else NULL,
    if (obj$trend) .trend(cN) else NULL
  )
  progress.bar <- txtProgressBar(max = iter, style = 3)
  progress <- function(n) setTxtProgressBar(progress.bar, n)
  cores <- detectCores()
  cluster <- makeCluster(max(cores - 1, 1), type = "SOCK")
  registerDoSNOW(cluster)
  tmp.stats <- foreach(
    i = 1:iter,
    .combine = c,
    .inorder = FALSE,
    .errorhandling = "remove",
    .packages = c("breaktest"),
    .options.snow = list(progress = progress)
  ) %dopar%
    {
      u <- rep(0, cLag + cN)
      eps <- sample(vEps, cN, replace = TRUE)
      if (cLag > 0) {
        for (s in 1:cN) {
          u[cLag + s] <- u[(cLag + s - 1):s] %*% vCoefs + eps[s]
        }
        u <- u[-(1:cLag)]
      } else {
        for (s in 1:cN) {
          u[s] <- eps[s]
        }
      }
      tmp.y <- as.matrix(cumsum(u))
      if (obj$recursive) {
        tmp.y <- detrend.recursively(
          tmp.y,
          mDeter,
          obj$recursive.params$cc,
          obj$recursive.params$gamma,
          obj$recursive.params$trim
        )
      }
      tmp.res <- OLS.reg(
        tmp.y,
        obj$model$exog
      )
      tmp.res$t.stats[1]
    }
  stopCluster(cluster)
  sum(tmp.stats < obj$t.alpha) / iter
}


#' @rdname bootstrap
#' @param obj An object of class `bt_kpss`.
#' @param iter Number of bootstrap iterations.
#' @param boot.type Type of bootstrapping:
#' * `"sample"`: sampling from residuals with replacement,
#' * `"Cavaliere-Taylor"`: multiplying residuals by \eqn{N(0, 1)}-distributed
#' variable,
#' * `"Rademacher"`: multiplying residuals by Rademacher-distributed variable.
#'
#' @import doSNOW
#' @import foreach
#' @import parallel
#' @importFrom stats rnorm
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#'
#' @keywords internal
#' @export
bootstrap.bt_kpss <- function(obj, iter = 999, boot.type = "sample", ...) {
  xreg <- obj$exog

  u <- obj$residuals

  max.lag <- obj$lr.var.max.lag

  kernel <- obj$lr.var.kernel

  cores <- detectCores()
  .progress <- txtProgressBar(max = iter, style = 3)
  progress <- function(n) setTxtProgressBar(.progress, n)
  cluster <- makeCluster(max(cores - 1, 1))
  registerDoSNOW(cluster)
  result <- foreach(
    i = 1:iter,
    .combine = c,
    .inorder = FALSE,
    .errorhandling = "remove",
    .packages = c("breaktest"),
    .options.snow = list(progress = progress)
  ) %dopar%
    {
      y.loop <- switch(
        boot.type,
        "sample" = sample(u, length(u), replace = TRUE),
        "Cavaliere-Taylor" = rnorm(length(u)) * u,
        "Rademacher" = sample(c(-1, 1), length(u), replace = TRUE) * u,
        stop(
          "ERROR! bootstsrap.bt_kpss: Unknown bootstrap type '",
          boot.type,
          "'"
        )
      )
      resids <- OLS.reg(y.loop, xreg)$residuals

      if (is.null(kernel)) {
        .kpss.statistic(resids, .lr.var.kurozumi(resids))
      } else {
        .kpss.statistic(resids, .lr.var.spc(resids, max.lag, kernel))
      }
    }
  stopCluster(cluster)
  (1 / iter) * sum(obj$statistic <= result)
}


#' @rdname bootstrap
#' @description
#' `SADF.bootstrap.test` is a wild bootstrapping procedure for estimating
#' critical and \eqn{p}-values for [uroot.SADF].
#'
#' `GSADF.bootstrap.test` is the same procedure but for `GSADF.test`.
#'
#' @details
#' Refactored original code by Kurozumi et al.
#'
#' @param y A time series of interest.
#' @param trim A trimming parameter to determine the lower and upper bounds for
#' a possible break point.
#' @param const Whether the constant needs to be included.
#' @param alpha The significance level of interest.
#' @param iter The number of iterations.
#' @param seed The seed parameter for the random number generator.
#'
#' @references
#' Kurozumi, Eiji, Anton Skrobotov, and Alexey Tsarev.
#' “Time-Transformed Test for Bubbles under Non-Stationary Volatility.”
#' Journal of Financial Econometrics, April 23, 2022.
#' https://doi.org/10.1093/jjfinec/nbac004.
#'
#' @import doSNOW
#' @import foreach
#' @import parallel
#' @importFrom stats quantile
#' @importFrom stats rnorm
#' @importFrom stats sd
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#'
#' @export
bootstrap.bt_SADF <- function(obj, iter = 999) {
  y <- obj$y

  trim <- obj$trim

  const <- obj$const

  N <- length(y)
  ## Find SADF.value.

  model <- uroot.SADF(y, trim, const)
  SADF.value <- model$SADF.value

  ## Do parallel.

  cores <- detectCores()
  progress.bar <- txtProgressBar(max = iter, style = 3)
  progress <- function(n) setTxtProgressBar(progress.bar, n)
  cluster <- makeCluster(max(cores - 1, 1))
  registerDoSNOW(cluster)
  SADF.bootstrap.values <- foreach(
    step = 1:iter,
    .combine = c,
    .options.snow = list(progress = progress)
  ) %dopar%
    {
      y.star <- cumsum(rnorm(N - 1) * .diffn(y, na = 0))
      uroot.SADF(y.star, trim, const)$SADF.value
    }
  stopCluster(cluster)
  sum(SADF.bootstrap.values > SADF.value) / iter
}


#' @rdname bootstrap
#'
#' @import doSNOW
#' @import foreach
#' @import parallel
#' @importFrom stats rnorm
#' @importFrom stats sd
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#'
#' @export
bootstrap.bt_GSADF <- function(obj, iter = 999) {
  y <- obj$y

  trim <- obj$trim

  const <- obj$const

  N <- length(y)
  ## Find GSADF.value.

  model <- uroot.GSADF(y, trim, const)
  GSADF.value <- model$GSADF.value

  ## Do parallel.

  cores <- detectCores()
  progress.bar <- txtProgressBar(max = iter, style = 3)
  progress <- function(n) setTxtProgressBar(progress.bar, n)
  cluster <- makeCluster(max(cores - 1, 1))
  clusterExport(cluster, c("GSADF.test", ".diffn"))
  registerDoSNOW(cluster)
  GSADF.bootstsrap.values <- foreach(
    step = 1:iter,
    .combine = c,
    .options.snow = list(progress = progress)
  ) %dopar%
    {
      y.star <- cumsum(rnorm(N - 1) * .diffn(y, na = 0))
      uroot.GSADF(y.star, trim, const)$GSADF.value
    }
  stopCluster(cluster)
  sum(GSADF.bootstsrap.values > GSADF.value) / iter
}


#' @rdname bootstrap
#'
#' @import doSNOW
#' @import foreach
#' @import parallel
#' @importFrom stats rnorm
#' @importFrom stats sd
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#'
#' @export
bootstrap.bt_mdfCHLT <- function(obj, iter = 999, y, ...) {
  N <- length(y)
  dy <- diff(y)
  trim <- obj$params$trim

  tb_dy <- obj$params$tb

  tau_lam_MZ <- obj$params$MZ$tau

  cbar_tau_lam_MZ <- obj$params$MZ$cbar

  tau_lam_ADF <- obj$params$ADF$tau

  cbar_tau_lam_ADF <- obj$params$ADF$cbar

  ## Bootstrap

  r <- c(
    0,
    OLS.reg(
      dy,
      cbind(.const(N), .du(tb_dy, N))[-1, ]
    )$residuals
  )
  cores <- detectCores()
  progress.bar <- txtProgressBar(max = iter, style = 3)
  progress <- function(n) setTxtProgressBar(progress.bar, n)
  cluster <- makeCluster(max(cores - 1, 1))
  registerDoSNOW(cluster)
  tmp.result <- foreach(
    i = 1:iter,
    .combine = rbind,
    .options.snow = list(progress = progress)
  ) %dopar%
    {
      z <- rnorm(N)
      y_wb <- cumsum(r * z)
      if (tau_lam_MZ < trim) {
        r_GLS_t_wb <- GLS.reg(
          y_wb,
          cbind(.const(N), .trend(N)),
          -13.5
        )$residuals

        MZ_wb <- .mz.statistics(r_GLS_t_wb, 0)
        MZa_wb <- MZ_wb$mza

        MSB_wb <- MZ_wb$msb

        MZt_wb <- MZ_wb$mzt

        rm(MZ_wb)
      } else {
        resid.wb <- GLS.bt(y, tau_lam_MZ, cbar_tau_lam_MZ)$residuals

        MZ_wb <- .mz.statistics(resid.wb, 0)
        MZa_wb <- MZ_wb$mza

        MSB_wb <- MZ_wb$msb

        MZt_wb <- MZ_wb$mzt

        rm(MZ_wb)
      }
      if (tau_lam_ADF < trim) {
        r_GLS_t_wb <- GLS.reg(
          y_wb,
          cbind(.const(N), .trend(N)),
          -13.5
        )$residuals

        ers_ADF_wb <- uroot.ADF(
          r_GLS_t_wb,
          const = FALSE,
          trend = FALSE,
          max.lag = 0,
          criterion = NULL
        )$t.alpha
      } else {
        r_GLS_wb <- GLS.bt(y, tau_lam_ADF, cbar_tau_lam_ADF)$residuals

        ers_ADF_wb <- uroot.ADF(
          r_GLS_wb,
          const = FALSE,
          trend = FALSE,
          max.lag = 0,
          criterion = NULL
        )$t.alpha
      }
      c(MZa_wb, MSB_wb, MZt_wb, ers_ADF_wb)
    }
  stopCluster(cluster)
  list(
    MZa = sort(tmp.result[, 1])[trunc(0.05 * iter)],
    MSB = sort(tmp.result[, 2])[trunc(0.05 * iter)],
    MZt = sort(tmp.result[, 3])[trunc(0.05 * iter)],
    ADF = sort(tmp.result[, 4])[trunc(0.05 * iter)]
  )
}

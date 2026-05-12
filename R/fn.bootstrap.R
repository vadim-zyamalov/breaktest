#' @title
#' Functions to receive bootstrapped p-values.
#'
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
#' @importFrom stats quantile
#' @importFrom stats rnorm
#' @importFrom stats sd
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#'
#' @keywords internal
#' @export
bootstrap <- function(obj, ...) UseMethod("bootstrap")


#' @rdname bootstrap
#' @keywords internal
#' @exportS3Method
bootstrap.bt_adf <- function(obj, iter = 999, ...) {
  vCoefs <- obj$model$coefficients[-1]
  vEps <- obj$model$residuals

  cLag <- obj$lag

  cN <- length(vEps)
  mDeter <- cbind(
    if (obj$const) .const(cN) else NULL,
    if (obj$trend) .trend(cN) else NULL
  )
  # progress.bar <- txtProgressBar(max = iter, style = 3)
  # progress <- function(n) setTxtProgressBar(progress.bar, n)
  # cores <- detectCores()
  # cluster <- makeCluster(max(cores - 1, 1), type = "SOCK")
  # registerDoSNOW(cluster)
  # tmp.stats <- foreach(
  #  i = 1:iter,
  #  .combine = c,
  #  .inorder = FALSE,
  #  .errorhandling = "remove",
  #  .packages = c("breaktest"),
  #  .options.snow = list(progress = progress)
  # ) %dopar%
  result <- NULL
  for (i in 1:iter) {
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
    result <- c(result, tmp.res$t.stats[1])
  }
  # stopCluster(cluster)
  sum(result < obj$t.alpha) / iter
}


#' @rdname bootstrap
#' @keywords internal
#' @exportS3Method
bootstrap.bt_kpss <- function(obj, iter = 999, boot.type = "sample", ...) {
  xreg <- obj$exog
  u <- obj$residuals
  max.lag <- obj$lr.var.max.lag
  kernel <- obj$lr.var.kernel

  result <- NULL
  for (i in 1:iter) {
    y.loop <- switch(boot.type,
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

    result <- c(
      result,
      if (kernel == "Kurozumi") {
        .kpss.statistic(resids, .lr.var.kurozumi(resids))
      } else {
        .kpss.statistic(resids, .lr.var.spc(resids, max.lag, kernel))
      }
    )
  }

  (1 / iter) * sum(obj$statistic <= result)
}


#' @rdname bootstrap
#' @keywords internal
#' @exportS3Method
bootstrap.bt_SADF <- function(obj, iter = 999) {
  y <- obj$y

  trim <- obj$trim

  const <- obj$const

  N <- length(y)
  ## Find SADF.value.

  model <- uroot.SADF(y, trim, const)
  SADF.value <- model$SADF.value

  ## Do parallel.

  # cores <- detectCores()
  # progress.bar <- txtProgressBar(max = iter, style = 3)
  # progress <- function(n) setTxtProgressBar(progress.bar, n)
  # cluster <- makeCluster(max(cores - 1, 1))
  # registerDoSNOW(cluster)
  # result <- foreach(
  #  step = 1:iter,
  #  .combine = c,
  #  .options.snow = list(progress = progress)
  # ) %dopar%
  result <- NULL
  for (i in 1:iter) {
    y.star <- cumsum(rnorm(N - 1) * .diffn(y, na = 0))
    result <- c(result, uroot.SADF(y.star, trim, const)$SADF.value)
  }
  # stopCluster(cluster)
  sum(result > SADF.value) / iter
}


#' @rdname bootstrap
#' @keywords internal
#' @exportS3Method
bootstrap.bt_GSADF <- function(obj, iter = 999) {
  y <- obj$y

  trim <- obj$trim

  const <- obj$const

  N <- length(y)
  ## Find GSADF.value.

  model <- uroot.GSADF(y, trim, const)
  GSADF.value <- model$GSADF.value

  # cores <- detectCores()
  # progress.bar <- txtProgressBar(max = iter, style = 3)
  # progress <- function(n) setTxtProgressBar(progress.bar, n)
  # cluster <- makeCluster(max(cores - 1, 1))
  # clusterExport(cluster, c("GSADF.test", ".diffn"))
  # registerDoSNOW(cluster)
  # result <- foreach(
  #  step = 1:iter,
  #  .combine = c,
  #  .options.snow = list(progress = progress)
  # ) %dopar%
  result <- NULL
  for (i in 1:iter) {
    y.star <- cumsum(rnorm(N - 1) * .diffn(y, na = 0))
    result <- c(result, uroot.GSADF(y.star, trim, const)$GSADF.value)
  }
  # stopCluster(cluster)
  sum(result > GSADF.value) / iter
}


#' @rdname bootstrap
#' @keywords internal
#' @exportS3Method
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
  # cores <- detectCores()
  # progress.bar <- txtProgressBar(max = iter, style = 3)
  # progress <- function(n) setTxtProgressBar(progress.bar, n)
  # cluster <- makeCluster(max(cores - 1, 1))
  # registerDoSNOW(cluster)
  # tmp.result <- foreach(
  #  i = 1:iter,
  #  .combine = rbind,
  #  .options.snow = list(progress = progress)
  # ) %dopar%
  result <- NULL
  for (i in 1:iter) {
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
    result <- rbind(result, c(MZa_wb, MSB_wb, MZt_wb, ers_ADF_wb))
  }
  # stopCluster(cluster)
  list(
    MZa = sort(result[, 1])[trunc(0.05 * iter)],
    MSB = sort(result[, 2])[trunc(0.05 * iter)],
    MZt = sort(result[, 3])[trunc(0.05 * iter)],
    ADF = sort(result[, 4])[trunc(0.05 * iter)]
  )
}

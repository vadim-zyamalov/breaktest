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
  ) %dopar% {
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
  ) %dopar% {
    y.loop <- switch(boot.type,
      "sample" = sample(u, length(u), replace = TRUE),
      "Cavaliere-Taylor" = rnorm(length(u)) * u,
      "Rademacher" = sample(c(-1, 1), length(u), replace = TRUE) * u,
      stop("ERROR! bootstsrap.bt_kpss: Unknown bootstrap type '", boot.type, "'")
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

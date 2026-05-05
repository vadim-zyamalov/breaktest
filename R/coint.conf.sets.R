#' @title
#' Confidence sets for the break date in cointegrating regressions
#'
#' @description
#' This procedure is to construct a confidence set for the change point in
#' cointegrating regressions.
#'
#' @param y A time series of interest.
#' @param trend Whether the trend is to be included.
#' @param zb I(1) regressors with break.
#' @param zf I(1) regressors without break.
#' @param nF,nL Number of leads and lags of `z` regressors.
#' If any is NULL then both are estimated using informational `criterion`.
#' @param conf.level Confidence level to obtain appropriate critical values.
#' @param trim The trimming parameter to find the lower and upper bounds of
#' possible break date.
#' @param criterion A criterion for lead and lag number estimation.
#'
#' @returns
#' A list of confidence sets.
#'
#' @references
#' Kurozumi, Eiji, and Anton Skrobotov.
#' “Confidence Sets for the Break Date in Cointegrating Regressions.”
#' Oxford Bulletin of Economics and Statistics 80, no. 3 (2018): 514–35.
#' https://doi.org/10.1111/obes.12223.
#'
#' @export
coint.conf.sets <- function(
  y,
  trend = FALSE,
  zb = NULL,
  zf = NULL,
  nF = NULL,
  nL = NULL,
  conf.level = 0.9,
  trim = 0.05,
  criterion = "bic"
) {
  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }
  if (!is.null(zb) && !is.matrix(zb)) {
    zb <- as.matrix(zb)
  }
  if (!is.null(zf) && !is.matrix(zf)) {
    zf <- as.matrix(zf)
  }

  if (is.null(nF) || is.null(nL)) {
    estLF <- segments.KS(y, trend, zb, zf, trim, criterion)
    nF <- estLF$lead
    nL <- estLF$lag
  }

  N <- nrow(y)
  td <- .trend(N)
  rows <- (nL + 2):(N - nF)

  wb <- cbind(
    .const(N),
    if (trend) .trend(N) / N else NULL,
    zb
  )
  z <- cbind(zb, zf)
  dz <- .diffn(z)

  Nzb <- if (!is.null(zb)) ncol(zb) else 0
  Nzf <- if (!is.null(zf)) ncol(zf) else 0

  wf <- cbind(zf, dz)
  for (f in seq_len(nF)) {
    wf <- cbind(wf, .lagn(dz, -f))
  }
  for (l in seq_len(nL)) {
    wf <- cbind(wf, .lagn(dz, l))
  }

  td <- td[rows]
  y <- y[rows, , drop = FALSE]
  wb <- wb[rows, , drop = FALSE]
  wf <- wf[rows, , drop = FALSE]

  N <- nrow(y)
  tb_L1 <- trunc(2 * trim * N)
  tb_U1 <- trunc((1 - 2 * trim) * N)
  tb_L2 <- trunc(trim * N)
  tb_U2 <- trunc((1 - trim) * N)

  csetSUP <- numeric(N)
  csetAVG <- numeric(N)
  csetEXP <- numeric(N)
  csetBLS <- numeric(N)

  w <- cbind(wb, wf)
  uhat <- OLS.reg(y, w)$residuals
  minSSR <- sum(uhat^2)
  Tb <- N

  for (tb in tb_L1:tb_U1) {
    wb1 <- rbind(
      matrix(0, tb, ncol(wb)),
      as.matrix(wb[(tb + 1):N, ])
    )
    w <- cbind(wb, wb1, wf)
    uhat <- OLS.reg(y, w)$residuals
    loopSSR <- sum(uhat^2)
    if (loopSSR < minSSR) {
      minSSR <- loopSSR
      Tb <- tb
    }
  }

  wb1e <- rbind(
    matrix(0, Tb, ncol(wb)),
    as.matrix(wb[(Tb + 1):N, ])
  )
  td[Tb] <- -1

  w <- cbind(wb1e, wb, wf)
  model <- OLS.reg(y, w)
  bhat <- model$coefficients
  uhat <- model$residuals

  lrvU <- .lr.var.kurozumi(uhat)

  l.hat <- (wb[Tb, ] %*% bhat[seq_len(ncol(wb))])^2 / lrvU
  c.bls <- if (conf.level == 0.9) {
    7.686962
  } else if (conf.level == 0.95) {
    11.03281
  }

  bdd <- trunc(c.bls / l.hat)
  blsL <- max(1, Tb - bdd - 1)
  blsU <- min(Tb + bdd + 1, N)

  csetBLS[blsL:blsU] <- 1

  for (tb in tb_L1:tb_U1) {
    lmb1 <- tb / N

    wb1 <- rbind(
      matrix(0, tb, ncol(wb)),
      as.matrix(wb[(tb + 1):N, ])
    )

    w <- cbind(wb, wb1, wf)

    yhat <- OLS.reg(y, w)$residuals

    we <- cbind(w, if (abs(tb - Tb) > ncol(wb)) wb1e else NULL)
    lrvU2 <- .lr.var.kurozumi(OLS.reg(y, we)$residuals)

    SUPstat <- 0
    AVGstat <- 0
    EXPstat <- 0

    nbreak <- 0
    dbreak <- 0

    for (tb2 in tb_L2:tb_U2) {
      lmb2 <- tb2 / N

      if (abs(lmb2 - lmb1) <= 0.05) {
        dbreak <- dbreak + 1
      } else {
        wb2 <- rbind(
          matrix(0, tb2, ncol(wb)),
          as.matrix(wb[(tb2 + 1):N, ])
        )

        r <- (wb2 - wb1)
        rhat <- r - w %*% solve(crossprod(w), crossprod(w, r))

        g <- crossprod(rhat, as.matrix(yhat))
        h <- crossprod(rhat)
        ghg <- c(t(g) %*% spdinv(h) %*% g) / lrvU2

        SUPstat <- min(SUPstat, ghg)
        AVGstat <- AVGstat + ghg
        EXPstat <- EXPstat + exp(ghg / 2)

        nbreak <- nbreak + 1
      }
    }

    AVGstat <- AVGstat / (nbreak - dbreak)
    EXPstat <- log(EXPstat / (nbreak - dbreak))

    cv <- get.cv.coint.conf.sets(
      lmb1,
      trend,
      conf.level,
      Nzb,
      Nzf
    )
    if (SUPstat <= cv$cval_sup) {
      csetSUP[tb] <- 1
    }
    if (AVGstat <= cv$cval_avg) {
      csetAVG[tb] <- 1
    }
    if (EXPstat <= cv$cval_exp) {
      csetEXP[tb] <- 1
    }
  }

  result <- list(
    td = c(rep(0, nL + 1), td, rep(0, nF)),
    cset.sup = c(rep(0, nL + 1), csetSUP, rep(0, nF)),
    cset.avg = c(rep(0, nL + 1), csetAVG, rep(0, nF)),
    cset.exp = c(rep(0, nL + 1), csetEXP, rep(0, nF)),
    cset.bls = c(rep(0, nL + 1), csetBLS, rep(0, nF))
  )
  class(result) <- "bt_confSet"

  result
}


#' @title
#' Estimating optimal number of leads and lags
#'
#' @details
#' The function is not intended to be used directly so it's not exported.
#'
#' @param y A time series of interest.
#' @param trend Whether the trend is to be included.
#' @param zb I(1) regressors with break.
#' @param zf I(1) regressors without break.
#' @param trim The trimming parameter to find the lower and upper bounds of
#' possible break date.
#' @param criterion A criterion for lead and lag number estimation.
#'
#' @return A list of estimated values of leads and lags.
#'
#' @references
#' Kurozumi, Eiji, and Anton Skrobotov.
#' “Confidence Sets for the Break Date in Cointegrating Regressions.”
#' Oxford Bulletin of Economics and Statistics 80, no. 3 (2018): 514–35.
#' https://doi.org/10.1111/obes.12223.
#'
#' @keywords internal
segments.KS <- function(
  y,
  trend = TRUE,
  zb = NULL,
  zf = NULL,
  trim = 0.05,
  criterion = "bic"
) {
  if (!criterion %in% c("bic", "aic", "hq", "lwz")) {
    stop("ERROR! Unknown criterion")
  }

  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }
  if (!is.null(zb) && !is.matrix(zb)) {
    zb <- as.matrix(zb)
  }
  if (!is.null(zf) && !is.matrix(zf)) {
    zf <- as.matrix(zf)
  }

  N <- nrow(y)

  tb_L <- trunc(2 * trim * N)
  tb_U <- trunc((1 - 2 * trim) * N)

  wb <- cbind(
    .const(N),
    if (trend) .trend(N) / N else NULL,
    zb
  )

  z <- cbind(zb, zf)
  w <- cbind(wb, zf)
  uhat <- OLS.reg(y, w)$residuals
  minSSR <- sum(uhat^2)
  Tb <- N
  wbb <- NULL

  for (tb in tb_L:tb_U) {
    loopWb <- rbind(matrix(0, tb, ncol(wb)), wb[(tb + 1):N, , drop = FALSE])
    w <- cbind(wb, loopWb, zf)
    uhat <- OLS.reg(y, w)$residuals
    loopSSR <- sum(uhat^2)

    if (loopSSR < minSSR) {
      minSSR <- loopSSR
      Tb <- tb
      wbb <- loopWb
    }
  }

  dz <- .diffn(z)
  wf <- cbind(zf, dz)

  y <- y[2:N, , drop = FALSE]
  wb <- wb[2:N, , drop = FALSE]
  wf <- wf[2:N, , drop = FALSE]
  wbb <- wbb[2:N, , drop = FALSE]
  dz <- dz[2:N, , drop = FALSE]

  N <- nrow(y)

  maxLF <- trunc(4 * (N / 100)^(1 / 4))
  if ((Tb - maxLF - 1) <= ncol(wb)) {
    maxLF <- Tb - ncol(wb) - 2
  } else if ((N - maxLF - Tb) <= ncol(wb)) {
    maxLF <- N - Tb - ncol(wb) - 1
  }

  y0 <- y[(maxLF + 1):(N - maxLF), ]
  w0 <- cbind(
    wb[(maxLF + 1):(N - maxLF), ],
    wbb[(maxLF + 1):(N - maxLF), ],
    wf[(maxLF + 1):(N - maxLF), ]
  )

  uhat <- OLS.reg(y0, w0)$residuals
  minIC <- info.criterions(uhat, ncol(w0))[[criterion]]
  estL <- 0
  estF <- 0

  for (loopL in 1:maxLF) {
    for (loopF in 1:maxLF) {
      loopW <- w0
      for (k in 1:loopL) {
        loopW <- cbind(
          loopW,
          as.matrix(
            dz[(maxLF + 1 - k):(N - maxLF - k), ]
          )
        )
      }
      for (k in 1:loopF) {
        loopW <- cbind(
          loopW,
          as.matrix(
            dz[(maxLF + 1 + k):(N - maxLF + k), ]
          )
        )
      }

      uhat <- OLS.reg(y0, loopW)$residuals

      loopIC <- info.criterions(uhat, ncol(loopW))[[criterion]]
      if (loopIC < minIC) {
        estL <- loopL
        estF <- loopF
        minIC <- loopIC
      }
    }
  }

  list(
    bp = Tb,
    lead = estF,
    lag = estL
  )
}

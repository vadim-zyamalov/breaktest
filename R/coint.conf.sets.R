#' @title
#' Confidence sets for the break date in cointegrating regressions
#'
#' @order 1
#'
#' @description
#' Procedure to construct a confidence set for the change point in
#' cointegrating regressions.
#'
#' @details
#' The function constructs confidence sets for a break date in cointegrating
#' regressions by inverting a test for the break location, which is obtained by maximizing the
#' weighted average of power. It is found in Kurozumi and Skrobotov (2018)
#' that the limiting distribution of the test depends on the number of I(1) regressors
#' whose coefficients sustain structural change and the number of I(1) regressors
#' whose coefficients are fixed throughout the sample.
#' By Monte Carlo simulations, the authors of the original paper showed
#' that compared with a confidence interval developed by using the existing method
#' based on the limiting distribution of the break point estimator under
#' the assumption of the shrinking shift, the confidence set proposed in the present paper has
#' a more accurate coverage rate, while the length of the confidence set is comparable.
#'
#' The model underlying the confidence set construction is
#' \deqn{
#' y_t = \beta_{b,c} + \beta_{b,\tau} t + z_{b,t}'\beta_{b,z}
#'   + \mathrm{1}(t > [\lambda_0 T])\left(\delta_{b,c} + \delta_{b,\tau} t + z_{b,t}'\delta_{b,z}\right)
#'   + z_{f,t}'\beta_{f,z} + u_t,
#' }
#' where \eqn{z_{b,t}} (`zb`) is the vector of \eqn{I(1)} regressors whose coefficients
#' sustain a structural break at the unknown date \eqn{T_0 = [\lambda_0 T]}, and
#' \eqn{z_{f,t}} (`zf`) is the vector of \eqn{I(1)} regressors whose coefficients are
#' fixed throughout the sample. The constant and, if `trend = TRUE`, the linear trend
#' are allowed to break jointly with \eqn{z_{b,t}}, corresponding to Model II-c of
#' Kurozumi and Skrobotov (2018); setting `trend = FALSE` reduces the specification
#' to the corresponding Model I-c.
#'
#' Since \eqn{z_{b,t}} and \eqn{z_{f,t}} are \eqn{I(1)} and endogenous, the regression
#' is augmented with leads and lags of their first differences (dynamic OLS), so that
#' the estimated equation becomes
#' \deqn{\begin{aligned}
#' y_t & = \beta_{b,c} + \beta_{b,\tau} t + z_{b,t}'\beta_{b,z}
#'   + \mathrm{1}(t > [\lambda_0 T])\left(\delta_{b,c} + \delta_{b,\tau} t + z_{b,t}'\delta_{b,z}\right)
#'   + z_{f,t}'\beta_{f,z} + {}\\
#'   & \quad + \sum_{j=-nL}^{nF} \pi_{b,j}'\Delta z_{b,t-j}
#'   + \sum_{j=-nL}^{nF} \pi_{f,j}'\Delta z_{f,t-j} + u_t,
#' \end{aligned}
#' }
#' where `nF` and `nL` are the numbers of leads and lags, respectively. If either is
#' `NULL`, both are selected by the information `criterion`.
#'
#' The confidence set for the break date \eqn{T_0} is obtained by inverting a test for
#' the break location that maximizes the weighted average power over the magnitude and
#' location of the break, at the confidence level `conf.level`. The lower and upper
#' bounds of admissible break dates are determined by the trimming parameter `trim`.
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
#' @return Named list of five \eqn{T_y \times 1} numeric matrices:
#'   \item{td}{Time index; \code{td[i] == -1} marks the estimated break date
#'     (all other entries equal their own time index, or 0 in the padded
#'     leading/trailing region lost to leads/lags/differencing).}
#'   \item{cset_sup}{Confidence set from the sup-test (1 = point included).}
#'   \item{cset_avg}{Confidence set from the avg-test.}
#'   \item{cset_exp}{Confidence set from the exp-test.}
#'   \item{cset_bls}{Confidence set from the Bai, Lumsdaine and Stock (1998) method.}
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
  y <- .msub(y, rows)
  wb <- .msub(wb, rows)
  wf <- .msub(wf, rows)

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

  # lrvU <- .lr.var.kurozumi(uhat)
  lrvU <- .lr.var.quad(uhat, 0.97)

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
    # lrvU2 <- .lr.var.kurozumi(OLS.reg(y, we)$residuals)
    lrvU2 <- .lr.var.quad(OLS.reg(y, we)$residuals, 0.97)

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
        ghg <- as.numeric(t(g) %*% spdinv(h) %*% g) / lrvU2

        SUPstat <- min(SUPstat, ghg)
        AVGstat <- AVGstat + ghg
        EXPstat <- EXPstat + exp(ghg / 2)
      }
      nbreak <- nbreak + 1
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

  for (tb in tb_L:tb_U) {
    loopWb <- rbind(matrix(0, tb, ncol(wb)), .msub(wb, (tb + 1):N))
    w <- cbind(wb, loopWb, zf)
    uhat <- OLS.reg(y, w)$residuals
    loopSSR <- sum(uhat^2)

    if (loopSSR < minSSR) {
      minSSR <- loopSSR
      Tb <- tb
    }
  }

  cat("SSR: ", minSSR, Tb, "\n")

  wbb <- rbind(matrix(0, Tb, ncol(wb)), .msub(wb, (Tb + 1):N))

  dz <- .diffn(z)
  wf <- cbind(zf, dz)

  y <- .msub(y, 2:N)
  wb <- .msub(wb, 2:N)
  wf <- .msub(wf, 2:N)
  wbb <- .msub(wbb, 2:N)
  dz <- .msub(dz, 2:N)

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

  minIC <- info.criterions(OLS.reg(y0, w0), criterion)
  estL <- 0
  estF <- 0

  for (loopL in .seqi(1, maxLF)) {
    for (loopF in .seqi(1, maxLF)) {
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

      loopIC <- info.criterions(OLS.reg(y0, loopW), criterion)
      if (loopIC < minIC) {
        estL <- loopL
        estF <- loopF
        minIC <- loopIC
      }
    }
  }

  cat("z_ld =", estF, " z_lg =", estL, "\n")

  list(
    bp = Tb,
    lead = estF,
    lag = estL
  )
}

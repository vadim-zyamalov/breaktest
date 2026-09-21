#' @rdname coint.conf.sets
#' @order 2
#'
#' @param centering Boolean. Whether the data should be detrended during the procedure.
#'
#' @export
cset_break_coint <- function(
  y,
  trend,
  zb = NULL,
  zf = NULL,
  nF = NULL,
  nL = NULL,
  conf.level = c(0.90, 0.95),
  trim = 0.05,
  centering = FALSE,
  criterion = "bic"
) {
  c_bls <- if (conf.level == 0.90) {
    7.686962
  } else if (conf.level == 0.95) {
    11.03281
  } else {
    stop("c_level must be 0.90 or 0.95")
  }

  y <- as.matrix(y)
  zb <- if (!is.null(zb)) as.matrix(zb)
  zf <- if (!is.null(zf)) as.matrix(zf)

  if (is.null(nF) || is.null(nL)) {
    ll <- .select.lead.lag(y, trend, zb, zf)
    nF <- ll$leads
    nL <- ll$lags
  }

  N <- nrow(y)
  td <- matrix(1:N, ncol = 1)

  wb <- cbind(
    .ones(N, 1),
    if (trend) .trend(N) / N else NULL,
    zb
  )

  z <- cbind(zb, zf)

  pzb <- if (is.null(zb)) 0 else ncol(zb)
  pzf <- if (is.null(zf)) 0 else ncol(zf)

  td <- .msub(td, .seqi(2, N))
  y <- .msub(y, .seqi(2, N))
  wb <- .msub(wb, .seqi(2, N))
  xf <- .msub(z, .seqi(2, N)) - .msub(z, .seqi(1, N - 1))

  wf <- cbind(
    .msub(zf, .seqi(2, N)),
    xf
  )
  N1 <- nrow(y)

  rng <- .seqi(nL + 1, N1 - nF)
  td <- .msub(td, rng)
  y <- .msub(y, rng)
  wb <- .msub(wb, rng)
  wf <- .msub(wf, rng)

  for (k1 in .seqi(1, nF)) {
    wf <- cbind(
      wf,
      .msub(xf, .seqi(nL + 1 + k1, N1 - nF + k1))
    )
  }
  for (k2 in .seqi(1, nL)) {
    wf <- cbind(
      wf,
      .msub(xf, .seqi(nL + 1 - k2, N1 - nF - k2))
    )
  }

  N2 <- nrow(y)

  st_t1 <- trunc(2 * trim * N2)
  ed_t1 <- trunc((1 - 2 * trim) * N2)
  st_t2 <- trunc(trim * N2)
  ed_t2 <- trunc((1 - trim) * N2)

  cset_sup <- .zeros(N2, 1)
  cset_avg <- .zeros(N2, 1)
  cset_exp <- .zeros(N2, 1)
  cset_bls <- .zeros(N2, 1)

  w <- cbind(wb, wf)
  ssr0 <- as.numeric(crossprod(OLS.reg(y, w)$residuals))
  est_date <- N2

  for (T1 in .seqi(st_t1, ed_t1)) {
    wb1 <- rbind(
      .zeros(T1, ncol(wb)),
      .msub(wb, .seqi(T1 + 1, N2))
    )
    w <- cbind(wb, wb1, wf)
    ssr1 <- as.numeric(crossprod(OLS.reg(y, w)$residuals))
    if (ssr1 < ssr0) {
      ssr0 <- ssr1
      est_date <- T1
    }
  }

  wb1e <- rbind(
    .zeros(est_date, ncol(wb)),
    .msub(wb, .seqi(est_date + 1, N2))
  )
  td[est_date, 1] <- -1

  if (centering) {
    if (trend) {
      level <- .ones(nrow(xf), 1)
      level <- rbind(
        .zeros(est_date + nL + 1, 1),
        .msub(level, .seqi(est_date + nL + 1 + 1, N - 1))
      )
      xf <- OLS.reg(xf, level)$residuals
    }

    wf <- cbind(
      .msub(zf, .seqi(2, N)),
      xf
    )
    wf <- .msub(wf, .seqi(nL + 1, N1 - nF))

    for (k1 in .seqi(1, nF)) {
      wf <- cbind(wf, .msub(xf, .seqi(nL + 1 + k1, N1 - nF + k1)))
    }
    for (k2 in .seqi(1, nL)) {
      wf <- cbind(wf, .msub(xf, .seqi(nL + 1 - k2, N1 - nF - k2)))
    }
    N2 <- nrow(y)
  }

  w <- cbind(wb1e, wb, wf)
  b_hat <- solve(crossprod(w), crossprod(w, y))
  u_hat <- y - w %*% b_hat

  lv <- .lr.matr.quad(u_hat, 0.97)
  lrv_u <- as.numeric(lv$Omega)

  L_hat <- as.numeric(
    (.msub(wb, est_date) %*% .msub(b_hat, seq_len(ncol(wb))))^2 / lrv_u
  )

  bdd <- trunc(c_bls / L_hat)
  bls_l <- max(est_date - bdd - 1, 1)
  bls_u <- min(est_date + bdd + 1, N2)
  cset_bls[bls_l:bls_u, 1] <- 1

  for (T1 in .seqi(st_t1, ed_t1)) {
    L1 <- T1 / N2
    wb1 <- rbind(.zeros(T1, ncol(wb)), .msub(wb, .seqi(T1 + 1, N2)))
    w <- cbind(wb, wb1, wf)

    b_hat <- solve(crossprod(w), crossprod(w, y))
    y_hat <- y - w %*% b_hat

    if (abs(T1 - est_date) > ncol(wb)) {
      we <- cbind(w, wb1e)
    } else {
      we <- w
    }
    be_hat <- solve(crossprod(we), crossprod(we, y))
    u_hat <- y - we %*% be_hat

    lv <- .lr.matr.quad(u_hat, 0.97)
    lrv_u <- as.numeric(lv$Omega)

    sup_stat <- 0
    avg_stat <- 0
    exp_stat <- 0
    nbreak <- 0
    dbreak <- 0

    for (T2 in .seqi(st_t2, ed_t2)) {
      L2 <- T2 / N2
      if (abs(L2 - L1) <= 0.05) {
        dbreak <- dbreak + 1
      } else {
        wb2 <- rbind(.zeros(T2, ncol(wb)), .msub(wb, .seqi(T2 + 1, N2)))
        r <- wb2 - wb1
        r_hat <- r - w %*% solve(crossprod(w), crossprod(w, r))

        g <- crossprod(r_hat, y_hat)
        h <- crossprod(r_hat)
        ghg <- as.numeric(t(g) %*% solve(h) %*% g / lrv_u)

        sup_stat <- min(sup_stat, ghg)
        avg_stat <- avg_stat + ghg
        exp_stat <- exp_stat + exp(ghg / 2)
      }
      nbreak <- nbreak + 1
    }

    sup_stat <- sup_stat
    avg_stat <- avg_stat / (nbreak - dbreak)
    exp_stat <- log(exp_stat / (nbreak - dbreak))

    cv <- get.cv.coint.conf.sets(L1, trend, conf.level, pzb, pzf)

    if (sup_stat <= cv$cval_sup) {
      cset_sup[T1, 1] <- 1
    }
    if (avg_stat <= cv$cval_avg) {
      cset_avg[T1, 1] <- 1
    }
    if (exp_stat <= cv$cval_exp) cset_exp[T1, 1] <- 1
  }

  pad0 <- .zeros(nL + 1, 1)
  pad1 <- .zeros(nF, 1)
  result <- list(
    td = rbind(pad0, td, pad1),
    cset_sup = rbind(pad0, cset_sup, pad1),
    cset_avg = rbind(pad0, cset_avg, pad1),
    cset_exp = rbind(pad0, cset_exp, pad1),
    cset_bls = rbind(pad0, cset_bls, pad1)
  )
  class(result) <- "bt_confSet"

  result
}


#' Select the lead/lag length for the first-differenced I(1) regressors by BIC.
#'
#' @param y         \eqn{T_y \times 1} numeric matrix or vector, the dependent time series.
#' @param trend     Boolean. FALSE = only a constant is included;
#'                  TRUE = a constant and a linear trend are included.
#' @param zb        \eqn{T_y \times pzb} matrix of I(1) regressors whose coefficient
#'                  sustains a structural change; \code{NULL} if there are none.
#' @param zf        \eqn{T_y \times pzf} matrix of I(1) regressors whose coefficient is
#'                  fixed (does not change at the break); \code{NULL} if there
#'                  are none.
#' @param trim      Trimming share defining lower and upper bounds of possible break moments.
#' @param maxLF     Maximum number of leads and lags. If NULL then
#'                  \eqn{\left[4 \left(\frac{T}{100}\right)^{\frac14}\right]} is used.
#' @param criterion Information criterion to select the number of leads and lags.
#'
#' @return Named list:
#'   \item{est_ldg1}{Selected lead length.}
#'   \item{est_ldg2}{Selected lag length.}
#'
#' @keywords internal
.select.lead.lag <- function(
  y,
  trend,
  zb = NULL,
  zf = NULL,
  trim = 0.05,
  maxLF = NULL,
  criterion = "bic"
) {
  y <- as.matrix(y)
  zb <- if (!is.null(zb)) as.matrix(zb)
  zf <- if (!is.null(zf)) as.matrix(zf)

  N <- nrow(y)

  wb <- cbind(
    .ones(N),
    if (trend) .trend(N) / N else NULL,
    zb
  )

  z <- cbind(zb, zf)

  st_t1 <- trunc(2 * trim * N)
  ed_t1 <- trunc((1 - 2 * trim) * N)

  w <- cbind(wb, zf)
  ssr0 <- as.numeric(crossprod(OLS.reg(y, w)$residuals))
  est_date <- N

  for (T1 in .seqi(st_t1, ed_t1)) {
    wb1 <- rbind(
      .zeros(T1, ncol(wb)),
      .msub(wb, .seqi(T1 + 1, N))
    )
    w <- cbind(wb, wb1, zf)
    ssr1 <- as.numeric(crossprod(OLS.reg(y, w)$residuals))
    if (ssr1 < ssr0) {
      ssr0 <- ssr1
      est_date <- T1
    }
  }

  wb1e <- rbind(
    .zeros(est_date, ncol(wb)),
    .msub(wb, .seqi(est_date + 1, N))
  )

  y <- .msub(y, .seqi(2, N))
  wb <- .msub(wb, .seqi(2, N))
  wb1e <- .msub(wb1e, .seqi(2, N))
  xf <- .msub(z, .seqi(2, N)) - .msub(z, .seqi(1, N - 1))

  zf <- .msub(zf, .seqi(2, N))
  wf <- cbind(zf, xf)

  N <- nrow(y)

  if (is.null(maxLF)) {
    maxLF <- trunc(4 * (N / 100)^0.25)
  }

  if (est_date - maxLF - 1 <= ncol(wb)) {
    maxLF <- est_date - ncol(wb) - 2
  } else if (N - maxLF - est_date <= ncol(wb)) {
    maxLF <- N - est_date - ncol(wb) - 1
  }

  idx0 <- .seqi(maxLF + 1, N - maxLF)
  y0 <- .msub(y, idx0)
  w0 <- cbind(
    .msub(wb, idx0),
    .msub(wb1e, idx0),
    .msub(wf, idx0)
  )

  ic0 <- info.criterions(OLS.reg(y0, w0), criterion)
  est_leads <- 0
  est_lags <- 0

  for (ld1 in .seqi(1, maxLF)) {
    for (lg1 in .seqi(1, maxLF)) {
      w1 <- w0
      for (k1 in 1:ld1) {
        w1 <- cbind(
          w1,
          .msub(xf, .seqi(maxLF + 1 - k1, N - maxLF - k1))
        )
      }
      for (k2 in 1:lg1) {
        w1 <- cbind(
          w1,
          .msub(xf, .seqi(maxLF + 1 + k2, N - maxLF + k2))
        )
      }
      ic1 <- info.criterions(OLS.reg(y0, w1), criterion)
      if (ic1 < ic0) {
        est_leads <- ld1
        est_lags <- lg1
        ic0 <- ic1
      }
    }
  }

  list(leads = est_leads, lags = est_lags)
}


#' Confidence sets for a break date based on optimal tests (Kurozumi & Yamamoto, 2015).
#'
#' @param y         \eqn{T_y \times 1} numeric matrix or vector, the dependent time series.
#' @param zb        \eqn{T_y \times pzb} matrix of I(1) regressors whose coefficient
#'                  sustains a structural change; \code{NULL} if there are none.
#' @param zf        \eqn{T_y \times pzf} matrix of I(1) regressors whose coefficient is
#'                  fixed (does not change at the break); \code{NULL} if there
#'                  are none.
#' @param trim      Trimming share defining lower and upper bounds of possible break moments.
#' @param c_level   Numeric, 0.90 or 0.95. Confidence level.
#'
#' @return Named list of five \eqn{T_y \times 1} numeric matrices (1 = point included
#'   in the confidence set):
#'   \item{cset_em}{Confidence set from the Elliott and Muller (EM) test.}
#'   \item{cset_mem}{Confidence set from the modified EM test.}
#'   \item{cset_sup}{Confidence set from the sup-test.}
#'   \item{cset_avg}{Confidence set from the avg-test.}
#'   \item{cset_exp}{Confidence set from the exp-test.}
coint.conf.sets.KY <- function(
  y,
  zb,
  zf = NULL,
  trim = 0.05,
  c_level = c(0.90, 0.95)
) {
  y <- as.matrix(y)
  zb <- if (!is.null(zb)) as.matrix(zb)
  zf <- if (!is.null(zf)) as.matrix(zf)

  N <- nrow(y)
  Nc <- ncol(zb)

  st1 <- trunc(trim * N)
  ed1 <- trunc((1 - trim) * N)

  cset_em <- .zeros(N, 1)
  cset_mem <- .zeros(N, 1)
  cset_sup <- .zeros(N, 1)
  cset_avg <- .zeros(N, 1)
  cset_exp <- .zeros(N, 1)

  w <- cbind(zb, zf)

  ssr0 <- as.numeric(crossprod(OLS.reg(y, w)$residuals))
  est_date <- N

  for (T1 in .seqi(st1, ed1)) {
    x2 <- rbind(
      .zeros(T1, Nc),
      .msub(zb, .seqi(T1 + 1, N))
    )

    w <- cbind(zb, x2, zf)
    ssr1 <- as.numeric(crossprod(OLS.reg(y, w)$residuals))
    if (ssr1 < ssr0) {
      ssr0 <- ssr1
      est_date <- T1
    }
  }

  x2b <- rbind(
    .zeros(est_date, Nc),
    .msub(zb, .seqi(est_date + 1, N))
  )

  for (T1 in .seqi(st1, ed1)) {
    L1 <- T1 / N
    x2 <- rbind(
      .zeros(T1, Nc),
      .msub(zb, .seqi(T1 + 1, N))
    )

    w <- cbind(zb, x2, zf)
    u_hat <- y - w %*% solve(crossprod(w), crossprod(w, y))
    xu <- zb * as.vector(u_hat)
    xu_em <- xu

    xu1 <- .msub(xu_em, .seqi(1, T1))
    xu2 <- .msub(xu_em, .seqi(T1 + 1, N))

    lrv1 <- .lr.matr.quad(xu1, 0.97)$Omega
    lrv2 <- .lr.matr.quad(xu2, 0.97)$Omega
    lrv1h <- .sym_mat_pow(lrv1, -0.5)
    lrv2h <- .sym_mat_pow(lrv2, -0.5)

    xu_em[.seqi(1, T1), ] <- .msub(xu_em, .seqi(1, T1)) %*% lrv1h
    xu_em[.seqi(T1 + 1, N), ] <- .msub(xu_em, .seqi(T1 + 1, N)) %*% lrv2h
    Sxu_em <- .cumsumc(xu_em)

    w <- cbind(w, if (abs(T1 - est_date) >= Nc) x2b else NULL)

    u_hat <- y - w %*% solve(crossprod(w), crossprod(w, y))
    xu_add <- zb * as.vector(u_hat)

    xu1 <- .msub(xu_add, .seqi(1, T1))
    xu2 <- .msub(xu_add, .seqi(T1 + 1, N))

    lrv1 <- .lr.matr.quad(xu1, 0.97)$Omega
    lrv2 <- .lr.matr.quad(xu2, 0.97)$Omega
    lrv1h <- .sym_mat_pow(lrv1, -0.5)
    lrv2h <- .sym_mat_pow(lrv2, -0.5)

    xu[.seqi(1, T1), ] <- .msub(xu, .seqi(1, T1)) %*% lrv1h
    xu[.seqi(T1 + 1, N), ] <- .msub(xu, .seqi(T1 + 1, N)) %*% lrv2h
    Sxu <- .cumsumc(xu)

    em_stat <- 0
    mem_stat <- 0
    sup_stat <- 0
    avg_stat <- 0
    exp_stat <- 0
    nbreak <- 0
    dbreak <- 0

    for (T2 in .seqi(st1, ed1)) {
      L2 <- T2 / N
      Gv_em <- as.vector(Sxu_em[T2, ])
      Gv <- as.vector(Sxu[T2, ])

      if (T1 >= T2) {
        em_stat <- em_stat + sum(Gv_em^2) / T1^2
        om <- L2 * (L1 - L2) / L1
        GG <- sum(Gv^2)
        GGT <- GG / N
        GGTo <- GGT / om
        mem_stat <- mem_stat + GG / T1^2
        if (sup_stat < GGT) {
          sup_stat <- GGT
        }
        avg_stat <- avg_stat + GGT
        if (abs(L2 - L1) > 0.05) {
          exp_stat <- exp_stat + om^(-Nc / 2) * exp(GGTo / 2)
        }
      } else {
        em_stat <- em_stat + sum(Gv_em^2) / (N - T1)^2
        om <- (L2 - L1) * (1 - L2) / (1 - L1)
        GG <- sum(Gv^2)
        GGT <- GG / N
        GGTo <- GGT / om
        mem_stat <- mem_stat + GG / (N - T1)^2
        if (sup_stat < GGT) {
          sup_stat <- GGT
        }
        avg_stat <- avg_stat + GGT
        if (abs(L2 - L1) > 0.05) {
          exp_stat <- exp_stat + om^(-Nc / 2) * exp(GGTo / 2)
        }
      }

      if (abs(L2 - L1) <= 0.05) {
        dbreak <- dbreak + 1
      }
      nbreak <- nbreak + 1
    }

    avg_stat <- avg_stat / (nbreak - 1)
    exp_stat <- log(exp_stat / (nbreak - dbreak))

    cv <- get.cv.KY(L1, c_level, Nc)

    if (em_stat <= cv$cval_em) {
      cset_em[T1, 1] <- 1
    }
    if (mem_stat <= cv$cval_em) {
      cset_mem[T1, 1] <- 1
    }
    if (sup_stat <= cv$cval_sup) {
      cset_sup[T1, 1] <- 1
    }
    if (avg_stat <= cv$cval_avg) {
      cset_avg[T1, 1] <- 1
    }
    if (exp_stat <= cv$cval_exp) cset_exp[T1, 1] <- 1
  }

  result <- list(
    cset_em = cset_em,
    cset_mem = cset_mem,
    cset_sup = cset_sup,
    cset_avg = cset_avg,
    cset_exp = cset_exp
  )
  class(result) <- c("bt_confSetOpt", "bt_confSet")

  result
}

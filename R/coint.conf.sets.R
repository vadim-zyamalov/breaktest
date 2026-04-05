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
#' @param z.lead,z.lag Number of leads and lags of `z` regressors.
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
coint.conf.sets <- function(y,
                            trend = FALSE,
                            zb = NULL,
                            zf = NULL,
                            z.lead = NULL,
                            z.lag = NULL,
                            conf.level = 0.9,
                            trim = 0.05,
                            criterion = "bic") {
  if (!is.matrix(y)) y <- as.matrix(y)
  if (!is.null(zb) && !is.matrix(zb)) zb <- as.matrix(zb)
  if (!is.null(zf) && !is.matrix(zf)) zf <- as.matrix(zf)

  if (is.null(z.lead) || is.null(z.lag)) {
    ll.est <- select.lead.lag.KS(y, trend, zb, zf, trim, criterion)
    z.lead <- ll.est$lead
    z.lag <- ll.est$lag
  }

  N <- nrow(y)
  td <- 1:N
  rows <- (z.lag + 2):(N - z.lead)

  wb <- cbind(
    .const(N),
    if (trend) .trend(N) / N else NULL,
    zb
  )
  z <- cbind(zb, zf)
  d.z <- .diffn(z)

  p.zb <- if (!is.null(zb)) ncol(zb) else 0
  p.zf <- if (!is.null(zf)) ncol(zf) else 0

  wf <- cbind(zf, d.z)

  if (z.lead > 0) {
    wf <- cbind(
      wf,
      apply(as.array(1:z.lead), 1, function(l) .lagn(d.z, -l))
    )
  }
  if (z.lag > 0) {
    wf <- cbind(
      wf,
      apply(as.array(1:z.lag), 1, function(l) .lagn(d.z, l))
    )
  }

  td <- td[rows, , drop = FALSE]
  y <- y[rows, , drop = FALSE]
  wb <- wb[rows, , drop = FALSE]
  wf <- wf[rows, , drop = FALSE]

  N2 <- nrow(y)
  first.break1 <- trunc(2 * trim * N2)
  last.break1 <- trunc((1 - 2 * trim) * N2)
  first.break2 <- trunc(trim * N2)
  last.break2 <- trunc((1 - trim) * N2)

  cset.sup <- numeric(N2)
  cset.avg <- numeric(N2)
  cset.exp <- numeric(N2)
  cset.bls <- numeric(N2)

  w <- cbind(wb, wf)
  u.hat <- .OLS(y, w)$residuals
  ssr.0 <- c(t(u.hat) %*% u.hat)
  est.date <- N2

  for (tb in first.break1:last.break1) {
    wb1 <- rbind(
      matrix(0, tb, ncol(wb)),
      as.matrix(wb[(tb + 1):N2, ])
    )
    w <- cbind(wb, wb1, wf)
    ww.inv <- solve(t(w) %*% w)
    u.hat <- .OLS(y, w)$residuals
    ssr.1 <- c(t(u.hat) %*% u.hat)
    if (ssr.1 < ssr.0) {
      ssr.0 <- ssr.1
      est.date <- tb
    }
  }

  wb1e <- rbind(
    matrix(0, est.date, ncol(wb)),
    as.matrix(wb[(est.date + 1):N2, ])
  )
  td[est.date, 1] <- -1

  w <- cbind(wb1e, wb, wf)
  b.hat <- solve(t(w) %*% w) %*% t(w) %*% y
  u.hat <- y - w %*% b.hat

  lrv.u <- lr.var(
    u.hat,
    demean = FALSE,
    kernel = "quadratic",
    limit.lags = TRUE,
    limit.selector = "Andrews"
  )

  l.hat <- (wb[est.date, ] %*% b.hat[seq_len(ncol(wb))])^2 / lrv.u
  if (conf.level == 0.9) {
    c.bls <- 7.686962
  } else if (conf.level == 0.95) {
    c.bls <- 11.03281
  }

  bdd <- trunc(c.bls / l.hat)
  bls.l <- est.date - bdd - 1
  bls.u <- est.date + bdd + 1
  if (bls.l < 1) {
    bls.l <- 1
  }
  if (bls.u > N2) {
    bls.u <- N2
  }

  cset.bls[bls.l:bls.u] <- 1

  for (tb in first.break1:last.break1) {
    lambda.1 <- tb / N2

    wb1 <- rbind(
      matrix(0, tb, ncol(wb)),
      as.matrix(wb[(tb + 1):N2, ])
    )

    w <- cbind(wb, wb1, wf)

    y.hat <- .OLS(y, w)$residuals

    if (abs(tb - est.date) > ncol(wb)) {
      we <- cbind(w, wb1e)
    } else {
      we <- w
    }
    be.hat <- solve(t(we) %*% we) %*% t(we) %*% y
    u.hat <- y - we %*% be.hat
    lrv.u2 <- lr.var(
      u.hat,
      demean = FALSE,
      kernel = "quadratic",
      limit.lags = TRUE,
      limit.selector = "Andrews"
    )

    sup.stat <- 0
    avg.stat <- 0
    exp.stat <- 0

    nbreak <- 0
    dbreak <- 0

    for (tb2 in first.break2:last.break2) {
      lambda.2 <- tb2 / N2

      if (abs(lambda.2 - lambda.1) <= 0.05) {
        dbreak <- dbreak + 1
      } else {
        wb2 <- rbind(
          matrix(0, tb2, ncol(wb)),
          as.matrix(wb[(tb2 + 1):N2, ])
        )

        r <- wb2 - wb1
        br.hat <- ww.inv %*% t(w) %*% r
        r.hat <- r - w %*% br.hat

        g <- t(r.hat) %*% y.hat
        h <- t(r.hat) %*% r.hat
        ghg <- c(t(g) %*% solve(h) %*% g) / lrv.u2

        if (sup.stat < ghg) {
          sup.stat <- ghg
        }

        avg.stat <- avg.stat + ghg
        exp.stat <- exp.stat + exp(ghg / 2)

        nbreak <- nbreak + 1
      }
    }

    avg.stat <- avg.stat / (nbreak - dbreak)
    exp.stat <- log(exp.stat / (nbreak - dbreak))

    cv <- get.cv.coint.conf.sets(
      lambda.1,
      trend,
      conf.level,
      p.zb,
      p.zf
    )
    if (sup.stat <= cv$cval_sup) {
      cset.sup[tb] <- 1
    }
    if (avg.stat <= cv$cval_avg) {
      cset.avg[tb] <- 1
    }
    if (exp.stat <= cv$cval_exp) {
      cset.exp[tb] <- 1
    }
  }

  if (z.lead == 0) {
    td <- rbind(
      matrix(0, (z.lag + 1), 1),
      td
    )
    cset.sup <- rbind(
      matrix(0, (z.lag + 1), 1),
      matrix(cset.sup, N2, 1)
    )
    cset.avg <- rbind(
      matrix(0, (z.lag + 1), 1),
      matrix(cset.avg, N2, 1)
    )
    cset.exp <- rbind(
      matrix(0, (z.lag + 1), 1),
      matrix(cset.exp, N2, 1)
    )
    cset.bls <- rbind(
      matrix(0, (z.lag + 1), 1),
      matrix(cset.bls, N2, 1)
    )
  } else {
    td <- rbind(
      matrix(0, (z.lag + 1), 1),
      td,
      matrix(0, z.lead, 1)
    )
    cset.sup <- rbind(
      matrix(0, (z.lag + 1), 1),
      matrix(cset.sup, N2, 1),
      matrix(0, z.lead, 1)
    )
    cset.avg <- rbind(
      matrix(0, (z.lag + 1), 1),
      matrix(cset.avg, N2, 1),
      matrix(0, z.lead, 1)
    )
    cset.exp <- rbind(
      matrix(0, (z.lag + 1), 1),
      matrix(cset.exp, N2, 1),
      matrix(0, z.lead, 1)
    )
    cset.bls <- rbind(
      matrix(0, (z.lag + 1), 1),
      matrix(cset.bls, N2, 1),
      matrix(0, z.lead, 1)
    )
  }

  list(
    td = td,
    cset.sup = cset.sup,
    cset.avg = cset.avg,
    cset.exp = cset.exp,
    cset.bls = cset.bls
  )
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
select.lead.lag.KS <- function(y,
                               trend = TRUE,
                               zb = NULL,
                               zf = NULL,
                               trim = 0.05,
                               criterion = "bic") {
  if (!criterion %in% c("bic", "aic", "hq", "lwz")) {
    stop("ERROR! Unknown criterion")
  }

  if (!is.matrix(y)) y <- as.matrix(y)
  if (!is.null(zb) && !is.matrix(zb)) zb <- as.matrix(zb)
  if (!is.null(zf) && !is.matrix(zf)) zf <- as.matrix(zf)

  N <- nrow(y)

  first.break <- trunc(2 * trim * N)
  last.break <- trunc((1 - 2 * trim) * N)

  wb <- cbind(
    .const(N),
    if (trend) .trend(N) / N else NULL,
    zb
  )

  z <- cbind(zb, zf)
  w <- cbind(wb, zf)
  d.z <- .diffn(z)

  wf <- cbind(zf, d.z)

  u.hat <- .OLS(y, w)$residuals
  ssr.0 <- drop(t(u.hat) %*% u.hat)
  est.date <- N

  for (t in first.break:last.break) {
    wb1 <- rbind(0, wb[(t + 1):N, , drop = FALSE])
    w <- cbind(wb, wb1, zf)
    u.hat <- .OLS(y, w)$residuals
    ssr.1 <- drop(t(u.hat) %*% u.hat)

    if (ssr.1 < ssr.0) {
      ssr.0 <- ssr.1
      est.date <- t
    }
  }

  est.dt <- wb[(est.date + 1):N, , drop = FALSE]
  y <- y[2:N, , drop = FALSE]
  wb <- wb[2:N, , drop = FALSE]
  wf <- wf[2:N, , drop = FALSE]

  N <- nrow(y)

  max.lead.lag <- trunc(4 * (N / 100)^(1 / 4))
  if ((est.date - max.lead.lag - 1) <= ncol(wb)) {
    max.lead.lag <- est.date - ncol(wb) - 2
  } else if ((N - max.lead.lag - est.date) <= ncol(wb)) {
    max.lead.lag <- N - est.date - ncol(wb) - 1
  }

  y.0 <- as.matrix(y[(max.lead.lag + 1):(N - max.lead.lag), ])
  w.0 <- cbind(
    as.matrix(wb[(max.lead.lag + 1):(N - max.lead.lag), ]),
    as.matrix(est.dt[(max.lead.lag + 1):(N - max.lead.lag), ]),
    as.matrix(wf[(max.lead.lag + 1):(N - max.lead.lag), ])
  )

  u.hat <- .OLS(y.0, w.0)$residuals
  min.ic <- info.criterion(u.hat, ncol(w.0))[[criterion]]
  est.lead <- 0
  est.lag <- 0

  for (cur.lead in 1:max.lead.lag) {
    for (cur.lag in 1:max.lead.lag) {
      w.1 <- w.0
      for (k in 1:cur.lead) {
        w.1 <- cbind(
          w.1,
          as.matrix(
            d.z[(max.lead.lag + 1 - k):(N - max.lead.lag - k), ]
          )
        )
      }
      for (k in 1:cur.lag) {
        w.1 <- cbind(
          w.1,
          as.matrix(
            d.z[(max.lead.lag + 1 + k):(N - max.lead.lag + k), ]
          )
        )
      }

      u.hat <- .OLS(y.0, w.1)$residuals

      cur.ic <- info.criterion(u.hat, ncol(w.1))[[criterion]]
      if (cur.ic < min.ic) {
        est.lead <- cur.lead
        est.lag <- cur.lag
        min.ic <- cur.ic
      }
    }
  }

  list(
    lead = est.lead,
    lag = est.lag
  )
}

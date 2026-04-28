#' @title
#' Perron-Yabu (2009) statistic for break at unknown date.
#'
#' @details
#' The code provided is the original Ox code by Skrobotov (2018)
#' ported to R.
#'
#' @param y A time series of interest.
#' @param const,trend Allowing the break in constant or trend.
#' @param criterion Needed information criterion: aic, bic, hq or lwz.
#' @param trim A trimming parameter to determine the lower and upper bounds for
#' a possible break point.
#' @param max.lag The maximum possible lag in the model.
#'
#' @return A list of the estimated Wald statistic as well as its c.v.
#'
#' @references
#' Perron, Pierre, and Tomoyoshi Yabu.
#' “Testing for Shifts in Trend With an Integrated or
#' Stationary Noise Component.”
#' Journal of Business & Economic Statistics 27, no. 3 (July 2009): 369–96.
#' https://doi.org/10.1198/jbes.2009.07268.
#'
#' @export
PY.statistic <- function(
  y,
  const = FALSE,
  trend = FALSE,
  criterion = "aic",
  trim = 0.15,
  max.lag = trunc(12 * (length(y) / 100)^(1 / 4))
) {
  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }
  if (!trim %in% c(0.01, 0.05, 0.10, 0.15, 0.25)) {
    stop("ERROR! PY.single: Illegal trim value")
  }

  model <- paste0(
    ifelse(const, "c", ""),
    ifelse(trend, "t", "")
  )
  if (model == "") {
    stop("ERROR! PY.single: Illegal model")
  }

  VR <- .cval_PY_single[[model]]$VR
  v.t <- .cval_PY_single[[model]]$v.t
  c.v <- .cval_PY_single[[model]]$c.v

  N <- nrow(y)
  first.break <- max(trunc(trim * N), max.lag + 2)
  last.break <- trunc((1 - trim) * N)

  vect1 <- NULL

  for (tb in first.break:last.break) {
    lambda <- (tb - 1) / N

    reg <- cbind(
      .const(N),
      .trend(N),
      if (const) .du(tb, N) else NULL,
      if (trend) .dt(tb, N) else NULL
    )

    # Estimation of alpha
    khat <- max(1, AR.reg(y, reg, max.lag, criterion)$lag)
    u <- OLS.reg(y, reg)$residuals
    du <- .diffn(u, na = 0)
    depu <- u[khat:N]
    regu <- .lagn(u, 1, na = 0)
    for (i in seq_len(khat - 1)) {
      regu <- cbind(regu, .lagn(du, i, na = 0))
    }
    regu <- regu[khat:N, , drop = FALSE]

    tmp.OLS <- OLS.reg(depu, regu)
    b <- tmp.OLS$coefficients
    ehat <- tmp.OLS$residuals

    VCV <- solve(t(regu) %*% regu) * sum(ehat^2) / length(ehat)

    ahat <- b[1]
    vahat <- VCV[1, 1]
    tau1 <- (ahat - 1) / sqrt(vahat)

    # Upper Biased Estimator
    tau05 <- v.t[ceiling(lambda * 10)]

    IP <- trunc((khat + 1) / 2)
    r <- ncol(reg)
    k <- 10

    c1 <- sqrt((1 + r) * N)
    c2 <- ((1 + r) * N - tau05^2 * (IP + N)) / (tau05 * (tau05 + k) * (IP + N))

    if (tau1 > tau05) {
      ctau <- -tau1
    } else if (tau1 <= tau05 && tau1 > -k) {
      ctau <- IP * tau1 / N - (r + 1) / (tau1 + c2 * (tau1 + k))
    } else if (tau1 <= -k && tau1 > -c1) {
      ctau <- IP * tau1 / N - (r + 1) / tau1
    } else if (tau1 <= -c1) {
      ctau <- 0
    }

    amus <- ahat + ctau * sqrt(vahat)
    if (amus >= 1) {
      amus <- 1
    } else if (amus <= -1) {
      amus <- -0.99
    }

    CR <- sqrt(N) * abs(amus - 1)
    if (CR <= 1) {
      amus <- 1
    }

    gdep <- rbind(
      y[1, , drop = FALSE],
      y[2:N, , drop = FALSE] -
        amus * y[1:(N - 1), , drop = FALSE]
    )
    greg <- rbind(
      reg[1, , drop = FALSE],
      reg[2:N, , drop = FALSE] -
        amus * reg[1:(N - 1), , drop = FALSE]
    )

    tmp.OLS <- OLS.reg(gdep, greg)
    b <- tmp.OLS$coefficients
    v <- tmp.OLS$residuals

    if (khat == 1) {
      h0 <- sum(v^2) / length(v)
    } else {
      if (amus == 1) {
        regv <- NULL
        for (ki in seq_len(khat - 1)) {
          regv <- cbind(regv, .lagn(v, ki, na = 0))
        }

        depv <- v[(khat - 1):N]
        regv <- regv[(khat - 1):N, , drop = FALSE]

        tmp.OLS <- OLS.reg(depv, regv)
        beta <- tmp.OLS$coefficients
        e <- tmp.OLS$residuals

        if (const && !trend) {
          vbeta <- matrix(0, nrow = khat - 1, ncol = 3)
          for (ki in seq_len(khat - 1)) {
            DUki <- .du(tb - ki, N)
            regki <- cbind(.const(N), .trend(N), DUki)
            gdepki <- rbind(
              y[1, , drop = FALSE],
              y[2:N, , drop = FALSE] - amus * y[1:(N - 1), , drop = FALSE]
            )
            gregki <- rbind(
              reg[1, ],
              regki[2:N, ] - amus * regki[1:(N - 1), ]
            )
            vbeta[ki, ] <- drop(OLS.reg(gdepki, gregki)$coefficients)
          }
          b[3] <- b[3] - drop(t(vbeta[, 3]) %*% beta)
          h0 <- sum(e^2) / (N - khat)
        }

        if (!const && trend) {
          h0 <- (sum(e^2) / (N - khat)) /
            ((1 - sum(beta))^2)
        }

        if (const && trend) {
          vbeta <- matrix(0, nrow = khat - 1, ncol = 4)
          for (ki in seq_len(khat - 1)) {
            regki <- cbind(
              .const(N),
              .trend(N),
              .du(tb - ki, N),
              .du(tb - ki, N) * (.trend(N) - tb)
            )
            gdepki <- rbind(
              y[1, , drop = FALSE],
              y[2:N, , drop = FALSE] - amus * y[1:(N - 1), , drop = FALSE]
            )
            gregki <- rbind(
              reg[1, ],
              regki[2:N, ] - amus * regki[1:(N - 1), ]
            )
            vbeta[ki, ] <- drop(OLS.reg(gdepki, gregki)$coefficients)
          }

          sige <- sum(e^2) / (N - khat)

          h0 <- sige / ((1 - sum(beta))^2)
          b[3] <- (sqrt(h0) / sqrt(sige)) *
            (b[3] - drop(t(vbeta[, 3]) %*% beta))
        }
      }

      if (abs(amus) < 1) {
        h0 <- .lr.var.quad(v)
      }
    }

    VCV <- h0 * solve(t(greg) %*% greg)
    vect1 <- c(
      vect1,
      t(VR %*% b) %*% solve(VR %*% VCV %*% t(VR)) %*% (VR %*% b)
    )
  }

  wald <- log(sum(exp(vect1 / 2)) / N)

  if (trim == 0.01) {
    cv <- c.v[1, ]
  }
  if (trim == 0.05) {
    cv <- c.v[2, ]
  }
  if (trim == 0.10) {
    cv <- c.v[3, ]
  }
  if (trim == 0.15) {
    cv <- c.v[4, ]
  }
  if (trim == 0.25) {
    cv <- c.v[5, ]
  }

  result <- list(
    statistic = wald,
    critical.value = cv
  )
  class(result) <- "bt_PY"

  result
}

#' @title
#' Sequential Perron-Yabu (2009) statistic for breaks at unknown date.
#'
#' @details
#' The code provided is the original Ox code by Skrobotov (2018)
#' ported to R.
#'
#' @param y A time series of interest.
#' @param const Allowing the break in constant.
#' @param breaks A number of breaks.
#' @param criterion Needed information criterion: aic, bic, hq or lwz.
#' @param trim A trimming value for a possible break date bounds.
#' @param max.lag The maximum possible lag in the model.
#'
#' @return An estimated Wald statistic.
#'
#' @references
#' Kejriwal, Mohitosh, and Pierre Perron.
#' “A Sequential Procedure to Determine the Number of Breaks in Trend
#' with an Integrated or Stationary Noise Component:
#' Determination of Number of Breaks in Trend.”
#' Journal of Time Series Analysis 31, no. 5 (September 2010): 305–28.
#' https://doi.org/10.1111/j.1467-9892.2010.00666.x.
#'
#' @export
KP.seq.statistic <- function(
  y,
  const = FALSE,
  breaks = 1,
  criterion = "aic",
  trim = 0.15,
  max.lag = trunc(12 * (length(y) / 100)^(1 / 4))
) {
  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }

  model <- ifelse(const, 2, 1)
  R <- .cval_PY_sequential[[model]]$R
  v.t <- .cval_PY_sequential[[model]]$v.t

  N <- nrow(y)
  h <- trunc(trim * N)

  if (breaks == 0) {
    datevec <- c(1, N + 1)
  } else {
    SSR.data <- SSR.matrix(y, cbind(.const(N), .trend(N)), h)
    dates <- segments.BP(
      y,
      cbind(.const(N), .trend(N)),
      breaks,
      h,
      SSR.data
    )
    datevec <- c(1, sort(drop(dates$break.point)), N + 1)
  }
  wald <- NULL

  for (i in seq_len(breaks + 1)) {
    T_i <- datevec[i + 1] - datevec[i]
    vect1 <- NULL

    t.low <- max(trunc(datevec[i] + T_i * trim), max.lag + 2)
    t.high <- trunc(datevec[i + 1] - T_i * trim)

    if (t.low < t.high - 1) {
      for (tb in t.low:t.high) {
        lam1 <- (tb - 1) / (datevec[i + 1] - 1)

        reg <- cbind(
          .const(N),
          if (const) .du(tb - 1, N) else NULL,
          .trend(N) - (datevec[i] - 1),
          .dt(tb - 1, N)
        )

        y_i <- y[datevec[i]:(datevec[i + 1] - 1), , drop = FALSE]
        reg_i <- reg[datevec[i]:(datevec[i + 1] - 1), , drop = FALSE]

        khat <- max(1, AR.reg(y_i, reg_i, max.lag, criterion)$lag)

        u <- OLS.reg(y_i, reg_i)$residuals
        du <- .diffn(u, na = 0)

        depu <- u[khat:length(u)]
        regu <- .lagn(u, 1, na = 0)
        for (l in seq_len(khat - 1)) {
          regu <- cbind(regu, .lagn(du, l, na = 0))
        }
        regu <- regu[khat:length(u), , drop = FALSE]

        tmp.OLS <- OLS.reg(depu, regu)
        b <- tmp.OLS$coefficients
        ehat <- tmp.OLS$residuals

        VCV <- solve(t(regu) %*% regu) * sum(ehat^2) / length(ehat)

        ahat <- b[1]
        vahat <- VCV[1, 1]
        tau1 <- (ahat - 1) / sqrt(vahat)

        # Upper Biased Estimator
        t05 <- v.t[ceiling(lam1 * 10)]

        IP <- trunc((khat + 1) / 2)
        r <- ncol(reg)
        k <- 10

        c1 <- sqrt((1 + r) * T_i)
        c2 <- ((1 + r) * T_i - t05^2 * (IP + T_i)) /
          (t05 * (t05 + k) * (IP + T_i))

        if (tau1 > t05) {
          c.tau <- -tau1
        }
        if (tau1 <= t05 && tau1 > -k) {
          c.tau <- IP * tau1 / N - (r + 1) / (tau1 + c2 * (tau1 + k))
        }
        if (tau1 <= -k && tau1 > -c1) {
          c.tau <- IP * tau1 / N - (r + 1) / tau1
        }
        if (tau1 <= -c1) {
          c.tau <- 0
        }

        amus <- ahat + c.tau * sqrt(vahat)
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
          y[datevec[i], , drop = FALSE],
          y[(datevec[i] + 1):(datevec[i + 1] - 1), , drop = FALSE] -
            amus * y[datevec[i]:(datevec[i + 1] - 2), , drop = FALSE] # nolint
        )
        greg <- rbind(
          reg[datevec[i], , drop = FALSE],
          reg[(datevec[i] + 1):(datevec[i + 1] - 1), , drop = FALSE] -
            amus * reg[datevec[i]:(datevec[i + 1] - 2), , drop = FALSE] # nolint
        )

        tmp.OLS <- OLS.reg(gdep, greg)
        b <- tmp.OLS$coefficients
        v <- tmp.OLS$residuals

        if (khat == 1) {
          h0 <- sum(v^2) / length(v)
        } else {
          if (amus == 1) {
            regv <- NULL
            for (l in seq_len(khat - 1)) {
              regv <- cbind(regv, .lagn(v, l, na = 0))
            }

            depv <- v[(khat - 1):length(v)]
            regv <- regv[(khat - 1):length(v), , drop = FALSE]

            tmp.OLS <- OLS.reg(depv, regv)
            beta <- tmp.OLS$coefficients
            e <- tmp.OLS$residuals

            if (!const) {
              h0 <- (sum(e^2) / (T_i - khat)) / ((1 - sum(beta))^2)
            }
            if (const) {
              vbeta <- matrix(0, nrow = khat - 1, ncol = 4)
              for (ki in seq_len(khat - 1)) {
                regki <- cbind(
                  .const(N),
                  .du(tb - ki - 1, N),
                  .trend(N),
                  .du(tb - ki - 1, N) * (.trend(N) - (tb - 1))
                )
                gdepki <- rbind(
                  y[datevec[i], , drop = FALSE],
                  y[(datevec[i] + 1):(datevec[i + 1] - 1), , drop = FALSE] - # nolint
                    amus * y[(datevec[i]):(datevec[i + 1] - 2), , drop = FALSE] # nolint
                )
                gregki <- rbind(
                  reg[datevec[i], ],
                  regki[(datevec[i] + 1):(datevec[i + 1] - 1), ] - # nolint
                    amus * regki[datevec[i]:(datevec[i + 1] - 2), ] # nolint
                )
                vbeta[ki, ] <- drop(OLS.reg(gdepki, gregki)$coefficients)
              }
              sige <- sum(e^2) / (T_i - khat)
              h0 <- sige / ((1 - sum(beta))^2)
              b[2] <- (sqrt(h0) / sqrt(sige)) *
                (b[2] - drop(t(vbeta[, 2]) %*% beta))
            }
          }

          if (abs(amus) < 1) {
            h0 <- .lr.var.quad(v)
          }
        }
        VCV <- h0 * qr.solve(t(greg) %*% greg)
        vect1 <- c(
          vect1,
          t(R %*% b) %*%
            qr.solve(R %*% VCV %*% t(R)) %*%
            (R %*% b)
        )
      }

      wald <- c(wald, log(sum(exp(vect1 / 2)) / T_i))
    }
  }

  max(wald)
}

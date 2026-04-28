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
    date.vec <- c(1, N + 1)
  } else {
    SSR.data <- SSR.matrix(y, cbind(.const(N), .trend(N)), h)
    dates <- segments.BP(
      y,
      cbind(.const(N), .trend(N)),
      breaks,
      h,
      SSR.data
    )
    date.vec <- c(1, drop(dates$break.point) + 1, N + 1)
  }

  wald <- numeric(breaks + 1)

  for (i in 1:(breaks + 1)) {
    N.i <- date.vec[i + 1] - date.vec[i]
    vect1 <- numeric(N)

    t.low <- max(trunc(date.vec[i] + N.i * trim - 1), max.lag + 2)
    t.high <- trunc(date.vec[i + 1] - N.i * trim - 1)

    if (t.low < t.high - 1) {
      for (tb in t.low:t.high) {
        lambda <- (tb - 1) / (date.vec[i + 1] - 1)

        x <- cbind(
          .const(N),
          if (const) .du(tb, N) else NULL,
          .trend(N) - date.vec[i] + 1,
          .dt(tb, N)
        )

        y.i <- y[date.vec[i]:(date.vec[i + 1] - 1), , drop = FALSE]
        x.i <- x[date.vec[i]:(date.vec[i + 1] - 1), , drop = FALSE]

        k.hat <- max(1, AR.reg(y.i, x.i, max.lag, criterion)$lag)

        resids <- OLS.reg(y.i, x.i)$residuals
        d.resid <- .diffn(resids, na = 0)

        y.u <- resids[k.hat:length(resids)]
        x.u <- .lagn(resids, 1, na = 0)
        if (k.hat > 1) {
          x.u <- cbind(
            x.u,
            apply(
              as.array(1:(k.hat - 1)),
              1,
              function(i) .lagn(d.resid, i, na = 0)
            )
          )
        }
        x.u <- x.u[k.hat:nrow(x.u), , drop = FALSE]

        tmp.OLS <- OLS.reg(y.u, x.u)
        beta.u <- tmp.OLS$coefficients
        u.resid <- tmp.OLS$residuals
        rm(tmp.OLS)

        VCV <- qr.solve(t(x.u) %*% x.u) * sum(u.resid^2) / nrow(u.resid)

        a.hat <- beta.u[1]
        var.a.hat <- VCV[1, 1]
        tau <- (a.hat - 1) / sqrt(var.a.hat)

        tau05 <- v.t[ceiling(lambda * 10)]

        IP <- trunc((k.hat + 1) / 2)

        k <- 10
        k.x <- ncol(x)

        c1 <- sqrt((1 + k.x) * N.i)
        c2 <- ((1 + k.x) * N.i - tau05^2 * (IP + N.i)) /
          (tau05 * (tau05 + k) * (IP + N.i))

        if (tau > tau05) {
          c.tau <- -tau
        }
        if (tau <= tau05 && tau > -k) {
          c.tau <- IP * tau / N - (k.x + 1) / (tau + c2 * (tau + k))
        }
        if (tau <= -k && tau > -c1) {
          c.tau <- IP * tau / N - (k.x + 1) / tau
        }
        if (tau <= -c1) {
          c.tau <- 0
        }

        a.hat.M <- a.hat + c.tau * sqrt(var.a.hat)
        if (a.hat.M >= 1) {
          a.hat.M <- 1
        } else if (a.hat.M <= -1) {
          a.hat.M <- -0.99
        }

        CR <- sqrt(N) * abs(a.hat.M - 1)
        if (CR <= 1) {
          a.hat.M <- 1
        }

        y.g <- rbind(
          y[date.vec[i], , drop = FALSE],
          y[(date.vec[i] + 1):(date.vec[i + 1] - 1), , drop = FALSE] -
            a.hat.M * y[date.vec[i]:(date.vec[i + 1] - 2), , drop = FALSE] # nolint
        )
        x.g <- rbind(
          x[date.vec[i], , drop = FALSE],
          x[(date.vec[i] + 1):(date.vec[i + 1] - 1), , drop = FALSE] -
            a.hat.M * x[date.vec[i]:(date.vec[i + 1] - 2), , drop = FALSE] # nolint
        )

        tmp.OLS <- OLS.reg(y.g, x.g)
        beta.g <- tmp.OLS$coefficients
        g.resid <- tmp.OLS$residuals
        rm(tmp.OLS)

        if (k.hat == 1) {
          h0 <- sum(g.resid^2) / length(g.resid)
        } else {
          if (a.hat.M == 1) {
            x.v <- apply(
              as.array(1:(k.hat - 1)),
              1,
              function(i) .lagn(g.resid, i, na = 0)
            )

            y.v <- g.resid[(k.hat - 1):length(g.resid)]
            x.v <- x.v[(k.hat - 1):length(g.resid), , drop = FALSE]

            tmp.OLS <- OLS.reg(y.v, x.v)
            beta.v <- tmp.OLS$coefficients
            v.resid <- tmp.OLS$residuals
            rm(tmp.OLS)

            if (!const) {
              h0 <- (sum(v.resid^2) / (N.i - k.hat)) / ((1 - sum(beta.v))^2)
            }
            if (const) {
              BETAS <- matrix(0, nrow = k.hat - 1, ncol = 4)
              for (k.i in 1:(k.hat - 1)) {
                x.ki <- cbind(
                  .const(N),
                  .du(tb - k.i, N),
                  .trend(N),
                  .dt(tb - k.i, N)
                )
                x.g.ki <- rbind(
                  x.ki[date.vec[i] + 1, ],
                  x.ki[(date.vec[i] + 2):(date.vec[i + 1] - 1), ] - # nolint
                    a.hat.M * x.ki[(date.vec[i] + 1):(date.vec[i + 1] - 2), ] # nolint
                )
                beta.ki <- OLS.reg(y.g, x.g.ki)$coefficients
                BETAS[k.i, ] <- drop(beta.ki)
                sig.e <- sum(v.resid^2) / (N.i - k.hat)
                beta.g[2] <- (sqrt(h0) / sqrt(sig.e)) *
                  (beta.g[2] - drop(t(BETAS[, 2]) %*% beta.v))
                h0 <- sig.e / ((1 - sum(beta.v))^2)
              }
            }
          }

          if (abs(a.hat.M) < 1) {
            h0 <- .lr.var.quad(g.resid)
          }
        }

        VCV <- h0 * qr.solve(t(x.g) %*% x.g)
        vect1[tb] <- t(R %*% beta.g) %*%
          qr.solve(R %*% VCV %*% t(R)) %*%
          (R %*% beta.g)
      }

      vect1 <- vect1[t.low:t.high]
      wald[i] <- log(sum(exp(vect1 / 2)) / (date.vec[i + 1] - date.vec[i])) # nolint
    }
  }

  max(wald)
}

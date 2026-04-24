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
PY.single <- function(y,
                      const = FALSE,
                      trend = FALSE,
                      criterion = "aic",
                      trim = 0.15,
                      max.lag) {
  if (!is.matrix(y)) y <- as.matrix(y)
  if (!trim %in% c(0.01, 0.05, 0.10, 0.15, 0.25)) {
    stop("ERROR! PY.single: Illegal trim value")
  }

  model <- paste0(
    ifelse(const, "c", ""),
    ifelse(trend, "t", "")
  )
  if (model == "") stop("ERROR! PY.single: Illegal model")

  VR <- .cval_PY_single[[model]]$VR
  v.t <- .cval_PY_single[[model]]$v.t
  c.v <- .cval_PY_single[[model]]$c.v

  N <- nrow(y)
  first.break <- max(trunc(trim * N), max.lag + 2) + 1
  last.break <- trunc((1 - trim) * N) + 1

  vect1 <- matrix(0, nrow = trunc((1 - 2 * trim) * N) + 2, ncol = 1)

  for (tb in first.break:last.break) {
    lambda <- tb / N

    x <- cbind(
      .const(N),
      .trend(N),
      if (const) .du(tb, N) else NULL,
      if (trend) .dt(tb, N) else NULL
    )

    k.hat <- max(1, AR.reg(y, x, max.lag, criterion)$lag)
    resids <- OLS.reg(y, x)$residuals
    d.resid <- .diffn(resids)

    y.u <- resids[k.hat:N, , drop = FALSE]
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
    x.u <- x.u[k.hat:N, , drop = FALSE]

    tmp.OLS <- OLS.reg(y.u, x.u)
    beta.u <- tmp.OLS$coefficients
    u.resid <- tmp.OLS$residuals
    rm(tmp.OLS)

    VCV <- qr.solve(t(x.u) %*% x.u) *
      drop(t(u.resid) %*% u.resid) / nrow(u.resid)

    a.hat <- beta.u[1]
    var.a.hat <- VCV[1, 1]
    tau <- (a.hat - 1) / sqrt(var.a.hat)

    tau05 <- v.t[ceiling(lambda * 10)]

    IP <- trunc((k.hat + 1) / 2)

    k <- 10
    k.x <- ncol(x)

    c1 <- sqrt((1 + k.x) * N)
    c2 <- ((1 + k.x) * N - tau05^2 * (IP + N)) /
      (tau05 * (tau05 + k) * (IP + N))

    if (tau > tau05) {
      c.tau <- -tau
    } else if (tau <= tau05 && tau > -k) {
      c.tau <- IP * tau / N - (k.x + 1) / (tau + c2 * (tau + k))
    } else if (tau <= -k && tau > -c1) {
      c.tau <- IP * tau / N - (k.x + 1) / tau
    } else if (tau <= -c1) {
      c.tau <- 0
    }

    a.hat.M <- a.hat + c.tau * sqrt(var.a.hat)
    if (a.hat.M >= 1) {
      a.hat.M <- 1
    } else if (a.hat.M <= -1) {
      a.hat.M <- -0.99
    }

    CR <- sqrt(N) * abs(a.hat.M - 1)
    if (CR <= 1) a.hat.M <- 1

    y.g <- rbind(
      y[1, , drop = FALSE],
      y[2:N, , drop = FALSE] -
        a.hat.M * y[1:(N - 1), , drop = FALSE]
    )
    x.g <- rbind(
      x[1, , drop = FALSE],
      x[2:N, , drop = FALSE] -
        a.hat.M * x[1:(N - 1), , drop = FALSE]
    )

    tmp.OLS <- OLS.reg(y.g, x.g)
    beta.g <- tmp.OLS$coefficients
    g.resid <- tmp.OLS$residuals
    rm(tmp.OLS)

    if (k.hat == 1) {
      h0 <- drop(t(g.resid) %*% g.resid) / nrow(g.resid)
    } else {
      if (a.hat.M == 1) {
        x.v <- NULL
        for (k.i in 1:(k.hat - 1)) {
          x.v <- cbind(x.v, .lagn(g.resid, k.i, na = 0))
        }

        y.v <- g.resid[(k.hat - 1):nrow(g.resid), , drop = FALSE]
        x.v <- x.v[(k.hat - 1):nrow(g.resid), , drop = FALSE]

        tmp.OLS <- OLS.reg(y.v, x.v)
        beta.v <- tmp.OLS$coefficients
        v.resid <- tmp.OLS$residuals
        rm(tmp.OLS)

        if (const && !trend) {
          BETAS <- matrix(0, nrow = k.hat - 1, ncol = 3)
          for (k.i in 1:(k.hat - 1)) {
            DU.ki <- .dt(tb - k.i, N)
            x.ki <- cbind(.const(N), DU.ki, .trend(N))
            x.g.ki <- rbind(
              x.ki[1, ],
              x.ki[2:N, ] - a.hat.M * x.ki[1:(N - 1), ]
            )
            beta.ki <- OLS.reg(y.g, x.g.ki)$coefficients
            BETAS[k.i, ] <- drop(beta.ki)
          }
          beta.g[2] <- beta.g[2] - drop(t(BETAS[, 2]) %*% beta.v)
          h0 <- drop(t(v.resid) %*% v.resid) / (N - k.hat)
        }

        if (!const && trend) {
          h0 <- (drop(t(v.resid) %*% v.resid) / (N - k.hat)) /
            ((1 - sum(beta.v))^2)
        }

        if (const && trend) {
          BETAS <- matrix(0, nrow = k.hat - 1, ncol = 4)
          for (k.i in 1:(k.hat - 1)) {
            x.ki <- cbind(
              .const(N),
              .trend(N),
              .du(tb - k.i, N),
              .dt(tb - k.i, N)
            )
            x.g.ki <- rbind(
              x.ki[1, ],
              x.ki[2:N, ] - a.hat.M * x.ki[1:(N - 1), ]
            )
            beta.ki <- OLS.reg(y.g, x.g.ki)$coefficients
            BETAS[k.i, ] <- drop(beta.ki)
          }

          sig.e <- drop(t(v.resid) %*% v.resid) / (N - k.hat)

          beta.g[2] <- (sqrt(h0) / sqrt(sig.e)) *
            (beta.g[2] - drop(t(BETAS[, 2]) %*% beta.v))
          h0 <- sig.e / ((1 - sum(beta.v))^2)
        }
      }

      if (abs(a.hat.M) < 1) {
        h0 <- .lr.var.quad(g.resid)
      }
    }

    VCV <- h0 * qr.solve(t(x.g) %*% x.g)
    vect1[tb - first.break + 1] <- t(VR %*% beta.g) %*%
      qr.solve(VR %*% VCV %*% t(VR)) %*% (VR %*% beta.g)
  }

  wald <- log(sum(exp(vect1 / 2)) / N)

  if (trim == 0.01) cv <- c.v[1, ]
  if (trim == 0.05) cv <- c.v[2, ]
  if (trim == 0.10) cv <- c.v[3, ]
  if (trim == 0.15) cv <- c.v[4, ]
  if (trim == 0.25) cv <- c.v[5, ]

  list(
    statistic = wald,
    critical.value = cv
  )
}

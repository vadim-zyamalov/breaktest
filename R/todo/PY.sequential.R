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
PY.sequential <- function(y,
                          const = FALSE,
                          breaks = 1,
                          criterion = "aic",
                          trim = 0.15,
                          max.lag = 1) {
    if (!is.matrix(y)) y <- as.matrix(y)

    if (!const) {
        R <- matrix(c(0, 0, 1), nrow = 1, ncol = 3, byrow = TRUE)
        v.t <- as.matrix(c(
            -4.27, -4.41, -4.51, -4.55, -4.56,
            -4.57, -4.51, -4.38, -4.26, -4.26
        ))
    } else if (const) {
        R <- matrix(c(
            0, 1, 0, 0,
            0, 0, 0, 1
        ), nrow = 2, ncol = 4, byrow = TRUE)
        v.t <- as.matrix(c(
            -4.38, -4.65, -4.78, -4.81, -4.90,
            -4.88, -4.75, -4.70, -4.41, -4.41
        ))
    }

    n.obs <- nrow(y)
    x.const <- rep(1, n.obs)
    x.trend <- 1:n.obs

    h <- trunc(trim * n.obs)

    if (breaks == 0) {
        date.vec <- c(1, n.obs + 1)
    } else {
        SSR.data <- SSR.matrix(y, cbind(x.const, x.trend), h)
        dates <- segments.OLS.N.breaks(
            y,
            cbind(x.const, x.trend),
            breaks,
            h,
            SSR.data
        )
        date.vec <- c(1, drop(dates$break.point) + 1, n.obs + 1)
    }

    wald <- rep(0, breaks + 1)

    for (i in 1:(breaks + 1)) {
        n.obs.i <- date.vec[i + 1] - date.vec[i]

        vect1 <- rep(0, n.obs)

        t.low <- max(trunc(date.vec[i] + n.obs.i * trim - 1), max.lag + 2)
        t.high <- trunc(date.vec[i + 1] - n.obs.i * trim - 1)

        if (t.low < t.high - 1) {
            for (tb in t.low:t.high) {
                lambda <- (tb - 1) / (date.vec[i + 1] - 1)

                DU <- as.numeric(x.trend > tb)
                DT <- DU * (x.trend - tb)

                x <- cbind(
                    x.const,
                    if (const) DU else NULL,
                    x.trend - date.vec[i] + 1,
                    DT
                )

                y.i <- y[date.vec[i]:(date.vec[i + 1] - 1), , drop = FALSE]
                x.i <- x[date.vec[i]:(date.vec[i + 1] - 1), , drop = FALSE]

                k.hat <- max(1, .estimate_ar(y.i, x.i, max.lag, criterion)$lag)

                resids <- .estimate_ols(y.i, x.i)$residuals

                d.resid <- as.matrix(c(0, diff(resids)))

                y.u <- resids[k.hat:nrow(resids), , drop = FALSE]
                x.u <- lagn(resids, 1, na = 0)
                if (k.hat > 1) {
                    for (l in 1:(k.hat - 1)) {
                        x.u <- cbind(x.u, lagn(d.resid, l, na = 0))
                    }
                }
                x.u <- x.u[k.hat:nrow(x.u), , drop = FALSE]

                tmp.OLS <- .estimate_ols(y.u, x.u)
                beta.u <- tmp.OLS$beta
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

                c1 <- sqrt((1 + k.x) * n.obs.i)
                c2 <- ((1 + k.x) * n.obs.i - tau05^2 * (IP + n.obs.i)) /
                    (tau05 * (tau05 + k) * (IP + n.obs.i))

                if (tau > tau05)
                    c.tau <- -tau
                if (tau <= tau05 && tau > -k) {
                    c.tau <- IP * tau / n.obs -
                        (k.x + 1) / (tau + c2 * (tau + k))
                }
                if (tau <= -k && tau > -c1) {
                    c.tau <- IP * tau / n.obs - (k.x + 1) / tau
                }
                if (tau <= -c1) {
                    c.tau <- 0
                }

                a.hat.M <- a.hat + c.tau * sqrt(var.a.hat)
                if (a.hat.M >= 1) {
                    a.hat.M <- 1
                } else if (abs(a.hat.M) < 1) {
                    a.hat.M <- a.hat.M
                } else {
                    a.hat.M <- -0.99
                }

                CR <- sqrt(n.obs) * abs(a.hat.M - 1)
                if (CR <= 1) a.hat.M <- 1

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

                tmp.OLS <- .estimate_ols(y.g, x.g)
                beta.g <- tmp.OLS$beta
                g.resid <- tmp.OLS$residuals
                rm(tmp.OLS)

                if (k.hat == 1) {
                    h0 <- drop(t(g.resid) %*% g.resid) / nrow(g.resid)
                } else {
                    if (a.hat.M == 1) {
                        x.v <- NULL
                        for (k.i in 1:(k.hat - 1))
                            x.v <- cbind(x.v, lagn(g.resid, k.i, na = 0))

                        y.v <- g.resid[(k.hat - 1):nrow(g.resid)]
                        x.v <- x.v[(k.hat - 1):nrow(g.resid), , drop = FALSE]

                        tmp.OLS <- .estimate_ols(y.v, x.v)
                        beta.v <- tmp.OLS$beta
                        v.resid <- tmp.OLS$residuals
                        rm(tmp.OLS)

                        if (!const) {
                            h0 <- (drop(t(v.resid) %*% v.resid) / (n.obs.i - k.hat)) / # nolint
                                ((1 - sum(beta.v))^2)
                        }
                        if (const) {
                            BETAS <- matrix(0, nrow = k.hat - 1, ncol = 4)
                            for (k.i in 1:(k.hat - 1)) {
                                DU.ki <- as.numeric(x.trend > tb - k.i)
                                DT.ki <- DU.ki * (x.trend - tb)
                                x.ki <- cbind(
                                    x.const,
                                    DU.ki,
                                    x.trend,
                                    DT.ki
                                )
                                x.g.ki <- rbind(
                                    x.ki[date.vec[i] + 1, ],
                                    x.ki[(date.vec[i] + 2):(date.vec[i + 1] - 1), ] - # nolint
                                    a.hat.M * x.ki[(date.vec[i] + 1):(date.vec[i + 1] - 2), ] # nolint
                                )
                                beta.ki <- .estimate_ols(y.g, x.g.ki)$beta
                                BETAS[k.i, ] <- drop(beta.ki)
                                sig.e <- drop(t(v.resid) %*% v.resid) / (n.obs.i - k.hat) # nolint
                                h0 <- sig.e / ((1 - sum(beta.v))^2)
                                beta.g[2] <- (sqrt(h0) / sqrt(sig.e)) *
                                    (beta.g[2] - drop(t(BETAS[, 2]) %*% beta.v))
                            }
                        }
                    }

                    if (abs(a.hat.M) < 1)
                        h0 <- lr.var.quadratic(g.resid)
                }

                VCV <- h0 * qr.solve(t(x.g) %*% x.g)
                vect1[tb] <- t(R %*% beta.g) %*%
                    qr.solve(R %*% VCV %*% t(R)) %*% (R %*% beta.g)
            }

            vect1 <- vect1[t.low:t.high]
            wald[i] <- log(sum(exp(vect1 / 2)) / (date.vec[i + 1] - date.vec[i])) # nolint
        }
    }

    return(max(wald))
}

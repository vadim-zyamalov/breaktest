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
        stop("ERROR! Illegal trim value")
    }

    if (const && !trend) {
        VR <- matrix(c(0, 1, 0), nrow = 1, ncol = 3, byrow = TRUE)
        v.t <- as.matrix(c(
            -4.30, -4.39, -4.39, -4.34, -4.32,
            -4.45, -4.42, -4.33, -4.27, -4.27
        ))
        c.v <- matrix(c(
            1.60, 2.07, 3.33,
            1.52, 1.97, 3.24,
            1.41, 1.88, 3.05,
            1.26, 1.74, 3.12,
            0.91, 1.33, 2.83
        ), ncol = 3, byrow = TRUE)
    } else if (!const && trend) {
        VR <- matrix(c(0, 0, 1), nrow = 1, ncol = 3, byrow = TRUE)
        v.t <- as.matrix(c(
            -4.27, -4.41, -4.51, -4.55, -4.56,
            -4.57, -4.51, -4.38, -4.26, -4.26
        ))
        c.v <- matrix(c(
            1.52, 2.02, 3.37,
            1.40, 1.93, 3.27,
            1.28, 1.86, 3.20,
            1.13, 1.67, 3.06,
            0.74, 1.28, 2.61
        ), ncol = 3, byrow = TRUE)
    } else if (const && trend) {
        VR <- matrix(c(
            0, 1, 0, 0,
            0, 0, 0, 1
        ), nrow = 2, ncol = 4, byrow = TRUE)
        v.t <- as.matrix(c(
            -4.38, -4.65, -4.78, -4.81, -4.90,
            -4.88, -4.75, -4.70, -4.41, -4.41
        ))
        c.v <- matrix(c(
            2.96, 3.55, 5.02,
            2.82, 3.36, 4.78,
            2.65, 3.16, 4.59,
            2.48, 3.12, 4.47,
            2.15, 2.79, 4.57
        ), ncol = 3, byrow = TRUE)
    } else {
        stop("ERROR! Unknown model")
    }

    n.obs <- nrow(y)
    x.const <- rep(1, n.obs)
    x.trend <- 1:n.obs

    first.break <- max(trunc(trim * n.obs), max.lag + 2) + 1
    last.break <- trunc((1 - trim) * n.obs) + 1

    vect1 <- matrix(0, nrow = trunc((1 - 2 * trim) * n.obs) + 2, ncol = 1)

    for (tb in first.break:last.break) {
        lambda <- tb / n.obs

        DU <- as.numeric(x.trend > tb)
        DT <- DU * (x.trend - tb)

        x <- cbind(
            x.const,
            x.trend,
            if (const) DU else NULL,
            if (trend) DT else NULL
        )

        k.hat <- max(1, .estimate_ar(y, x, max.lag, criterion)$lag)

        resids <- .estimate_ols(y, x)$residuals

        d.resid <- c(0, diff(resids))

        y.u <- resids[k.hat:n.obs, , drop = FALSE]
        x.u <- lagn(resids, 1, na = 0)
        if (k.hat > 1) {
            for (l in 1:(k.hat - 1)) {
                x.u <- cbind(x.u, lagn(d.resid, l, na = 0))
            }
        }
        x.u <- x.u[k.hat:n.obs, , drop = FALSE]

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

        c1 <- sqrt((1 + k.x) * n.obs)
        c2 <- ((1 + k.x) * n.obs - tau05^2 * (IP + n.obs)) /
            (tau05 * (tau05 + k) * (IP + n.obs))

        if (tau > tau05)
            c.tau <- -tau
        else if (tau <= tau05 && tau > -k)
            c.tau <- IP * tau / n.obs - (k.x + 1) / (tau + c2 * (tau + k))
        else if (tau <= -k && tau > -c1)
            c.tau <- IP * tau / n.obs - (k.x + 1) / tau
        else if (tau <= -c1)
            c.tau <- 0

        a.hat.M <- a.hat + c.tau * sqrt(var.a.hat)
        if (a.hat.M >= 1)
            a.hat.M <- 1
        else if (abs(a.hat.M) < 1)
            a.hat.M <- a.hat.M
        else
            a.hat.M <- -0.99

        CR <- sqrt(n.obs) * abs(a.hat.M - 1)
        if (CR <= 1) a.hat.M <- 1

        y.g <- rbind(
            y[1, , drop = FALSE],
            y[2:n.obs, , drop = FALSE] -
                a.hat.M * y[1:(n.obs - 1), , drop = FALSE]
        )
        x.g <- rbind(
            x[1, , drop = FALSE],
            x[2:n.obs, , drop = FALSE] -
                a.hat.M * x[1:(n.obs - 1), , drop = FALSE]
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

                y.v <- g.resid[(k.hat - 1):nrow(g.resid), , drop = FALSE]
                x.v <- x.v[(k.hat - 1):nrow(g.resid), , drop = FALSE]

                tmp.OLS <- .estimate_ols(y.v, x.v)
                beta.v <- tmp.OLS$beta
                v.resid <- tmp.OLS$residuals
                rm(tmp.OLS)

                if (const && !trend) {
                    BETAS <- matrix(0, nrow = k.hat - 1, ncol = 3)
                    for (k.i in 1:(k.hat - 1)) {
                        DU.ki <- as.numeric(x.trend > tb - k.i)
                        x.ki <- cbind(
                            x.const,
                            DU.ki,
                            x.trend
                        )
                        x.g.ki <- rbind(
                            x.ki[1, ],
                            x.ki[2:n.obs, ] - a.hat.M * x.ki[1:(n.obs - 1), ]
                        )
                        beta.ki <- .estimate_ols(y.g, x.g.ki)$beta
                        BETAS[k.i, ] <- drop(beta.ki)
                    }
                    beta.g[2] <- beta.g[2] -
                        drop(t(BETAS[, 2]) %*% beta.v)
                    h0 <- drop(t(v.resid) %*% v.resid) / (n.obs - k.hat)
                } else if (!const && trend) {
                    h0 <- (drop(t(v.resid) %*% v.resid) / (n.obs - k.hat)) /
                        ((1 - sum(beta.v))^2)
                } else {
                    BETAS <- matrix(0, nrow = k.hat - 1, ncol = 4)
                    for (k.i in 1:(k.hat - 1)) {
                        DU.ki <- as.numeric(x.trend > tb - k.i)
                        DT.ki <- DU.ki * (x.trend - tb)
                        x.ki <- cbind(
                            rep(1, n.obs),
                            DU.ki,
                            1:n.obs,
                            DT.ki
                        )
                        x.g.ki <- rbind(
                            x.ki[1, ],
                            x.ki[2:n.obs, ] - a.hat.M * x.ki[1:(n.obs - 1), ]
                        )
                        beta.ki <- .estimate_ols(y.g, x.g.ki)$beta
                        BETAS[k.i, ] <- drop(beta.ki)
                        sig.e <- drop(t(v.resid) %*% v.resid) / (n.obs - k.hat)
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
        vect1[tb - first.break + 1] <- t(VR %*% beta.g) %*%
            qr.solve(VR %*% VCV %*% t(VR)) %*% (VR %*% beta.g)
    }

    wald <- log(sum(exp(vect1 / 2)) / n.obs)
    if (trim == 0.01) cv <- c.v[1, ]
    if (trim == 0.05) cv <- c.v[2, ]
    if (trim == 0.10) cv <- c.v[3, ]
    if (trim == 0.15) cv <- c.v[4, ]
    if (trim == 0.25) cv <- c.v[5, ]

    return(
        list(
            wald = wald,
            critical.value = cv
        )
    )
}

#' @title
#' A set of residual based tests for cointegration
#'
#' @details
#' The code provided is the original GAUSS code by Perron and Rodríguez
#' ported to R.
#'
#' @param y,x Variables of interest. `x` can be a matrix of several variables.
#' @param deter A value equal to
#' * 1: quasi-demeaned y and x,
#' * 2: quasi-detrended y and x,
#' * 3: quasi-demeaned y and quasi-detrended x.
#' @param min.lag A minimum number of lags to be used in the tests.
#'
#' @return A list of:
#' * 7x1-matrix of test statistics values,
#' * estimated number of lags.
#'
#' @references
#' Perron, Pierre, and Gabriel Rodríguez.
#' “Residuals‐based Tests for Cointegration with Generalized Least‐squares
#' Detrended Data.”
#' The Econometrics Journal 19, no. 1 (February 1, 2016): 84–111.
#' https://doi.org/10.1111/ectj.12056.
#'
#' @export
coint.test.PR <- function(y,
                          x,
                          deter,
                          min.lag = 0) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.matrix(x)) x <- as.matrix(x)

    n.obs <- nrow(y)
    n.var <- ncol(x)

    if (n.var == 1) {
        options <- c(-13.75, -20.50, -13.50)
    } else if (n.var == 2) {
        options <- c(-18.25, -23.75, -18.00)
    } else if (n.var == 3) {
        options <- c(-22.25, -27.25, -23.00)
    } else if (n.var == 4) {
        options <- c(-26.25, -30.75, -26.00)
    } else {
        options <- c(-30.00, -33.75, -29.75)
    }

    opt.cbar <- options[deter]

    max.lag <- round(4 * (n.obs / 100)^(1 / 4))

    zy <- if (deter == 1 || deter == 3) {
        as.matrix(rep(1, n.obs))
    } else if (deter == 2) {
        cbind(rep(1, n.obs), (1:n.obs))
    }

    zx <- if (deter == 1) {
        as.matrix(rep(1, n.obs))
    } else if (deter == 2 || deter == 3) {
        cbind(rep(1, n.obs), (1:n.obs))
    }

    y.d <- .estimate_gls(y, zy, opt.cbar)$residuals
    x.d <- .estimate_gls(x, zx, opt.cbar)$residuals

    model <- .estimate_ols(y.d, x.d)
    ud.hat <- cbind(model$residuals)

    result <- resid.tests.PR(ud.hat, min.lag, max.lag, opt.cbar, deter)

    return(result)
}


#' @title
#' Internal procedure for calculating test statistics from
#' Perron-Rodriguez (2016)
#'
#' @param ud A vector of residuals for testing.
#' @param min.lag,max.lag Minimum and maximum lag number.
#' @param c.bar A `c` parameter used for GLS detrending purposes.
#' @param deter A value equal to
#' * 1: quasi-demeaned y and x,
#' * 2: quasi-detrended y and x,
#' * 3: quasi-demeaned y and quasi-detrended x.
#'
#' @return A list of:
#' * 7x1-matrix of test statistics values,
#' * estimated number of lags.
#'
#' @keywords internal
resid.tests.PR <- function(ud,
                           min.lag,
                           max.lag,
                           c.bar,
                           deter) {
    n.obs <- nrow(ud)

    if (ncol(ud) > 1) {
        stop("ERROR:")
    }

    gls.tests <- matrix(0, 7, 1)

    lag.ud <- ud[1:(n.obs - 1), 1, drop = FALSE]
    d.ud <- diffn(ud, na = 0)
    sum.ud.sq <- drop(t(lag.ud) %*% lag.ud)

    model.1 <- .estimate_ols(
        ud[2:n.obs, 1, drop = FALSE],
        ud[1:(n.obs - 1), 1, drop = FALSE]
    )

    rho.hat <- as.matrix(model.1$beta)
    omega <- as.matrix(model.1$residuals)
    s2.ud <- c(t(omega) %*% omega) / (nrow(omega) - 1)
    t.rho <- (rho.hat - 1) / sqrt(s2.ud / sum.ud.sq)

    fin.bic <- Inf
    fin.lag <- min.lag
    lag.bic <- min.lag
    while (lag.bic <= max.lag) {
        tmp.reg <- lagn(ud, 1, na = 0)

        h <- 1
        while (h <= lag.bic) {
            tmp.reg <- cbind(
                tmp.reg,
                lagn(d.ud, h, na = 0)
            )
            h <- h + 1
        }

        tmp.reg <-
            tmp.reg[(lag.bic + 2):nrow(tmp.reg), , drop = FALSE]

        model.2 <- .estimate_ols(
            d.ud[(lag.bic + 2):nrow(d.ud), , drop = FALSE],
            tmp.reg
        )

        eta <- cbind(model.2$residuals)
        s2.eta <- c(t(eta) %*% eta) / (nrow(eta) - ncol(tmp.reg))
        xtx.inv <- solve(t(tmp.reg) %*% tmp.reg)

        if (lag.bic == 0) {
            sumb <- 0
        } else {
            sumb <- sum(model.2$beta[2:(lag.bic + 1)])
        }

        s2.adj <- s2.eta / ((1 - sumb)^2)

        cur.bic <- log(c(t(eta) %*% eta) / (n.obs - max.lag)) +
            log(n.obs - max.lag) * lag.bic / (n.obs - max.lag)
        if (cur.bic < fin.bic) {
            fin.bic <- cur.bic
            fin.lag <- lag.bic

            gls.tests[1, 1] <- (ud[n.obs, 1]^2 / n.obs - s2.adj) /
                (2 * sum.ud.sq / n.obs^2)
            gls.tests[2, 1] <- sqrt(2 * sum.ud.sq / (n.obs^2 * s2.adj))
            gls.tests[3, 1] <- gls.tests[1, 1] * gls.tests[2, 1]
            gls.tests[4, 1] <- model.2$beta[1] / sqrt(s2.eta * xtx.inv[1, 1])
            gls.tests[5, 1] <- (n.obs - 1) * (rho.hat - 1) -
                (s2.adj - s2.ud) / (2 * sum.ud.sq / n.obs^2)
            gls.tests[6, 1] <- sqrt(s2.ud / s2.adj) * t.rho -
                (s2.adj - s2.ud) / sqrt(4 * s2.adj * sum.ud.sq / n.obs^2)
        }

        lag.bic <- lag.bic + 1
    }

    if (deter == 1 || deter == 3) {
        gls.tests[7, 1] <- (c.bar^2 * sum.ud.sq / n.obs^2 -
            c.bar * ud[n.obs, 1]^2 / n.obs) / s2.adj
    } else if (deter == 2) {
        gls.tests[7, 1] <- (c.bar^2 * sum.ud.sq / n.obs^2 +
            (1 - c.bar) * ud[n.obs, 1]^2 / n.obs) / s2.adj
    } else {
        stop("ERROR! Unknown `det.comp` value")
    }

    rownames(gls.tests) <- c(
        "MZ(rho)",
        "MSB",
        "MZ(t.rho)",
        "ADF",
        "Z(rho)",
        "Z(t.rho)",
        if (deter == 1 || deter == 3) {
            "MP(T, demeaned)"
        } else if (deter == 2) {
            "MP(T, detrended)"
        }
    )

    result <- list()
    result$gls.tests <- gls.tests
    result$lag <- fin.lag

    return(result)
}

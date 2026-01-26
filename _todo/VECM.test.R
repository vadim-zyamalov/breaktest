#' @title
#' Test of the cointegration rank with a possible break in trend
#'
#' @description
#' This procedure is aimed on the problem of testing for the cointegration rank
#' of a vector autoregressive process in the case where a trend break may
#' potentially be present in the data.
#'
#' The test is based on estimating the quasi log likelihood for two situations,
#' with break, and without it. The one with the smallest value is considered to
#' be the result.
#'
#' @details
#' The code provided is the original GAUSS code by Harris et al.
#' ported to R.
#'
#' @param y A matrix of \eqn{n} VAR variables.
#' @param r The cointegration rank tested against the alternative of \eqn{n}.
#' @param max.lag The maximum number of lags.
#' @param trim The trimming parameter to find the lower and upper bounds of
#' possible break dates.
#'
#' @return A list of:
#' * the indicator of the rejection of null.
#' * the estimated break point.
#' * the estimated lag number.
#'
#' @references
#' Harris, David, Stephen J. Leybourne, and A. M. Robert Taylor.
#' “Tests of the Co-Integration Rank in VAR Models in the Presence
#' of a Possible Break in Trend at an Unknown Point.”
#' Journal of Econometrics, Innovations in Multiple Time Series Analysis,
#' 192, no. 2 (June 1, 2016): 451–67.
#' https://doi.org/10.1016/j.jeconom.2016.02.010.
#'
#' @export
VECM.test <- function(y,
                      r,
                      max.lag,
                      trim = 0.15) {
    n.obs <- nrow(y)
    n.var <- ncol(y)

    if (r >= n.var) {
        stop("ERROR: Invalid values r or y")
    }

    first.break <- trunc(trim * n.obs)
    last.break <- trunc((1 - trim) * n.obs)
    breaks.list <- as.matrix(first.break:last.break)

    cv.br.4 <- matrix(
        c(
            17.440, 18.395, 18.753, 18.711, 18.778, 18.327, 17.379,
            34.623, 36.270, 37.423, 37.661, 37.541, 36.357, 34.546,
            55.826, 58.317, 59.084, 59.467, 59.285, 57.853, 55.507,
            80.611, 83.863, 84.718, 85.169, 84.704, 83.404, 80.790
        ),
        nrow = 4, ncol = 7, byrow = TRUE
    )
    cv.trend.4 <- c(12.511, 25.984, 43.195, 64.240)

    loglp <- NULL
    loglbp <- NULL

    for (i in 1:max.lag) {
        loglp <- cbind(loglp, VECM.logl(y, i))
        loglbp <- cbind(loglbp, VECM.break.logl(y, i, breaks.list))
    }

    bhatrp.index <- t(apply(loglbp, 2, which.min))

    rc <- 0:(max.lag - 1) * (n.var + 1) + (r + 1)
    nc <- 0:(max.lag - 1) * (n.var + 1) + (n.var + 1)

    p <- which.min(
        -2 * t(loglp[1, nc]) + log(n.obs) * (n.var^2 * (1:max.lag))
    )
    pbr <- which.min(
        -2 * diag(loglbp[bhatrp.index[1, rc], nc]) +
        log(n.obs) * (n.var^2 * (1:max.lag))
    )

    brphat.index <- bhatrp.index[rc[pbr]]

    tr0.phat <- ifelse(
        2 * (loglp[nc[p]] - loglp[rc[p]]) > cv.trend.4[n.var - r],
        1, 0
    )
    tr1.VECM.phat <- ifelse(
        2 * (loglbp[brphat.index, nc[pbr]] - loglbp[brphat.index, rc[pbr]]) >
        cv.br.4[n.var - r, round(10 * breaks.list[brphat.index] / n.obs) - 1],
        1, 0
    )

    SC1r <- -2 * loglbp[brphat.index, rc[pbr]] +
        (n.var + r + 2 + (n.var^2) * pbr) * log(n.obs)
    SC0r <- -2 * loglp[rc[p]] + (n.var^2) * p * log(n.obs)

    result <- list()

    if (SC1r < SC0r) {
        result$tr1.SCVECM.phat <- tr1.VECM.phat
        result$b.SCVECM <- breaks.list[brphat.index]
        result$p.SCVECM <- pbr
    } else {
        result$tr1.SCVECM.phat <- tr0.phat
        result$b.SCVECM <- 0
        result$p.SCVECM <- p
    }

    return(result)
}


#' @title
#' Quasi log likelihood for VECM without and with breaks
#' @order 1
#'
#' @param y A matrix of \eqn{n} VAR variables.
#' @param p Number of lags.
#' @param breaks.list Vector of possible break dates.
#'
#' @return The vector or matrix with the values of quasi log likelihood
#' for all possible values of cointegration rank and different break dates,
#' if the function with breaks is called.
#'
#' @importFrom Rfast spdinv
#'
#' @keywords internal
VECM.logl <- function(y,
                      p) {
    if (!is.matrix(y)) y <- as.matrix(y)

    n.obs <- nrow(y)

    trend <- matrix(1:n.obs, n.obs, 1)

    d.y <- diff(y)

    z0 <- d.y[p:(n.obs - 1), , drop = FALSE]
    z1 <- cbind(
        y[p:(n.obs - 1), , drop = FALSE],
        trend[(p + 1):n.obs, , drop = FALSE]
    )

    Xp <- matrix(1, (n.obs - p), 1)
    j <- 1
    while (j <= (p - 1)) {
        Xp <- cbind(
            Xp,
            d.y[(p - j):(n.obs - 1 - j), ]
        )
        j <- j + 1
    }

    r0 <- z0 - Xp %*% solve(t(Xp) %*% Xp) %*% t(Xp) %*% z0
    r1 <- z1 - Xp %*% solve(t(Xp) %*% Xp) %*% t(Xp) %*% z1
    Li <- solve(t(chol(t(r1) %*% r1)))

    lam <- cbind(
        c(
            0,
            rev(sort(eigen(
                Li %*% t(r1) %*% r0 %*%
                Rfast::spdinv(t(r0) %*% r0) %*%
                t(r0) %*% r1 %*% t(Li)
            )$values))[2:nrow(Li)]
        )
    )
    logL <- -(n.obs - p) / 2 * log(det((t(r0) %*% r0) / (n.obs - p))) +
        cumsum(log(1 - lam))

    rownames(logL) <- NULL
    colnames(logL) <- NULL
    return(t(logL))
}


#' @rdname VECM.logl
#' @order 2
#' @importFrom Rfast spdinv
#' @keywords internal
VECM.break.logl <- function(y,
                            p,
                            breaks.list) {
    if (!is.matrix(y)) y <- as.matrix(y)
    n.obs <- nrow(y)

    trend <- matrix(1:n.obs, n.obs, 1)

    d.y <- diff(y)
    z0 <- as.matrix(d.y[p:(n.obs - 1), ])

    logL <- NULL

    bc <- 1
    while (bc <= nrow(breaks.list)) {
        b <- breaks.list[bc, 1]

        D <- NULL
        j <- 1
        while (j <= p) {
            D <- cbind(D, ifelse(trend == (b + j), 1, 0))
            j <- j + 1
        }

        E1 <- ifelse(trend <= b, 1, 0)
        E2 <- ifelse(trend > b, 1, 0)
        tE <- apply(cbind(E1, E2), 2, cumsum)

        z1 <- cbind(y[p:(n.obs - 1), ], tE[(p + 1):nrow(tE), ])
        Xp <- cbind(E1, E2, D)[(p + 1):nrow(E1), ]

        j <- 1
        while (j <= (p - 1)) {
            Xp <- cbind(Xp, d.y[(p - j):(n.obs - 1 - j)])
            j <- j + 1
        }

        r0 <- z0 - Xp %*% solve(t(Xp) %*% Xp) %*% t(Xp) %*% z0
        r1 <- z1 - Xp %*% solve(t(Xp) %*% Xp) %*% t(Xp) %*% z1
        Li <- solve(t(chol(t(r1) %*% r1)))
        lam <- cbind(
            c(
                0,
                rev(sort(eigen(
                    Li %*% t(r1) %*% r0 %*%
                    Rfast::spdinv(t(r0) %*% r0) %*%
                    t(r0) %*% r1 %*% t(Li)
                )$values))[2:nrow(Li)]
            )
        )
        logLb <- -(n.obs - p) / 2 * log(det((t(r0) %*% r0) / (n.obs - p))) +
            cumsum(log(1 - lam))

        logL <- rbind(logL, logLb)

        bc <- bc + 1
    }

    rownames(logL) <- NULL
    colnames(logL) <- NULL

    return(logL)
}

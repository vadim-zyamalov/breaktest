#' @order 1
#' @title
#' Calculating long-run variance or covariance matrix
#'
#' @details
#' The code provided is based on the original code by Kurozumi, Sul et al.
#' ported to R.
#'
#' @param y A time series of interest.
#' @param demean Whether the demeaning is needed.
#' @param kernel A kernel to be used:
#' * `truncated`: \eqn{\left\{\begin{array}{ll}
#' 1 & |x| \leq 1 \\
#' 0 & \textrm{otherwize}
#' \end{array}\right.}
#' * `bartlett`: \eqn{\left\{\begin{array}{ll}
#' 1 - |x| & |x| \leq 1 \\
#' 0 & \textrm{otherwize}
#' \end{array}\right.}
#' * `parzen`: \eqn{\left\{\begin{array}{ll}
#' 1 - 6 x^2 + 6 {|x|}^3 & |x| \leq 1/2 \\
#' 2 (1 - |x|)^3 & 1/2 \leq |x| \leq 1 \\
#' 0 & \textrm{otherwize}
#' \end{array}\right.}
#' * `tukey-hanning`: \eqn{\left\{\begin{array}{ll}
#' (1 + \cos(\pi x))/2 & |x| \leq 1 \\
#' 0 & \textrm{otherwize}
#' \end{array}\right.}
#' * `quadratic`: \eqn{
#' \frac{25}{12 \pi^2 x^2}
#' \left(\frac{\sin(6 \pi x / 5)}{6 \pi x / 5} - \cos(6 \pi x / 5)\right)}
#' @param limit.lags Whether all lags shoult be used in formulae.
#' @param limit.selector Way of limit selection:
#' * `kpss-q`: \eqn{4 (T / 100)^{1 / 4}}.
#' * `kpss-m`: \eqn{12 (T / 100)^{1 / 4}}.
#' * `Andrews`: kernel-specific formula from Andrews (1991).
#' * `Kurozumi`: kernel-specific formula from Andrews (1991)
#' with Kurozumi (2002) proposal.
#' @param upper.rho.limit The upper limit for the value or AR-coefficient.
#' @param upper.lag.limit The value used to calculate the upper limit
#' for Kurozumi (2002) proposal.
#' @param recolor Whether the correction by Sul et al. (2005) should be used.
#' This option resets `limit.lags` to `TRUE`, and `limit.selector` to `Andrews`.
#' @param max.lag Maximum number of lags used in AR regresion during
#' recolorization. Otherwize ignored.
#' @param criterion The information crietreion: bic, aic or lwz.
#'
#' @references
#' Andrews, Donald W. K.
#' “Heteroskedasticity and Autocorrelation Consistent
#' Covariance Matrix Estimation.”
#' Econometrica 59, no. 3 (1991): 817–58.
#' https://doi.org/10.2307/2938229.
#'
#' Kurozumi, Eiji.
#' “Testing for Stationarity with a Break.”
#' Journal of Econometrics 108, no. 1 (May 1, 2002): 63–99.
#' https://doi.org/10.1016/S0304-4076(01)00106-3.
#'
#' Sul, Donggyu, Peter C. B. Phillips, and Chi-Young Choi.
#' “Prewhitening Bias in HAC Estimation.”
#' Oxford Bulletin of Economics and Statistics 67, no. 4 (August 2005): 517–46.
#' https://doi.org/10.1111/j.1468-0084.2005.00130.x.
#'
#' @importFrom stats na.omit
#'
#' @keywords internal
lr.var <- function(y,
                   demean = TRUE,
                   kernel = "bartlett",
                   limit.lags = FALSE,
                   limit.selector = "kpss-q",
                   upper.rho.limit = 0.97,
                   upper.lag.limit = 0.8,
                   recolor = FALSE,
                   max.lag = 0,
                   criterion = "bic") {
    if (!is.matrix(y)) y <- as.matrix(y)

    n.var <- ncol(y)

    if (!kernel %in% c("truncated",
                       "bartlett",
                       "parzen",
                       "tukey-hanning",
                       "quadratic")) {
        stop("ERROR! Unknown kernel")
    }
    if (!limit.selector %in% c("kpss-q", "kpss-m", "Andrews", "Kurozumi")) {
        stop("ERROR! Unknown limit selector")
    }
    if (limit.selector == "Kurozumi" && is.null(upper.lag.limit)) {
        stop("ERROR! Upper limit is needed for Kurozumi proposal")
    }
    if (limit.selector == "Kurozumi" && n.var > 1) {
        stop("ERROR! Kurozumi proposal is for a single variable case")
    }
    if (recolor && n.var > 1) {
        stop("ERROR! Recolorization is for a single variable case")
    }
    if (!criterion %in% c("bic", "aic", "lwz")) {
        stop("ERROR! Unknown criterion")
    }
    if (recolor) {
        limit.lags <- TRUE
        limit.selector <- "Andrews"
    }

    n.obs <- nrow(y)

    if (n.var == 1) {
        funcs <- lr.var.kernel(kernel, .lr.var.alpha.single, n.obs)
    } else {
        funcs <- lr.var.kernel(kernel, .lr.var.alpha.multi, n.obs)
    }

    if (recolor) {
        info.crit.min <- log(drop(t(y) %*% y) / (n.obs - max.lag))

        tmp.AR <- .estimate_ar(y, NULL, max.lag, criterion)
        info.crit <- info.criterion(tmp.AR$residuals, tmp.AR$lag)[[criterion]]

        if (info.crit.min < info.crit) {
            rho <- 0
            k <- 0
        } else {
            rho <- tmp.AR$beta
            k <- tmp.AR$lag
            y <- na.omit(tmp.AR$residuals)
            n.obs <- nrow(y)
        }
    } else {
        k <- 1
    }

    if (demean) {
        for (i in 1:n.var) {
            mean.y <- mean(y[, i])
            y[, i] <- y[, i] - mean.y
        }
    }

    if (!limit.lags) {
        limit <- n.obs - 1
    } else if (limit.selector == "kpss-q") {
        limit <- 4 * ((n.obs / 100)^(1 / 4))
    } else if (limit.selector == "kpss-m") {
        limit <- 12 * ((n.obs / 100)^(1 / 4))
    } else {
        if (k > 0) {
            limit <- funcs$limit(y, upper.rho.limit)

            if (limit.selector == "Kurozumi") {
                upper.lag.limit <- funcs$limit(upper.lag.limit, upper.rho.limit)
                limit <- min(limit, upper.lag.limit)
            }
        } else {
            limit <- 0
        }
    }
    limit <- trunc(limit)
    limit <- min(limit, n.obs - 1)

    lrv <- (t(y) %*% y) / n.obs
    if (k > 0) {
        for (i in 1:limit) {
            if (i == 0) break
            if (i < n.obs - 1) {
                lrv <- lrv + funcs$weight(i, limit) * (
                    (t(y[1:(n.obs - i), ]) %*% y[(1 + i):n.obs, ]) / n.obs +
                    t(t(y[1:(n.obs - i), ]) %*% y[(1 + i):n.obs, ]) / n.obs
                )
            } else {
                lrv <- lrv + funcs$weight(i, limit) * as.vector(
                    (t(y[1:(n.obs - i), ]) %*% y[(1 + i):n.obs, ]) / n.obs +
                    t(t(y[1:(n.obs - i), ]) %*% y[(1 + i):n.obs, ]) / n.obs
                )
            }
        }
    }

    if (recolor) {
        print(rho)
        lrv.recolored <- lrv / (1 - sum(rho))^2
        lrv <- min(lrv.recolored, n.obs * 0.15 * lrv)
    }

    return(drop(lrv))
}

#' @rdname lr.var
#' @order 2
lr.var.bartlett <- function(y) {
    return(lr.var(
        y,
        limit.lags = TRUE,
        limit.selector = "kpss-q"
    ))
}

#' @rdname lr.var
#' @order 3
lr.var.quadratic <- function(y) {
    return(lr.var(
        y,
        kernel = "quadratic",
        limit.lags = TRUE,
        limit.selector = "Andrews"
    ))
}

#' @rdname lr.var
#' @order 4
lr.var.bartlett.AK <- function(y) {
    return(lr.var(
        y,
        kernel = "bartlett",
        limit.lags = TRUE,
        limit.selector = "Kurozumi"
    ))
}

#' @rdname lr.var
#' @order 5
lr.var.SPC <- function(y,
                       max.lag = 0,
                       kernel = "bartlett",
                       criterion = "bic") {
    return(lr.var(
        y,
        max.lag = max.lag,
        kernel = kernel,
        criterion = criterion,
        recolor = TRUE
    ))
}


lr.var.kernel <- function(kernel,
                          alpha,
                          n.obs) {
    if (kernel == "truncated") {
        f.limit <- function(y, l) {
            return(0.6611 * (n.obs * alpha(y, l)$q2)^(1 / 5))
        }
        f.weight <- function(i, l) {
            if (abs(i / (l + 1)) <= 1) {
                return(1)
            }
            return(0)
        }
    } else if (kernel == "bartlett") {
        f.limit <- function(y, l) {
            return(1.1447 * (n.obs * alpha(y, l)$q1)^(1 / 3))
        }
        f.weight <- function(i, l) {
            x <- i / (l + 1)
            if (abs(x) <= 1) {
                return(1 - abs(x))
            }
            return(0)
        }
    } else if (kernel == "parzen") {
        f.limit <- function(y, l) {
            return(2.6614 * (n.obs * alpha(y, l)$q2)^(1 / 5))
        }
        f.weight <- function(i, l) {
            x <- i / (l + 1)
            if (abs(x) <= 0.5) {
                return(1 - 6 * x^2 + 6 * abs(x)^3)
            } else if (abs(x) <= 1) {
                return(2 * (1 - abs(x))^3)
            }
            return(0)
        }
    } else if (kernel == "tukey-hanning") {
        f.limit <- function(y, l) {
            return(1.7462 * (n.obs * alpha(y, l)$q2)^(1 / 5))
        }
        f.weight <- function(i, l) {
            x <- i / (l + 1)
            if (abs(x) <= 1) {
                return((1 + cos(pi * x)) / 2)
            }
            return(0)
        }
    } else if (kernel == "quadratic") {
        f.limit <- function(y, l) {
            return(1.3221 * (n.obs * alpha(y, l)$q2)^(1 / 5))
        }
        f.weight <- function(i, l) {
            x <- i / (l + 1)
            delta <- 6 * pi * x / 5
            return((3 / delta^2) * (sin(delta) / delta - cos(delta)))
        }
    }

    return(
        list(
            limit = f.limit,
            weight = f.weight
        )
    )
}

.lr.var.alpha.single <- function(y,
                                 upper.rho.limit) {
    n.obs <- nrow(y)

    if (!is.null(n.obs) && n.obs > 1) {
        r <- drop(
            (t(y[1:(n.obs - 1), 1]) %*% y[2:n.obs, 1]) /
            (t(y[1:(n.obs - 1), 1]) %*% y[1:(n.obs - 1), 1])
        )

        if (r > upper.rho.limit) {
            r <- upper.rho.limit
        } else if (r < -upper.rho.limit) {
            r <- -upper.rho.limit
        }
    } else {
        r <- y
    }

    return(
        list(
            q1 = 4 * r^2 / (1 + r)^2 / (1 - r)^2,
            q2 = 4 * r^2 / (1 - r)^4
        )
    )
}

.lr.var.alpha.multi <- function(y,
                                upper.rho.limit) {
    n.obs <- nrow(y)
    n.var <- ncol(y)

    nominator_1 <- 0
    nominator_2 <- 0
    denominator <- 0

    for (i in 1:n.var) {
        r <- (t(y[1:(n.obs - 1), i]) %*% y[2:n.obs, i]) /
            (t(y[1:(n.obs - 1), i]) %*% y[1:(n.obs - 1), i])
        r <- drop(r)

        if (r > upper.rho.limit) {
            r <- upper.rho.limit
        } else if (r < -upper.rho.limit) {
            r <- -upper.rho.limit
        }

        resids <- y[2:n.obs, i] - y[1:(n.obs - 1), i] * r
        s2 <- mean(resids^2)

        nominator_1 <- nominator_1 + 4 * r^2 * s2^2 / (1 - r)^6 / (1 + r)^2
        nominator_2 <- nominator_2 + 4 * r^2 * s2^2 / (1 - r)^8
        denominator <- denominator + s2^2 / (1 - r)^4
    }

    return(
        list(
            q1 = nominator_1 / denominator,
            q2 = nominator_2 / denominator
        )
    )
}

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
#' @param limit_lags Whether all lags shoult be used in formulae.
#' @param limit_selector Way of limit selection:
#' * `kpss-q`: \eqn{4 (T / 100)^{1 / 4}}.
#' * `kpss-m`: \eqn{12 (T / 100)^{1 / 4}}.
#' * `Andrews`: kernel-specific formula from Andrews (1991).
#' * `Kurozumi`: kernel-specific formula from Andrews (1991)
#' with Kurozumi (2002) proposal.
#' @param upper.rho.limit The upper limit for the value or AR-coefficient.
#' @param upper.lag.limit The value used to calculate the upper limit
#' for Kurozumi (2002) proposal.
#' @param recolor Whether the correction by Sul et al. (2005) should be used.
#' This option resets `limit_lags` to `TRUE`, and `limit_selector` to `Andrews`.
#' @param max_lag Maximum number of lags used in AR regresion during
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
.lr.variance <- function(
  y,
  demean = TRUE,
  kernel = "bartlett",
  limit.lags = FALSE,
  limit.selector = "kpss-q",
  rho.upper.bound = 0.97,
  lag.upper.bound = 0.8,
  recolor = FALSE,
  max.lag = 0,
  criterion = "bic"
) {
  if (!is.matrix(y)) y <- as.matrix(y)

  Nc <- ncol(y)

  if (!kernel %in% c(
    "truncated",
    "bartlett",
    "parzen",
    "tukey-hanning",
    "quadratic"
  )) {
    stop("ERROR! lr.variance: Unknown kernel")
  }
  if (!limit.selector %in% c("kpss-q", "kpss-m", "Andrews", "Kurozumi")) {
    stop("ERROR! lr.variance: Unknown limit selector")
  }
  if (limit.selector == "Kurozumi" && is.null(lag.upper.bound)) {
    stop("ERROR! lr.variance: Upper limit is needed for Kurozumi proposal")
  }
  if (limit.selector == "Kurozumi" && Nc > 1) {
    stop("ERROR! lr.variance: Kurozumi proposal is for a single variable case")
  }
  if (recolor && Nc > 1) {
    stop("ERROR! lr.variance: Recolorization is for a single variable case")
  }
  if (!criterion %in% c("bic", "aic", "lwz")) {
    stop("ERROR! lr.variance: Unknown criterion")
  }

  if (recolor) {
    limit.lags <- TRUE
    limit.selector <- "Andrews"
  }

  N <- nrow(y)

  funcs <- if (Nc == 1) {
    .lr.var.kernel(kernel, .lr.alpha.single, N)
  } else {
    .lr.var.kernel(kernel, .lr.alpha.multi, N)
  }

  if (recolor) {
    min.ic <- log(drop(t(y) %*% y) / (N - max.lag))

    model.est <- .AR(y, NULL, max.lag, criterion)
    ic.values <- .ic.values(model.est$residuals, model.est$lag)[[criterion]]

    if (min.ic < ic.values) {
      rho <- 0
      k <- 0
    } else {
      rho <- model.est$beta
      k <- model.est$lag
      y <- na.omit(model.est$residuals)
      N <- nrow(y)
    }
  } else {
    k <- 1
  }

  if (demean) sweep(y, 2, colMeans(y), FUN = "-")

  if (!limit.lags) {
    limit <- N - 1
  } else if (limit.selector == "kpss-q") {
    limit <- 4 * ((N / 100)^(1 / 4))
  } else if (limit.selector == "kpss-m") {
    limit <- 12 * ((N / 100)^(1 / 4))
  } else {
    if (k > 0) {
      limit <- funcs$limit(y, rho.upper.bound)

      if (limit.selector == "Kurozumi") {
        lag.upper.bound <- funcs$limit(lag.upper.bound, rho.upper.bound)
        limit <- min(limit, lag.upper.bound)
      }
    } else {
      limit <- 0
    }
  }
  limit <- min(trunc(limit), N - 1)

  lrv <- (t(y) %*% y) / N
  if (k > 0) {
    for (i in 1:limit) {
      if (i == 0) break
      if (i < N - 1) {
        lrv <- lrv + funcs$weight(i, limit) * (
          (t(y[1:(N - i), ]) %*% y[(1 + i):N, ]) / N +
            t(t(y[1:(N - i), ]) %*% y[(1 + i):N, ]) / N
        )
      } else {
        lrv <- lrv + funcs$weight(i, limit) * as.vector(
          (t(y[1:(N - i), ]) %*% y[(1 + i):N, ]) / N +
            t(t(y[1:(N - i), ]) %*% y[(1 + i):N, ]) / N
        )
      }
    }
  }

  if (recolor) {
    lrv.recolored <- lrv / (1 - sum(rho))^2
    lrv <- min(lrv.recolored, N * 0.15 * lrv)
  }

  drop(lrv)
}

#' @rdname dot-lr.variance
#' @order 2
.lr.var.bartlett <- function(y) {
  .lr.variance(
    y,
    limit.lags = TRUE,
    limit.selector = "kpss-q"
  )
}

#' @rdname dot-lr.variance
#' @order 3
.lr.var.quad <- function(y) {
  .lr.variance(
    y,
    kernel = "quadratic",
    limit.lags = TRUE,
    limit.selector = "Andrews"
  )
}

#' @rdname dot-lr.variance
#' @order 4
.lr.var.kurozumi <- function(y) {
  .lr.variance(
    y,
    kernel = "bartlett",
    limit.lags = TRUE,
    limit.selector = "Kurozumi"
  )
}

#' @rdname dot-lr.variance
#' @order 5
.lr.var.spc <- function(
  y,
  max.lag = 0,
  kernel = "bartlett",
  criterion = "bic"
) {
  .lr.variance(
    y,
    max.lag = max.lag,
    kernel = kernel,
    criterion = criterion,
    recolor = TRUE
  )
}


.lr.var.kernel <- function(
  kernel,
  alpha,
  N
) {
  limit.func <- switch(kernel,
    truncated = function(y, l) {
      0.6611 * (N * alpha(y, l)$q2)^(1 / 5)
    },
    bartlett = function(y, l) {
      1.1447 * (N * alpha(y, l)$q1)^(1 / 3)
    },
    parzen = function(y, l) {
      2.6614 * (N * alpha(y, l)$q2)^(1 / 5)
    },
    "tukey-hanning" = function(y, l) {
      1.7462 * (N * alpha(y, l)$q2)^(1 / 5)
    },
    quadratic = function(y, l) {
      1.3221 * (N * alpha(y, l)$q2)^(1 / 5)
    },
    stop("Unknown kernel!")
  )

  weight.func <- switch(kernel,
    truncated = function(i, l) {
      if (abs(i / (l + 1)) <= 1) {
        return(1)
      }
      0
    },
    bartlett = function(i, l) {
      x <- i / (l + 1)
      if (abs(x) <= 1) {
        return(1 - abs(x))
      }
      0
    },
    parzen = function(i, l) {
      x <- i / (l + 1)
      if (abs(x) <= 0.5) {
        return(1 - 6 * x^2 + 6 * abs(x)^3)
      } else if (abs(x) <= 1) {
        return(2 * (1 - abs(x))^3)
      }
      0
    },
    "tukey-hanning" = function(i, l) {
      x <- i / (l + 1)
      if (abs(x) <= 1) {
        return((1 + cos(pi * x)) / 2)
      }
      0
    },
    quadratic = function(i, l) {
      x <- i / (l + 1)
      delta <- 6 * pi * x / 5
      (3 / delta^2) * (sin(delta) / delta - cos(delta))
    },
    stop("Unknown kernel!")
  )

  list(
    limit  = limit.func,
    weight = weight.func
  )
}

.lr.alpha.single <- function(y, rho.upper.bound) {
  N <- nrow(y)

  if (!is.null(N) && N > 1) {
    r <- drop(
      (t(y[1:(N - 1), 1]) %*% y[2:N, 1]) /
        (t(y[1:(N - 1), 1]) %*% y[1:(N - 1), 1])
    )

    if (r > rho.upper.bound) {
      r <- rho.upper.bound
    } else if (r < -rho.upper.bound) {
      r <- -rho.upper.bound
    }
  } else {
    r <- y
  }

  list(
    q1 = 4 * r^2 / (1 + r)^2 / (1 - r)^2,
    q2 = 4 * r^2 / (1 - r)^4
  )
}

.lr.alpha.multi <- function(y, rho.upper.bound) {
  N <- nrow(y)
  Nc <- ncol(y)

  nom.1 <- 0
  nom.2 <- 0
  denom <- 0

  for (i in 1:Nc) {
    r <- (t(y[1:(N - 1), i]) %*% y[2:N, i]) /
      (t(y[1:(N - 1), i]) %*% y[1:(N - 1), i])
    r <- drop(r)

    if (r > rho.upper.bound) {
      r <- rho.upper.bound
    } else if (r < -rho.upper.bound) {
      r <- -rho.upper.bound
    }

    resids <- y[2:N, i] - y[1:(N - 1), i] * r
    s2 <- mean(resids^2)

    nom.1 <- nom.1 + 4 * r^2 * s2^2 / (1 - r)^6 / (1 + r)^2
    nom.2 <- nom.2 + 4 * r^2 * s2^2 / (1 - r)^8
    denom <- denom + s2^2 / (1 - r)^4
  }

  list(
    q1 = nom.1 / denom,
    q2 = nom.2 / denom
  )
}

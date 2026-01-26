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
.variance_lr <- function(
  y,
  demean = TRUE,
  kernel = "bartlett",
  limit_lags = FALSE,
  limit_selector = "kpss-q",
  rho_upper_bound = 0.97,
  lag_upper_bound = 0.8,
  recolor = FALSE,
  max_lag = 0,
  criterion = "bic"
) {
  if (!is.matrix(y)) y <- as.matrix(y)

  n_var <- ncol(y)

  if (!kernel %in% c(
    "truncated",
    "bartlett",
    "parzen",
    "tukey-hanning",
    "quadratic"
  )) {
    stop("ERROR! Unknown kernel")
  }
  if (!limit_selector %in% c("kpss-q", "kpss-m", "Andrews", "Kurozumi")) {
    stop("ERROR! Unknown limit selector")
  }
  if (limit_selector == "Kurozumi" && is.null(lag_upper_bound)) {
    stop("ERROR! Upper limit is needed for Kurozumi proposal")
  }
  if (limit_selector == "Kurozumi" && n_var > 1) {
    stop("ERROR! Kurozumi proposal is for a single variable case")
  }
  if (recolor && n_var > 1) {
    stop("ERROR! Recolorization is for a single variable case")
  }
  if (!criterion %in% c("bic", "aic", "lwz")) {
    stop("ERROR! Unknown criterion")
  }

  if (recolor) {
    limit_lags <- TRUE
    limit_selector <- "Andrews"
  }

  n_obs <- nrow(y)

  if (n_var == 1) {
    funcs <- .kernel_lr(kernel, .alpha_lr_single, n_obs)
  } else {
    funcs <- .kernel_lr(kernel, .alpha_lr_multi, n_obs)
  }

  if (recolor) {
    .ic_min <- log(drop(t(y) %*% y) / (n_obs - max_lag))

    .model <- .estimate_ar(y, NULL, max_lag, criterion)
    .ic <- .info_criterion(.model$residuals, .model$lag)[[criterion]]

    if (.ic_min < .ic) {
      rho <- 0
      k <- 0
    } else {
      rho <- .model$beta
      k <- .model$lag
      y <- na.omit(.model$residuals)
      n_obs <- nrow(y)
    }
  } else {
    k <- 1
  }

  if (demean) sweep(y, 2, colMeans(y), FUN = "-")

  if (!limit_lags) {
    limit <- n_obs - 1
  } else if (limit_selector == "kpss-q") {
    limit <- 4 * ((n_obs / 100)^(1 / 4))
  } else if (limit_selector == "kpss-m") {
    limit <- 12 * ((n_obs / 100)^(1 / 4))
  } else {
    if (k > 0) {
      limit <- funcs$limit(y, rho_upper_bound)

      if (limit_selector == "Kurozumi") {
        lag_upper_bound <- funcs$limit(lag_upper_bound, rho_upper_bound)
        limit <- min(limit, lag_upper_bound)
      }
    } else {
      limit <- 0
    }
  }
  limit <- trunc(limit)
  limit <- min(limit, n_obs - 1)

  lrv <- (t(y) %*% y) / n_obs
  if (k > 0) {
    for (i in 1:limit) {
      if (i == 0) break
      if (i < n_obs - 1) {
        lrv <- lrv + funcs$weight(i, limit) * (
          (t(y[1:(n_obs - i), ]) %*% y[(1 + i):n_obs, ]) / n_obs +
            t(t(y[1:(n_obs - i), ]) %*% y[(1 + i):n_obs, ]) / n_obs
        )
      } else {
        lrv <- lrv + funcs$weight(i, limit) * as.vector(
          (t(y[1:(n_obs - i), ]) %*% y[(1 + i):n_obs, ]) / n_obs +
            t(t(y[1:(n_obs - i), ]) %*% y[(1 + i):n_obs, ]) / n_obs
        )
      }
    }
  }

  if (recolor) {
    print(rho)
    lrv_recolored <- lrv / (1 - sum(rho))^2
    lrv <- min(lrv_recolored, n_obs * 0.15 * lrv)
  }

  drop(lrv)
}

#' @rdname lr.var
#' @order 2
.variance_lr_bartlett <- function(y) {
  .variance_lr(
    y,
    limit_lags = TRUE,
    limit_selector = "kpss-q"
  )
}

#' @rdname lr.var
#' @order 3
.variance_lr_quad <- function(y) {
  .variance_lr(
    y,
    kernel = "quadratic",
    limit_lags = TRUE,
    limit_selector = "Andrews"
  )
}

#' @rdname lr.var
#' @order 4
.variance_lr_kurozumi <- function(y) {
  .variance_lr(
    y,
    kernel = "bartlett",
    limit_lags = TRUE,
    limit_selector = "Kurozumi"
  )
}

#' @rdname lr.var
#' @order 5
.variance_lr_spc <- function(
  y,
  max_lag = 0,
  kernel = "bartlett",
  criterion = "bic"
) {
  .variance_lr(
    y,
    max_lag = max_lag,
    kernel = kernel,
    criterion = criterion,
    recolor = TRUE
  )
}


.kernel_lr <- function(
  kernel,
  alpha,
  n_obs
) {
  limit_func <- switch(kernel,
    truncated = function(y, l) {
      0.6611 * (n_obs * alpha(y, l)$q2)^(1 / 5)
    },
    bartlett = function(y, l) {
      1.1447 * (n_obs * alpha(y, l)$q1)^(1 / 3)
    },
    parzen = function(y, l) {
      2.6614 * (n_obs * alpha(y, l)$q2)^(1 / 5)
    },
    "tukey-hanning" = function(y, l) {
      1.7462 * (n_obs * alpha(y, l)$q2)^(1 / 5)
    },
    quadratic = function(y, l) {
      1.3221 * (n_obs * alpha(y, l)$q2)^(1 / 5)
    },
    stop("Unknown kernel!")
  )

  weight_func <- switch(kernel,
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
    limit  = limit_func,
    weight = weight_func
  )
}

.alpha_lr_single <- function(y, rho_upper_bound) {
  n_obs <- nrow(y)

  if (!is.null(n_obs) && n_obs > 1) {
    r <- drop(
      (t(y[1:(n_obs - 1), 1]) %*% y[2:n_obs, 1]) /
        (t(y[1:(n_obs - 1), 1]) %*% y[1:(n_obs - 1), 1])
    )

    if (r > rho_upper_bound) {
      r <- rho_upper_bound
    } else if (r < -rho_upper_bound) {
      r <- -rho_upper_bound
    }
  } else {
    r <- y
  }

  list(
    q1 = 4 * r^2 / (1 + r)^2 / (1 - r)^2,
    q2 = 4 * r^2 / (1 - r)^4
  )
}

.alpha_lr_multi <- function(y, rho_upper_bound) {
  n_obs <- nrow(y)
  n_var <- ncol(y)

  nom_1 <- 0
  nom_2 <- 0
  denom <- 0

  for (i in 1:n_var) {
    r <- (t(y[1:(n_obs - 1), i]) %*% y[2:n_obs, i]) /
      (t(y[1:(n_obs - 1), i]) %*% y[1:(n_obs - 1), i])
    r <- drop(r)

    if (r > rho_upper_bound) {
      r <- rho_upper_bound
    } else if (r < -rho_upper_bound) {
      r <- -rho_upper_bound
    }

    resids <- y[2:n_obs, i] - y[1:(n_obs - 1), i] * r
    s2 <- mean(resids^2)

    nom_1 <- nom_1 + 4 * r^2 * s2^2 / (1 - r)^6 / (1 + r)^2
    nom_2 <- nom_2 + 4 * r^2 * s2^2 / (1 - r)^8
    denom <- denom + s2^2 / (1 - r)^4
  }

  list(
    q1 = nom_1 / denom,
    q2 = nom_2 / denom
  )
}

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
LR.variance.single <- function(
  u,
  kernel = "bartlett",
  bw.selector = "kpss-q",
  bw.limit = NULL,
  lag.selector = NULL,
  recolor = FALSE,
  recolor.lag = 0,
  criterion = "bic"
) {
  N <- length(u)

  if (recolor) {
    minIC <- log(sum(u^2) / (N - recolor.lag))

    arModel <- AR.reg(u, NULL, recolor.lag, criterion)
    arIC <- info.criterions(arModel$residuals, arModel$lag)[[
      criterion
    ]]

    if (minIC < arIC) {
      # Sul, Phillips and Choi (2003)
      return((sum(u^2) / N) * min(1, N * 0.15))
    } else {
      arCoefs <- arModel$coefficients
      u <- arModel$residuals
      N <- length(u)
    }
  }

  u <- u - mean(u)

  rho <- sum(u[1:(N - 1)] * u[2:N]) / sum(u[2:N]^2)

  lmfL <- .lr.lag.limit(lag.selector)
  lmtF <- .lr.bandwidth(bw.selector, .alpha.single, N)
  wgtF <- .lr.weight(kernel)

  bw <- lmtF(rho, N)
  if (!is.null(bw.limit)) {
    bw <- min(bw, lmtF(bw.limit, N))
  }
  bw <- trunc(bw)

  lrv <- sum(u^2) / N
  for (i in 1:lmfL(N)) {
    lrv <- lrv + 2 * sum(u[1:(N - i)] * u[(1 + i):N]) * wgtF(i, bw) / N
  }

  if (recolor) {
    lrvR <- lrv / (1 - sum(arCoefs))^2
    return(min(lrvR, N * 0.15 * lrv))
  }

  lrv
}


#' @rdname LR.variance.single
#' @order 2
.lr.var.bartlett <- function(y) {
  LR.variance.single(
    y,
    kernel = "Bartlett",
    lag.selector = "kpss-q"
  )
}

#' @rdname LR.variance.single
#' @order 3
.lr.var.quad <- function(y) {
  LR.variance.single(
    y,
    kernel = "Quadratic",
    bw.selector = "Quadratic"
  )
}

#' @rdname LR.variance.single
#' @order 4
.lr.var.kurozumi <- function(y) {
  LR.variance.single(
    y,
    kernel = "Bartlett",
    bw.selector = "Bartlett",
    bw.limit = 0.8
  )
}

#' @rdname LR.variance.single
#' @order 5
.lr.var.spc <- function(
  y,
  max.lag = 0,
  kernel = "Bartlett",
  criterion = "bic"
) {
  LR.variance.single(
    y,
    max.lag = max.lag,
    kernel = kernel,
    criterion = criterion,
    recolor = TRUE
  )
}


.lr.lag.limit <- function(selector) {
  switch(
    selector,
    "kpss-q" = function(N) {
      4 * (N / 100)^(1 / 4)
    },
    "kpss-m" = function(N) {
      12 * (N / 100)^(1 / 4)
    },
    function(N) {
      N - 1
    }
  )
}


.lr.bandwidth <- function(selector, alpha, N) {
  switch(
    selector,
    "Bartlett" = function(r, N) {
      1.1447 * (N * alpha(r)$q1)^(1 / 3)
    },
    "Parzen" = function(r, N) {
      2.6614 * (N * alpha(r)$q2)^(1 / 5)
    },
    "Tuckey-Hanning" = function(r, N) {
      1.7462 * (N * alpha(r)$q2)^(1 / 5)
    },
    "Quadratic" = function(r, N) {
      1.3221 * (N * alpha(r)$q2)^(1 / 5)
    },
    stop("LR.variance: Unknown banwidth selector!")
  )
}

.lr.weight <- function(kernel) {
  switch(
    kernel,
    "truncated" = function(i, l) {
      if (abs(i / (l + 1)) <= 1) {
        return(1)
      }
      0
    },
    "Bartlett" = function(i, l) {
      x <- i / (l + 1)
      if (abs(x) < 1) {
        return(1 - abs(x))
      }
      0
    },
    "Parzen" = function(i, l) {
      x <- i / (l + 1)
      if (abs(x) <= 0.5) {
        return(1 - 6 * x^2 + 6 * abs(x)^3)
      } else if (abs(x) <= 1) {
        return(2 * (1 - abs(x))^3)
      }
      0
    },
    "Tukey-Hanning" = function(i, l) {
      x <- i / (l + 1)
      if (abs(x) <= 1) {
        return((1 + cos(pi * x)) / 2)
      }
      0
    },
    "Quadratic" = function(i, l) {
      delta <- 6 * pi * i / (5 * l)
      (3 / delta^2) * (sin(delta) / delta - cos(delta))
    },
    stop("Unknown kernel!")
  )
}

.alpha.single <- function(r) {
  list(
    q1 = 4 * r^2 / (1 + r)^2 / (1 - r)^2,
    q2 = 4 * r^2 / (1 - r)^4
  )
}

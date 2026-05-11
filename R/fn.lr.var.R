#' Calculating long-run variance or covariance matrix
#' @name LR.variance
#'
#' @param y A time series of interest.
#' @param kernel A kernel to be used:
#' * `Truncated`: \eqn{\left\{\begin{array}{ll}
#' 1 & |x| \leq 1 \\
#' 0 & \textrm{otherwize}
#' \end{array}\right.}
#' * `Bartlett`: \eqn{\left\{\begin{array}{ll}
#' 1 - |x| & |x| \leq 1 \\
#' 0 & \textrm{otherwize}
#' \end{array}\right.},
#' * `Parzen`: \eqn{\left\{\begin{array}{ll}
#' 1 - 6 x^2 + 6 {|x|}^3 & |x| \leq 1/2 \\
#' 2 (1 - |x|)^3 & 1/2 \leq |x| \leq 1 \\
#' 0 & \textrm{otherwize}
#' \end{array}\right.},
#' * `Tukey-Hanning`: \eqn{\left\{\begin{array}{ll}
#' (1 + \cos(\pi x))/2 & |x| \leq 1 \\
#' 0 & \textrm{otherwize}
#' \end{array}\right.},
#' * `Quadratic`: \eqn{
#' \frac{25}{12 \pi^2 x^2}
#' \left(\frac{\sin(6 \pi x / 5)}{6 \pi x / 5} - \cos(6 \pi x / 5)\right)}.
#' @param k A limiting parameter for Kurozumi's proposal.
#' @param kmax A maximum number of lars for recoloring procedure from Sul et al. (2005).
#' @param criterion An information criterion for recoloring lag selection.
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
NULL


#' @rdname LR.variance
#' @order 1
.lr.var.bartlett <- function(y) {
  N <- length(y)
  m <- round(.lr.bandwidth(NULL, N, "kpss-q"))
  wgtF <- .lr.weight("Bartlett")

  lrv <- sum(y^2) / N
  for (i in 1:(N - 1)) {
    lrv <- lrv + 2 * sum(y[1:(N - i)] * y[(1 + i):N]) * wgtF(i, m) / N
  }

  lrv
}

#' @rdname LR.variance
#' @order 2
.lr.var.quad <- function(y) {
  N <- length(y)
  a <- sum(y[1:(N - 1)] * y[2:N]) / sum(y[2:N]^2)
  a <- .alpha.single(a)$q2
  m <- .lr.bandwidth(a, N, "Quadratic")
  wgtF <- .lr.weight("Quadratic")

  lrv <- sum(y^2) / N
  for (i in 1:(N - 1)) {
    lrv <- lrv + 2 * sum(y[1:(N - i)] * y[(1 + i):N]) * wgtF(i - 1, m) / N
  }

  lrv
}

#' @rdname LR.variance
#' @order 3
.lr.var.kurozumi <- function(y, k = 0.8) {
  N <- length(y)

  a <- sum(y[1:(N - 1)] * y[2:N]) / sum(y[2:N]^2)
  a <- .alpha.single(a)$q1
  k <- .alpha.single(k)$q1

  m <- min(
    .lr.bandwidth(a, N, "Bartlett"),
    .lr.bandwidth(k, N, "Bartlett")
  )
  m <- trunc(m)

  wgtF <- .lr.weight("Bartlett")

  lrv <- sum(y^2) / N
  for (i in 1:m) {
    lrv <- lrv + 2 * sum(y[1:(N - i)] * y[(1 + i):N]) * wgtF(i, m) / N
  }

  lrv
}

#' @rdname LR.variance
#' @order 4
.lr.var.spc <- function(
  y,
  kmax = NULL,
  kernel = "Bartlett",
  criterion = "bic"
) {
  N <- length(y)

  if (is.null(kmax) || kmax < 0) kmax <- .lr.bandwidth(a, N, kernel)
  kmax <- max(kmax, 0)

  min_IC <- log(sum(y^2) / (N - kmax))

  arModel <- AR.reg(y, NULL, kmax, criterion)
  IC <- info.criterions(arModel$residuals, arModel$lag)[[criterion]]

  # Sul, Phillips and Choi (2003)
  if (min_IC < IC) {
    return((sum(y^2) / N) * min(1, N * 0.15))
  }

  rho <- arModel$coefficients
  res <- arModel$residuals
  N <- length(res)

  a <- drop(sum(y[1:(N - 1)] * y[2:N]) / sum(y[2:N]^2))
  a <- .alpha.single(a)$q1
  m <- trunc(.lr.bandwidth(a, N, kernel))
  wgtF <- .lr.weight(kernel)

  lrv <- sum(y^2) / N
  for (i in 1:m) {
    lrv <- lrv + 2 * sum(y[1:(N - i)] * y[(1 + i):N]) * wgtF(i, m) / N
  }

  lrv_recolored <- lrv / (1 - sum(rho))^2

  min(lrv_recolored, lrv * N * 0.15)
}


.lr.bandwidth <- function(alpha, N, selector = "Bartlett") {
  switch(selector,
    "Bartlett" = 1.1447 * (N * alpha)^(1 / 3),
    "Parzen" = 2.6614 * (N * alpha)^(1 / 5),
    "Tuckey-Hanning" = 1.7462 * (N * alpha)^(1 / 5),
    "Quadratic" = 1.3221 * (N * alpha)^(1 / 5),
    "kpss-q" = 4 * (N / 100)^(1 / 4),
    "kpss-m" = 12 * (N / 100)^(1 / 4),
    "full" = N - 1,
    stop("LR.variance: Unknown banwidth selector!")
  )
}


.lr.weight <- function(kernel) {
  switch(kernel,
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

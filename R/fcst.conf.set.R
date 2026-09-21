#' @title
#' Critical value for the Eo-Morley likelihood ratio statistic
#'
#' @param level A confidence level in (0, 1).
#' @param omega Two positive scale factors from Proposition 1.
#'
#' @return The quantile of equation (6), on the scale of twice the difference
#' in log likelihoods. With both scales equal to one, the 90%, 95% and 99%
#' critical values are approximately 5.94, 7.35 and 10.59.
#'
#' @keywords internal
get.cv.EoMorley <- function(level = 0.95, omega = c(1, 1)) {
  if (!is.numeric(level) || length(level) != 1 || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("ERROR! get.cv.EoMorley: level should be in (0, 1)")
  }
  if (!is.numeric(omega) || length(omega) != 2 ||
      any(!is.finite(omega)) || any(omega <= 0)) {
    stop("ERROR! get.cv.EoMorley: omega should contain two positive finite scales")
  }

  # Equivalent to -2*log(1-sqrt(level)), without cancellation near one.
  upper <- 2 * (log1p(sqrt(level)) - log1p(-level))
  scale <- max(omega)
  if (omega[1] == omega[2]) {
    return(unname(scale * upper))
  }
  root <- uniroot(function(z) .eo.morley.cdf(z * upper, omega / scale) - level,
                  interval = c(0, 1), tol = 1e-12)$root
  unname(scale * upper * root)
}


.eo.morley.cdf <- function(x, omega) {
  result <- numeric(length(x))
  positive <- x > 0
  result[positive] <- (-expm1(-x[positive] / (2 * omega[1]))) *
    (-expm1(-x[positive] / (2 * omega[2])))
  result
}


#' @title
#' Eo-Morley confidence set for the date of a structural break
#'
#' @description
#' Inverts the likelihood ratio tests in Eo and Morley (2015), equations
#' (4)-(6). This implementation considers a univariate regression with one
#' break in all coefficients, with either a constant or a regime-specific
#' error variance. The accepted dates need not form a contiguous interval.
#'
#' @param y A dependent variable. Incomplete rows of `y` and `x` are removed
#' jointly; all returned dates refer to the complete-case sample.
#' @param x Explanatory variables. If `NULL` a single constant is used.
#' @param level A confidence level in (0, 1).
#' @param trim A trimming fraction in (0, 0.5). Each regime must contain at
#' least `max(ceiling(trim * N), k + 1)` observations and have full column rank.
#' @param variance `regime` allows the variance to change at the coefficient
#' break; `constant` imposes a common variance at every candidate date.
#' @param vcov.type How to estimate the limit-distribution scales:
#' * `iid`: errors independent of regressors and serially independent within
#'   each regime; uses the model-based score covariance and empirical fourth
#'   moments. With constant variance this gives Corollary 1, `omega = c(1, 1)`;
#' * `HC`: empirical score covariances, allowing conditional heteroskedasticity;
#' * `HAC`: within-regime long-run score covariances, allowing serial correlation.
#' @param bandwidth An optional positive integer passed to [HAC.variance]
#' for `vcov.type = "HAC"`. Must not exceed the shorter fitted regime length.
#' If `NULL`, `max(1, floor(n_j^(1/4)))` is used separately in each regime.
#' @param time.index Optional original integer-period positions for HAC lag
#' products, with one entry per input row. Does not change break-date indexing.
#' The estimation sample must be contiguous; internal missing rows are rejected.
#' @param kernel A HAC kernel: `Bartlett`, `Parzen` or `Tukey-Hanning`.
#'
#' @details
#' For every admissible date \eqn{b}, coefficients are reestimated by OLS.
#' The profile Gaussian quasi-log-likelihood uses ML variances, `RSS_j/n_j`
#' for `regime` and `(RSS_1+RSS_2)/N` for `constant`. The break estimator
#' maximises that same profile. The confidence set contains exactly the dates
#' satisfying \eqn{LR(b)=2[\ell(\hat b)-\ell(b)]\leq c_{level}}.
#'
#' Two [SSR.recursive] passes from `breaktest` compute all candidate residual
#' sums of squares. [OLS.reg] supplies the fitted regimes at the maximum.
#' The scales \eqn{\omega_i=\Gamma_i^2/\Psi_i} are estimated once there,
#' following Proposition 1, and held fixed when testing candidate dates.
#' The CDF is \eqn{(1-\exp(-c/(2\omega_1)))(1-\exp(-c/(2\omega_2)))}.
#' HAC calculations reuse [HAC.variance] and the package's kernel functions.
#'
#' In the univariate notation, let \eqn{v_i=\hat\sigma_i^2},
#' \eqn{d=\hat\beta_2-\hat\beta_1} and \eqn{M_i=X_i'X_i/n_i}.
#' Then \eqn{B_1=(v_2-v_1)/v_2}, \eqn{B_2=(v_2-v_1)/v_1},
#' \eqn{Q_1=M_1/v_2}, \eqn{Q_2=M_2/v_1} and
#' \eqn{\Psi_i=B_i^2/2+d'Q_i d}.
#' The score sequences are \eqn{X_1 u_1/v_2} and \eqn{X_2 u_2/v_1}; their
#' covariance estimates give \eqn{\Pi_i}. For `iid`, these simplify to
#' \eqn{\Pi_1=v_1 M_1/v_2^2} and \eqn{\Pi_2=v_2 M_2/v_1^2}.
#' Covariances of \eqn{u_i^2/v_i-1} give \eqn{\Omega_i}, and
#' \eqn{\Gamma_i^2=B_i^2\Omega_i/4+d'\Pi_i d}.
#'
#' Validity requires the assumptions of the paper, including stabilising
#' moments within regimes (no unit roots or deterministic trends), an
#' identifiable break and the zero third-moment condition used in Proposition
#' 1. This is inference on timing conditional on a break, not a test for its
#' existence. Multiple breaks, multivariate systems, partial coefficient
#' restrictions and a variance-only restricted model are not implemented here.
#' For inference about a coefficient break allowing a variance break, use
#' `variance = "regime"`; this does not constrain the coefficients to be equal.
#'
#' @return An object of class `bt_fcstConfset`, a list of:
#' * `bp`: the profile-likelihood maximiser (earliest in an exact tie),
#' * `set`: every accepted last-pre-break date,
#' * `lower`, `upper`: the envelope of the set, not an interval replacing it,
#' * `intervals`: contiguous components with columns `lower` and `upper`,
#' * `profile`: all candidate dates, log likelihoods, LR values and acceptance,
#' * `critical.value`, `level`, `omega`: the cutoff and its calibration,
#' * `scale`: the estimated terms from Proposition 1,
#' * `fit`: the two pilot [OLS.reg] regressions,
#' * `variance`, `vcov.type`, `bandwidth`, `kernel`, `sigma`, `N`, `k`, `trim`:
#'   model settings, pilot standard deviations and sample information.
#'
#' @references
#' Eo, Yunjong, and James Morley.
#' "Likelihood-Ratio-Based Confidence Sets for the Timing of Structural Breaks."
#' Quantitative Economics 6, no. 2 (2015): 463-497.
#' https://doi.org/10.3982/QE186.
#'
#' @export
fcst.conf.set.EoMorley <- function(
  y,
  x = NULL,
  level = 0.95,
  trim = 0.15,
  variance = c("regime", "constant"),
  vcov.type = c("iid", "HC", "HAC"),
  bandwidth = NULL,
  kernel = c("Bartlett", "Parzen", "Tukey-Hanning"),
  time.index = NULL
) {
  variance <- match.arg(variance)
  vcov.type <- match.arg(vcov.type)
  kernel <- match.arg(kernel)
  get.cv.EoMorley(level)  # Validate the level before fitting any model.
  if (!is.numeric(trim) || length(trim) != 1 || !is.finite(trim) ||
      trim <= 0 || trim >= 0.5) {
    stop("ERROR! fcst.conf.set.EoMorley: trim should be in (0, 0.5)")
  }
  if (!is.null(bandwidth) &&
      (!is.numeric(bandwidth) || length(bandwidth) != 1 ||
       !is.finite(bandwidth) || bandwidth < 1 || bandwidth != floor(bandwidth))) {
    stop("ERROR! fcst.conf.set.EoMorley: bandwidth should be a positive integer")
  }
  if (!is.null(bandwidth) && vcov.type != "HAC") {
    stop("ERROR! fcst.conf.set.EoMorley: bandwidth is used only with HAC")
  }

  time.index <- .fcst.time.index(time.index, NROW(y))
  .d <- .fcst.data(y, x)
  time.index <- time.index[.d$rows]
  .fcst.require.regular(time.index, "fcst.conf.set.EoMorley")
  y <- .d$y
  x <- .d$x
  N <- .d$N
  k <- .d$k
  if (!is.numeric(y) || !is.numeric(x) || any(!is.finite(y)) ||
      any(!is.finite(x)) || k == 0) {
    stop("ERROR! fcst.conf.set.EoMorley: finite numeric data required")
  }
  width <- max(ceiling(trim * N), k + 1)
  if (2 * width > N) {
    stop("ERROR! fcst.conf.set.EoMorley: trimming leaves no admissible break dates")
  }
  if (qr(x[seq_len(width), , drop = FALSE])$rank < k ||
      qr(x[seq.int(N - width + 1, N), , drop = FALSE])$rank < k) {
    stop("ERROR! fcst.conf.set.EoMorley: candidate regimes must have full column rank")
  }

  dates <- seq.int(width, N - width)
  ssr.pre <- .fcst.ssr.recursive(y, x, 1, N, width)[dates]
  reverse <- rev(seq_len(N))
  ssr.post <- .fcst.ssr.recursive(y[reverse, , drop = FALSE],
                             x[reverse, , drop = FALSE], 1, N, width)[N - dates]
  invalid.variance <- if (variance == "constant") {
    any(ssr.pre + ssr.post <= 0)
  } else {
    any(ssr.pre <= 0) || any(ssr.post <= 0)
  }
  if (any(!is.finite(c(ssr.pre, ssr.post))) || invalid.variance) {
    stop("ERROR! fcst.conf.set.EoMorley: candidate regimes need positive residual variances")
  }
  loglik <- if (variance == "constant") {
    -0.5 * N * (log(2 * pi) + 1 + log((ssr.pre + ssr.post) / N))
  } else {
    -0.5 * (N * (log(2 * pi) + 1) + dates * log(ssr.pre / dates) +
              (N - dates) * log(ssr.post / (N - dates)))
  }
  best <- which.max(loglik)
  bp <- dates[best]
  pre <- seq_len(bp)
  post <- seq.int(bp + 1, N)
  fit <- list(pre = OLS.reg(y[pre, , drop = FALSE], x[pre, , drop = FALSE]),
              post = OLS.reg(y[post, , drop = FALSE], x[post, , drop = FALSE]))
  sigma.sq <- c(pre = mean(fit$pre$residuals^2),
                post = mean(fit$post$residuals^2))
  if (variance == "constant") {
    sigma.sq[] <- (bp * sigma.sq[1] + (N - bp) * sigma.sq[2]) / N
  }
  if (vcov.type == "HAC") {
    if (!is.null(bandwidth) && bandwidth > min(bp, N - bp)) {
      stop("ERROR! fcst.conf.set.EoMorley: bandwidth exceeds the shorter fitted regime")
    }
    bandwidth <- if (is.null(bandwidth)) {
      pmax(1, floor(c(bp, N - bp)^(1 / 4)))
    } else {
      rep(bandwidth, 2)
    }
    names(bandwidth) <- c("pre", "post")
  }
  scale <- .eo.morley.scale(fit, sigma.sq, vcov.type, bandwidth, kernel,
                            list(time.index[pre], time.index[post]))
  critical.value <- get.cv.EoMorley(level, scale$omega)
  LR <- pmax(0, 2 * (loglik[best] - loglik))
  accepted <- LR <= critical.value
  set <- dates[accepted]
  groups <- cumsum(c(1, diff(set) != 1))
  intervals <- data.frame(lower = as.integer(tapply(set, groups, min)),
                          upper = as.integer(tapply(set, groups, max)))

  result <- list(
    bp = bp,
    set = set,
    lower = min(set),
    upper = max(set),
    intervals = intervals,
    profile = data.frame(bp = dates, loglik = loglik, LR = LR, accepted = accepted),
    critical.value = critical.value,
    level = level,
    omega = scale$omega,
    scale = scale,
    fit = fit,
    variance = variance,
    vcov.type = vcov.type,
    bandwidth = bandwidth,
    kernel = if (vcov.type == "HAC") kernel else NULL,
    sigma = sqrt(sigma.sq),
    N = N,
    k = k,
    trim = trim
  )
  class(result) <- "bt_fcstConfset"
  result
}


# Univariate specialisation of the quantities preceding Proposition 1.
# All moments are taken within the fitted regimes, not at candidate dates.
.eo.morley.scale <- function(fit, sigma.sq, vcov.type, bandwidth, kernel,
                              time.index = list(NULL, NULL)) {
  delta <- fit$post$coefficients - fit$pre$coefficients
  B <- diff(unname(sigma.sq)) / rev(sigma.sq)
  Q <- Pi <- vector("list", 2)
  Omega <- Psi <- Gamma.sq <- numeric(2)
  for (j in 1:2) {
    x <- fit[[j]]$exog
    u <- fit[[j]]$residuals
    n <- length(u)
    opposite <- sigma.sq[3 - j]
    Q[[j]] <- crossprod(x) / (n * opposite)
    eta.sq <- u^2 / sigma.sq[j] - 1
    if (vcov.type == "iid") {
      Pi[[j]] <- Q[[j]] * sigma.sq[j] / opposite
      Omega[j] <- mean(eta.sq^2)
    } else {
      bw <- if (vcov.type == "HC") 1 else bandwidth[j]
      Pi[[j]] <- as.matrix(HAC.variance(x * (u / opposite), bw, kernel,
                                       time.index = time.index[[j]]))
      Omega[j] <- HAC.variance(eta.sq, bw, kernel, time.index = time.index[[j]])
    }
    Psi[j] <- B[j]^2 / 2 + drop(crossprod(delta, Q[[j]] %*% delta))
    Gamma.sq[j] <- B[j]^2 * Omega[j] / 4 +
      drop(crossprod(delta, Pi[[j]] %*% delta))
  }
  # Corollary 1 also avoids an artificial 0/0 at an exactly zero fitted shift.
  omega <- if (vcov.type == "iid" && sigma.sq[1] == sigma.sq[2]) {
    c(1, 1)
  } else {
    Gamma.sq / Psi
  }
  if (any(!is.finite(omega)) || any(omega <= 0)) {
    stop("ERROR! fcst.conf.set.EoMorley: scale factors are not positive; break is not identified")
  }
  names(omega) <- names(B) <- names(Q) <- names(Pi) <-
    names(Omega) <- names(Psi) <- names(Gamma.sq) <- c("pre", "post")
  list(omega = omega, B = B, Q = Q, Pi = Pi, Omega = Omega,
        Psi = Psi, Gamma.sq = Gamma.sq, delta = delta)
}

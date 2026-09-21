#' @title
#' Observation weights for forecasting under structural breaks
#' @name fcst.weights
#'
#' @description
#' A family of generators returning weights summing to unity. Multiple-break
#' PPP optimal weights may be negative. The weights are passed to [WLS.reg] or to
#' the corresponding named strategy of [fcst.forecast].
#'
#' The family contains baseline weights and the source-specific PPP optimal
#' and robust rules. Under the discrete-break assumptions of Pesaran et al.
#' (2013), oracle weights are piecewise constant across regimes. This does not
#' make every member of the family optimal for an arbitrary regression.
#'
#' @references
#' Pesaran, M. Hashem, Andreas Pick, and Mikhail Pranovich.
#' “Optimal Forecasts in the Presence of Structural Breaks.”
#' Journal of Econometrics 177, no. 2 (2013): 134–52.
#' https://doi.org/10.1016/j.jeconom.2013.04.002.
#'
#' Pesaran, M. Hashem, and Andreas Pick.
#' “Forecast Combination Across Estimation Windows.”
#' Journal of Business and Economic Statistics 29, no. 2 (2011): 307–18.
#'
#' @keywords internal
NULL


#' @rdname fcst.weights
#' @order 1
#'
#' @param N A number of observations.
#'
#' @details
#' `weights.equal` reproduces the full-sample OLS estimator and serves as the
#' benchmark all the other schemes are compared against.
#' @exportS3Method NULL
weights.equal <- function(N) {
  .fcst.weights.N(N)
  rep(1 / N, N)
}


# Validate the common sample-size contract of all weight generators.
.fcst.weights.N <- function(N) {
  if (!is.numeric(N) || length(N) != 1L || !is.finite(N) ||
      N < 1 || N != floor(N)) {
    stop("ERROR! fcst.weights: N must be a positive integer")
  }
  invisible(N)
}


#' @rdname fcst.weights
#' @order 2
#'
#' @param R A window size, the number of the most recent observations used.
#'
#' @details
#' `weights.rolling` gives equal weight to the last `R` observations and zero
#' weight to the rest. With `R = N - bp` it produces the post-break forecast.
#' @exportS3Method NULL
weights.rolling <- function(N, R) {
  .fcst.weights.N(N)
  if (!is.numeric(R) || length(R) != 1L || !is.finite(R) ||
      R != floor(R) || R < 1 || R > N) {
    stop("ERROR! weights.rolling: R must be an integer in [1, N]")
  }

  w <- numeric(N)
  w[(N - R + 1):N] <- 1 / R
  w
}


#' @rdname fcst.weights
#' @order 3
#'
#' @param gamma A down-weighting parameter, \eqn{0 \leq \gamma \leq 1} for
#' exponential weights. For `weights.AA.kernel`, any finite nonnegative
#' pre/post relative weight is allowed; see Details.
#' @param time.index Increasing original integer positions for exponential
#' weights. Defaults to 1:N. Missing periods retain their elapsed time.
#'
#' @details
#' `weights.exponential` implements the exponential smoothing (constant gain
#' least squares) scheme \eqn{w_t \propto \gamma^{N - t}}. Pesaran et al. (2013)
#' derive their geometric approximation for a scalar random-walk level model.
#' Their finite-sample optimum (eq. 3) includes boundary corrections; eq. 10
#' gives this approximation, with
#' \eqn{\gamma = 1 + \delta^2 / 2 - \delta (1 + \delta^2 / 4)^{1/2}} and
#' \eqn{\delta^2 = \sigma^2_v / \sigma^2_\varepsilon}.
#' gamma=0 places all mass on the last observation; a multiple-regressor
#' model is generally unidentified there. Adaptive rules are in [fcst.ExpS].
#' @exportS3Method NULL
weights.exponential <- function(N, gamma = 0.95, time.index = NULL) {
  if (!is.numeric(N) || length(N) != 1L || !is.finite(N) || N < 1 || N != floor(N)) {
    stop("ERROR! weights.exponential: N must be a positive integer")
  }
  if (!is.numeric(gamma) || length(gamma) != 1L || !is.finite(gamma) || gamma < 0 || gamma > 1) {
    stop("ERROR! weights.exponential: gamma should be in [0, 1]")
  }
  time.index <- .fcst.time.index(time.index, N)
  # gamma=0 is the last-observation limit (0^0=1); an arbitrary regression
  # need not be identified at this boundary. Gaps count as elapsed periods.
  w <- gamma^(max(time.index) - time.index)
  w / sum(w)
}


#' @rdname fcst.weights
#' @order 4
#'
#' @param delta A ratio \eqn{\sigma_v / \sigma_\varepsilon} of the standard
#' deviation of the coefficient innovations to that of the equation error.
#'
#' @details
#' `weights.ExpS.delta` is a convenience wrapper translating the signal-to-noise
#' ratio of the random walk coefficient model into the down-weighting parameter
#' of the geometric approximation, not the exact finite-sample oracle weights.
#' @exportS3Method NULL
weights.ExpS.delta <- function(N, delta) {
  if (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta) || delta < 0) {
    stop("ERROR! weights.ExpS.delta: delta must be finite and nonnegative")
  }

  # Rationalisation avoids subtracting nearly equal large terms. Scale first
  # for large delta, so neither delta^2 nor the denominator can overflow.
  root.gamma <- if (delta <= 1) {
    2 / (sqrt(delta^2 + 4) + delta)
  } else {
    (2 / delta) / (sqrt(1 + (2 / delta)^2) + 1)
  }
  weights.exponential(N, root.gamma^2)
}


#' @rdname fcst.weights
#' @order 5
#'
#' @param bp A break point (the last observation of the first regime), or a
#' vector of break points for the multiple break case.
#' @param phi A normalised break magnitude
#' \eqn{\phi = x_{T+1}' \lambda / (x_{T+1}' \Omega_{xx}^{-1} x_{T+1})^{1/2}}
#' with \eqn{\lambda = (\beta_{(1)} - \beta_{(2)}) / \sigma_{(2)}}, or a vector
#' of such magnitudes for the multiple break case. See [fcst.break.size].
#' @param q A ratio \eqn{\sigma_{(1)} / \sigma_{(2)}} of the pre- to post-break
#' error standard deviations. Must equal one for multiple breaks, where the
#' original derivation assumes a common error variance.
#'
#' @details
#' `weights.PPP.optimal` implements the MSFE-optimal weights of Pesaran et al.
#' (2013), eqs. (31)-(32) for a single break and eqs. (38)-(39) for `n` breaks.
#' The weights are constant within a regime and differ across regimes.
#'
#' Under multiple breaks the weights need not increase towards the end of the
#' sample and may even be negative: biases of the opposite sign coming from
#' different regimes can offset each other, so an early regime may receive the
#' largest weight.
#'
#' The oracle scheme requires the break dates and magnitudes to be known.
#' Estimating them introduces uncertainty; [weights.PPP.robust] avoids these
#' inputs through a first-order asymptotic approximation.
#' @exportS3Method NULL
weights.PPP.optimal <- function(N, bp, phi, q = 1) {
  .fcst.weights.N(N)
  n <- length(bp)

  if (!is.numeric(phi) || length(phi) != n || !n || any(!is.finite(phi))) {
    stop("ERROR! weights.PPP.optimal: finite phi and nonempty bp must have equal length")
  }
  if (!is.numeric(bp) || any(!is.finite(bp)) || any(bp != floor(bp)) ||
      any(bp < 1) || any(bp >= N) || is.unsorted(bp, strictly = TRUE)) {
    stop("ERROR! weights.PPP.optimal: bp should be an increasing sequence in [1, N)")
  }
  if (!is.numeric(q) || length(q) != 1L || !is.finite(q) || q <= 0) {
    stop("ERROR! weights.PPP.optimal: q must be finite and positive")
  }
  if (n > 1L && q != 1) {
    stop("ERROR! weights.PPP.optimal: multiple breaks require the common-variance assumption q = 1")
  }

  b <- bp / N

  if (n == 1) {
    # Evaluate the pre/post relative weight without overflowing squared inputs.
    scale <- max(1, q, abs(phi))
    relative <- c((1 / scale)^2, (q / scale)^2 + bp * (phi / scale)^2)
    w <- c(rep(relative[1], bp), rep(relative[2], N - bp))
    return(w / sum(w))
  } else {
    # Multiple breaks, eqs. (38)-(39). The error variance is assumed stable.
    b.aug <- c(0, b, 1)
    len <- diff(b.aug)            # regime lengths, n + 1 of them
    # Algebraically identical to eqs. (38)-(39): a.n=1+N*Var(phi).
    # Centred moments and common scaling avoid cancellation/overflow in a.n.
    scale <- max(1, abs(phi))
    magnitude <- c(phi / scale, 0)
    average <- sum(len * magnitude)
    a.n <- (1 / scale)^2 + N * sum(len * (magnitude - average)^2)
    w.reg <- vapply(magnitude, function(value) {
      ((1 / scale)^2 + N * sum(len * magnitude * (magnitude - value))) / (N * a.n)
    }, 0.0)
  }

  if (any(!is.finite(w.reg))) {
    stop("ERROR! weights.PPP.optimal: magnitudes exceed the numerical range of the multiple-break formula")
  }
  regime <- findInterval(seq_len(N), bp + 1) + 1
  w.reg[regime]
}


#' @rdname fcst.weights
#' @order 6
#'
#' @param b.lo,b.hi The support of the uniform prior for the break fraction.
#' When both are NULL, the full-support approximation of eqs. (46)-(48) is
#' used. Explicit bounds use the finite-support expression of eq. (44), so
#' supplying `1 / N, 1 - 1 / N` need not reproduce the default weights.
#'
#' @details
#' `weights.PPP.robust` integrates the optimal weights over a uniformly
#' distributed break date. The first order term of the resulting expansion does
#' not depend on the break magnitude at all, and a break in the error variance
#' enters only at order \eqn{T^{-1}}. With the full support this collapses to
#' \deqn{w_t^* = \frac{-\log(1 - t / N)}{N - 1}, \quad t = 1, \dots, N - 1,
#' \qquad w_N^* = \frac{\log N}{N - 1},}
#' rescaled to sum to unity.
#'
#' This asymptotic robust rule needs neither the break date nor its magnitude.
#' Its forecast performance depends on the design and break size; it does not
#' uniformly dominate equal weights, which minimise variance with no break.
#' @exportS3Method NULL
weights.PPP.robust <- function(N, b.lo = NULL, b.hi = NULL) {
  .fcst.weights.N(N)
  if (N == 1L) {
    if (!is.null(b.lo) || !is.null(b.hi)) {
      stop("ERROR! weights.PPP.robust: a break prior requires at least two observations")
    }
    return(1)
  }
  a <- seq_len(N) / N

  if (is.null(b.lo) && is.null(b.hi)) {
    w <- c(-log(1 - a[-N]), log(N)) / (N - 1)
  } else {
    b.lo <- if (is.null(b.lo)) 1 / N else b.lo
    b.hi <- if (is.null(b.hi)) 1 - 1 / N else b.hi

    if (!is.numeric(b.lo) || length(b.lo) != 1L || !is.finite(b.lo) ||
        !is.numeric(b.hi) || length(b.hi) != 1L || !is.finite(b.hi) ||
        b.lo <= 0 || b.hi >= 1 || b.lo >= b.hi) {
      stop("ERROR! weights.PPP.robust: 0 < b.lo < b.hi < 1 is required")
    }

    w <- numeric(N)
    w[a >= b.lo & a <= b.hi] <-
      -log((1 - a[a >= b.lo & a <= b.hi]) / (1 - b.lo)) / (b.hi - b.lo)
    w[a > b.hi] <- -log((1 - b.hi) / (1 - b.lo)) / (b.hi - b.lo)
    w <- w / N
  }

  w / sum(w)
}


#' @rdname fcst.weights
#' @order 7
#'
#' @details
#' `weights.AA.kernel` assigns relative weight `gamma` before the break and
#' unity after it, then normalises the weights. On `[0, 1]` this is the discrete
#' Aitchison-Aitken kernel used by Lee et al. (2022), including its zero limit.
#' With `gamma = 1` it gives full-sample OLS, with zero post-break OLS.
#' Values above one are also accepted to represent the effective WGLS weight
#' \eqn{\gamma^* = \gamma_{kernel}/q^2}; the kernel parameter itself remains
#' bounded by one. Here \eqn{q} is the pre/post standard deviation ratio.
#'
#' Pesaran et al. (2013) weights are a special case with
#' \eqn{\gamma = w_{(1)} / w_{(2)}}, so estimating `gamma` by cross-validation
#' (see [fcst.WGLS.cv]) avoids having to specify the break magnitude and the
#' variance ratio separately.
#'
#' @references
#' Lee, Tae-Hwy, Shahnaz Parsaeian, and Aman Ullah.
#' “Forecasting under Structural Breaks Using Improved Weighted Estimation.”
#' Oxford Bulletin of Economics and Statistics, 2022.
#' @exportS3Method NULL
weights.AA.kernel <- function(N, bp, gamma) {
  if (!is.numeric(N) || length(N) != 1 || !is.finite(N) ||
      N < 2 || N != floor(N)) {
    stop("ERROR! weights.AA.kernel: N should be an integer of at least two")
  }
  if (!is.numeric(gamma) || length(gamma) != 1 || !is.finite(gamma) || gamma < 0) {
    stop("ERROR! weights.AA.kernel: gamma should be finite and nonnegative")
  }
  if (!is.numeric(bp) || length(bp) != 1 || !is.finite(bp) ||
      bp != floor(bp) || bp < 1 || bp >= N) {
    stop("ERROR! weights.AA.kernel: bp is out of range")
  }

  # Scaling first avoids overflow when effective weights exceed one.
  scale <- max(1, gamma)
  w <- c(rep(gamma / scale, bp), rep(1 / scale, N - bp))
  w / sum(w)
}


#' @title
#' Estimate the normalised break magnitude and the variance ratio
#'
#' @description
#' Feasible counterparts of the quantities \eqn{\phi} and \eqn{q} entering the
#' optimal weights of Pesaran et al. (2013).
#'
#' @param y A dependent variable.
#' @param x Explanatory variables.
#' @param bp A break point or a vector of break points.
#' @param x.new A (1 x k) vector of regressor values for the forecast origin.
#' @param time.index Original increasing integer positions; the retained sample
#' must be regular without internal missing periods for the source PPP pilots.
#'
#' @return A list of:
#' * `phi`: normalised break magnitudes, one per break,
#' * `q`: the ratio of the pre- to post-break error standard deviation
#'   (only meaningful for a single break),
#' * `lambda`: raw coefficient differences scaled by the last regime sigma,
#' * `sigma`: per-regime residual standard deviations.
#'
#' @details
#' \eqn{\phi_{(i)} = x_{T+1}' \lambda_{(i)} /
#' (x_{T+1}' \Omega_{xx}^{-1} x_{T+1})^{1/2}} with
#' \eqn{\lambda_{(i)} = (\beta_{(i)} - \beta_{(n+1)}) / \sigma_{(n+1)}} and
#' \eqn{\Omega_{xx}} estimated by the full-sample second moment matrix.
#'
#' @keywords internal
fcst.break.size <- function(y, x, bp, x.new = NULL, time.index = NULL) {
  time.index <- .fcst.time.index(time.index, NROW(y))
  .d <- .fcst.data(y, x, x.new)
  .fcst.require.regular(time.index[.d$rows], "fcst.break.size")
  y <- .d$y
  x <- .d$x
  x.new <- .d$x.new
  N <- .d$N

  if (!is.numeric(bp) || !length(bp) || any(!is.finite(bp)) ||
      any(bp != floor(bp)) || any(bp < 1 | bp >= N) || anyDuplicated(bp)) {
    stop("ERROR! fcst.break.size: bp must contain distinct integer dates in [1, N)")
  }
  if (nrow(x.new) != 1L || any(!is.finite(x.new))) {
    stop("ERROR! fcst.break.size: one finite forecast row is required")
  }
  bp <- sort(bp)
  edges <- c(0, bp, N)
  n <- length(bp)

  betas <- vector("list", n + 1)
  sigmas <- numeric(n + 1)

  for (i in seq_len(n + 1)) {
    idx <- (edges[i] + 1):edges[i + 1]
    if (length(idx) <= .d$k || qr(x[idx, , drop = FALSE])$rank < .d$k) {
      stop("ERROR! fcst.break.size: every regime needs a full-rank design and positive residual degrees of freedom")
    }
    .m <- OLS.reg(y[idx, , drop = FALSE], x[idx, , drop = FALSE])
    betas[[i]] <- .m$coefficients
    sigmas[i] <- sqrt(sum(.m$residuals^2, na.rm = TRUE) / length(idx))
  }

  omega.xx <- crossprod(x) / N
  scale <- drop(sqrt(x.new %*% solve(omega.xx) %*% t(x.new)))
  if (!is.finite(scale) || scale <= 0 || any(!is.finite(sigmas)) ||
      any(sigmas <= 0)) {
    stop("ERROR! fcst.break.size: positive forecast scale and regime residual variances required")
  }

  lambda <- lapply(
    seq_len(n),
    function(i) (betas[[i]] - betas[[n + 1]]) / sigmas[n + 1]
  )
  phi <- sapply(lambda, function(l) drop(x.new %*% l) / scale)

  list(
    phi = phi,
    q = if (n == 1) sigmas[1] / sigmas[2] else NA_real_,
    lambda = lambda,
    sigma = sigmas
  )
}

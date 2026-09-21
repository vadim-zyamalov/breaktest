#' @title
#' Forecast combination across estimation windows
#'
#' @description
#' Averages the forecasts produced by the same model over a set of nested
#' estimation windows,
#' \deqn{\hat y_{N+1} = \frac{1}{m} \sum_{i=1}^m \hat y_{N+1}(w_i),
#' \qquad w_i = w_{min} + \frac{i - 1}{m - 1} (1 - w_{min}),}
#' where \eqn{w_i} is the window length as a fraction of the sample.
#'
#' Pesaran and Pick (2011) establish bias and MSFE comparisons under their
#' break-model assumptions and study their robustness in simulations. These
#' are not uniform guarantees for arbitrary regressions or break sizes.
#' Averaging requires no estimate of the break date or magnitude, but the
#' chosen range of windows still affects the forecast.
#'
#' @param y A dependent variable.
#' @param x Explanatory variables. If `NULL` a single constant is used.
#' @param x.new A (1 x k) vector of regressor values for the forecast origin.
#' @param w.min The shortest window as a fraction of the sample.
#' @param w.max The longest window as a fraction of the sample.
#' @param m A number of windows. If `NULL` every admissible window is used,
#' i.e. the windows are one observation apart.
#' @param weights Optional finite nonnegative relative combination weights,
#' with at least one positive entry. They are normalised to sum to one.
#' Defaults to equal weights as in the source paper.
#'
#' @return An object of class `bt_fcstCombine`, a list of the combined
#' `forecast`, the vector of individual `forecasts`, the window lengths used and
#' the implied per-observation averaging weights. For an intercept-only model
#' these reproduce the forecast directly; for a general regression they
#' describe window membership and are not an equivalent WLS fit.
#'
#' @references
#' Pesaran, M. Hashem, and Andreas Pick.
#' “Forecast Combination Across Estimation Windows.”
#' Journal of Business and Economic Statistics 29, no. 2 (2011): 307–18.
#'
#' @export
fcst.AveW <- function(
  y,
  x = NULL,
  x.new = NULL,
  w.min = 0.1,
  w.max = 1,
  m = NULL,
  weights = NULL
) {
  .d <- .fcst.data(y, x, x.new)
  y <- .d$y
  x <- .d$x
  x.new <- .d$x.new
  N <- .d$N
  k <- .d$k
  if (!is.numeric(w.min) || length(w.min) != 1L || !is.finite(w.min) ||
      !is.numeric(w.max) || length(w.max) != 1L || !is.finite(w.max) ||
      w.min <= 0 || w.min > w.max || w.max > 1) {
    stop("ERROR! fcst.AveW: require 0 < w.min <= w.max <= 1")
  }
  if (!is.null(m) && (!is.numeric(m) || length(m) != 1L || !is.finite(m) ||
      m < 1 || m != floor(m))) stop("ERROR! fcst.AveW: m must be a positive integer")

  R.min <- max(round(w.min * N), k + 1)
  R.max <- min(round(w.max * N), N)

  if (R.min > R.max) {
    stop("ERROR! fcst.AveW: empty set of estimation windows")
  }

  R.all <- if (is.null(m)) {
    R.min:R.max
  } else {
    unique(round(seq(R.min, R.max, length.out = m)))
  }

  if (is.null(weights)) {
    weights <- rep(1 / length(R.all), length(R.all))
  }
  if (!is.numeric(weights) || length(weights) != length(R.all) ||
      any(!is.finite(weights)) || any(weights < 0) || !any(weights > 0)) {
    stop("ERROR! fcst.AveW: weights must be finite, nonnegative, nonzero and conformable with the windows")
  }
  weights <- weights / max(weights)
  weights <- weights / sum(weights)

  fc <- rep(NA_real_, length(R.all))
  coefficients <- numeric(k)
  obs.w <- numeric(N)

  for (i in seq_along(R.all)) {
    idx <- (N - R.all[i] + 1):N
    if (qr(x[idx, , drop = FALSE])$rank < k) {
      stop("ERROR! fcst.AveW: singular design in an estimation window")
    }
    .m <- OLS.reg(y[idx, , drop = FALSE], x[idx, , drop = FALSE])
    fc[i] <- .fcst.predict(.m, x.new)
    coefficients <- coefficients + weights[i] * .m$coefficients
    obs.w[idx] <- obs.w[idx] + weights[i] / R.all[i]
  }

  result <- list(
    forecast = drop(crossprod(weights, fc)),
    coefficients = coefficients,
    forecasts = fc,
    R.all = R.all,
    comb.weights = weights,
    obs.weights = obs.w,
    N = N
  )

  class(result) <- "bt_fcstCombine"
  result
}


#' @title
#' Forecast averaging over a quasi-likelihood ratio region
#'
#' @description
#' The AveLR scheme of Koo and Seo (2015). Instead of conditioning on a point
#' estimate of the break date, the forecast is averaged over every date that is
#' not rejected by the rescaled quasi-likelihood ratio statistic:
#' \deqn{\Phi = \left\{\gamma : \log S_N(\gamma) - \log S_N(\hat\gamma)
#' \leq c \, N^{-2/3}\right\}.}
#'
#' The motivation is that when the single break model is only an approximation
#' to a continuously changing data generating process ("strong
#' mis-specification"), the oracle property fails and the convergence rate of
#' the break date estimator drops to \eqn{N^{-1/3}} or slower. A point estimate
#' therefore cannot be trusted and the rescaling by \eqn{N^{-2/3}} rather than
#' \eqn{N^{-1}} reflects the slower rate.
#'
#' @param y A dependent variable.
#' @param x Explanatory variables. If `NULL` a single constant is used.
#' @param x.new A (1 x k) vector of regressor values for the forecast origin.
#' @param trim A trimming fraction bounding the admissible break dates.
#' @param c A tuning constant. Koo and Seo (2015) suggest 5 as a rule of thumb
#' after extensive experimentation and report that the results are not very
#' sensitive to it.
#' @param break.vars `NULL` or all columns of `x`. The published AveLR
#' uses unrestricted coefficients in both regimes; partial changes are rejected.
#' @param time.index Original integer-period positions. Internal gaps in the
#' estimation block are not supported.
#'
#' @return An object of class `bt_fcstCombine` with the averaged `forecast`,
#' the LR averaging region `set` and the individual `forecasts`. The tuning
#' constant does not define a confidence set with a nominal coverage level.
#'
#' @references
#' Koo, Bonsoo, and Myung Hwan Seo.
#' “Structural-Break Models under Mis-Specification:
#' Implications for Forecasting.”
#' Journal of Econometrics 188, no. 1 (2015): 166–81.
#'
#' @export
fcst.AveLR <- function(
  y,
  x = NULL,
  x.new = NULL,
  trim = 0.15,
  c = 5,
  break.vars = NULL,
  time.index = NULL
) {
  time.index <- .fcst.time.index(time.index, NROW(y))
  .d <- .fcst.data(y, x, x.new)
  .fcst.inference.data(.d)
  .fcst.require.regular(time.index[.d$rows], "fcst.AveLR")
  y <- .d$y
  x <- .d$x
  x.new <- .d$x.new
  N <- .d$N
  k <- .d$k
  .fcst.window.trim(trim)
  if (!is.numeric(c) || length(c) != 1L || !is.finite(c) || c < 0) {
    stop("ERROR! fcst.AveLR: c must be finite and nonnegative")
  }
  if (!is.null(break.vars) &&
      !identical(.fcst.break.columns(break.vars, k), seq_len(k))) {
    stop("ERROR! fcst.AveLR: only unrestricted changes in all coefficients are supported")
  }
  break.vars <- NULL

  bp.min <- max(floor(trim * N), k + 1)
  bp.max <- min(ceiling((1 - trim) * N), N - k - 1)

  if (bp.min > bp.max) {
    stop("ERROR! fcst.AveLR: trimming leaves no admissible break dates")
  }

  bp.all <- bp.min:bp.max

  # Unrestricted changes: use the package-wide recursive SSR matrix.
  SSR.data <- .fcst.ssr.matrix(y, x, k + 1)
  ssr <- SSR.data[1, bp.all] + SSR.data[cbind(bp.all + 1, N)]

  if (!any(is.finite(ssr))) stop("ERROR! fcst.AveLR: no identified candidate break")
  best <- min(ssr)
  lr <- if (best == 0) ifelse(ssr == 0, 0, Inf) else log(ssr) - log(best)
  keep <- bp.all[lr <= c * N^(-2 / 3)]

  if (length(keep) == 0) {
    keep <- bp.all[which.min(ssr)]
  }

  coefs <- vapply(keep, function(bp) {
    .s <- .fcst.split.design(x, bp, break.vars)
    .m <- OLS.reg(y, .s$xs)
    as.numeric(.m$coefficients[.s$map.2])
  }, numeric(k))
  coefs <- matrix(coefs, nrow = k)
  fc <- drop(x.new %*% coefs)

  result <- list(
    forecast = mean(fc),
    coefficients = rowMeans(coefs),
    forecasts = fc,
    set = keep,
    bp = bp.all[which.min(ssr)],
    LR = lr,
    bp.all = bp.all,
    N = N
  )

  class(result) <- "bt_fcstCombine"
  result
}


#' @title
#' Limiting distribution of the break date estimator of Bai (1997)
#' @name fcst.Bai.cdf
#'
#' @description
#' The distribution function of \eqn{\arg\max_s \{W(s) - |s| / 2\}} generalised
#' to different pre- and post-break second moments and error variances. A port
#' of the `bai_qnt` procedure of the GAUSS replication code
#' `forecast_breaks_R1.prg`.
#'
#' @param x A vector of arguments.
#' @param xi,phi Post/pre ratios of projected regressor second moments and
#' projected long-run score variances, respectively. With iid errors the
#' latter includes both the error-variance and second-moment ratios.
#' Both equal to one gives the classical symmetric case.
#'
#' @return `.bai.cdf` returns a vector of probabilities.
#'
#' @references
#' Bai, Jushan.
#' “Estimation of a Change Point in Multiple Regression Models.”
#' The Review of Economics and Statistics 79, no. 4 (1997): 551–63.
#'
#' @importFrom stats pnorm
#'
#' @keywords internal
.bai.cdf <- function(x, xi = 1, phi = 1) {
  if (!is.numeric(x) || anyNA(x) || !is.numeric(xi) || length(xi) != 1L ||
      !is.finite(xi) || xi <= 0 || !is.numeric(phi) || length(phi) != 1L ||
      !is.finite(phi) || phi <= 0) {
    stop("ERROR! .bai.cdf: x must be numeric and xi, phi finite and positive")
  }
  res <- numeric(length(x))

  neg <- is.finite(x) & x < 0
  if (any(neg)) {
    xn <- x[neg]
    xp <- xi / phi
    a <- xp * (1 + xp) / 2
    b <- 0.5 + xp
    cc <- phi * (phi + 2 * xi) / ((phi + xi) * xi)
    d <- (phi + 2 * xi)^2 / (xi * (phi + xi))

    res[neg] <- -sqrt(-xn / (2 * pi)) * exp(xn / 8) -
      cc * exp(-xn * a + pnorm(-b * sqrt(-xn), log.p = TRUE)) +
      (d - 2 - xn / 2) * pnorm(-sqrt(-xn) / 2)
  }

  pos <- is.finite(x) & x > 0
  if (any(pos)) {
    xq <- x[pos]
    a <- (phi + xi) / 2
    b <- (2 * phi + xi) / (2 * sqrt(phi))
    cc <- xi * (2 * phi + xi) / (phi * (phi + xi))
    d <- (2 * phi + xi)^2 / (phi * (phi + xi))

    res[pos] <- 1 + xi * sqrt(xq / (2 * pi * phi)) *
      exp(-xq * xi^2 / (8 * phi)) +
      cc * exp(xq * a + pnorm(-sqrt(xq) * b, log.p = TRUE)) +
      (-d + 2 - xq * xi^2 / (2 * phi)) * pnorm(-sqrt(xq / phi) * xi / 2)
  }

  if (any(x == 0)) {
    # Probability that the maximum occurs on the left half-line.
    res[x == 0] <- xi / (xi + phi)
  }
  res[is.infinite(x) & x < 0] <- 0
  res[is.infinite(x) & x > 0] <- 1

  pmin(pmax(res, 0), 1)
}


#' @rdname fcst.Bai.cdf
#'
#' @param level A confidence level.
#' @param range Initial bracketing half-width; doubled until both requested
#' tails are covered. It is not a truncation of the limiting distribution.
#' @param step Absolute root-finding tolerance (retained argument name).
#'
#' @return `.bai.quantiles` returns the pair of quantiles bounding the
#' confidence interval.
#'
#' @keywords internal
.bai.quantiles <- function(level = 0.95, xi = 1, phi = 1, range = 50, step = 0.01) {
  .fcst.test.level(level)
  if (!is.numeric(range) || length(range) != 1L || !is.finite(range) || range <= 0 ||
      !is.numeric(step) || length(step) != 1L || !is.finite(step) || step <= 0) {
    stop("ERROR! .bai.quantiles: range and step must be finite and positive")
  }
  probabilities <- c(lower = (1 - level) / 2, upper = (1 + level) / 2)
  endpoints <- c(-range, range)
  for (i in 1:2) {
    for (j in seq_len(100)) {
      probability <- .bai.cdf(endpoints[i], xi, phi)
      if (!is.finite(probability)) stop("ERROR! .bai.quantiles: nonfinite limiting CDF")
      if ((i == 1 && probability <= probabilities[1]) ||
          (i == 2 && probability >= probabilities[2])) break
      endpoints[i] <- endpoints[i] * 2
      if (j == 100) stop("ERROR! .bai.quantiles: could not bracket both limiting tails")
    }
  }
  vapply(probabilities, function(p) {
    stats::uniroot(function(z) .bai.cdf(z, xi, phi) - p, endpoints,
                   tol = step, check.conv = TRUE)$root
  }, 0.0)
}


#' @title
#' Confidence interval for the break date of Bai (1997)
#'
#' @param y A dependent variable.
#' @param x Explanatory variables.
#' @param bp An LS estimate of the break date, or an estimate with established
#' asymptotic equivalence. If `NULL`, use [fcst.window.LS]. An arbitrary known
#' date does not have this estimator's sampling distribution. The published
#' Altansukh--Osborn forecasting algorithm may supply its own estimated date
#' while constructing this interval on the original data.
#' @param level A confidence level.
#' @param het Whether the pre- and post-break second moments and variances are
#' allowed to differ. With `FALSE` common moments are estimated by pooling
#' regime-residual contributions and the symmetric distribution is imposed.
#' @param trim Trimming used to estimate an unknown break date. Interval bounds
#' themselves are clipped only to `[1, N-1]`.
#' @param vcov.type "iid" assumes homoskedastic uncorrelated errors within
#' each regime and uses sigma_j^2 Q_j. "HC" estimates score second moments;
#' "HAC" estimates long-run score covariances for serially correlated errors.
#' @param bandwidth,kernel,time.index Settings of [HAC.variance]. Original
#' integer time positions must form a contiguous estimation block.
#'
#' @return A list of the break date `bp`, the interval bounds `lower`, `upper`
#' and the scaling factor `L`.
#'
#' @details
#' The interval is based on
#' \deqn{\frac{(\Delta' Q \Delta)^2}{\Delta' \Omega \Delta}
#' (\hat T_b - T_b^0) \Rightarrow \arg\max_s \{W(s) - |s| / 2\},}
#' with \eqn{\Delta = \hat\beta_2 - \hat\beta_1}.
#' For heterogeneous regimes, xi=(Delta'Q2 Delta)/(Delta'Q1 Delta) and
#' phi=(Delta'Omega2 Delta)/(Delta'Omega1 Delta). HC and HAC use the
#' projected scores (x_t'Delta)*u_t, rather than substituting an error-only
#' long-run variance. These estimators require the corresponding moment,
#' weak-dependence and shrinking-break assumptions of Bai (1997). A zero
#' estimated break or nonpositive projected score variance is unidentified
#' for this approximation and produces an explicit error.
#' This implementation uses stabilising moments within regimes. It does not
#' implement Bai's separate local normalisation for deterministic trends or
#' provide coverage after arbitrary same-sample regressor selection. A supplied
#' date's provenance and these model assumptions cannot be verified from its
#' integer value; passing `bp` asserts the stated estimator contract.
#'
#' Note that the coverage of this interval is known to fall well below the
#' nominal level when the break is small, which is exactly why Altansukh and
#' Osborn (2022) prefer the likelihood-ratio confidence set of Eo and Morley
#' (2015) for forecast averaging; see [fcst.conf.set.EoMorley].
#'
#' @export
fcst.conf.int.Bai <- function(
  y,
  x = NULL,
  bp = NULL,
  level = 0.95,
  het = TRUE,
  trim = 0.15,
  vcov.type = c("iid", "HC", "HAC"),
  bandwidth = NULL,
  kernel = "Bartlett",
  time.index = NULL
) {
  vcov.type <- match.arg(vcov.type)
  .fcst.test.level(level)
  .fcst.window.trim(trim)
  if (!is.logical(het) || length(het) != 1L || is.na(het)) {
    stop("ERROR! fcst.conf.int.Bai: het must be TRUE or FALSE")
  }
  time.index <- .fcst.time.index(time.index, NROW(y))
  .d <- .fcst.data(y, x)
  y <- .d$y
  x <- .d$x
  N <- .d$N
  k <- .d$k
  time.index <- time.index[.d$rows]
  .fcst.require.regular(time.index, "fcst.conf.int.Bai")
  .fcst.inference.data(.d)

  if (is.null(bp)) {
    bp <- fcst.window.LS(y, x, trim)$bp
  }
  if (!is.numeric(bp) || length(bp) != 1L || !is.finite(bp) ||
      bp != floor(bp) || bp <= k || N - bp <= k) {
    stop("ERROR! fcst.conf.int.Bai: both regimes need more than k observations")
  }

  i1 <- 1:bp
  i2 <- (bp + 1):N
  if (qr(x[i1, , drop = FALSE])$rank < k || qr(x[i2, , drop = FALSE])$rank < k) {
    stop("ERROR! fcst.conf.int.Bai: both regime designs must have full rank")
  }

  m1 <- OLS.reg(y[i1, , drop = FALSE], x[i1, , drop = FALSE])
  m2 <- OLS.reg(y[i2, , drop = FALSE], x[i2, , drop = FALSE])

  delta <- m2$coefficients - m1$coefficients

  q1 <- crossprod(x[i1, , drop = FALSE]) / length(i1)
  q2 <- crossprod(x[i2, , drop = FALSE]) / length(i2)

  s1 <- mean(m1$residuals^2, na.rm = TRUE)
  s2 <- mean(m2$residuals^2, na.rm = TRUE)

  dq1 <- drop(delta %*% q1 %*% delta)
  dq2 <- drop(delta %*% q2 %*% delta)
  score <- c(drop(x[i1, , drop = FALSE] %*% delta) * as.numeric(m1$residuals),
             drop(x[i2, , drop = FALSE] %*% delta) * as.numeric(m2$residuals))
  projected.Q <- c(pre = dq1, post = dq2)
  if (vcov.type == "iid") {
    projected.Omega <- c(pre = s1 * dq1, post = s2 * dq2)
  } else if (vcov.type == "HC") {
    projected.Omega <- c(pre = mean(score[i1]^2), post = mean(score[i2]^2))
  } else {
    if (is.null(bandwidth)) bandwidth <- max(1L, floor(min(bp, N - bp)^(1 / 4)))
    projected.Omega <- c(
      pre = HAC.variance(score[i1], bandwidth, kernel, demean = FALSE, time.index = time.index[i1]),
      post = HAC.variance(score[i2], bandwidth, kernel, demean = FALSE, time.index = time.index[i2])
    )
  }
  if (!het) {
    fraction <- c(bp, N - bp) / N
    projected.Q[] <- sum(fraction * projected.Q)
    # Under common homoskedastic errors estimate sigma^2 and Q separately.
    common <- if (vcov.type == "iid") {
      sum(fraction * c(s1, s2)) * projected.Q[1]
    } else sum(fraction * projected.Omega)
    projected.Omega[] <- common
  }
  if (any(!is.finite(projected.Q)) || any(projected.Q <= 0) ||
      any(!is.finite(projected.Omega)) || any(projected.Omega <= 0)) {
    stop("ERROR! fcst.conf.int.Bai: nonzero break and positive projected score variances required")
  }

  xi <- unname(projected.Q[2] / projected.Q[1])
  phi <- unname(projected.Omega[2] / projected.Omega[1])

  # Scaling of the argmax distribution
  L <- unname(projected.Q[1]^2 / projected.Omega[1])
  if (!is.finite(L) || L <= 0) stop("ERROR! fcst.conf.int.Bai: invalid distribution scale")

  qs <- .bai.quantiles(level, xi, phi)

  list(
    bp = bp,
    lower = max(1, bp - ceiling(qs["upper"] / L)),
    upper = min(N - 1, bp - floor(qs["lower"] / L)),
    L = L,
    xi = xi,
    phi = phi,
    quantiles = qs,
    projected.Q = projected.Q,
    projected.Omega = projected.Omega,
    vcov.type = vcov.type,
    het = het,
    bandwidth = bandwidth,
    kernel = kernel,
    time.index = time.index,
    rows = .d$rows,
    level = level
  )
}


#' @title
#' Forecast averaging over a confidence set for the break date
#'
#' @description
#' The confidence set forecast of Altansukh and Osborn (2022):
#' \deqn{\hat y_{N+1} = |C|^{-1}\sum_{b\in C}
#' \hat\beta_{b+1:N}' x_{N+1}.}
#' Every accepted break date is used once; gaps in a disconnected set are
#' preserved. A date `b` is the last pre-break observation, so its estimation
#' window starts at `b + 1`.
#'
#' The scheme is automatically adaptive. A large break gives a narrow interval
#' and the forecast is close to the post-break one; a small break gives a wide
#' interval and the forecast approaches plain window averaging.
#'
#' @param y A dependent variable.
#' @param x Explanatory variables. If `NULL` a single constant is used.
#' @param x.new A (1 x k) vector of regressor values for the forecast origin.
#' @param interval An optional two-element vector with break-date interval
#' bounds. Retained for compatibility; do not also supply `set`.
#' @param level A confidence level for the set or interval.
#' @param trim A trimming fraction bounding the admissible window starts.
#' @param conf.method `eo.morley` (default) uses [fcst.conf.set.EoMorley];
#' `bai` retains the previous [fcst.conf.int.Bai] method. Ignored for supplied
#' `interval` or `set`.
#' @param set An optional vector of accepted break dates, possibly with gaps.
#' Dates index the complete-case sample. Following the source rule for window
#' starts, require `b + 1 <= N * (1 - trim)` and more than k post-break rows.
#' An empty remaining set is an error.
#' The temporal dispatcher may attach a named `multiplicity` attribute when
#' distinct original dates map to the same complete-case split. Its positive
#' integer counts preserve equal weights over the original dates.
#' @param ... Arguments for the selected confidence-set routine, e.g.
#' `variance`, `vcov.type` and `bandwidth` for Eo-Morley or `het` for Bai.
#' @param time.index Original increasing integer time positions, filtered
#' together with incomplete rows before confidence-set inference.
#' Inferred sets require a contiguous estimation block. A supplied set only
#' defines a combination: its confidence interpretation belongs to its source.
#'
#' @return An object of class `bt_fcstCombine`. `break.set` contains the dates
#' used, `set` the corresponding window starts (as in earlier versions),
#' `interval` their envelope and `confidence.set` the inference result when
#' computed internally. `comb.weights` weights the distinct split forecasts;
#' these are equal unless original-date multiplicities were supplied by the
#' temporal dispatcher. `set.multiplicity` records those counts.
#' `source.break.set` retains the original dates; `excluded.break.set` records
#' the dates omitted by the deterministic window rule. The filtered collection
#' is not asserted to be a new confidence set with the original coverage.
#'
#' @references
#' Altansukh, Gantungalag, and Denise R. Osborn.
#' “Using Structural Break Inference for Forecasting Time Series.”
#' Empirical Economics 63 (2022): 1–41.
#'
#' @export
fcst.CIavg <- function(
  y,
  x = NULL,
  x.new = NULL,
  interval = NULL,
  level = 0.95,
  trim = 0.1,
  conf.method = c("eo.morley", "bai"),
  set = NULL,
  time.index = NULL,
  ...
) {
  conf.method <- match.arg(conf.method)
  if (!is.null(interval) && !is.null(set)) {
    stop("ERROR! fcst.CIavg: supply either interval or set, not both")
  }
  if (!is.numeric(trim) || length(trim) != 1 || !is.finite(trim) ||
      trim <= 0 || trim >= 0.5) {
    stop("ERROR! fcst.CIavg: trim should be in (0, 0.5)")
  }
  time.index <- .fcst.time.index(time.index, NROW(y))
  .d <- .fcst.data(y, x, x.new)
  y <- .d$y
  x <- .d$x
  x.new <- .d$x.new
  N <- .d$N
  k <- .d$k
  time.index <- time.index[.d$rows]
  if (!is.numeric(y) || !is.numeric(x) || !is.numeric(x.new) ||
      any(!is.finite(y)) || any(!is.finite(x)) || any(!is.finite(x.new)) ||
      nrow(x.new) != 1) {
    stop("ERROR! fcst.CIavg: finite numeric data and one row of x.new required")
  }

  .ci <- NULL
  if (is.null(interval) && is.null(set)) {
    if (conf.method == "eo.morley") {
      .ci <- fcst.conf.set.EoMorley(y, x, level = level, trim = trim,
                                   time.index = time.index, ...)
      set <- .ci$set
    } else {
      .ci <- fcst.conf.int.Bai(y, x, level = level, trim = trim,
                              time.index = time.index, ...)
      interval <- c(.ci$lower, .ci$upper)
    }
  }
  if (!is.null(interval)) {
    if (!is.numeric(interval) || length(interval) != 2 ||
        any(!is.finite(interval)) || interval[1] > interval[2]) {
      stop("ERROR! fcst.CIavg: invalid interval bounds")
    }
    lo <- max(1, round(interval[1]))
    hi <- min(N - 1, round(interval[2]))
    set <- if (lo <= hi) seq.int(lo, hi) else integer()
  } else if (!is.numeric(set) || any(!is.finite(set)) ||
             any(set != floor(set)) || any(set < 1) || any(set >= N)) {
    stop("ERROR! fcst.CIavg: set must contain integer break dates in [1, N)")
  }
  multiplicity <- attr(set, "multiplicity", exact = TRUE)
  if (!is.null(multiplicity)) {
    labels <- names(multiplicity)
    if (!is.numeric(multiplicity) || length(multiplicity) != length(unique(set)) ||
        any(!is.finite(multiplicity)) || any(multiplicity <= 0) ||
        any(multiplicity != floor(multiplicity)) || is.null(labels) ||
        anyNA(labels) || anyDuplicated(labels) ||
        !setequal(labels, as.character(unique(set)))) {
      stop("ERROR! fcst.CIavg: multiplicity must contain positive integer counts named by each distinct set date")
    }
  }
  source.dates <- sort(unique(as.numeric(set)))
  dates <- source.dates[source.dates + 1 <= N * (1 - trim) & source.dates < N - k]
  if (length(dates) == 0) {
    stop("ERROR! fcst.CIavg: no admissible dates remain in the confidence set")
  }
  multiplicity <- if (is.null(multiplicity)) {
    setNames(rep(1, length(dates)), as.character(dates))
  } else multiplicity[as.character(dates)]
  weights <- multiplicity / max(multiplicity)
  weights <- weights / sum(weights)

  starts <- dates + 1
  coefs <- vapply(starts, function(s) {
    idx <- s:N
    if (qr(x[idx, , drop = FALSE])$rank < k) {
      stop("ERROR! fcst.CIavg: unidentified forecast component; no dates are silently discarded")
    }
    as.numeric(OLS.reg(y[idx, , drop = FALSE], x[idx, , drop = FALSE])$coefficients)
  }, numeric(k))
  coefs <- matrix(coefs, nrow = k)
  fc <- drop(x.new %*% coefs)

  result <- list(
    forecast = sum(weights * fc),
    coefficients = drop(coefs %*% weights),
    forecasts = fc,
    comb.weights = unname(weights),
    set.multiplicity = multiplicity,
    set = starts,
    break.set = dates,
    source.break.set = source.dates,
    excluded.break.set = setdiff(source.dates, dates),
    interval = range(dates),
    confidence.set = .ci,
    conf.method = if (is.null(.ci)) "supplied" else conf.method,
    N = N
  )

  class(result) <- "bt_fcstCombine"
  result
}

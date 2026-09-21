#' @title
#' Forecast evaluation in the presence of instabilities
#' @name fcst.evaluation
#'
#' @description
#' Standard average measures of predictive ability are misleading when the
#' forecasting performance itself is unstable: positive and negative loss
#' differentials in different sub-samples cancel out. Rossi (2021) shows a
#' striking example where the classical Mincer-Zarnowitz test does not reject
#' rationality of the Federal Reserve and SPF inflation forecasts
#' (\eqn{p = 0.81} and \eqn{0.80}) purely because the systematic
#' under-prediction of the 1970s offsets the over-prediction after 1980.
#'
#' It is worth stressing that breaks in the *parameters* are neither necessary
#' nor sufficient for instability in the *forecasting performance*, so one
#' should test the latter directly rather than the former.
#'
#' @references
#' Diebold, Francis X., and Roberto S. Mariano.
#' “Comparing Predictive Accuracy.”
#' Journal of Business and Economic Statistics 13, no. 3 (1995): 253–63.
#'
#' Clark, Todd E., and Kenneth D. West.
#' “Approximately Normal Tests for Equal Predictive Accuracy
#' in Nested Models.”
#' Journal of Econometrics 138, no. 1 (2007): 291–311.
#'
#' Giacomini, Raffaella, and Barbara Rossi.
#' “Forecast Comparisons in Unstable Environments.”
#' Journal of Applied Econometrics 25, no. 4 (2010): 595–620.
#'
#' Rossi, Barbara, and Tatevik Sekhposyan.
#' “Forecast Rationality Tests in the Presence of Instabilities,
#' with Applications to Federal Reserve and Survey Forecasts.”
#' Journal of Applied Econometrics 31, no. 3 (2016): 507–32.
#'
#' Rossi, Barbara (2021). "Forecasting in the Presence of Instabilities:
#' How We Know Whether Models Predict Well and How to Improve Them."
#' Journal of Economic Literature 59, no. 4: 1135–90.
#' https://doi.org/10.1257/jel.20201479.
#'
#' Harvey, D., S. Leybourne and P. Newbold (1997). "Testing the equality of
#' prediction mean squared errors." International Journal of Forecasting
#' 13, 281-291. https://doi.org/10.1016/S0169-2070(96)00719-4.
#' Harvey, D., S. Leybourne and E. Whitehouse (2017). "Forecast evaluation
#' tests and negative long-run variance estimates in small samples."
#' International Journal of Forecasting 33, 833-847.
#' https://doi.org/10.1016/j.ijforecast.2017.05.001.
#'
#' @keywords internal
NULL


#' @rdname fcst.evaluation
#' @order 1
#'
#' @param loss.1,loss.2 Sequences of forecast error losses of the two competing
#' methods, e.g. squared forecast errors.
#' @param bandwidth Bartlett bandwidth B: lags 1 through B-1 receive weights
#' 1-j/B. If `NULL`, max(h, floor(P^(1/4))) is used, with P the variance
#' estimation sample size. For HLN, only the rectangular estimator through
#' h-1 is allowed; bandwidth must be NULL or h.
#' @param hln Whether the Harvey-Leybourne-Newbold small sample correction
#' should be applied to the rectangular long-run variance through h-1,
#' with Student's t reference distribution with P-1 degrees of freedom.
#' This assumes the loss differential is (h-1)-dependent. The correction
#' does not apply to Bartlett HAC (Harvey et al., 2017). Nonpositive
#' rectangular variance estimates are refused, without a substitute test.
#'
#' @return `fcst.test.DM` returns an object of class `bt_fcstDM` with the
#' statistic, its one- and two-sided p-values and the mean loss differential.
#'
#' @details
#' `fcst.test.DM` is the Diebold-Mariano statistic
#' \eqn{\sqrt{P} \bar{\Delta L} / \sqrt{\hat V}}, asymptotically standard
#' normal under equal predictive accuracy, a suitable central limit theorem
#' for the loss differential and a consistent long-run variance estimate.
#' These assumptions are not verified from the supplied forecasts. In
#' particular, arbitrary nesting, model selection and parameter instability
#' do not automatically satisfy them.
#'
#' @param h Positive integer forecast horizon. The default HAC bandwidth is
#' at least h, so lag products through h-1 are included.
#' @param time.index Original integer period positions for HAC lag products.
#' Every test requires complete, consecutive evaluation dates. Missing
#' outcomes or forecast failures must not be silently selected out of a
#' planned evaluation sample when interpreting its nominal calibration.
#'
#' @importFrom stats pnorm pt
#'
#' @export
fcst.test.DM <- function(loss.1, loss.2, bandwidth = NULL, hln = FALSE,
                         h = 1, time.index = NULL) {
  if (!is.numeric(loss.1) || !is.numeric(loss.2) ||
      length(loss.1) != length(loss.2)) {
    stop("ERROR! fcst.test.DM: losses must be numeric and have equal lengths")
  }
  if (!is.logical(hln) || length(hln) != 1L || is.na(hln)) {
    stop("ERROR! fcst.test.DM: hln must be TRUE or FALSE")
  }
  h <- .fcst.horizons(h)
  if (length(h) != 1L) stop("ERROR! fcst.test.DM: h must be a single horizon")
  time.index <- .fcst.time.index(time.index, length(loss.1))
  d <- as.numeric(loss.1) - as.numeric(loss.2)
  if (any(!is.finite(loss.1)) || any(!is.finite(loss.2)) || any(!is.finite(d))) {
    stop("ERROR! fcst.test.DM: inference requires complete finite losses")
  }
  .fcst.require.regular(time.index, "fcst.test.DM")
  P <- length(d)

  if (P < 3 || P <= h) {
    stop("ERROR! fcst.test.DM: not enough observations")
  }

  if (hln) {
    if (!is.null(bandwidth) && (!is.numeric(bandwidth) || length(bandwidth) != 1L ||
        !is.finite(bandwidth) || bandwidth != h)) {
      stop("ERROR! fcst.test.DM: HLN requires rectangular lags through h-1; bandwidth must be NULL or h")
    }
    bandwidth <- h
    centred <- d - mean(d)
    V <- sum(centred^2) / P
    if (h > 1L) {
      V <- V + 2 * sum(vapply(seq_len(h - 1L), function(j) {
        sum(centred[(j + 1L):P] * centred[seq_len(P - j)]) / P
      }, 0))
    }
  } else {
    bandwidth <- .fcst.evaluation.bandwidth(P, h, bandwidth)
    V <- HAC.variance(d, bandwidth, time.index = time.index)
  }
  if (!is.finite(V) || V <= 0) {
    stop("ERROR! fcst.test.DM: loss differential has no positive long-run variance")
  }
  stat <- sqrt(P) * mean(d) / sqrt(V)

  if (hln) {
    stat <- stat * sqrt((P + 1 - 2 * h + h * (h - 1) / P) / P)
  }

  result <- list(
    statistic = stat,
    p.value = if (hln) 2 * pt(-abs(stat), df = P - 1) else 2 * pnorm(-abs(stat)),
    p.value.one.sided = if (hln) pt(stat, df = P - 1) else pnorm(stat),
    mean.diff = mean(d),
    n.obs = P,
    bandwidth = bandwidth,
    h = h,
    hln = hln,
    variance = V,
    kernel = if (hln) "rectangular" else "Bartlett",
    test = "Diebold-Mariano",
    alternative = "less",
    time.index = time.index,
    assumptions.verified = FALSE,
    inference = if (hln) "HLN approximation under (h-1)-dependence of the loss differential"
      else "asymptotic normal calibration conditional on the DM central limit theorem and consistent HAC"
  )

  class(result) <- "bt_fcstDM"
  result
}


#' @rdname fcst.evaluation
#' @order 2
#'
#' @param forecast A sequence of forecasts of the larger (nesting) model.
#' @param actual A sequence of realised values.
#' @param benchmark Forecasts from the smaller nested model. A scalar is
#' repeated; the default zero retains the martingale-difference benchmark.
#'
#' @return `fcst.test.CW` returns a `bt_fcstCW` (also inheriting `bt_fcstDM`).
#' Its `p.value` is one-sided, against improvement by the larger model.
#'
#' @details
#' `fcst.test.CW` implements the Clark-West adjustment for nested models,
#' using the adjusted differential
#' \eqn{(y-f_0)^2-(y-f_1)^2+(f_0-f_1)^2}. Positive values favour the larger
#' model. This removes the null bias induced by estimation of its extra
#' parameters. The approximate normal calibration requires nested models
#' and the estimation and regularity conditions of Clark and West (2007);
#' it is not a general test for arbitrary competing forecasts. The paper's
#' formal critical-value argument covers one-step conditionally
#' homoskedastic forecasts, and its multistep/heteroskedastic extension has
#' one additional parameter. Broader configurations have simulation support,
#' not a universal normal-limit guarantee. Nested selection procedures need
#' a separate justification. The zero
#' benchmark special case is implemented in `CW_test_nanJAE.m`.
#'
#' @export
fcst.test.CW <- function(forecast, actual, bandwidth = NULL, h = 1,
                         time.index = NULL, benchmark = 0) {
  if (!is.numeric(forecast) || !is.numeric(actual) ||
      length(forecast) != length(actual)) {
    stop("ERROR! fcst.test.CW: forecast and actual must be numeric and have equal lengths")
  }
  if (!is.numeric(benchmark) || !length(benchmark) ||
      !(length(benchmark) %in% c(1L, length(actual)))) {
    stop("ERROR! fcst.test.CW: benchmark must be numeric, scalar or the length of actual")
  }
  h <- .fcst.horizons(h)
  if (length(h) != 1L) stop("ERROR! fcst.test.CW: h must be a single horizon")
  forecast <- as.numeric(forecast)
  actual <- as.numeric(actual)
  benchmark <- rep_len(as.numeric(benchmark), length(actual))
  time.index <- .fcst.time.index(time.index, length(actual))

  ok <- is.finite(forecast) & is.finite(actual) & is.finite(benchmark)
  if (!all(ok)) stop("ERROR! fcst.test.CW: inference requires complete finite forecasts and actuals")
  .fcst.require.regular(time.index, "fcst.test.CW")

  P <- length(actual)
  if (P < 3 || P <= h) {
    stop("ERROR! fcst.test.CW: not enough observations")
  }

  # Algebraically identical to the adjusted squared-error differential,
  # while avoiding cancellation between three potentially large squares.
  d <- 2 * (actual - benchmark) * (forecast - benchmark)
  bandwidth <- .fcst.evaluation.bandwidth(P, h, bandwidth)
  if (any(!is.finite(d))) stop("ERROR! fcst.test.CW: nonfinite adjusted loss differential")
  V <- HAC.variance(d, bandwidth, time.index = time.index)
  if (!is.finite(V) || V <= 0) {
    stop("ERROR! fcst.test.CW: adjusted loss differential has no positive long-run variance")
  }
  stat <- sqrt(P) * mean(d) / sqrt(V)

  result <- list(
    statistic = stat,
    p.value = pnorm(stat, lower.tail = FALSE),
    p.value.one.sided = pnorm(stat, lower.tail = FALSE),
    mean.diff = mean(d),
    n.obs = P,
    bandwidth = bandwidth,
    h = h,
    hln = FALSE,
    test = "Clark-West",
    alternative = "greater",
    time.index = time.index,
    benchmark = benchmark,
    assumptions.verified = FALSE,
    inference = "one-sided approximate normal calibration conditional on the Clark-West nesting and regularity conditions"
  )

  class(result) <- c("bt_fcstCW", "bt_fcstDM")
  result
}


#' @rdname fcst.evaluation
#' @order 3
#'
#' @param window A size of the rolling window the local statistic is computed
#' over.
#' @param level A nominal size, either 0.05 or 0.10.
#'
#' @return `fcst.test.fluctuation` returns an object of class `bt_fcstFluct`
#' with the whole path of the local statistic, its supremum, the critical value
#' and the test decision.
#'
#' @details
#' `fcst.test.fluctuation` is the Fluctuation test of Giacomini and Rossi
#' (2010): rolling sums of loss differentials standardised by a common
#' full-evaluation-sample long-run variance,
#' \eqn{\sup_t |F_{t,m}|}. It is designed for smooth and persistent changes in
#' relative performance and is what uncovers the "pockets of predictability"
#' that average tests miss. Unlike the local-variance variant in
#' `FluctuationL.m`, the common denominator follows the published test.
#' Critical values are provided by [get.cv.fluctuation].
#' Tabulated calibration uses m/N, where N is the full evaluation sample,
#' and requires complete, equally spaced observations and 0.1 <= m/N <= 0.9.
#' Missing observations, gaps and invalid local variances are refused rather
#' than silently changing the set of windows. Critical values are linearly
#' interpolated Monte Carlo approximations to the asymptotic distribution.
#' Their validity requires the source paper's functional central limit
#' theorem and stable long-run variance under the null. Arbitrary variance
#' breaks or arbitrary expanding/adaptive model estimates are not covered
#' merely because the loss sequences are finite.
#'
#' @export
fcst.test.fluctuation <- function(
  loss.1,
  loss.2,
  window,
  level = 0.05,
  bandwidth = NULL,
  h = 1,
  time.index = NULL
) {
  if (!is.numeric(loss.1) || !is.numeric(loss.2) ||
      length(loss.1) != length(loss.2)) {
    stop("ERROR! fcst.test.fluctuation: losses must be numeric and have equal lengths")
  }
  d <- as.numeric(loss.1) - as.numeric(loss.2)
  N <- length(d)
  time.index <- .fcst.time.index(time.index, N)

  .fcst.evaluation.window(window, N, time.index, minimum = 3L)
  if (any(!is.finite(d))) {
    stop("ERROR! fcst.test.fluctuation: tabulated calibration requires complete finite losses")
  }
  cv <- get.cv.fluctuation(window / N, level)
  bandwidth <- .fcst.evaluation.bandwidth(N, h, bandwidth)
  V <- HAC.variance(d, bandwidth, time.index = time.index)
  if (!is.finite(V) || V <= 0) {
    stop("ERROR! fcst.test.fluctuation: no positive full-sample long-run variance")
  }
  path <- sapply(window:N, function(j) {
    idx <- (j - window + 1):j
    dj <- d[idx]
    sqrt(window) * mean(dj) / sqrt(V)
  })
  if (any(!is.finite(path))) {
    stop("ERROR! fcst.test.fluctuation: nonfinite window statistic")
  }

  i.max <- which.max(abs(path))

  result <- list(
    statistic = path[i.max],
    at = window + i.max - 1,
    path = path,
    index = window:N,
    cv = cv,
    reject = abs(path[i.max]) > cv,
    window = window,
    level = level,
    mu = window / N,
    h = h,
    bandwidth = bandwidth,
    n.obs = N,
    time.index = time.index,
    calibrated = TRUE,
    assumptions.verified = FALSE,
    variance = V,
    variance.scope = "full evaluation sample",
    inference = "approximate tabulated calibration conditional on the source FCLT and stable null long-run variance",
    test = "Giacomini-Rossi fluctuation"
  )

  class(result) <- "bt_fcstFluct"
  result
}


#' @rdname fcst.evaluation
#' @order 4
#'
#' @param fc.error A sequence of forecast errors, actual minus forecast.
#' @param x1 A matrix of regressors whose coefficients are tested to be zero.
#' For the Mincer-Zarnowitz regression this is a constant and the forecast.
#' @param x2 An optional matrix of the remaining regressors, not tested.
#'
#' @return `fcst.test.rationality` returns an object of class `bt_fcstFluct`.
#'
#' @details
#' `fcst.test.rationality` is the Fluctuation rationality test of Rossi and
#' Sekhposyan (2016): a Wald test of \eqn{H_0: \beta = 0} in
#' \eqn{e_{t+h} = x_{1t}' \beta + x_{2t}' \gamma + u_t} recomputed in rolling
#' windows, with a HAC covariance matrix, and the supremum taken over time.
#' Port of `FluctuationRationalityRecNW.m` together with `olsWaldtestNW.m`.
#'
#' This implementation uses the survey/model-free calibration in Corollary
#' 9 and Table A.1c of the authors' November 2014 paper and appendix, at 5%
#' or 10%. Forecasts and instruments must satisfy the paper's regularity
#' conditions with negligible forecast-parameter estimation error (F=0).
#' The population-model correction for non-negligible parameter estimation
#' error is not implemented. A HAC covariance of the auxiliary regression
#' alone does not provide that correction.
#' Instruments must be available at the forecast origin. Any nuisance
#' coefficients on x2 are constant under the source specification. The
#' required score functional central limit theorem, exogeneity and F=0
#' conditions are caller assumptions, not properties checked by this code.
#'
#' Complete, equally spaced data are required, with 0.1 <= window/N <= 0.9.
#' Every local regression and tested covariance must have full rank; skipping
#' unidentified windows would invalidate the tabulated supremum calibration.
#' Original Monte Carlo critical values are linearly interpolated, so nominal
#' significance is approximate and asymptotic. The 5% default preserves the
#' threshold used by earlier versions, which had incorrectly labelled it 10%.
#'
#' @param inference "model.free" uses the survey/model-free calibration.
#' "population" is explicitly unsupported and raises an error.
#'
#' @export
fcst.test.rationality <- function(fc.error, x1, x2 = NULL, window,
                                  h = 1, bandwidth = NULL, time.index = NULL,
                                  level = 0.05, inference = c("model.free", "population")) {
  inference <- match.arg(inference)
  if (inference == "population") {
    stop("ERROR! fcst.test.rationality: population-model parameter-estimation correction is not implemented")
  }
  fc.error <- as.matrix(fc.error)
  x1 <- as.matrix(x1)
  if (!is.null(x2)) {
    x2 <- as.matrix(x2)
  }

  P <- nrow(fc.error)
  p <- ncol(x1)
  time.index <- .fcst.time.index(time.index, P)
  if (!is.numeric(fc.error) || !is.numeric(x1) || !p ||
      ncol(fc.error) != 1L || nrow(x1) != P ||
      (!is.null(x2) && (!is.numeric(x2) || nrow(x2) != P))) {
    stop("ERROR! fcst.test.rationality: incompatible regression dimensions")
  }

  x <- if (is.null(x2)) x1 else cbind(x1, x2)
  .fcst.evaluation.window(window, P, time.index, minimum = ncol(x) + 2L)
  if (any(!is.finite(fc.error)) || any(!is.finite(x))) {
    stop("ERROR! fcst.test.rationality: tabulated calibration requires complete finite regression data")
  }
  cv <- get.cv.rationality(p, window / P, level)
  bandwidth <- .fcst.evaluation.bandwidth(window, h, bandwidth)

  path <- rep(NA_real_, P - window + 1)
  coefs <- matrix(NA_real_, P - window + 1, ncol(x))

  for (i in seq_along(path)) {
    idx <- i:(i + window - 1)

    xx <- x[idx, , drop = FALSE]
    yy <- fc.error[idx, , drop = FALSE]

    if (qr(xx)$rank < ncol(xx) || rcond(crossprod(xx)) < .Machine$double.eps) {
      stop(sprintf("ERROR! fcst.test.rationality: rank-deficient regression in window ending at %s", max(idx)))
    }

    .m <- OLS.reg(yy, xx)
    cf <- matrix(.m$coefficients, ncol = 1)
    resid <- as.numeric(.m$residuals)

    vcov <- .fcst.vcov(xx, resid, "HAC", bandwidth, time.index[idx])
    tested <- vcov[seq_len(p), seq_len(p), drop = FALSE]
    root <- tryCatch(chol(tested), error = function(e) NULL)
    if (any(!is.finite(tested)) || is.null(root) ||
        rcond(tested) < .Machine$double.eps) {
      stop(sprintf("ERROR! fcst.test.rationality: singular tested covariance in window ending at %s", max(idx)))
    }
    standardised <- forwardsolve(t(root), cf[seq_len(p), , drop = FALSE])
    path[i] <- sum(standardised^2)
    coefs[i, ] <- drop(cf)
  }

  if (any(!is.finite(path))) {
    stop("ERROR! fcst.test.rationality: nonfinite window statistic")
  }
  i.max <- which.max(path)

  result <- list(
    statistic = path[i.max],
    at = window + i.max - 1,
    path = path,
    index = window:P,
    cv = cv,
    reject = path[i.max] > cv,
    coefficients = coefs,
    window = window,
    level = level,
    mu = window / P,
    h = h,
    bandwidth = bandwidth,
    n.obs = P,
    n.restrictions = p,
    time.index = time.index,
    calibrated = TRUE,
    assumptions.verified = FALSE,
    inference = "model.free: approximate asymptotic calibration with negligible forecast-parameter estimation error",
    test = "Rossi-Sekhposyan fluctuation rationality"
  )

  class(result) <- "bt_fcstFluct"
  result
}


.fcst.evaluation.window <- function(window, n, time.index, minimum) {
  if (!is.numeric(window) || length(window) != 1L || !is.finite(window) ||
      window != floor(window) || window < minimum || window > n) {
    stop(sprintf("ERROR! fcst.evaluation: window must be a finite integer between %s and %s", minimum, n))
  }
  .fcst.require.regular(time.index, "fcst.evaluation")
  invisible(NULL)
}


.fcst.evaluation.bandwidth <- function(n, h, bandwidth) {
  horizons <- .fcst.horizons(h)
  if (length(horizons) != 1L || n <= h) {
    stop("ERROR! fcst.evaluation: a single horizon below the evaluation sample size is required")
  }
  if (is.null(bandwidth)) bandwidth <- max(h, floor(n^(1 / 4)))
  if (!is.numeric(bandwidth) || length(bandwidth) != 1L ||
      !is.finite(bandwidth) || bandwidth != floor(bandwidth) || bandwidth < h) {
    stop("ERROR! fcst.evaluation: bandwidth must include at least h-1 lags")
  }
  bandwidth
}

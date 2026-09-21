#' @title
#' Two-step break inference and forecasting of Altansukh and Osborn
#'
#' @description
#' Detects a coefficient break and a possibly distinct variance break, then
#' repeats coefficient inference after feasible GLS transformation. Uses
#' [fcst.test.Wald] and the package's [OLS.reg].
#'
#' @param y,x,x.new An aligned regression and one forecast row.
#' @param trim Symmetric trimming; .1 is the paper's setting.
#' @param level Nominal asymptotic level of the sup-Wald tests.
#' @param inference "two.step" (FGLS) or the paper's "HC" comparator.
#' @param stepwise Test individual coefficient changes after detecting a break.
#' @param step.level Nominal level for the individual asymptotic two-sided t
#' tests. Remove the least significant change and refit until all remaining
#' changes are significant. These are conditional selection decisions, not
#' multiplicity-adjusted or post-selection p-values.
#' @param hc.type HC covariance for the pilot and HC comparator: HC0 or HC3.
#' @param average Average forecasts over a confidence set for the coefficient
#' break after the final stability test rejects. The variance transformation is
#' held fixed and coefficient selection is repeated at every accepted date.
#' @param conf.method Confidence set method: "eo.morley" or "bai".
#' @param conf.level Confidence probability, distinct from the test size `level`.
#'
#' @details
#' The pilot is a robust coefficient test. Its least-squares date is used only
#' if the test rejects. Fit that mean model and test sqrt(pi/2)*abs(residual)
#' for a constant-mean break using homoskedastic covariance. If detected, the
#' two regime means estimate standard deviations (not variances). Divide y
#' and every x column by these deviations, and repeat the coefficient test
#' with homoskedastic covariance. There is no iteration of the variance step.
#' The forecast uses the final regime's coefficients, with stable coefficients
#' estimated over the full sample when stepwise=TRUE. Without a detected
#' coefficient break it uses full-sample OLS or FGLS as appropriate.
#'
#' Dates index complete regression rows and denote the last pre-break response.
#' With `average=TRUE`, the confidence set is computed on the original data,
#' allowing variance and regressor moments to change across coefficient regimes
#' (Section 2.3). It is not computed on the FGLS-transformed sample. Each date
#' in this set gives a separate OLS/FGLS and optional stepwise fit; forecasts
#' receive equal weights, preserving disconnected sets. The last window start
#' cannot exceed N-trim*N, as in Section 2.1. Failure of a component fit is an
#' error, rather than a silent change to the averaging measure.
#' [fcst.CIavg] remains the separate, unconditional raw-data OLS combination.
#' The source assumes serially uncorrelated innovations; a variance break
#' does not substitute for a serial-correlation model.
#' The estimation sample must be contiguous after trimming leading or trailing
#' missing observations. Internal missing rows and external selection of the
#' regressor set are not covered by this procedure. The published stepwise
#' restrictions on changes of given coefficients remain available.
#'
#' @param time.index Original integer-period positions of observations.
#' Internal gaps in the estimation block are not supported.
#' @return A bt_fcstFGLS with forecast, coefficients, bp, variance.bp,
#' transformation scales, all three tests, and the stepwise trace. Averaging
#' adds the confidence set, accepted dates, component coefficients, forecasts,
#' and stepwise traces. `fit` describes the single estimated-date fit; it is
#' not a regression representation of the averaged forecast.
#' @references Altansukh and Osborn (2022), Empirical Economics 63, 1-41,
#' Sections 2.1-2.3. https://doi.org/10.1007/s00181-021-02137-w.
#' @importFrom stats pnorm
#' @export
fcst.FGLS <- function(y, x = NULL, x.new = NULL, trim = .1, level = .05,
                      inference = c("two.step", "HC"), stepwise = FALSE,
                      step.level = level, hc.type = c("HC0", "HC3"),
                      average = FALSE, conf.method = c("eo.morley", "bai"),
                      conf.level = .95, time.index = NULL) {
  inference <- match.arg(inference)
  hc.type <- match.arg(hc.type)
  conf.method <- match.arg(conf.method)
  .fcst.test.level(level)
  .fcst.test.level(step.level)
  .fcst.test.level(conf.level)
  if (!is.logical(average) || length(average) != 1L || is.na(average)) {
    stop("ERROR! fcst.FGLS: average must be TRUE or FALSE")
  }
  if (!is.logical(stepwise) || length(stepwise) != 1L || is.na(stepwise)) {
    stop("ERROR! fcst.FGLS: stepwise must be TRUE or FALSE")
  }
  time.index <- .fcst.time.index(time.index, NROW(y))
  d <- .fcst.data(y, x, x.new)
  .fcst.inference.data(d)
  time.index <- time.index[d$rows]
  .fcst.require.regular(time.index, "fcst.FGLS")
  pilot <- fcst.test.Wald(d$y, d$x, trim = trim, vcov.type = hc.type, level = level)
  pilot.bp <- if (pilot$reject["sup"]) pilot$bp.LS else 0L
  scale <- rep(1, d$N)
  sigma <- NULL
  variance <- NULL
  variance.bp <- 0L
  final <- pilot
  if (inference == "two.step") {
    mean.fit <- .fcst.fgls.fit(d$y, d$x, pilot.bp, scale, seq_len(d$k), hc.type)
    absolute <- sqrt(pi / 2) * abs(mean.fit$fit$residuals)
    variance <- fcst.test.Wald(absolute, trim = trim, vcov.type = "const", level = level)
    variance.bp <- if (variance$reject["sup"]) variance$bp.LS else 0L
    sigma <- if (variance.bp > 0L) {
      c(pre = mean(absolute[seq_len(variance.bp)]),
        post = mean(absolute[seq.int(variance.bp + 1L, d$N)]))
    } else c(full = mean(absolute))
    if (any(!is.finite(sigma)) || any(sigma <= 0)) {
      stop("ERROR! fcst.FGLS: estimated standard deviations must be positive")
    }
    if (variance.bp > 0L) scale <- rep(sigma, c(variance.bp, d$N - variance.bp))
    final <- fcst.test.Wald(d$y / scale, d$x / scale, trim = trim,
                            vcov.type = "const", level = level)
  }
  bp <- if (final$reject["sup"]) final$bp.LS else 0L
  covariance <- if (inference == "two.step") "const" else hc.type
  fitted <- .fcst.fgls.stepwise(d$y, d$x, bp, scale, covariance, stepwise, step.level)
  combination <- if (average && bp > 0L) {
    .fcst.fgls.average(d, bp, scale, covariance, stepwise, step.level,
                       trim, conf.method, conf.level)
  } else NULL
  coefficients <- if (is.null(combination)) fitted$coefficients else combination$coefficients
  structure(list(forecast = drop(d$x.new %*% coefficients),
                  coefficients = coefficients, bp = bp, variance.bp = variance.bp,
                  break.vars = fitted$break.vars, stepwise = stepwise,
                  step.trace = fitted$trace, scale = scale, sigma = sigma,
                  average = average, averaged = !is.null(combination),
                  confidence.set = combination$confidence.set,
                  break.set = combination$break.set,
                  component.coefficients = combination$component.coefficients,
                  forecasts = combination$forecasts,
                  component.weights = combination$weights,
                  component.break.vars = combination$break.vars,
                  component.traces = combination$traces,
                  conf.method = if (average) conf.method else NULL,
                  conf.level = if (average) conf.level else NULL,
                  weights = 1 / scale^2, tests = list(pilot = pilot, variance = variance,
                                                     final = final),
                  pilot.bp = pilot.bp, fit = fitted$fit, N = d$N, k = d$k,
                  rows = d$rows, time.index = time.index, inference = inference, level = level,
                  step.level = step.level, hc.type = hc.type,
                  inference.note = "nominal asymptotic decisions under the source assumptions; stepwise tests condition on an estimated date"),
            class = "bt_fcstFGLS")
}


.fcst.fgls.average <- function(d, bp, scale, vcov.type, stepwise, step.level,
                                trim, method, level) {
  confidence <- if (method == "eo.morley") {
    fcst.conf.set.EoMorley(d$y, d$x, level = level, trim = trim,
                           variance = "regime", vcov.type = "HC")
  } else {
    fcst.conf.int.Bai(d$y, d$x, bp = bp, level = level, het = TRUE,
                       trim = trim, vcov.type = "HC")
  }
  dates <- if (method == "eo.morley") confidence$set else {
    seq.int(ceiling(confidence$lower), floor(confidence$upper))
  }
  dates <- dates[dates > d$k & dates < d$N - d$k &
                   dates + 1 <= d$N * (1 - trim)]
  if (!length(dates)) {
    stop("ERROR! fcst.FGLS: no confidence-set dates allow an admissible forecast window")
  }
  fits <- lapply(dates, function(date) {
    .fcst.fgls.stepwise(d$y, d$x, date, scale, vcov.type, stepwise, step.level)
  })
  coefficients <- do.call(rbind, lapply(fits, function(fit) fit$coefficients))
  rownames(coefficients) <- dates
  forecasts <- drop(coefficients %*% as.numeric(d$x.new))
  weights <- rep(1 / length(dates), length(dates))
  list(coefficients = colMeans(coefficients), confidence.set = confidence,
       break.set = dates, component.coefficients = coefficients,
       forecasts = forecasts, weights = weights,
       break.vars = lapply(fits, function(fit) fit$break.vars),
       traces = lapply(fits, function(fit) fit$trace))
}


.fcst.fgls.fit <- function(y, x, bp, scale, columns, vcov.type) {
  if (bp == 0L || !length(columns)) {
    design <- x
    map <- seq_len(ncol(x))
    split <- NULL
  } else {
    split <- .fcst.split.design(x, bp, columns)
    design <- split$xs
    map <- split$map.2
  }
  design <- design / scale
  if (nrow(design) <= ncol(design) || qr(design)$rank < ncol(design)) {
    stop("ERROR! fcst.FGLS: unidentified regime design")
  }
  fit <- OLS.reg(y / scale, design)
  list(fit = fit, coefficients = setNames(as.numeric(fit$coefficients[map]), colnames(x)),
       split = split, vcov = .fcst.vcov(design, fit$residuals, vcov.type))
}


.fcst.fgls.stepwise <- function(y, x, bp, scale, vcov.type, stepwise, level) {
  columns <- if (bp > 0L) seq_len(ncol(x)) else integer()
  trace <- data.frame(iteration = integer(), column = integer(), statistic = numeric(),
                       p.value = numeric(), dropped = logical())
  iteration <- 0L
  repeat {
    fitted <- .fcst.fgls.fit(y, x, bp, scale, columns, vcov.type)
    if (!stepwise || !length(columns)) break
    iteration <- iteration + 1L
    split <- fitted$split
    difference <- fitted$fit$coefficients[split$idx.2] - fitted$fit$coefficients[split$idx.1]
    V <- fitted$vcov
    variance <- diag(V)[split$idx.1] + diag(V)[split$idx.2] -
      2 * V[cbind(split$idx.1, split$idx.2)]
    if (any(!is.finite(variance)) || any(variance <= 0)) {
      stop("ERROR! fcst.FGLS: singular covariance in stepwise inference")
    }
    statistic <- difference / sqrt(variance)
    p <- 2 * pnorm(-abs(statistic))
    drop <- which.max(p)
    remove <- p[drop] >= level
    trace <- rbind(trace, data.frame(iteration = iteration, column = columns,
                                     statistic = as.numeric(statistic), p.value = as.numeric(p),
                                     dropped = seq_along(columns) == drop & remove,
                                     row.names = NULL))
    if (!remove) break
    columns <- columns[-drop]
  }
  fitted$break.vars <- columns
  fitted$trace <- trace
  fitted
}

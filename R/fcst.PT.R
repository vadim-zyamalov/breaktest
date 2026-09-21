#' @title
#' Pesaran-Timmermann weighted and pooled window forecasts
#'
#' @description
#' Average expanding-window forecasts, either with weights inverse to past
#' MSFE (WA) or equal weights (pooled). An optional single coefficient-break pretest
#' chooses between the full sample and break-restricted averaging.
#'
#' @param y,x,x.new Regression inputs, as in [fcst.forecast].
#' @param bp Known last pre-break response row on the original grid. Windows
#' start no later than bp+1. NULL with pretest=FALSE averages all candidates.
#' @param trim Trimming for an unknown-date pretest and default minimum window.
#' @param min.window Minimum estimation length is strictly greater than this
#' integer (omega in the paper). Default max(k,round(trim*N)).
#' @param omega First inner-CV origin (T minus omega-tilde in the paper),
#' default round(.75*T). Only WA uses this setting. Each validation target
#' must already have been observed at the outer forecast origin.
#' @param pretest Run [fcst.test.Wald] on the outer training sample. Failure
#' to reject gives the full-sample forecast; rejection restricts starts to
#' at most bp.LS+1. Cannot be combined with a supplied bp.
#' @param level,vcov.type,bandwidth Pretest settings; NULL covariance chooses
#' HC0, or HAC for direct h>1. This does not establish post-selection validity.
#' @param ... Horizon, lags, selection and other common [fcst.forecast] settings.
#'
#' @details
#' WA computes an ordinary recursive forecast error for each start and each
#' common validation origin. Weights are proportional to 1/MSFE. The prose of
#' Pesaran-Timmermann Section 3.4 and Altansukh-Osborn Appendix A specify inverse
#' weights; the printed numerators in PT (22)-(23) omit the inverse. Here the
#' inverse is used. Zero-error candidates share all mass; a candidate failing
#' any usable fold has zero weight and its failure is recorded. The paper's
#' bound start<=T-omega-omega-tilde ensures all validation fits exceed omega.
#'
#' Pooled averaging uses equal weights and the bound start<=T-omega. Its
#' break-restricted variant is the forecast average in PT Section 3.5, not
#' a single regression fitted to pooled observations. Observation weights are
#' descriptive inclusion weights, not equivalent regression weights.
#'
#' Predictors must be specified before fitting; nested OCMT selection is not
#' part of this source procedure. The optional outer pretest conditions the
#' candidate set on information at the forecast origin,
#' as in the paper; inner errors are not an unbiased evaluation of the entire
#' pretest strategy. [fcst.compare] supplies that outer evaluation.
#' For iterated h>1, each window's recursive path is computed separately and
#' the paths are averaged with fixed weights. The returned coefficient vector
#' is a weighted summary and need not generate that averaged recursive path.
#'
#' @return A bt_fcst with window forecasts, weights, validation errors,
#' selection diagnostics and pretest in details.
#' @references Pesaran and Timmermann (2007), Journal of Econometrics 137,
#' 134-161, Sections 3.3-3.6. https://doi.org/10.1016/j.jeconom.2006.03.010.
#' Altansukh and Osborn (2022), Appendix A.
#' @export
fcst.WA <- function(y, x = NULL, x.new = NULL, bp = NULL, trim = .15,
                    min.window = NULL, omega = NULL, pretest = FALSE,
                    level = .05, vcov.type = NULL, bandwidth = NULL, ...) {
  fcst.forecast(y, x, x.new, method = "wa", bp = bp, trim = trim,
                min.window = min.window, omega = omega, pretest = pretest,
                level = level, vcov.type = vcov.type, bandwidth = bandwidth, ...)
}


#' @rdname fcst.WA
#' @export
fcst.Pooled <- function(y, x = NULL, x.new = NULL, bp = NULL, trim = .15,
                        min.window = NULL, pretest = FALSE, level = .05,
                        vcov.type = NULL, bandwidth = NULL, ...) {
  fcst.forecast(y, x, x.new, method = "pooled", bp = bp, trim = trim,
                min.window = min.window, pretest = pretest,
                level = level, vcov.type = vcov.type, bandwidth = bandwidth, ...)
}


.fcst.pt.horizon <- function(history, data, row, h, type, bp, trim, method,
                              x.future, future, min.window = NULL, omega = NULL,
                              pretest = FALSE, level = .05, vcov.type = NULL,
                              bandwidth = NULL) {
  if (!is.null(history$selection)) {
    stop("ERROR! fcst.WA: nested predictor selection is not supported by the source procedure")
  }
  .fcst.require.regular(data$map$target, "fcst.WA")
  if (!is.logical(pretest) || length(pretest) != 1L || is.na(pretest)) {
    stop("ERROR! fcst.WA: pretest must be TRUE or FALSE")
  }
  .fcst.test.level(level)
  if (is.null(vcov.type)) vcov.type <- if (type == "direct" && h > 1L) "HAC" else "HC0"
  if (pretest && !is.null(bp)) stop("ERROR! fcst.WA: supply bp or request a pretest, not both")
  if (is.null(min.window)) min.window <- max(data$k, round(trim * data$N))
  if (!is.numeric(min.window) || length(min.window) != 1L || !is.finite(min.window) ||
      min.window != floor(min.window) || min.window < data$k || min.window >= data$N) {
    stop("ERROR! fcst.WA: min.window must be an integer from k to N-1")
  }
  mapped.bp <- .fcst.break.row(bp, data, history$N)
  test <- NULL
  full.only <- FALSE
  if (pretest) {
    if (type == "direct" && h > 1L && vcov.type == "HAC" && is.null(bandwidth)) {
      bandwidth <- max(h, floor(data$N^(1 / 4)))
    }
    test <- fcst.test.Wald(data$y, data$x, trim = trim, vcov.type = vcov.type,
                            bandwidth = bandwidth, time.index = data$map$target, level = level)
    full.only <- !test$reject["sup"]
    mapped.bp <- if (full.only) 0L else test$bp.LS
  }
  cv <- NULL
  max.start <- data$N - min.window
  if (!is.null(mapped.bp)) max.start <- min(max.start, mapped.bp + 1L)
  if (method == "wa" && !full.only) {
    if (is.null(omega)) omega <- min(round(.75 * history$N), history$N - h)
    cv <- .fcst.cv.folds(history, h, type, omega, x.future, allow.iterated = TRUE)
    for (fold in cv$folds) .fcst.require.regular(fold$data$map$target, "fcst.WA")
    last.starts <- vapply(cv$folds, function(fold) {
      n <- fold$data$N - max(min.window, fold$data$k)
      if (n < 1L) return(0L)
      sum(data$map$target <= fold$data$map$target[n])
    }, 0L)
    max.start <- min(max.start, last.starts)
  }
  if (max.start < 1L) stop("ERROR! fcst.WA: no window is long enough in every validation fold")
  starts <- if (full.only) 1L else seq_len(max.start)
  failures <- data.frame(origin = integer(), start = integer(), reason = character())
  errors <- NULL
  msfe <- NULL
  if (!is.null(cv)) {
    errors <- matrix(NA_real_, length(cv$folds), length(starts),
                       dimnames = list(cv$origins, starts))
    for (j in seq_along(starts)) {
      first <- data$map$target[starts[j]]
      for (i in seq_along(cv$folds)) {
        fold <- cv$folds[[i]]
        idx <- which(fold$data$map$target >= first)
        attempt <- tryCatch({
          coef <- .fcst.pt.coefficients(fold$data$y, fold$data$x, idx)
          error <- fold$actual - .fcst.fold.predict(coef, fold, history, h, type)
          if (!is.finite(error)) stop("ERROR! fcst.WA: nonfinite validation error")
          error
        }, error = function(e) e)
        if (inherits(attempt, "error")) {
          failures <- rbind(failures, data.frame(origin = fold$origin, start = starts[j],
                                                  reason = conditionMessage(attempt)))
        } else errors[i, j] <- attempt
      }
    }
    msfe <- .fcst.cv.score(errors, average = TRUE)
    weights <- .fcst.inverse.msfe(msfe)
  } else weights <- rep(1 / length(starts), length(starts))
  beta <- matrix(0, data$k, length(starts))
  paths <- matrix(0, if (type == "iterated") h else 1L, length(starts))
  obs.weights <- numeric(data$N)
  for (j in seq_along(starts)) {
    if (weights[j] == 0) next
    idx <- seq.int(starts[j], data$N)
    beta[, j] <- .fcst.pt.coefficients(data$y, data$x, idx)
    paths[, j] <- if (type == "iterated" && h > 1L) {
      .fcst.iterate(.fcst.expand.coef(beta[, j], data$columns, history$k),
                     history, history$N, h, future)
    } else drop(row %*% beta[, j])
    obs.weights[idx] <- obs.weights[idx] + weights[j] / length(idx)
  }
  path <- as.numeric(paths %*% weights)
  coefficients <- drop(beta %*% weights)
  details <- list(method = method, starts = starts, start.rows = data$map$target[starts],
                    window.length = data$N - starts + 1L, min.window = min.window,
                    comb.weights = weights, coefficients = coefficients,
                    window.coefficients = beta, window.paths = paths,
                    forecasts = paths[nrow(paths), ], msfe = msfe,
                    obs.weights = obs.weights, pretest = test, full.only = full.only,
                    bp = mapped.bp, cv.origins = cv$origins, cv.targets = cv$targets,
                    cv.errors = errors, cv.excluded = cv$excluded, failures = failures,
                    cv.selection = lapply(cv$folds, function(f) f$data$selection))
  # A failed candidate has zero mass, but zero is not its estimated forecast.
  details$window.coefficients[, weights == 0] <- NA_real_
  details$window.paths[, weights == 0] <- NA_real_
  details$forecasts[weights == 0] <- NA_real_
  structure(list(forecast = tail(path, 1L), coefficients = coefficients,
                   path = if (type == "iterated") path else NULL,
                   method = method, weights = obs.weights, bp = mapped.bp,
                   N = data$N, details = details), class = "bt_fcst")
}


.fcst.pt.coefficients <- function(y, x, idx) {
  xx <- x[idx, , drop = FALSE]
  if (length(idx) <= ncol(x) || qr(xx)$rank < ncol(x)) {
    stop("ERROR! fcst.WA: short or singular window")
  }
  coef <- OLS.reg(y[idx, , drop = FALSE], xx)$coefficients
  if (any(!is.finite(coef))) stop("ERROR! fcst.WA: nonfinite window coefficients")
  as.numeric(coef)
}


.fcst.inverse.msfe <- function(msfe) {
  if (any(msfe == 0)) return(as.numeric(msfe == 0) / sum(msfe == 0))
  weights <- min(msfe) / msfe
  weights / sum(weights)
}

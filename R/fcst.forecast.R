#' @title
#' Direct and iterated forecasts under structural breaks
#'
#' @description
#' A common interface to all estimation-window, weighting and combination
#' strategies. Fits use the existing regression and break-dating routines.
#'
#' @param y A univariate series on a regular observation grid. Keep missing
#' periods as NA rows; lags and horizons count original rows, not complete cases.
#' @param x With lags = NULL, an already aligned design: row t explains `y[t]`
#' and must be available at t-1. Include any constant explicitly. With explicit
#' lags, external predictors observed at their row's date; the generated
#' one-step design uses `x[t-1]`. NULL means no external predictors.
#' @param x.new One forecast-design row available at the current origin, for
#' the already aligned input. Required for h > 1 when x is supplied. The
#' legacy h = 1 default is the last complete design row. Do not supply with lags.
#' @param method One of full, postbreak, window, exps, ppp, ppp.robust, wgls,
#' avew, avelr, ciavg, stein, bp.combine, fgls, wa, pooled.
#' @param window.method A rule in [fcst.window].
#' @param bp Last pre-break response row. With explicit lags or h > 1, indexes
#' the original y; the same applies with predictor selection enabled.
#' It is mapped to the horizon-specific training sample. Legacy h = 1
#' matrix calls retain complete-case indexing. NULL estimates it.
#' FGLS, WA and pooled always use original-row dates.
#' FGLS estimates its own dates, so bp must be NULL for that method.
#' @param trim A break-date trimming fraction in (0, 0.5).
#' @param ... Method-specific arguments. WGLS uses the published post-break
#' LOO criterion at every horizon; cv.method="loo" is its only supported value.
#' ExpS accepts gamma.method="fixed" (default), "cv", "cvl", "ml" or "tradeoff"
#' and the corresponding [fcst.ExpS] options.
#' Supplied CIavg set/interval use the same date convention as bp.
#' @param h Positive integer horizon(s). A scalar returns a single forecast;
#' a vector returns results for each distinct horizon.
#' @param type "direct" fits a separate h-step regression; "iterated" estimates
#' a one-step rule and recursively applies its fixed coefficients.
#' @param lags NULL preserves the supplied-design interface. An explicit
#' vector, e.g. 1:2, constructs `y[t-1]`, `y[t-2]` automatically. integer(0)
#' selects the generated design without AR terms. Required for iterated
#' h > 1 with nonconstant x.
#' @param constant Add a constant in the generated design. Ignored when
#' lags is NULL.
#' @param x.future External x at T+1,...,T+h-1 for iteration. Supply a matrix
#' or a function function(origin, h, y, x) returning those rows. The function
#' receives only histories through origin and is required for inner CV and
#' [fcst.compare] with external predictors. Direct forecasts do not use it.
#' @param index Optional original row labels, defaulting to time(y) for ts
#' and row numbers otherwise. Labels do not insert missing periods.
#' @param selection NULL disables predictor selection. A named list of
#' [fcst.select] options (always, p, delta, vcov.type, bandwidth) enables
#' unweighted single-stage OCMT on each training history. Supported downstream
#' strategies are full, fixed-gamma exps, ppp and ppp.robust. Other combinations
#' and predictor selection inside tuning CV are not implemented. Returned
#' coefficients retain all original columns, with zeros for unselected ones.
#'
#' @details
#' The generated direct regression pairs `y[t+h]` with a constant, `y[t+1-lags]`
#' and `x[t]`, using only pairs with t+h <= T. An already aligned matrix instead
#' pairs design row r with `y[r+h-1]`. No future realised y enters an iteration.
#' Window/weight selection is frozen at the real origin. WA, pooled and FGLS
#' confidence-set averaging use individual recursive paths. The current
#' AveW/AveLR/CIavg and Boot-Pick combinations are direct-only for h > 1.
#'
#' IJR and adaptive window/ExpS CV/CVL support direct h-step regressions,
#' not nonlinear iterated tuning. FGLS supports h=1 and iterated forecasts;
#' its residual-variance procedure is not applied to overlapping direct errors.
#' Unsupported method/horizon/selection combinations fail explicitly.
#'
#' Stein uses its published iid/FGLS coefficient estimate, never a HAC-modified
#' shrinkage rule. Direct Eo-Morley/Bai and Boot-Pick use HAC by default at h>1.
#' Source assumptions remain necessary for inference; no blanket risk dominance
#' is claimed. Simple recursion of a fixed estimate is a plug-in baseline.
#' Statistical methods require contiguous complete training observations.
#' Full OLS, fixed ExpS and a supplied post-break window retain descriptive
#' complete-case forecasts. These are point forecasts, not prediction intervals.
#'
#' @return A bt_fcst with forecast, coefficients, weights, bp, details, h,
#' type, origin and target. Generated/multi-step calls also return training
#' (original origin/target row map), bp.row (internal complete-case break
#' index), and the recursive path for iteration.
#' With multiple horizons a bt_fcstMulti contains forecasts and by.horizon,
#' retaining all individual results.
#'
#' @references
#' Marcellino, M., J. H. Stock and M. W. Watson (2006).
#' "A Comparison of Direct and Iterated Multistep AR Methods for Forecasting
#' Macroeconomic Time Series." Journal of Econometrics 135, 499-526.
#' https://doi.org/10.1016/j.jeconom.2005.07.020.
#'
#' @importFrom stats time is.ts
#' @importFrom utils tail
#'
#' @export
fcst.forecast <- function(
  y,
  x = NULL,
  x.new = NULL,
  method = c("full", "postbreak", "window", "exps", "ppp", "ppp.robust",
             "avew", "avelr", "ciavg", "stein", "bp.combine", "wgls",
             "fgls", "wa", "pooled"),
  window.method = "cv",
  bp = NULL,
  trim = 0.15,
  ...,
  h = 1,
  type = c("direct", "iterated"),
  lags = NULL,
  constant = TRUE,
  x.future = NULL,
  index = NULL,
  selection = NULL
) {
  method <- match.arg(method)
  type <- match.arg(type)
  horizons <- .fcst.horizons(h)
  if (!is.numeric(trim) || length(trim) != 1 || !is.finite(trim) ||
      trim <= 0 || trim >= 0.5) stop("ERROR! fcst.forecast: trim must be in (0, 0.5)")
  options <- list(...)
  .fcst.check.options(options, .fcst.method.options(method, window.method))
  if (method == "exps" && !is.null(options$gamma.method)) {
    options$gamma.method <- match.arg(options$gamma.method,
                                      c("fixed", "cv", "cvl", "ml", "tradeoff"))
  }
  selection <- .fcst.selection.options(selection)
  for (hh in horizons) .fcst.source.options(method, window.method, hh, type, options, selection)
  history <- .fcst.history(y, x, lags, constant, index, selection)
  by.horizon <- lapply(horizons, function(hh) {
    legacy <- .fcst.legacy(horizons, lags, method, options, selection)
    if (legacy) {
      args <- options
      args$cv.method <- NULL
      f <- do.call(.fcst.forecast.one,
                   c(list(y = y, x = x, x.new = x.new, method = method,
                          window.method = window.method, bp = bp, trim = trim), args))
      f$h <- 1L
      f$type <- type
      f$origin <- history$N
      f$target <- history$N + 1L
      f$path <- if (type == "iterated") f$forecast else NULL
    } else {
      f <- .fcst.forecast.h(history, hh, type, x.new, x.future, method,
                             window.method, bp, trim, options)
    }
    if (length(f$forecast) != 1L || !is.finite(f$forecast)) {
      stop("ERROR! fcst.forecast: the method did not return a finite scalar forecast")
    }
    if (!is.null(selection) || .fcst.adaptive(method, options, bp)) {
      f <- .fcst.suppress.inference(f)
    }
    f
  })
  names(by.horizon) <- as.character(horizons)
  if (length(horizons) == 1L) return(by.horizon[[1]])
  structure(list(h = horizons, type = type, method = method,
                 forecasts = vapply(by.horizon, function(f) f$forecast, 0.0),
                 by.horizon = by.horizon), class = "bt_fcstMulti")
}


# Check supported algorithms before fitting, including fcst.compare's direct
# calls to the temporal dispatcher. A warning cannot validate a new estimator.
.fcst.source.options <- function(method, window.method, h, type, options, selection) {
  if (!is.null(options$cv.method) &&
      (method != "wgls" || !identical(options$cv.method, "loo"))) {
    stop("ERROR! fcst.forecast: WGLS supports only cv.method = 'loo'")
  }
  adaptive.exps <- method == "exps" && !is.null(options$gamma.method) &&
    options$gamma.method != "fixed"
  if (!is.null(selection) &&
      (!method %in% c("full", "exps", "ppp", "ppp.robust") || adaptive.exps)) {
    stop("ERROR! fcst.forecast: selection is supported only with full, fixed ExpS, PPP and robust PPP")
  }
  if (h > 1L && type == "direct" && method == "fgls") {
    stop("ERROR! fcst.forecast: FGLS supports h = 1 or iterated forecasts, not overlapping direct errors")
  }
  if (h > 1L && type == "iterated") {
    if (method %in% c("avew", "avelr", "ciavg", "bp.combine")) {
      stop("ERROR! fcst.forecast: this combination supports only direct forecasts for h > 1")
    }
    if ((method == "window" && window.method %in% c("cv", "cv.pre", "cvl", "cvl.pre", "ijr")) ||
        (adaptive.exps && options$gamma.method %in% c("cv", "cvl"))) {
      stop("ERROR! fcst.forecast: adaptive window/ExpS tuning requires direct forecasts for h > 1")
    }
  }
  invisible(NULL)
}


.fcst.gaps.allowed <- function(method, options, bp) {
  method == "full" || (method == "postbreak" && !is.null(bp)) ||
    (method == "exps" && (is.null(options$gamma.method) || options$gamma.method == "fixed"))
}


.fcst.adaptive <- function(method, options, bp) {
  !method %in% c("full", "ppp.robust", "avew") &&
    !(method == "postbreak" && !is.null(bp)) &&
    !(method == "exps" && (is.null(options$gamma.method) || options$gamma.method == "fixed"))
}


# Only unchanged fixed-model, one-step calls use complete-case break dates.
.fcst.legacy <- function(horizons, lags, method, options, selection) {
  length(horizons) == 1L && horizons == 1L && is.null(lags) && is.null(selection) &&
    !method %in% c("fgls", "wa", "pooled") &&
    (is.null(options$cv.method) || options$cv.method == "loo") &&
    !(method == "exps" && !is.null(options$gamma.method) && options$gamma.method != "fixed")
}


# Explicit routing prevents misspelled settings from being silently ignored.
.fcst.method.options <- function(method, window.method = "cv") {
  switch(method,
    full = character(),
    postbreak = c("max.breaks", "criterion"),
    window = switch(window.method,
      ls = character(), bp = c("max.breaks", "criterion"),
      tradeoff = "n.grid", cv = "omega", cv.pre = "omega",
      cvl = "omega", cvl.pre = "omega",
      ijr = c("R0", "R.range", "c.max", "cv.pretest", "level", "vcov.type",
              "bandwidth", "omega", "pretest"),
      stop("ERROR! fcst.forecast: unknown window.method")),
    exps = c("gamma", "gamma.method", "gamma.grid", "eta.grid", "omega", "prior", "sigma.sq"),
    ppp = character(), ppp.robust = c("b.lo", "b.hi"),
    wgls = c("gamma.grid", "gamma.star.grid", "cv.method"),
    avew = c("w.min", "w.max", "m", "weights"),
    avelr = c("c", "break.vars"),
    ciavg = c("interval", "level", "conf.method", "set", "variance",
              "vcov.type", "bandwidth", "kernel", "het"),
    stein = c("tau", "gls", "vcov.type", "bandwidth", "kernel"),
    bp.combine = c("level", "target", "break.vars", "vcov.type", "cv.table", "bandwidth"),
    fgls = c("level", "inference", "stepwise", "step.level", "hc.type",
              "average", "conf.method", "conf.level"),
    wa = c("min.window", "omega", "pretest", "level", "vcov.type", "bandwidth"),
    pooled = c("min.window", "pretest", "level", "vcov.type", "bandwidth"),
    stop("ERROR! fcst.forecast: unknown method"))
}


.fcst.check.options <- function(options, allowed) {
  nms <- names(options)
  if (length(options) && (is.null(nms) || any(nms == "") || anyDuplicated(nms))) {
    stop("ERROR! fcst.forecast: method options must have distinct names")
  }
  unknown <- setdiff(nms, allowed)
  if (length(unknown)) {
    stop("ERROR! fcst.forecast: unsupported option(s): ", paste(unknown, collapse = ", "))
  }
}


# Fit a strategy to a regression that has already been aligned.
.fcst.forecast.one <- function(
  y,
  x = NULL,
  x.new = NULL,
  method = c(
    "full", "postbreak", "window", "exps", "ppp", "ppp.robust",
    "avew", "avelr", "ciavg", "stein", "bp.combine", "wgls", "fgls"
  ),
  window.method = "cv",
  bp = NULL,
  trim = 0.15,
  time.index = NULL,
  ...
) {
  method <- match.arg(method)

  time.index <- .fcst.time.index(time.index, NROW(y))
  # The forecast date belongs to the original clock, including a missing tail.
  forecast.time <- max(time.index) + 1L
  .d <- .fcst.data(y, x, x.new)
  time.index <- time.index[.d$rows]
  if (!.fcst.gaps.allowed(method, list(...), bp)) {
    .fcst.require.regular(time.index, "fcst.forecast")
  }
  y <- .d$y
  x <- .d$x
  x.new <- .d$x.new
  N <- .d$N

  w <- NULL
  details <- NULL

  fc <- switch(
    method,
    full = {
      w <- weights.equal(N)
      details <- WLS.reg(y, x, w)
      .fcst.predict(details, x.new)
    },
    postbreak = {
      if (is.null(bp)) {
        bp <- fcst.window.BP(y, x, trim, time.index = time.index, ...)$bp
      }
      w <- weights.rolling(N, N - bp)
      details <- WLS.reg(y, x, w)
      .fcst.predict(details, x.new)
    },
    window = {
      args <- list(y = y, x = x, method = window.method, trim = trim,
                   time.index = time.index, ...)
      if (window.method == "ijr") {
        args$x.new <- x.new
        args$time.index <- time.index
        if (is.null(args$forecast.time)) args$forecast.time <- forecast.time
      }
      if (!is.null(bp) && window.method %in%
          c("tradeoff", "cv", "cv.pre", "cvl", "cvl.pre")) args$bp <- bp
      details <- do.call(fcst.window, args)
      w <- weights.rolling(N, details$R)
      details$fit <- WLS.reg(y, x, w)
      bp <- details$details$bp
      .fcst.predict(details$fit, x.new)
    },
    exps = {
      args <- list(...)
      if (!is.null(args$gamma.method) && args$gamma.method != "fixed") {
        stop("ERROR! fcst.forecast: adaptive ExpS requires the temporal dispatcher")
      }
      args$gamma.method <- NULL
      .fcst.check.options(args, c("gamma", "time.index"))
      args$time.index <- time.index
      w <- do.call(weights.exponential, c(list(N = N), args))
      details <- WLS.reg(y, x, w)
      .fcst.predict(details, x.new)
    },
    ppp = {
      if (is.null(bp)) {
        bp <- fcst.window.LS(y, x, trim, time.index = time.index)$bp
      }
      .bs <- fcst.break.size(y, x, bp, x.new, time.index = time.index)
      w <- weights.PPP.optimal(N, bp, .bs$phi, if (length(bp) == 1L) .bs$q else 1)
      details <- list(fit = WLS.reg(y, x, w), break.size = .bs)
      .fcst.predict(details$fit, x.new)
    },
    ppp.robust = {
      w <- weights.PPP.robust(N, ...)
      details <- WLS.reg(y, x, w)
      .fcst.predict(details, x.new)
    },
    wgls = {
      details <- fcst.WGLS.cv(y, x, x.new, bp = bp, trim = trim, time.index = time.index, ...)
      w <- details$weights
      bp <- details$bp
      details$forecast
    },
    avew = {
      details <- fcst.AveW(y, x, x.new, ...)
      w <- details$obs.weights
      details$forecast
    },
    avelr = {
      details <- fcst.AveLR(y, x, x.new, trim = trim, time.index = time.index, ...)
      bp <- details$bp
      details$forecast
    },
    ciavg = {
      details <- fcst.CIavg(y, x, x.new, trim = trim, time.index = time.index, ...)
      bp <- details$confidence.set$bp
      details$forecast
    },
    stein = {
      details <- fcst.stein(y, x, x.new, bp = bp, trim = trim, time.index = time.index, ...)
      bp <- details$bp
      details$forecast
    },
    bp.combine = {
      details <- fcst.test.BP(y, x, x.new, trim = trim, time.index = time.index, ...)
      bp <- details$bp
      details$forecast["combined"]
    },
    fgls = {
      if (!is.null(bp)) stop("ERROR! fcst.FGLS: the two-step procedure estimates its own dates")
      details <- fcst.FGLS(y, x, x.new, trim = trim, time.index = time.index, ...)
      bp <- details$bp
      details$forecast
    }
  )

  result <- list(
    forecast = drop(unname(fc)),
    coefficients = if (!is.null(details$coefficients)) details$coefficients
      else details$fit$coefficients,
    method = method,
    weights = w,
    bp = bp,
    N = N,
    details = details
  )

  class(result) <- "bt_fcst"
  result
}

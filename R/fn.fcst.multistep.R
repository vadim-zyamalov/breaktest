# Shared time alignment and recursion for fcst.forecast and fcst.compare.
# Rows are periods on the original grid; incomplete rows are removed only
# after the response has been shifted and the lagged design has been built.

.fcst.horizons <- function(h) {
  if (!is.numeric(h) || !is.null(dim(h)) || !length(h) ||
      any(!is.finite(h)) || any(h < 1) || any(h != floor(h))) {
    stop("ERROR! fcst.forecast: h must contain positive integer horizons")
  }
  sort(unique(as.integer(h)))
}


.fcst.history <- function(y, x, lags, constant, index = NULL, selection = NULL) {
  if (is.null(index)) index <- if (is.ts(y)) as.numeric(time(y)) else seq_len(NROW(y))
  y <- as.matrix(y)
  if (!is.numeric(y) || ncol(y) != 1 || any(is.infinite(y))) {
    stop("ERROR! fcst.forecast: y must be a numeric univariate series")
  }
  N <- nrow(y)
  if (length(index) != N || anyNA(index) || anyDuplicated(index)) {
    stop("ERROR! fcst.forecast: index must uniquely label every original row")
  }
  if (!is.null(x)) {
    x <- as.matrix(x)
    if (!is.numeric(x) || nrow(x) != N || !ncol(x) || any(is.infinite(x))) {
      stop("ERROR! fcst.forecast: x must be a numeric matrix with one row per period")
    }
  }
  generated <- !is.null(lags)
  if (generated) {
    if (!is.numeric(lags) || !is.null(dim(lags)) ||
        any(!is.finite(lags)) || any(lags < 1) || any(lags != floor(lags)) ||
        anyDuplicated(lags) || any(lags >= N)) {
      stop("ERROR! fcst.forecast: lags must be distinct positive integers below length(y)")
    }
    if (!is.logical(constant) || length(constant) != 1 || is.na(constant)) {
      stop("ERROR! fcst.forecast: constant must be TRUE or FALSE")
    }
    lags <- as.integer(lags)
    z <- matrix(numeric(), N, 0)
    if (constant) z <- cbind(z, const = .const(N))
    for (lag in lags) {
      z <- cbind(z, .lagn(y, lag))
      colnames(z)[ncol(z)] <- paste0("y.l", lag)
    }
    if (!is.null(x)) z <- cbind(z, .lagn(x, 1))
    if (!ncol(z)) stop("ERROR! fcst.forecast: the design has no regressors")
  } else {
    z <- if (is.null(x)) matrix(.const(N), N, 1, dimnames = list(NULL, "const")) else x
  }
  labels <- colnames(z)
  if (is.null(labels)) labels <- rep("", ncol(z))
  unnamed <- is.na(labels) | labels == ""
  labels[unnamed] <- paste0("x", which(unnamed))
  colnames(z) <- labels
  list(y = y, x = x, z = z, N = N, k = ncol(z), index = index,
       lags = lags, constant = constant, generated = generated, selection = selection)
}


.fcst.align <- function(history, origin, h) {
  last <- origin - h + 1L
  if (last < 1) stop("ERROR! fcst.forecast: horizon leaves no training pairs")
  predictor <- seq_len(last)
  target <- predictor + h - 1L
  y <- history$y[target, , drop = FALSE]
  x <- history$z[predictor, , drop = FALSE]
  keep <- is.finite(y[, 1]) & rowSums(!is.finite(x)) == 0
  minimum <- if (is.null(history$selection)) ncol(x) + 1L else 2L
  if (sum(keep) <= minimum) {
    stop("ERROR! fcst.forecast: not enough complete training pairs")
  }
  list(y = y[keep, , drop = FALSE], x = x[keep, , drop = FALSE],
       N = sum(keep), k = ncol(x),
       map = data.frame(row = seq_len(sum(keep)),
                        origin = predictor[keep] - 1L,
                        target = target[keep]))
}


# x is observed at the forecast origin: the one-step regression uses x[t-1].
# Only lagged y columns are updated with predictions during recursion.
.fcst.regressor <- function(history, values, origin, step = 1L, future = NULL) {
  row <- numeric()
  if (history$constant) row <- c(row, 1)
  if (length(history$lags)) {
    positions <- origin + step - history$lags
    if (any(positions < 1)) stop("ERROR! fcst.forecast: insufficient lag history")
    row <- c(row, values[positions])
  }
  if (!is.null(history$x)) {
    row <- c(row, if (step == 1L) history$x[origin, ] else future[step - 1L, ])
  }
  matrix(row, nrow = 1)
}


.fcst.new <- function(history, origin, x.new = NULL) {
  if (history$generated) {
    if (!is.null(x.new)) {
      stop("ERROR! fcst.forecast: x.new is built automatically when lags is supplied")
    }
    row <- .fcst.regressor(history, history$y[, 1], origin)
  } else if (!is.null(x.new)) {
    row <- if (is.matrix(x.new)) x.new else matrix(x.new, nrow = 1)
  } else if (is.null(history$x)) {
    row <- matrix(1, 1, 1)
  } else if (origin < history$N) {
    # Aligned-design contract: x[s+1] is available at origin s.
    row <- history$z[origin + 1L, , drop = FALSE]
  } else {
    stop("ERROR! fcst.forecast: supply x.new for a multi-step aligned design")
  }
  if (!is.numeric(row) || nrow(row) != 1 || ncol(row) != history$k ||
      any(!is.finite(row))) {
    stop("ERROR! fcst.forecast: forecast regressors must be one finite conformable row")
  }
  row
}


.fcst.future <- function(history, origin, h, x.future) {
  if (h == 1L || is.null(history$x)) return(NULL)
  if (is.null(x.future)) {
    stop("ERROR! fcst.forecast: iterated forecasts with external x require x.future")
  }
  future <- if (is.function(x.future)) {
    x.future(origin = origin, h = h,
             y = history$y[seq_len(origin), , drop = FALSE],
             x = history$x[seq_len(origin), , drop = FALSE])
  } else {
    if (origin != history$N) {
      stop("ERROR! fcst.forecast: inner CV requires an x.future function for each origin")
    }
    x.future
  }
  future <- as.matrix(future)
  if (!is.numeric(future) || nrow(future) < h - 1L ||
      ncol(future) != ncol(history$x) ||
      any(!is.finite(future[seq_len(h - 1L), , drop = FALSE]))) {
    stop("ERROR! fcst.forecast: x.future must supply h-1 finite rows of external regressors")
  }
  future[seq_len(h - 1L), , drop = FALSE]
}


.fcst.iterate <- function(coefficients, history, origin, h, future) {
  if (!history$generated) {
    # The only aligned design with automatic recursion is a constant.
    if (!is.null(history$x)) {
      stop("ERROR! fcst.forecast: automatic recursion requires explicit lags")
    }
    return(rep(as.numeric(coefficients), h))
  }
  values <- c(history$y[seq_len(origin), 1], rep(NA_real_, h))
  for (step in seq_len(h)) {
    row <- .fcst.regressor(history, values, origin, step, future)
    if (any(!is.finite(row))) stop("ERROR! fcst.forecast: missing lag in forecast recursion")
    values[origin + step] <- drop(row %*% coefficients)
    if (!is.finite(values[origin + step])) {
      stop("ERROR! fcst.forecast: forecast recursion is not finite")
    }
  }
  values[origin + seq_len(h)]
}


# A supplied date always refers to the last pre-break RESPONSE on the
# original grid. An estimated direct-regression break is horizon-specific.
.fcst.break.row <- function(bp, data, origin) {
  if (is.null(bp)) return(NULL)
  if (!is.numeric(bp) || length(bp) != 1 || !is.finite(bp) ||
      bp != floor(bp) || bp < 0 || bp >= origin) {
    stop("ERROR! fcst.forecast: bp must be a last pre-break row below the forecast origin")
  }
  sum(data$map$target <= bp)
}


.fcst.map.options <- function(options, data, origin) {
  if (!is.null(options$breaks)) {
    options$breaks <- vapply(options$breaks, .fcst.break.row, 0L, data = data, origin = origin)
  }
  if (!is.null(options$set) || !is.null(options$interval)) {
    if (!is.null(options$set) && !is.null(options$interval)) {
      stop("ERROR! fcst.forecast: supply set or interval, not both")
    }
    dates <- options$set
    if (!is.null(options$interval)) {
      bounds <- options$interval
      if (!is.numeric(bounds) || length(bounds) != 2 || any(!is.finite(bounds)) ||
          any(bounds != floor(bounds)) || bounds[1] > bounds[2]) {
        stop("ERROR! fcst.forecast: invalid interval on the original row grid")
      }
      dates <- seq.int(bounds[1], bounds[2])
    }
    if (!is.numeric(dates) || any(!is.finite(dates)) ||
        any(dates != floor(dates)) || any(dates < 1) || any(dates >= origin)) {
      stop("ERROR! fcst.forecast: invalid set on the original row grid")
    }
    mapped <- vapply(unique(dates), .fcst.break.row, 0L, data = data, origin = origin)
    counts <- table(mapped[mapped > 0 & mapped < data$N])
    options$set <- as.integer(names(counts))
    attr(options$set, "multiplicity") <- setNames(as.numeric(counts), names(counts))
    options$interval <- NULL
  }
  options
}


.fcst.forecast.h <- function(history, h, type, x.new, x.future,
                             method, window.method, bp, trim, options) {
  if (type == "iterated" && h > 1L && !history$generated && !is.null(history$x)) {
    stop("ERROR! fcst.forecast: iterated h > 1 requires explicit lags (integer(0) for no AR terms)")
  }
  .fcst.source.options(method, window.method, h, type, options, history$selection)
  origin <- history$N
  fit.h <- if (type == "direct") h else 1L
  data <- .fcst.align(history, origin, fit.h)
  if (!is.null(history$selection) || !.fcst.gaps.allowed(method, options, bp)) {
    .fcst.require.regular(data$map$target, "fcst.forecast")
  }
  if (h == 1L && !history$generated && is.null(x.new)) x.new <- tail(data$x, 1L)
  row <- .fcst.new(history, origin, x.new)
  data <- .fcst.select.data(data, history$selection, h, type)
  row <- row[, data$columns, drop = FALSE]
  future <- if (type == "iterated") .fcst.future(history, origin, h, x.future) else NULL
  mapped.bp <- if (method == "ppp" && length(bp) > 1L) {
    vapply(bp, .fcst.break.row, 0L, data = data, origin = origin)
  } else .fcst.break.row(bp, data, origin)
  args <- .fcst.map.options(options, data, origin)
  if (!is.null(data$selection) && !is.null(args$break.vars)) {
    cols <- args$break.vars
    if (!is.numeric(cols) || !length(cols) || any(!is.finite(cols)) ||
        any(cols != floor(cols) | cols < 1 | cols > history$k) || anyDuplicated(cols)) {
      stop("ERROR! fcst.forecast: break.vars must index distinct original design columns")
    }
    args$break.vars <- match(intersect(cols, data$columns), data$columns)
    if (!length(args$break.vars)) stop("ERROR! fcst.forecast: selection retained no break.vars")
  }
  args$cv.method <- NULL
  predict.coef <- function(coef, at = origin, xn = row, xf = future) {
    if (type == "direct" || h == 1L) drop(xn %*% coef)
    else .fcst.iterate(.fcst.expand.coef(coef, data$columns, history$k), history, at, h, xf)[h]
  }
  if (method %in% c("wa", "pooled")) {
    if (type == "direct" && h > 1L && is.null(args$vcov.type)) args$vcov.type <- "HAC"
    result <- do.call(.fcst.pt.horizon,
                      c(list(history = history, data = data, row = row, h = h, type = type,
                             bp = bp, trim = trim, method = method, x.future = x.future,
                             future = future), args))
  } else if (method == "window" && window.method %in% c("cv", "cv.pre", "cvl", "cvl.pre")) {
    result <- .fcst.window.horizon(history, data, h, type, bp, trim,
                                    window.method, x.future, args)
  } else if (method == "exps" && !is.null(args$gamma.method) && args$gamma.method != "fixed") {
    result <- .fcst.exps.horizon(history, data, h, type, x.future, args, row)
  } else {
    # Direct h-step disturbances can be serially dependent. These covariance
    # defaults do not extend any risk-dominance or confidence-coverage theorem.
    if (type == "direct" && h > 1L && method %in% c("bp.combine", "ciavg")) {
      if (method == "bp.combine" || (is.null(args$set) && is.null(args$interval))) {
        if (is.null(args$vcov.type)) args$vcov.type <- "HAC"
        if (args$vcov.type == "HAC") {
          if (is.null(args$bandwidth)) args$bandwidth <- max(h, floor(data$N^(1 / 4)))
          args$time.index <- data$map$target
        }
      }
    }
    if (method %in% c("ciavg", "bp.combine")) args$time.index <- data$map$target
    if (method == "exps") args$time.index <- data$map$target
    if (method == "window" && window.method == "ijr") {
      args$time.index <- data$map$target
      args$forecast.time <- origin + fit.h
      if (is.null(args$bandwidth)) args$bandwidth <- max(fit.h, floor(data$N^(1 / 4)))
      if (is.null(args$R0)) {
        pilot.omega <- if (is.null(args$omega)) max(history$k + 1L, round(.75 * origin)) else args$omega
        args$pilot.fun <- function() {
          pilot <- .fcst.window.horizon(history, data, h, type, NULL, .1,
                                         "cv", x.future, list(omega = pilot.omega))
          pilot$details$details
        }
      }
    }
    result <- do.call(.fcst.forecast.one,
                     c(list(y = data$y, x = data$x, x.new = row, method = method,
                            window.method = window.method, bp = mapped.bp, trim = trim), args))
  }

  if (!method %in% c("wa", "pooled")) {
    result$forecast <- predict.coef(result$coefficients)
    result$path <- if (type == "iterated") {
      if (h == 1L) result$forecast else .fcst.iterate(
        .fcst.expand.coef(result$coefficients, data$columns, history$k), history, origin, h, future)
    } else NULL
  }
  if (method == "fgls" && isTRUE(result$details$averaged) && type == "iterated" && h > 1L) {
    coefficients <- result$details$component.coefficients
    paths <- vapply(seq_len(nrow(coefficients)), function(j) {
      .fcst.iterate(.fcst.expand.coef(coefficients[j, ], data$columns, history$k),
                     history, origin, h, future)
    }, numeric(h))
    paths <- matrix(paths, nrow = h)
    result$path <- drop(paths %*% result$details$component.weights)
    result$forecast <- result$path[h]
    result$details$component.paths <- paths
  }
  if (!is.null(data$selection)) {
    result$selection <- data$selection
    result$coefficients <- setNames(.fcst.expand.coef(result$coefficients, data$columns, history$k),
                                     colnames(history$z))
  }
  result$h <- h
  result$type <- type
  result$details$forecast.horizon <- fit.h
  result$details$selection.horizon <- h
  result$origin <- origin
  result$target <- origin + h
  result$origin.index <- history$index[origin]
  result$training <- data$map
  result$training$origin.index <- history$index[pmax(1L, data$map$origin)]
  result$training$origin.index[data$map$origin == 0] <- NA
  result$training$target.index <- history$index[data$map$target]
  result$bp.row <- result$bp
  if (!is.null(result$bp)) {
    result$bp <- if (!is.null(bp)) tail(bp, 1L) else if (result$bp == 0) 0L else data$map$target[result$bp]
    result$bp.index <- if (result$bp == 0) NA else history$index[result$bp]
  }
  if (method == "ciavg") {
    result$break.set <- data$map$target[result$details$break.set]
  }
  if (method == "fgls") {
    result$break.set <- data$map$target[result$details$break.set]
    result$variance.bp <- if (result$details$variance.bp == 0L) 0L else data$map$target[result$details$variance.bp]
  }
  result$forecast.scope <- if (type == "direct") {
    "linear direct regression; inference is conditional on the source-paper assumptions"
  } else {
    "fixed-estimate recursion; no general h-step risk dominance is asserted"
  }
  if (method == "stein" && h > 1L) {
    result$details$risk.eligible <- FALSE
    result$details$risk.conditions <- "the source risk result is not asserted for this h-step forecast"
  }
  if (!is.null(data$selection) || .fcst.adaptive(method, options, bp)) {
    result <- .fcst.suppress.inference(result)
  }
  result
}

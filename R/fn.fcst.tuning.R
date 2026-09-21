# Horizon-aware, rolling-origin validation. Every candidate is scored on the
# same folds, and every fold trains only on responses already observed there.

.fcst.cv.folds <- function(history, h, type, omega, x.future, allow.iterated = FALSE) {
  if (type == "iterated" && h > 1L && !is.null(history$x) && !is.function(x.future)) {
    stop("ERROR! fcst.forecast: inner CV with external x requires an x.future function")
  }
  if (!is.null(history$selection)) {
    stop("ERROR! fcst.forecast: predictor selection inside tuning CV is unsupported")
  }
  if (type == "iterated" && h > 1L && !allow.iterated) {
    stop("ERROR! fcst.forecast: adaptive CV requires a direct regression for h > 1")
  }
  last <- history$N - h
  if (is.null(omega)) omega <- min(round(0.9 * history$N), last)
  if (!is.numeric(omega) || length(omega) != 1 || !is.finite(omega) ||
      omega != floor(omega) || omega < 1 || omega > last) {
    stop("ERROR! fcst.forecast: omega must leave observed h-step validation targets")
  }
  folds <- list()
  excluded <- data.frame(origin = integer(), reason = character())
  for (s in seq.int(omega, last)) {
    fold <- tryCatch({
      actual <- history$y[s + h, 1]
      if (!is.finite(actual)) stop("missing validation response")
      data <- .fcst.align(history, s, if (type == "direct") h else 1L)
      .fcst.require.regular(data$map$target, "fcst.forecast CV")
      row <- .fcst.new(history, s)
      data$columns <- seq_len(data$k)
      future <- if (type == "iterated") .fcst.future(history, s, h, x.future) else NULL
      list(data = data, row = row, future = future, actual = actual, origin = s,
           target = s + h)
    }, error = function(e) e)
    if (inherits(fold, "error")) {
      excluded <- rbind(excluded, data.frame(origin = s, reason = conditionMessage(fold)))
    } else {
      folds[[length(folds) + 1L]] <- fold
    }
  }
  if (nrow(excluded)) {
    stop("ERROR! fcst.forecast: every scheduled CV fold must be usable: ", excluded$reason[1])
  }
  list(folds = folds, excluded = excluded,
       origins = vapply(folds, function(f) f$origin, 0L),
       targets = vapply(folds, function(f) f$target, 0L))
}


.fcst.fold.predict <- function(coef, fold, history, h, type) {
  if (type == "direct" || h == 1L) drop(fold$row %*% coef)
  else .fcst.iterate(.fcst.expand.coef(coef, fold$data$columns, history$k),
                      history, fold$origin, h, fold$future)[h]
}


.fcst.cv.score <- function(errors, average = FALSE) {
  valid <- colSums(!is.finite(errors)) == 0L
  criterion <- rep(Inf, ncol(errors))
  criterion[valid] <- colSums(errors[, valid, drop = FALSE]^2)
  if (average) criterion <- criterion / nrow(errors)
  if (!any(is.finite(criterion))) {
    stop("ERROR! fcst.forecast: no candidate is valid on every inner-CV fold")
  }
  criterion
}


.fcst.window.horizon <- function(history, data, h, type, bp, trim,
                                  method, x.future, options) {
  if (length(setdiff(names(options), "omega"))) {
    stop("ERROR! fcst.forecast: unsupported horizon-CV window option")
  }
  cv <- .fcst.cv.folds(history, h, type, options$omega, x.future)
  max.start <- max(1L, round((1 - trim) * data$N))
  mapped.bp <- .fcst.break.row(bp, data, history$N)
  if (method %in% c("cv.pre", "cvl.pre")) {
    if (is.null(mapped.bp)) mapped.bp <- fcst.window.LS(data$y, data$x, trim)$bp
    max.start <- min(max.start, mapped.bp + 1L)
  }
  errors <- matrix(NA_real_, length(cv$folds), max.start)
  for (j in seq_len(max.start)) {
    first.target <- data$map$target[j]
    for (i in seq_along(cv$folds)) {
      fold <- cv$folds[[i]]
      idx <- which(fold$data$map$target >= first.target)
      errors[i, j] <- tryCatch({
        if (length(idx) <= fold$data$k) stop("short window")
        x <- fold$data$x[idx, , drop = FALSE]
        if (qr(x)$rank < fold$data$k) stop("singular window")
        coef <- OLS.reg(fold$data$y[idx, , drop = FALSE], x)$coefficients
        fold$actual - .fcst.fold.predict(coef, fold, history, h, type)
      }, error = function(e) NA_real_)
    }
  }
  criterion <- .fcst.cv.score(errors)
  posterior <- NULL
  if (method %in% c("cvl", "cvl.pre")) {
    s.sq <- mean(OLS.reg(data$y, data$x)$residuals^2)
    if (!is.finite(s.sq) || s.sq < 0) {
      stop("ERROR! fcst.forecast: Laplace CV requires a finite nonnegative horizon-error scale")
    }
    delta <- criterion - min(criterion)
    posterior <- if (s.sq == 0) as.numeric(delta == 0) else exp(-0.5 * delta / s.sq)
    posterior[!is.finite(posterior)] <- 0
    posterior <- posterior / sum(posterior)
    start <- max(1L, round(sum(seq_len(max.start) * posterior)))
    if (!is.finite(criterion[start])) {
      stop("ERROR! fcst.forecast: the CVL mean selects an unidentifiable window")
    }
  } else {
    start <- which.min(criterion)
  }
  weights <- weights.rolling(data$N, data$N - start + 1L)
  fit <- WLS.reg(data$y, data$x, weights)
  inner <- list(start = start, crit = criterion, posterior = posterior,
                bp = mapped.bp, cv.method = "rolling", h = h, type = type,
                cv.origins = cv$origins, cv.targets = cv$targets,
                cv.errors = errors, cv.excluded = cv$excluded,
                cv.selection = lapply(cv$folds, function(f) f$data$selection))
  details <- list(start = start, R = data$N - start + 1L, N = data$N,
                  method = method, details = inner, fit = fit)
  class(details) <- "bt_fcstWindow"
  structure(list(method = "window", coefficients = fit$coefficients,
                 weights = weights, bp = mapped.bp, N = data$N, details = details),
            class = "bt_fcst")
}

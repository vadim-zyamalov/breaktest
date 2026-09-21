#' @title
#' Compare forecasts on common origins and targets
#'
#' @description
#' An expanding-origin exercise using only the history available at each
#' origin. Every method's RMSFE uses the same valid rows within a horizon.
#' Failures are reported. If any planned row is excluded, scores are labelled
#' conditional and nominal DM inference is suppressed.
#'
#' @param y,x Data as in [fcst.forecast]. No rows are removed before choosing
#' origins, constructing lags or locating targets.
#' @param methods Strategies in [fcst.forecast]. A named vector permits aliases,
#' e.g. c(AR.direct = "full", AR.iterated = "full").
#' @param n.init Original row of the first forecast origin. The default is
#' max(round(N/2), 10*k), using the largest design among the methods. With
#' selection enabled a method contributes k=1 to this heuristic, avoiding
#' an impossible starting point in high dimensions; set n.init explicitly
#' when the conditioning set or chosen model needs a longer training history.
#' @param benchmark Method name or alias used for relative RMSFE and DM.
#' Added to methods if absent.
#' @param verbose Report progress and the final number of failed fits.
#' @param ... Shared method arguments. Each receives only supported options.
#' Unknown options are errors. Prefer method.options for different settings.
#' @param h,type,lags,constant,x.future,index See [fcst.forecast]. Iterated
#' forecasts with external x require an x.future function; its inputs are
#' always restricted to the relevant origin's history.
#' @param method.options Named list of method-specific argument lists, keyed
#' by method name or alias. Overrides shared settings. Allows different types,
#' lags, window.method, grids, trimming, supplied break dates and other options.
#' Supplied break dates index the original response rows at every origin.
#' @param common.horizons Also use the intersection of valid origins across
#' horizons. Otherwise each horizon ends at N-h and has its own common mask.
#' @param selection Predictor-selection options as in [fcst.forecast]. Each
#' method can override these in method.options. selections retains each
#' origin/method's complete selection diagnostics.
#' @param dm.bandwidth Bartlett bandwidth for ordinary DM, at least h.
#' The default is max(h, floor(P^(1/4))). With dm.hln=TRUE this must be NULL
#' or h, selecting the rectangular estimator through h-1.
#' @param dm.hln Use the Harvey-Leybourne-Newbold correction, requiring an
#' (h-1)-dependent loss differential; see [fcst.test.DM].
#' @param n.cores Number of local SOCK workers. The default 1 is sequential.
#' Parallel tasks are (horizon, origin) pairs; methods and CV within each
#' task remain sequential.
#' @param seed Optional integer for stochastic x.future functions. Assigns
#' one L'Ecuyer-CMRG stream per task, independently of n.cores, and restores
#' the caller's RNG state. With NULL the sequential RNG behaviour is unchanged;
#' random callbacks are not reproducible across parallel runs.
#' @param parallel.export Named list of global objects needed by user
#' callbacks on workers. Self-contained closures need no explicit export.
#' @param parallel.packages Additional installed packages to load on workers.
#'
#' @return For one horizon, a bt_fcstCompare containing forecasts, errors,
#' RMSFE, rel.RMSFE, DM, origins, targets, common (the scoring mask), n.common,
#' n.available per method, failures (origin/target/method/reason), settings,
#' score.scope, DM.messages, and original origin.index/target.index labels.
#' No common valid rows gives
#' NA scores. With multiple horizons, a bt_fcstCompareMulti contains
#' by.horizon with these results and a long summary table.
#'
#' @details
#' A supplied aligned design asserts that `x[t+1]` is known at t; this cannot be
#' inferred from numeric values. Prefer explicit lags for time-series work.
#' Forecasts are fitted even when the realised target is missing, but that
#' target is excluded from every method's descriptive score. Any such
#' exclusion, including a forecast failure, prevents a nominal DM statistic
#' from being reported: conditioning on successful outcomes can alter its
#' null distribution even if the retained dates happen to be consecutive.
#' On a complete sample, DM still requires its loss-process central limit
#' theorem and variance assumptions. A generic expanding-origin comparison
#' does not establish them for every strategy or nesting relationship.
#'
#' @importFrom stats setNames
#'
#' @export
fcst.compare <- function(
  y, x = NULL,
  methods = c("full", "postbreak", "ppp.robust", "avew"),
  n.init = NULL, benchmark = "full", verbose = TRUE, ...,
  h = 1, type = c("direct", "iterated"), lags = NULL, constant = TRUE,
  x.future = NULL, index = NULL, method.options = list(),
  common.horizons = FALSE, dm.bandwidth = NULL, dm.hln = FALSE,
  n.cores = 1L, seed = NULL, parallel.export = list(),
  parallel.packages = character(), selection = NULL
) {
  n.cores <- .fcst.n.cores(n.cores)
  seed <- .fcst.seed(seed)
  .fcst.parallel.options(parallel.export, parallel.packages)
  horizons <- .fcst.horizons(h)
  type <- match.arg(type)
  if (!is.character(methods) || !length(methods) || anyNA(methods) ||
      !is.character(benchmark) || length(benchmark) != 1L || is.na(benchmark)) {
    stop("ERROR! fcst.compare: invalid methods or benchmark")
  }
  labels <- names(methods)
  if (is.null(labels)) labels <- methods
  if (anyNA(labels) || any(labels == "") || anyDuplicated(labels)) {
    stop("ERROR! fcst.compare: use unique names to compare variants of a method")
  }
  names(methods) <- labels
  if (!benchmark %in% labels) methods <- c(setNames(benchmark, benchmark), methods)
  labels <- names(methods)
  if (!is.list(method.options)) stop("ERROR! fcst.compare: method.options must be a list")
  .fcst.check.options(method.options, labels)
  shared <- list(...)
  base.names <- c("window.method", "bp", "trim", "type", "lags", "constant", "x.future", "selection")
  selection <- .fcst.selection.options(selection)
  defaults <- list(type = type, lags = lags, constant = constant, x.future = x.future, selection = selection)
  configs <- lapply(labels, function(label) {
    own <- method.options[[label]]
    if (is.null(own)) own <- list()
    if (!is.list(own)) stop("ERROR! fcst.compare: each method option must be a list")
    window.method <- if (!is.null(own$window.method)) own$window.method
      else if (!is.null(shared$window.method)) shared$window.method else "cv"
    allowed <- union(base.names, .fcst.method.options(methods[[label]], window.method))
    .fcst.check.options(own, allowed)
    cfg <- defaults
    for (nm in intersect(names(shared), allowed)) cfg[nm] <- shared[nm]
    for (nm in names(own)) cfg[nm] <- own[nm]
    for (hh in horizons) {
      .fcst.source.options(methods[[label]], window.method, hh,
                            match.arg(cfg$type, c("direct", "iterated")), cfg,
                            .fcst.selection.options(cfg$selection))
    }
    attr(cfg, "allowed") <- allowed
    cfg
  })
  names(configs) <- labels
  .fcst.check.options(shared, unique(unlist(lapply(configs, attr, "allowed"))))
  histories <- lapply(configs, function(cfg) {
    cfg$type <- match.arg(cfg$type, c("direct", "iterated"))
    if (cfg$type == "iterated" && max(horizons) > 1L && !is.null(x)) {
      if (is.null(cfg$lags) || !is.function(cfg$x.future)) {
        stop("ERROR! fcst.compare: iterated external-x models need explicit lags and an x.future function")
      }
    }
    .fcst.history(y, x, cfg$lags, cfg$constant, index, .fcst.selection.options(cfg$selection))
  })
  history <- histories[[1]]
  N <- history$N
  if (is.null(n.init)) {
    n.init <- max(round(0.5 * N), 10L * max(vapply(histories, function(z) {
      if (is.null(z$selection)) z$k else 1L
    }, 0L)))
  }
  if (!is.numeric(n.init) || length(n.init) != 1L || !is.finite(n.init) ||
      n.init != floor(n.init) || n.init < 1 || n.init > N - max(horizons)) {
    stop("ERROR! fcst.compare: n.init leaves no forecasts for the largest horizon")
  }
  if (!is.logical(common.horizons) || length(common.horizons) != 1L ||
      is.na(common.horizons)) stop("ERROR! fcst.compare: common.horizons must be logical")
  # One cluster covers all horizons. Stable task order also preserves the
  # sequential call order and assigns stochastic streams independently of
  # worker scheduling. Only origin histories are passed to forecast callbacks.
  tasks <- do.call(rbind, lapply(horizons, function(hh) {
    data.frame(h = hh, origin = seq.int(n.init,
               N - if (common.horizons) max(horizons) else hh))
  }))
  context <- list(history = history, histories = histories, configs = configs,
                  methods = methods, labels = labels)
  fitted <- .fcst.parallel.map(
    nrow(tasks), .fcst.compare.origin, context, n.cores = n.cores,
    make.task = function(i) tasks[i, ], seed = seed, verbose = verbose,
    label = "fcst.compare origins", export = parallel.export,
    packages = parallel.packages
  )
  by.horizon <- lapply(horizons, function(hh) {
    origins <- seq.int(n.init, N - if (common.horizons) max(horizons) else hh)
    targets <- origins + hh
    values <- fitted[tasks$h == hh]
    forecasts <- do.call(rbind, lapply(values, function(z) z$forecasts))
    dimnames(forecasts) <- list(as.character(origins), labels)
    failures <- do.call(rbind, lapply(values, function(z) z$failures))
    rownames(failures) <- NULL
    errors <- history$y[targets, 1] - forecasts
    result <- list(forecasts = forecasts, errors = errors,
                   common = rowSums(!is.finite(errors)) == 0L,
                   n.available = colSums(is.finite(errors)), failures = failures,
                   methods = labels, strategies = methods, benchmark = benchmark,
                   h = hh, origins = origins, targets = targets,
                   origin.index = history$index[origins], target.index = history$index[targets],
                   settings = configs,
                   selections = setNames(lapply(values, function(z) z$selections), as.character(origins)))
    result
  })
  if (common.horizons) {
    common <- Reduce(intersect, lapply(by.horizon, function(z) z$origins[z$common]))
    by.horizon <- lapply(by.horizon, function(z) {
      z$common <- z$origins %in% common
      z
    })
  }
  by.horizon <- lapply(by.horizon, .fcst.compare.score,
                       dm.bandwidth = dm.bandwidth, dm.hln = dm.hln)
  names(by.horizon) <- as.character(horizons)
  failed <- sum(vapply(by.horizon, function(z) nrow(z$failures), 0L))
  if (verbose && failed) message("fcst.compare: ", failed, " failed fits; see failures")
  if (length(horizons) == 1L) return(by.horizon[[1]])
  summary <- do.call(rbind, lapply(by.horizon, function(z) {
    data.frame(h = z$h, method = z$methods, RMSFE = unname(z$RMSFE),
               rel.RMSFE = unname(z$rel.RMSFE), DM = unname(z$DM),
               n.common = z$n.common, n.available = unname(z$n.available),
               score.scope = z$score.scope)
  }))
  rownames(summary) <- NULL
  structure(list(h = horizons, by.horizon = by.horizon, summary = summary,
                 methods = labels, benchmark = benchmark,
                 common.horizons = common.horizons), class = "bt_fcstCompareMulti")
}


# A complete origin is the unit of parallel work. Keeping methods together
# avoids nested clusters and leaves all horizon-aware inner CV unchanged.
.fcst.compare.origin <- function(task, context) {
  s <- task$origin
  hh <- task$h
  history <- context$history
  forecasts <- setNames(rep(NA_real_, length(context$labels)), context$labels)
  selections <- setNames(vector("list", length(context$labels)), context$labels)
  failures <- data.frame(origin = integer(), target = integer(),
                         method = character(), reason = character())
  yy <- history$y[seq_len(s), , drop = FALSE]
  xx <- if (!is.null(history$x)) history$x[seq_len(s), , drop = FALSE] else NULL
  for (label in context$labels) {
    cfg <- context$configs[[label]]
    f <- tryCatch({
      xn <- if (is.null(cfg$lags))
        context$histories[[label]]$z[s + 1L, , drop = FALSE] else NULL
      args <- cfg
      # The legacy matrix API uses complete-case break indexing, whereas
      # comparison accepts dates on the original response time grid.
      if (.fcst.legacy(hh, cfg$lags, context$methods[[label]], cfg, cfg$selection)) {
        data <- .fcst.align(context$histories[[label]], s, 1L)
        if (!is.null(args$bp)) args$bp <- .fcst.break.row(args$bp, data, s)
        args <- .fcst.map.options(args, data, s)
      }
      do.call(fcst.forecast,
              c(list(y = yy, x = xx, x.new = xn, method = context$methods[[label]],
                     h = hh, index = history$index[seq_len(s)]), args))
    }, error = function(e) e)
    if (inherits(f, "error")) {
      failures <- rbind(failures, data.frame(origin = s, target = s + hh,
                          method = label, reason = conditionMessage(f)))
    } else {
      forecasts[label] <- f$forecast
      selections[label] <- list(f$selection)
    }
  }
  list(forecasts = forecasts, failures = failures, selections = selections)
}


.fcst.compare.score <- function(result, dm.bandwidth, dm.hln) {
  err <- result$errors
  err[!result$common, ] <- NA_real_
  n <- sum(result$common)
  rmsfe <- if (n) sqrt(colSums(err^2, na.rm = TRUE) / n)
    else setNames(rep(NA_real_, ncol(err)), colnames(err))
  dm <- setNames(rep(NA_real_, ncol(err)), colnames(err))
  dm.messages <- setNames(rep(NA_character_, ncol(err)), colnames(err))
  for (method in setdiff(result$methods, result$benchmark)) {
    if (!all(result$common)) {
      dm.messages[method] <- "DM suppressed: planned evaluation rows were excluded; common-row scores are conditional on valid outcomes"
      next
    }
    test <- tryCatch(
      fcst.test.DM(err[, method]^2, err[, result$benchmark]^2,
                   bandwidth = dm.bandwidth, hln = dm.hln, h = result$h,
                   time.index = result$origins),
      error = function(e) e
    )
    if (inherits(test, "error")) dm.messages[method] <- conditionMessage(test)
    else dm[method] <- test$statistic
  }
  result$RMSFE <- rmsfe
  result$rel.RMSFE <- if (is.finite(rmsfe[result$benchmark]) &&
                         rmsfe[result$benchmark] > 0) rmsfe / rmsfe[result$benchmark]
    else rmsfe * NA_real_
  result$DM <- dm
  result$DM.messages <- dm.messages
  result$score.scope <- if (all(result$common)) "complete planned evaluation sample"
    else "conditional on common valid outcomes"
  result$DM.assumptions.verified <- FALSE
  result$n.common <- n
  class(result) <- "bt_fcstCompare"
  result
}

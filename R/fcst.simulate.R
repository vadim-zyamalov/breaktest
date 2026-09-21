#' @title
#' Reproducible forecast experiments under parameter instability
#'
#' @description
#' The ten mean-instability designs in Section 4 of Hirano-Wright, and an
#' end-of-sample comparison using [fcst.compare]. Forecasts see only their
#' training history. The future signal and observations are used for scoring.
#'
#' @param N Training sample length; a vector in fcst.simulate.
#' @param design Design numbers 1:10; a vector in fcst.simulate. 1: stable;
#' 2:4: a single jump of 10/sqrt(N) after N*(.25,.5,.75);
#' 5:7: random-walk signal innovations with SD (1,5,10)/N;
#' 8:10: compound Poisson jumps with rate (1,2,3)/N per period and
#' jump SD 3/sqrt(N). The pre-sample signal `b[0]` is zero in every design;
#' a stochastic increment may already occur at observation 1.
#' @param h Positive forecast horizon(s).
#' @param rho AR(1) coefficient of the observation disturbance, |rho|<1.
#' Its innovation SD is 1-rho, so the long-run variance is one. The initial
#' disturbance is drawn from its stationary distribution. A vector is
#' allowed in fcst.simulate.
#' @param seed Optional integer seed. Explicit seeds restore the caller's
#' RNG state; NULL consumes it normally. Simulation uses an independent
#' L'Ecuyer stream per (design,rho,N,replication) task, in fixed grid order.
#'
#' @details
#' fcst.DGP returns all N+max(h) observations, signal and disturbances,
#' actual jump sizes/counts and dates. Multiple Poisson events in one period
#' are accumulated, not discarded. The signal starts from `b[0]=0`, as
#' specified in the paper. AR disturbances implement the paper's
#' serial-correlation experiments, which the supplied mcprog.m omits.
#'
#' The future latent signal `b[N+h]` is the target of the paper's scaled risk.
#' It is distinct from both observed `y[N+h]` and `E[y[N+h]|b[N],u[N]]`, which is
#' `b[N]`+rho^h*`u[N]` for these end-of-sample designs. All three targets are
#' returned explicitly. Simulated latent states never enter forecast fitting.
#' The stationary AR initialisation makes an arbitrary burn-in unnecessary.
#' Changing N, rho or h defines a new experiment within these DGP formulas;
#' it does not claim to reproduce an unreported table or extend a forecasting
#' method's theoretical domain. Each method retains its own source restrictions.
#'
#' @return fcst.DGP returns a bt_fcstDGP; fcst.simulate returns a
#' bt_fcstSimulation with settings, per-replication losses, failures, breaks,
#' summary (MSFE, RMSFE, latent-signal and oracle-mean risk, MC standard
#' errors), conditional.summary, and optional generated data. Unconditional
#' summary risk and MCSE are returned only if every planned replication has
#' valid forecasts from every method at every requested horizon. Otherwise
#' these entries are NA; conditional.summary reports the selected common
#' successes with an explicit scope label. Those conditional risks are not
#' estimates of the original experiment's unconditional risk. Raw losses and
#' failures remain visible. MCSE is across replications, not time points.
#'
#' @references
#' Hirano, K. and J. H. Wright (2022). "Analyzing Cross-Validation for
#' Forecasting with Structural Instability." Journal of Econometrics 226,
#' 139-154. Section 4 and equation (4.1).
#' https://doi.org/10.1016/j.jeconom.2020.10.009.
#'
#' @importFrom stats rnorm rpois sd
#' @export
fcst.DGP <- function(N = 100L, design = 1L, h = 1L, rho = 0, seed = NULL) {
  .fcst.sim.arguments(N, design, rho)
  if (length(N) != 1L || length(design) != 1L || length(rho) != 1L) {
    stop("ERROR! fcst.DGP: N, design and rho must be scalars")
  }
  horizons <- .fcst.horizons(h)
  seed <- .fcst.seed(seed)
  if (!is.null(seed)) {
    state <- .fcst.rng.state()
    on.exit(.fcst.rng.restore(state), add = TRUE)
    set.seed(seed)
  }
  total <- N + max(horizons)
  signal <- jumps <- numeric(total)
  counts <- integer(total)
  if (design %in% 2:4) {
    bp <- floor(N * c(.25, .5, .75)[design - 1L])
    jumps[bp + 1L] <- 10 / sqrt(N)
    counts[bp + 1L] <- 1L
    signal <- cumsum(jumps)
  } else if (design %in% 5:7) {
    signal <- cumsum(rnorm(total, sd = c(1, 5, 10)[design - 4L] / N))
  } else if (design %in% 8:10) {
    counts <- rpois(total, lambda = (design - 7L) / N)
    # The sum of c independent N(0,s^2) jumps is N(0,c*s^2).
    jumps <- rnorm(total, sd = 3 * sqrt(counts / N))
    signal <- cumsum(jumps)
  }
  disturbance <- numeric(total)
  previous <- rnorm(1L, sd = (1 - rho) / sqrt(1 - rho^2))
  innovations <- rnorm(total, sd = 1 - rho)
  for (t in seq_len(total)) {
    disturbance[t] <- rho * previous + innovations[t]
    previous <- disturbance[t]
  }
  target <- N + horizons
  breaks <- which(counts > 0L) - 1L
  structure(list(y = signal + disturbance, signal = signal,
                  disturbance = disturbance, innovations = innovations,
                  jumps = jumps, jump.counts = counts, breaks = breaks,
                  N = N, design = design, rho = rho, h = horizons, seed = seed,
                  targets = target, future.signal = signal[target],
                  future.actual = (signal + disturbance)[target],
                  oracle.mean = signal[N] + rho^horizons * disturbance[N]),
            class = "bt_fcstDGP")
}


.fcst.sim.arguments <- function(N, design, rho) {
  if (!is.numeric(N) || !length(N) || any(!is.finite(N)) ||
      any(N < 6 | N != floor(N))) stop("ERROR! fcst.simulate: N must contain integers >= 6")
  if (!is.numeric(design) || !length(design) || any(!is.finite(design)) ||
      any(design != floor(design) | design < 1 | design > 10)) {
    stop("ERROR! fcst.simulate: design must contain integers 1:10")
  }
  if (!is.numeric(rho) || !length(rho) || any(!is.finite(rho)) || any(abs(rho) >= 1)) {
    stop("ERROR! fcst.simulate: rho must be strictly between -1 and 1")
  }
  invisible(NULL)
}


#' @rdname fcst.DGP
#' @param n.rep Number of independent replications per design cell.
#' @param methods,benchmark,method.options As in [fcst.compare].
#' @param forecast.options Named list of shared forecasting options, e.g.
#' list(lags=1:2, type="direct", trim=.2). No future latent states are passed.
#' @param n.cores Number of SOCK workers. Parallelism is across replications;
#' each inner forecast comparison remains sequential.
#' @param keep.data Retain complete generated paths; FALSE keeps only the
#' targets, losses, jump dates and settings needed for compact reporting.
#' @param verbose Print progress from the main process.
#' @export
fcst.simulate <- function(N = 100L, design = 1:10, n.rep = 100L, h = 1L,
                          rho = 0, methods = c("full", "postbreak", "exps"),
                          benchmark = "full", method.options = list(),
                          forecast.options = list(), seed = 123L, n.cores = 1L,
                          keep.data = FALSE, verbose = TRUE) {
  .fcst.sim.arguments(N, design, rho)
  h <- .fcst.horizons(h)
  if (!is.numeric(n.rep) || length(n.rep) != 1L || !is.finite(n.rep) ||
      n.rep < 1 || n.rep != floor(n.rep)) {
    stop("ERROR! fcst.simulate: n.rep must be a positive integer")
  }
  n.rep <- as.integer(n.rep)
  n.cores <- .fcst.n.cores(n.cores)
  seed <- .fcst.seed(seed)
  if (!is.logical(keep.data) || length(keep.data) != 1L || is.na(keep.data)) {
    stop("ERROR! fcst.simulate: keep.data must be logical")
  }
  if (!is.list(forecast.options) || !is.list(method.options)) {
    stop("ERROR! fcst.simulate: forecast.options and method.options must be lists")
  }
  .fcst.check.options(forecast.options, names(forecast.options))
  reserved <- c("y", "x", "h", "n.init", "common.horizons", "methods", "benchmark",
                "method.options", "n.cores", "seed", "verbose", "index")
  if (any(names(forecast.options) %in% reserved)) {
    stop("ERROR! fcst.simulate: forecast.options contains a simulation-controlled argument")
  }
  grid <- expand.grid(replication = seq_len(n.rep), N = unique(N),
                       design = unique(design), rho = unique(rho))
  context <- list(h = h, methods = methods, benchmark = benchmark,
                  method.options = method.options, options = forecast.options,
                  keep.data = keep.data)
  values <- .fcst.parallel.map(nrow(grid), .fcst.simulate.rep, context,
                                n.cores = n.cores, make.task = function(i) grid[i, ],
                                seed = seed, verbose = verbose, label = "fcst.simulate replications")
  losses <- do.call(rbind, lapply(values, function(z) z$losses))
  failures <- do.call(rbind, lapply(values, function(z) z$failures))
  rownames(losses) <- rownames(failures) <- NULL
  key <- interaction(losses$N, losses$design, losses$rho, losses$h, losses$method, drop = TRUE)
  conditional.summary <- do.call(rbind, lapply(unique(key), function(label) {
    rows <- losses[key == label, , drop = FALSE]
    good <- rows[rows$common, , drop = FALSE]
    average <- function(x) if (length(x)) mean(x) else NA_real_
    mcse <- function(x) if (length(x) > 1L) sd(x) / sqrt(length(x)) else NA_real_
    cbind(rows[1L, c("N", "design", "rho", "h", "method"), drop = FALSE],
          data.frame(n.rep = nrow(rows), n.available = sum(is.finite(rows$forecast)),
                     n.common = nrow(good), MSFE = average(good$observation.loss),
                     RMSFE = sqrt(average(good$observation.loss)),
                     MCSE.MSFE = mcse(good$observation.loss),
                     signal.risk = average(good$signal.loss),
                     scaled.risk = rows$N[1] * average(good$signal.loss),
                     MCSE.scaled.risk = rows$N[1] * mcse(good$signal.loss),
                     mean.risk = average(good$mean.loss), MCSE.mean.risk = mcse(good$mean.loss),
                     risk.scope = "conditional on joint forecast success"))
  }))
  rownames(conditional.summary) <- NULL
  summary <- conditional.summary
  complete <- summary$n.common == summary$n.rep
  risk.fields <- c("MSFE", "RMSFE", "MCSE.MSFE", "signal.risk", "scaled.risk",
                   "MCSE.scaled.risk", "mean.risk", "MCSE.mean.risk")
  summary[!complete, risk.fields] <- NA_real_
  summary$risk.scope <- ifelse(complete, "all planned replications",
                               "unconditional risk unavailable: incomplete replications")
  structure(list(summary = summary, conditional.summary = conditional.summary,
                  losses = losses, failures = failures,
                  breaks = lapply(values, function(z) z$breaks), grid = grid,
                  data = if (keep.data) lapply(values, function(z) z$data) else NULL,
                  settings = c(context, list(N = N, design = design, rho = rho,
                                             n.rep = n.rep, seed = seed, n.cores = n.cores))),
            class = "bt_fcstSimulation")
}


.fcst.simulate.rep <- function(task, context) {
  d <- fcst.DGP(task$N, task$design, context$h, task$rho)
  compare <- do.call(fcst.compare, c(list(y = d$y, h = context$h, n.init = task$N,
                                         common.horizons = TRUE, methods = context$methods,
                                         benchmark = context$benchmark,
                                         method.options = context$method.options,
                                         n.cores = 1L, verbose = FALSE), context$options))
  results <- if (length(context$h) == 1L) list(compare) else compare$by.horizon
  losses <- failures <- vector("list", length(results))
  for (j in seq_along(results)) {
    result <- results[[j]]
    f <- as.numeric(result$forecasts[1L, ])
    losses[[j]] <- data.frame(replication = task$replication, N = task$N,
                              design = task$design, rho = task$rho, h = result$h,
                              method = result$methods, forecast = f,
                              actual = d$future.actual[j], signal = d$future.signal[j],
                              oracle.mean = d$oracle.mean[j], common = result$common[1L],
                              observation.loss = (f - d$future.actual[j])^2,
                              signal.loss = (f - d$future.signal[j])^2,
                              mean.loss = (f - d$oracle.mean[j])^2)
    fail <- result$failures
    failures[[j]] <- cbind(task[rep(1L, nrow(fail)), , drop = FALSE], fail)
  }
  list(losses = do.call(rbind, losses), failures = do.call(rbind, failures),
       breaks = d$breaks, data = if (context$keep.data) d else NULL)
}

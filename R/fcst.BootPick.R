#' @title
#' Test whether modelling a structural break improves forecast accuracy
#'
#' @description
#' Boot and Pick (2020) point out that the classical break tests answer the
#' wrong question. Under the MSFE loss there is a bias-variance trade-off, so a
#' break that is statistically detectable may still be too small to be worth
#' modelling. Their null hypothesis is therefore *equal forecast accuracy*
#' between the post-break (or the combined) forecast and the full-sample
#' forecast, which maps into a non-zero critical break magnitude.
#'
#' Two features distinguish the resulting statistic from the sup-Wald of
#' Andrews (1993):
#' * the coefficient difference is weighted by the regressors of the forecast
#'   origin, so a large break in a coefficient can be irrelevant for the
#'   forecast while a small one can be decisive;
#' * under a known break date the statistic is \eqn{\chi^2(1, 1)} rather than
#'   \eqn{\chi^2(1)}, and the 5% critical value is 7.00 instead of 3.84.
#'
#' @param y A dependent variable.
#' @param x Explanatory variables. If `NULL` a single constant is used.
#' @param x.new A (1 x k) vector of regressor values for the forecast origin.
#' Defaults to the last row of `x`.
#' @param trim A trimming fraction bounding the admissible break dates.
#' @param level A nominal size, one of 0.01, 0.05 or 0.10.
#' @param target Which comparison is tested:
#' * `postbreak`: post-break forecast versus full sample,
#' * `combined`: the combined forecast versus full sample.
#' @param break.vars Indices of the columns of `x` subject to the break.
#' `NULL` means a pure structural change in all the coefficients. The empirical
#' part of the original paper allows a break in the intercept only, which
#' corresponds to `break.vars = 1`.
#' @param vcov.type One of `const`, `HC0`, `HC3` or `HAC`. The original code uses
#' `HC3` in its heteroskedasticity robust mode.
#' @param cv.table An optional replacement for the built-in critical values,
#' either the full output of [fcst.BP.simulate] or its target-specific branch.
#' Custom branches require metadata with `trim` and `target` fields.
#'
#' @param bandwidth Optional HAC bandwidth (products through bandwidth-1).
#' @param time.index Original integer-period positions. The estimation sample
#' must be a contiguous block; internal missing observations are unsupported.
#'
#' @return An object of class `bt_fcstBP`, a list of:
#' * `W`, `S`: the two test statistics evaluated at the estimated break date,
#' * `cv.W`, `cv.S`, `zeta.sqrt`: critical values and the equal performance
#'   boundary at the estimated break date,
#' * `reject.W`, `reject.S`: test decisions,
#' * `bp`, `tau`: the estimated break date and fraction,
#' * `forecast`: a named vector with the `full`, `postbreak` and `combined`
#'   point forecasts,
#' * `weight.full`: the weight \eqn{1 / (1 + W)} placed on the full-sample
#'   forecast in the combination,
#' * `W.all`: the whole profile of the statistic over the admissible dates.
#' * `calibrated`, `inference`: whether the table matches the search region
#'   and a diagnostic. Unsupported calibration retains point forecasts and
#'   W, but S, critical values and test decisions are NA.
#'
#' @details
#' The statistic is
#' \deqn{\sup_\tau W(\tau) = \sup_\tau
#' \frac{\left[x_{N+1}' (\hat\beta_1(\tau) - \hat\beta_2(\tau))\right]^2}
#' {x_{N+1}' \widehat{Var}(\hat\beta_1(\tau) - \hat\beta_2(\tau)) x_{N+1}},}
#' and the size-corrected version is
#' \eqn{S(\hat\tau) = \sqrt{W(\hat\tau)} - \zeta^{1/2}(\hat\tau)}, whose
#' critical values are almost independent of the break date.
#'
#' The combined forecast reproduces the optimally weighted forecast of Pesaran
#' et al. (2013) written as a convex combination:
#' \deqn{\hat y^c = \frac{1}{1 + W} \hat y^{full} +
#' \frac{W}{1 + W} \hat y^{post}.}
#' The entire requested trimming region must be identified. Singular split
#' designs, nonpositive projected variances or sample-size restrictions that
#' shorten that region produce an error, rather than a different search rule.
#' Partial changes and consistent HC/HAC covariance estimates are supported
#' under the paper's moment and functional-limit assumptions. These assumptions
#' are not established by the covariance option or the table compatibility flag.
#'
#' @references
#' Boot, Tom, and Andreas Pick.
#' “Does Modeling a Structural Break Improve Forecast Accuracy?”
#' Journal of Econometrics 215 (2020): 35-59. Local author version: 2019.
#' https://doi.org/10.1016/j.jeconom.2019.07.007.
#'
#' @export
fcst.test.BP <- function(
  y,
  x = NULL,
  x.new = NULL,
  trim = 0.15,
  level = 0.05,
  target = c("postbreak", "combined"),
  break.vars = NULL,
  vcov.type = "HC3",
  cv.table = NULL,
  bandwidth = NULL,
  time.index = NULL
) {
  target <- match.arg(target)
  .fcst.test.level(level)
  if (!is.numeric(trim) || length(trim) != 1L || !is.finite(trim) ||
      trim <= 0 || trim >= .5) {
    stop("ERROR! fcst.test.BP: trim must be in (0,.5)")
  }

  time.index <- .fcst.time.index(time.index, NROW(y))
  .d <- .fcst.data(y, x, x.new)
  .fcst.inference.data(.d)
  time.index <- time.index[.d$rows]
  .fcst.require.regular(time.index, "fcst.test.BP")
  y <- .d$y
  x <- .d$x
  x.new <- .d$x.new
  N <- .d$N
  k <- .d$k

  break.vars <- .fcst.break.columns(break.vars, k)
  nb <- length(break.vars)
  if (N <= k + nb) stop("ERROR! fcst.test.BP: no residual degrees of freedom")

  requested <- c(max(1L, ceiling(trim * N)), min(N - 1L, floor((1 - trim) * N)))
  bp.min <- max(requested[1], nb + 1L)
  bp.max <- min(requested[2], N - nb - 1L)

  if (bp.min > bp.max) {
    stop("ERROR! fcst.test.BP: trimming leaves no admissible break dates")
  }
  if (bp.min != requested[1] || bp.max != requested[2]) {
    stop("ERROR! fcst.test.BP: sample size does not identify the full requested search region")
  }

  bp.all <- bp.min:bp.max
  W.all <- rep(NA_real_, length(bp.all))
  fcast.post <- rep(NA_real_, length(bp.all))

  for (i in seq_along(bp.all)) {
    .s <- .fcst.split.design(x, bp.all[i], break.vars)

    xs <- .s$xs
    if (qr(xs)$rank < ncol(xs) || rcond(crossprod(xs)) < .Machine$double.eps) {
      stop("ERROR! fcst.test.BP: singular design inside the requested search region")
    }

    .m <- OLS.reg(y, xs)
    vcov <- .fcst.vcov(xs, .m$residuals, vcov.type, bandwidth, time.index)

    # R b = beta_1 - beta_2 for the breaking regressors
    R <- matrix(0, nrow = length(break.vars), ncol = ncol(xs))
    R[cbind(seq_along(.s$idx.1), .s$idx.1)] <- 1
    R[cbind(seq_along(.s$idx.2), .s$idx.2)] <- -1

    xb <- x.new[, break.vars, drop = FALSE]

    d <- drop(xb %*% R %*% .m$coefficients)
    v <- drop(xb %*% R %*% vcov %*% t(R) %*% t(xb))

    if (!is.finite(v) || v <= 0) {
      stop("ERROR! fcst.test.BP: nonpositive projected variance inside the requested search region")
    }

    W.all[i] <- d^2 / v
    fcast.post[i] <- drop(x.new %*% .m$coefficients[.s$map.2])
  }

  i.max <- which.max(W.all)
  W <- W.all[i.max]
  bp <- bp.all[i.max]
  tau <- bp / N

  .cv <- get.cv.BootPick(tau, level, target, cv.table, trim = trim)
  S <- sqrt(W) - .cv$zeta.sqrt

  full.fit <- OLS.reg(y, x)
  y.full <- .fcst.predict(full.fit, x.new)
  y.post <- fcast.post[i.max]
  w.full <- 1 / (1 + W)
  y.comb <- w.full * y.full + (1 - w.full) * y.post
  split <- .fcst.split.design(x, bp, break.vars)
  beta.post <- OLS.reg(y, split$xs)$coefficients[split$map.2]
  beta.full <- full.fit$coefficients

  result <- list(
    W = W,
    S = S,
    cv.W = .cv$cv.W,
    cv.S = .cv$cv.S,
    zeta.sqrt = .cv$zeta.sqrt,
    reject.W = W > .cv$cv.W,
    reject.S = S > .cv$cv.S,
    bp = bp,
    tau = tau,
    forecast = c(full = y.full, postbreak = y.post, combined = y.comb),
    weight.full = w.full,
    coefficients = w.full * beta.full + (1 - w.full) * beta.post,
    beta.full = beta.full,
    beta.post = beta.post,
    bp.all = bp.all,
    W.all = W.all,
    level = level,
    target = target,
    calibrated = .cv$calibrated,
    inference = .cv$inference,
    table.metadata = .cv$metadata,
    trim = trim,
    effective.trim = range(bp.all[is.finite(W.all)]) / N,
    rows = .d$rows,
    time.index = time.index,
    break.vars = break.vars,
    vcov.type = vcov.type,
    bandwidth = bandwidth
  )

  class(result) <- "bt_fcstBP"
  result
}


#' @title
#' Simulate the asymptotic distribution of the Boot-Pick test
#'
#' @description
#' Regenerates the equal performance boundary \eqn{\zeta^{1/2}(\tau_b)} and the
#' critical values of the `W` and `S` statistics for an arbitrary trimming and
#' grid. This is a port of `loopbatch.m`, `invar_calc.m`, `det_equalperf.m` and
#' `loopbatch_onlyboundary.m` and `calc_cvsimple.m` of the original MATLAB
#' package, condensed into a single routine. A first simulation estimates the
#' risk boundary; a second independent simulation evaluates the process at
#' that exact interpolated boundary, without rounding the break magnitude.
#'
#' The published tables shipped with the package (see [fcst.BP.tables]) were
#' produced with `n.rep = 50000`, `n.grid = 1000`, `trim = 0.15` and
#' `tau.step = 0.05`. Reproducing them takes hours; the defaults here are much
#' lighter and are meant for experimentation.
#'
#' @param trim A trimming fraction bounding the admissible break dates.
#' @param tau.step A step of the break fraction grid.
#' @param theta.max,theta.step The range and the step of the standardised break
#' magnitude grid.
#' @param n.grid A number of points the unit interval is divided into.
#' @param n.rep A number of simulated Brownian motions.
#' @param levels Nominal sizes the critical values are needed for.
#' @param seed An optional random seed. As in the original simulator, set.seed
#' is applied in the calling process and simulation advances its RNG stream.
#' Innovations for both passes are drawn there in replication order for every
#' n.cores. Each pass uses n.grid*n.rep normal draws; the second pass shares
#' innovations across all break fractions and both forecast targets.
#' @param verbose Whether to report the progress.
#' @param n.cores Number of local SOCK workers. The default 1 is sequential.
#' Work is divided into fixed blocks of 64 replications; at most one worker
#' per block is started. Critical values and risks agree across worker counts.
#'
#' @return A list of:
#' * `tau`: the break fraction grid,
#' * `theta`: the break magnitude grid,
#' * `postbreak`, `combined`: each contains `tau`, the equal-performance
#'   boundary `zeta.sqrt` and matrices `W` / `S` with critical values by `levels`,
#' * `risk`: a list of the three risk profiles over the `(tau, theta)` grid,
#' * `n.rep`: the number of replications per pass,
#' * `metadata`: settings needed to reproduce and validate the calibration.
#'
#' @details
#' Under the local-to-zero asymptotics the statistic converges to
#' \deqn{Q^*(\tau) = \left[Z(\tau) + \mu(\tau; \theta_{\tau_b})\right]^2,
#' \quad Z(\tau) = \frac{B(\tau) - \tau B(1)}{\sqrt{\tau (1 - \tau)}},}
#' \deqn{\mu(\tau; \theta) = \theta \left[
#' \sqrt{\tfrac{1 - \tau}{\tau}} \tau_b I(\tau_b < \tau) +
#' \sqrt{\tfrac{\tau}{1 - \tau}} (1 - \tau_b) I(\tau_b \geq \tau)\right].}
#' The standardised risks of the three forecasts are then
#' \eqn{r^{full} = B(1) + \tau_b \theta},
#' \eqn{r^{post} = (1 - \hat\tau)^{-1}[B(1) - B(\hat\tau) +
#' \theta (\tau_b - \hat\tau) I(\hat\tau < \tau_b)]} and
#' \eqn{r^{comb} = (1 + W)^{-1} r^{full} + W (1 + W)^{-1} r^{post}}.
#'
#' @importFrom stats rnorm quantile approx
#'
#' @export
fcst.BP.simulate <- function(
  trim = 0.15,
  tau.step = 0.05,
  theta.max = 40,
  theta.step = 1,
  n.grid = 1000,
  n.rep = 2000,
  levels = c(0.01, 0.05, 0.10),
  seed = NULL,
  verbose = TRUE,
  n.cores = 1L
) {
  n.cores <- .fcst.n.cores(n.cores)
  .fcst.seed(seed)
  for (value in list(n.grid, n.rep)) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 1 || value != floor(value)) {
      stop("ERROR! fcst.BP.simulate: n.grid and n.rep must be positive integers")
    }
  }
  for (value in list(tau.step, theta.max, theta.step)) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0) {
      stop("ERROR! fcst.BP.simulate: grid steps and theta.max must be positive")
    }
  }
  if (!is.numeric(trim) || length(trim) != 1L || !is.finite(trim) ||
      trim <= 0 || trim >= 0.5 || round(trim * n.grid) < 1 ||
      round((1 - trim) * n.grid) >= n.grid ||
      ceiling(trim * n.grid) > floor((1 - trim) * n.grid) || theta.step > theta.max) {
    stop("ERROR! fcst.BP.simulate: empty or invalid simulation grid")
  }
  if (!is.numeric(levels) || !length(levels) || any(!is.finite(levels)) ||
      any(levels <= 0 | levels >= 1) || anyDuplicated(levels)) {
    stop("ERROR! fcst.BP.simulate: levels must be probabilities in (0, 1)")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  tau.all <- seq(trim, 1 - trim, by = tau.step)
  # Both endpoints are needed to evaluate S at every possible estimated date.
  if (tail(tau.all, 1L) < 1 - trim - 1e-12) tau.all <- c(tau.all, 1 - trim)
  else tau.all[length(tau.all)] <- 1 - trim
  theta.all <- seq(0, theta.max, by = theta.step)

  n.tau <- length(tau.all)
  n.theta <- length(theta.all)

  # Grid of the scanning region
  t.idx <- ceiling(trim * n.grid):floor((1 - trim) * n.grid)
  t.val <- t.idx / n.grid
  n.t <- length(t.val)

  # Quantities invariant to the Brownian motion realisation (invar_calc.m)
  J1a <- 1 / sqrt(t.val * (1 - t.val))
  mu.arr <- array(0, dim = c(n.t, n.tau, n.theta))
  for (i.th in seq_len(n.theta)) {
    for (i.tb in seq_len(n.tau)) {
      tb <- tau.all[i.tb]
      pre <- t.val < tb
      mu.arr[pre, i.tb, i.th] <-
        sqrt(t.val[pre] / (1 - t.val[pre])) * (1 - tb) * theta.all[i.th]
      mu.arr[!pre, i.tb, i.th] <-
        sqrt((1 - t.val[!pre]) / t.val[!pre]) * tb * theta.all[i.th]
    }
  }

  r.full <- matrix(0, n.tau, n.theta)
  r.post <- matrix(0, n.tau, n.theta)
  r.comb <- matrix(0, n.tau, n.theta)

  # Fixed blocks keep reduction order independent of the worker count.
  # Innovations are drawn on the master, lazily, preserving the legacy RNG
  # sequence without allocating an n.grid by n.rep matrix.
  block.size <- 64L
  blocks <- ceiling(n.rep / block.size)
  context <- list(n.grid = n.grid, t.idx = t.idx, t.val = t.val, J1a = J1a,
                  mu.arr = mu.arr, tau.all = tau.all, theta.all = theta.all)
  make.task <- function(block) {
    first <- (block - 1L) * block.size + 1L
    size <- min(block.size, n.rep - first + 1L)
    list(first = first, innovations = matrix(rnorm(n.grid * size), n.grid, size))
  }
  receive <- function(value) {
    r.full <<- r.full + value$risk$full
    r.post <<- r.post + value$risk$postbreak
    r.comb <<- r.comb + value$risk$combined
    invisible(NULL)
  }
  .fcst.parallel.map(blocks, .fcst.BP.chunk, context, n.cores,
                      make.task = make.task, receive = receive, verbose = verbose,
                      label = "fcst.BP.simulate blocks")

  r.full <- r.full / n.rep
  r.post <- r.post / n.rep
  r.comb <- r.comb / n.rep

  # Break magnitude at which the two risks cross (det_equalperf.m)
  .crossing <- function(r.ref, r.alt) {
    sapply(seq_len(n.tau), function(j) {
      d <- r.ref[j, ] - r.alt[j, ]
      i <- which(d[-1] >= 0 & d[-n.theta] < 0)
      if (length(i) == 0) {
        return(NA_real_)
      }
      i <- i[1]
      phi <- abs(d[i]) / abs(d[i + 1] - d[i])
      (1 - phi) * theta.all[i] + phi * theta.all[i + 1]
    })
  }

  theta.post <- .crossing(r.full, r.post)
  theta.comb <- .crossing(r.full, r.comb)

  zeta.post <- sqrt(tau.all * (1 - tau.all)) * theta.post
  zeta.comb <- sqrt(tau.all * (1 - tau.all)) * theta.comb

  # Independent second pass, at the interpolated null boundary (the authors'
  # loopbatch_onlyboundary.m). Every target/date shares the same innovations.
  boundary <- cbind(postbreak = theta.post, combined = theta.comb)
  zeta <- cbind(postbreak = zeta.post, combined = zeta.comb)
  drift <- array(NA_real_, c(n.t, n.tau, 2L))
  for (target in 1:2) for (j in seq_len(n.tau)) {
    tb <- tau.all[j]
    drift[, j, target] <- boundary[j, target] * ifelse(t.val < tb,
      sqrt(t.val / (1 - t.val)) * (1 - tb), sqrt((1 - t.val) / t.val) * tb)
  }
  supW <- supS <- array(NA_real_, c(n.rep, n.tau, 2L))
  boundary.context <- list(n.grid = n.grid, t.idx = t.idx, t.val = t.val,
                           J1a = J1a, drift = drift, tau = tau.all, zeta = zeta)
  receive.boundary <- function(value) {
    idx <- seq.int(value$first, length.out = dim(value$W)[1L])
    supW[idx, , ] <<- value$W
    supS[idx, , ] <<- value$S
    invisible(NULL)
  }
  .fcst.parallel.map(blocks, .fcst.BP.boundary.chunk, boundary.context, n.cores,
                      make.task = make.task, receive = receive.boundary,
                      verbose = verbose, label = "fcst.BP.simulate boundary blocks")

  .cvals <- function(target) {
    cv.W <- matrix(NA_real_, length(levels), n.tau)
    cv.S <- matrix(NA_real_, length(levels), n.tau)
    rownames(cv.W) <- rownames(cv.S) <- vapply(levels, format, "", nsmall = 2)

    for (j in seq_len(n.tau)) {
      if (anyNA(supW[, j, target])) next
      cv.W[, j] <- quantile(supW[, j, target], probs = 1 - levels, names = FALSE)
      cv.S[, j] <- quantile(supS[, j, target], probs = 1 - levels, names = FALSE)
    }

    list(W = cv.W, S = cv.S)
  }

  cv.post <- .cvals(1L)
  cv.comb <- .cvals(2L)
  metadata <- list(trim = trim, tau.step = tau.step, theta.max = theta.max,
                   theta.step = theta.step, n.grid = n.grid, n.rep = n.rep,
                   levels = levels, seed = seed, rng.kind = RNGkind(),
                   search.endpoints = range(t.val), boundary.method = "exact-process",
                   simulation.passes = 2L, quantile.type = 7L)

  list(
    tau = tau.all,
    theta = theta.all,
    postbreak = list(
      tau = tau.all,
      zeta.sqrt = zeta.post,
      W = cv.post$W,
      S = cv.post$S,
      metadata = c(metadata, list(target = "postbreak"))
    ),
    combined = list(
      tau = tau.all,
      zeta.sqrt = zeta.comb,
      W = cv.comb$W,
      S = cv.comb$S,
      metadata = c(metadata, list(target = "combined"))
    ),
    risk = list(full = r.full, postbreak = r.post, combined = r.comb),
    n.rep = n.rep,
    metadata = metadata
  )
}


# One block of the original Brownian-motion experiment. No worker draws RNG:
# its columns of innovations were supplied in replication order by the master.
.fcst.BP.chunk <- function(task, context) {
  t.idx <- context$t.idx
  t.val <- context$t.val
  tau.all <- context$tau.all
  theta.all <- context$theta.all
  n.tau <- length(tau.all)
  n.theta <- length(theta.all)
  size <- ncol(task$innovations)
  r.full <- r.post <- r.comb <- matrix(0, n.tau, n.theta)
  for (m in seq_len(size)) {
    bm <- cumsum(task$innovations[, m]) / sqrt(context$n.grid)
    bridge <- (bm[t.idx] - t.val * bm[context$n.grid]) * context$J1a
    for (i.th in seq_len(n.theta)) {
      for (i.tb in seq_len(n.tau)) {
        q <- (bridge + context$mu.arr[, i.tb, i.th])^2
        i.loc <- which.max(q)
        W <- q[i.loc]
        loc <- t.val[i.loc]
        tb <- tau.all[i.tb]
        th <- theta.all[i.th]
        e.full <- bm[context$n.grid] + tb * th
        e.post <- (bm[context$n.grid] - bm[t.idx[i.loc]] +
          th * (tb - loc) * (loc < tb)) / (1 - loc)
        e.comb <- (e.full + W * e.post) / (1 + W)
        r.full[i.tb, i.th] <- r.full[i.tb, i.th] + e.full^2
        r.post[i.tb, i.th] <- r.post[i.tb, i.th] + e.post^2
        r.comb[i.tb, i.th] <- r.comb[i.tb, i.th] + e.comb^2
      }
    }
  }
  list(first = task$first, risk = list(full = r.full, postbreak = r.post, combined = r.comb))
}


# The second pass has no RNG or plug-in rounding on the workers. Incomplete
# boundaries leave that target uncalibrated instead of interpolating over NA.
.fcst.BP.boundary.chunk <- function(task, context) {
  n.tau <- length(context$tau)
  size <- ncol(task$innovations)
  W <- S <- array(NA_real_, c(size, n.tau, 2L))
  for (m in seq_len(size)) {
    bm <- cumsum(task$innovations[, m]) / sqrt(context$n.grid)
    bridge <- (bm[context$t.idx] - context$t.val * bm[context$n.grid]) * context$J1a
    for (target in 1:2) {
      if (anyNA(context$zeta[, target])) next
      for (j in seq_len(n.tau)) {
        q <- (bridge + context$drift[, j, target])^2
        location <- which.max(q)
        W[m, j, target] <- q[location]
        S[m, j, target] <- sqrt(q[location]) - .grid.interp(context$zeta[, target],
          context$tau, context$t.val[location])
      }
    }
  }
  list(first = task$first, W = W, S = S)
}

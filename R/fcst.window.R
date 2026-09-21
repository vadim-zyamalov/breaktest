#' @title
#' Selection of the estimation window for end-of-sample forecasting
#' @name fcst.window
#'
#' @description
#' A family of routines choosing the start date of the estimation window
#' \eqn{[h, N]} used to produce a forecast for period \eqn{N + 1}.
#'
#' The methods address the bias-variance cost of retaining pre-break data.
#' Under the source-model assumptions, a longer window can reduce estimation
#' variance while introducing break-related bias (Pesaran and Timmermann, 2007).
#' Neither unbiasedness nor an interior optimal window is universal.
#'
#' The following methods are available:
#' * `ls`: least squares (quasi-ML) estimator of the break date,
#' * `bp`: last break date of the Bai-Perron sequential procedure with the
#'   number of breaks selected by BIC,
#' * `tradeoff`: plug-in minimiser of the asymptotic risk,
#' * `cv`: cross-validation of Pesaran and Timmermann (2007),
#' * `cv.pre`: the same, restricted to windows starting no later than the
#'   first post-break observation,
#' * `cvl`: Laplace (quasi-Bayesian) cross-validation of Hirano and Wright
#'   (2022),
#' * `cvl.pre`: the same, restricted as in `cv.pre`,
#' * `ijr`: optimal window of Inoue, Jin and Rossi (2017).
#'
#' @references
#' Pesaran, M. Hashem, and Allan Timmermann.
#' “Selection of Estimation Window in the Presence of Breaks.”
#' Journal of Econometrics 137, no. 1 (2007): 134–61.
#'
#' Hirano, Keisuke, and Jonathan H. Wright.
#' “Analyzing Cross-Validation for Forecasting with Structural Instability.”
#' Journal of Econometrics 226 (2022): 139-154.
#' https://doi.org/10.1016/j.jeconom.2020.10.009.
#'
#' Inoue, Atsushi, Lu Jin, and Barbara Rossi.
#' “Rolling Window Selection for Out-of-Sample Forecasting
#' with Time-Varying Parameters.”
#' Journal of Econometrics 196, no. 1 (2017): 55–67.
#'
#' @keywords internal
NULL


#' @title
#' Cross-validation criterion over estimation window start dates
#'
#' @description
#' \deqn{C(h) = \sum_{t = \omega}^{N - 1}
#' \left(y_{t+1} - x_{t+1}' \hat\beta_{h:t}\right)^2 .}
#'
#' A single implementation shared by the plain cross-validation, its
#' pre-break-restricted version and both Laplace variants. In the replication
#' codes this loop is written out four times (`empirfunc.m`) and once more in a
#' slightly different parameterisation (`PTcv.m`).
#'
#' @param y A dependent variable.
#' @param x Explanatory variables.
#' @param h.max The largest admissible window start date.
#' @param omega The first observation of the pseudo out-of-sample evaluation
#' period. Defaults to `round(0.9 * N)`, the setting of Hirano and Wright
#' (2022). Pesaran and Timmermann (2007) reserve 25% of the sample, which
#' corresponds to `round(0.75 * N)`.
#'
#' @return A vector of length `h.max` with the criterion values. Window start
#' dates leaving too few observations for estimation get `Inf`.
#'
#' @keywords internal
.fcst.cv.criterion <- function(y, x, h.max, omega = NULL) {
  N <- nrow(y)
  k <- ncol(x)

  if (is.null(omega)) {
    omega <- round(0.9 * N)
  }
  if (!is.numeric(omega) || length(omega) != 1L || !is.finite(omega) ||
      omega != floor(omega) || omega <= k || omega >= N) {
    stop("ERROR! fcst.window: omega leaves no validation observations")
  }
  if (!is.numeric(h.max) || length(h.max) != 1L || !is.finite(h.max) ||
      h.max != floor(h.max) || h.max < 1 || h.max > N) {
    stop("ERROR! fcst.window: invalid largest window start")
  }

  crit <- rep(Inf, h.max)

  for (h in seq_len(h.max)) {
    if (omega - h + 1 <= k) {
      next
    }

    acc <- 0
    ok <- TRUE

    for (j in omega:(N - 1)) {
      xj <- x[h:j, , drop = FALSE]
      yj <- y[h:j, , drop = FALSE]

      xx <- crossprod(xj)
      if (rcond(xx) < .Machine$double.eps) {
        ok <- FALSE
        break
      }

      bh <- OLS.reg(yj, xj)$coefficients
      acc <- acc + (y[j + 1, 1] - drop(x[j + 1, , drop = FALSE] %*% bh))^2
    }

    if (ok) {
      crit[h] <- acc
    }
  }

  crit
}


#' @rdname fcst.window
#' @order 1
#'
#' @param y A dependent variable.
#' @param x Explanatory variables. If `NULL` a single constant is used.
#' @param trim A trimming fraction bounding the admissible break dates.
#' @param time.index Original increasing integer observation positions. The
#' retained sample must be equally spaced without internal missing periods.
#'
#' @return `fcst.window.LS` returns a list with the estimated break date `bp`,
#' the implied window start `start` and `Wald`, in chi-square units, of the
#' homoskedastic sup-Wald test. `F.stat` retains that same value for backward
#' compatibility; `F` is the conventional statistic divided by k.
#'
#' @details
#' `fcst.window.LS` minimises the SSR of the two-regime regression over the
#' admissible break dates. Under a local-to-zero break the estimator is not
#' consistent. Its local asymptotic distribution depends on break magnitude;
#' weak changes need not locate the break precisely (Hirano and Wright, 2022).
fcst.window.LS <- function(y, x = NULL, trim = 0.15, time.index = NULL) {
  time.index <- .fcst.time.index(time.index, NROW(y))
  .d <- .fcst.data(y, x)
  .fcst.require.regular(time.index[.d$rows], "fcst.window.LS")
  y <- .d$y
  x <- .d$x
  N <- .d$N
  k <- .d$k
  .fcst.window.trim(trim)

  h.min <- max(round(trim * N), k + 1)
  h.max <- min(round((1 - trim) * N), N - k - 1)

  if (h.min > h.max) {
    stop("ERROR! fcst.window.LS: trimming leaves no admissible break dates")
  }

  SSR.data <- .fcst.ssr.matrix(y, x, k + 1)

  # segments.CSS minimises SSR.data[beg, bp] + SSR.data[bp + 1, end] over bp,
  # which is exactly the concentrated two-regime SSR
  .seg <- segments.CSS(1, N, h.min, h.max, N, SSR.data)

  bp <- .seg$break.point
  ssr.u <- .seg$SSR
  ssr.r <- sum(OLS.reg(y, x)$residuals^2, na.rm = TRUE)
  if (!is.finite(ssr.u)) stop("ERROR! fcst.window.LS: no identified split regression")
  statistic <- if (ssr.u == 0) {
    if (ssr.r == 0) 0 else Inf
  } else max(0, (ssr.r - ssr.u) / (ssr.u / (N - 2 * k)))

  list(
    bp = bp,
    start = bp + 1,
    F.stat = statistic,
    Wald = statistic,
    F = statistic / k,
    SSR = ssr.u,
    SSR.data = SSR.data
  )
}


#' @rdname fcst.window
#' @order 2
#'
#' @param max.breaks The largest number of breaks considered.
#' @param criterion An information criterion name accepted by
#' [info.criterions] (`bic`, `aic`, `hq`, `lwz`). NULL uses the `2 * k` penalty
#' per break of Hirano and Wright's `empirfunc.m`. Otherwise the segmented
#' model is refitted and scored by the package-wide
#' criterion, which counts all \eqn{(m + 1) k} coefficients rather than
#' applying an ad hoc per-break penalty.
#'
#' @return `fcst.window.BP` returns a list with the selected number of breaks
#' `n.breaks`, the full vector of break dates `bp.all`, the last break date
#' `bp` and the implied window start `start`.
#'
#' @details
#' `fcst.window.BP` calls the package-wide [segments.BP] implementation of the
#' Bai-Perron sequential procedure, so that no separate dynamic-programming
#' routine is needed. The number of breaks is chosen by BIC and the window
#' starts right after the last detected break.
fcst.window.BP <- function(
  y,
  x = NULL,
  trim = 0.15,
  max.breaks = 5,
  criterion = NULL,
  time.index = NULL
) {
  time.index <- .fcst.time.index(time.index, NROW(y))
  .d <- .fcst.data(y, x)
  .fcst.require.regular(time.index[.d$rows], "fcst.window.BP")
  y <- .d$y
  x <- .d$x
  N <- .d$N
  k <- .d$k
  .fcst.window.trim(trim)
  if (!is.numeric(max.breaks) || length(max.breaks) != 1L ||
      !is.finite(max.breaks) || max.breaks < 0 || max.breaks != floor(max.breaks)) {
    stop("ERROR! fcst.window.BP: max.breaks must be a nonnegative integer")
  }

  width <- max(round(trim * N), k + 1)
  max.breaks <- max(0, min(max.breaks, floor(N / width) - 1))
  if (!is.null(criterion)) {
    if (!is.character(criterion) || length(criterion) != 1L || is.na(criterion) ||
        !tolower(criterion) %in% c("bic", "aic", "hq", "hqic", "lwz")) {
      stop("ERROR! fcst.window.BP: unknown information criterion")
    }
    criterion <- tolower(criterion)
  }
  # The core dispatch calls Hannan-Quinn "HQIC"; accept the forecast API's
  # historical "hq" spelling without modifying the shared implementation.
  core.criterion <- if (identical(criterion, "hq")) "hqic" else criterion

  m.full <- OLS.reg(y, x)
  ssr.r <- sum(m.full$residuals^2, na.rm = TRUE)
  SSR.data <- if (max.breaks > 0) .fcst.ssr.matrix(y, x, width) else NULL

  # Design matrix of the model segmented at an arbitrary set of break dates
  .segmented <- function(bps) {
    edges <- c(0, sort(bps), N)
    do.call(cbind, lapply(seq_along(edges[-1]), function(i) {
      d <- numeric(N)
      d[(edges[i] + 1):edges[i + 1]] <- 1
      x * d
    }))
  }

  bic <- if (is.null(criterion)) {
    log(ssr.r / N)
  } else {
    info.criterions(m.full, core.criterion)
  }

  segs <- vector("list", max.breaks)

  for (m in seq_len(max.breaks)) {
    segs[[m]] <- tryCatch(
      segments.BP(y, x, m, width, SSR.data),
      error = function(e) NULL
    )

    bic <- c(
      bic,
      if (is.null(segs[[m]])) {
        Inf
      } else if (is.null(criterion)) {
        log(segs[[m]]$SSR / N) + 2 * k * log(N) / N * m
      } else {
        tryCatch(
          info.criterions(
            OLS.reg(y, .segmented(segs[[m]]$break.point)),
            core.criterion
          ),
          error = function(e) Inf
        )
      }
    )
  }

  n.breaks <- which.min(bic) - 1

  if (n.breaks == 0) {
    return(
      list(n.breaks = 0, bp.all = integer(0), bp = 0, start = 1, BIC = bic,
           IC = bic, criterion = if (is.null(criterion)) "HW-BIC" else criterion)
    )
  }

  bp.all <- sort(segs[[n.breaks]]$break.point)

  list(
    n.breaks = n.breaks,
    bp.all = bp.all,
    bp = max(bp.all),
    start = max(bp.all) + 1,
    BIC = bic,
    IC = bic,
    criterion = if (is.null(criterion)) "HW-BIC" else criterion
  )
}


#' @title
#' Asymptotic forecasting risk under a single break
#'
#' @description
#' \deqn{R(\eta) = \mu' \mu \left[\max\left(0,
#' \frac{c - \eta}{1 - \eta}\right)\right]^2 + \frac{K}{1 - \eta},}
#' the local asymptotic risk of the rolling window forecast when the model has a
#' single break at fraction `c` of standardised magnitude \eqn{\mu}
#' (Hirano and Wright, 2022, eq. A.2). Port of `func1.m`.
#'
#' The first term is the squared bias from including pre-break data, the second
#' the estimation variance from a short window.
#'
#' @param eta A vector of candidate window start fractions.
#' @param c A break fraction.
#' @param mu A vector of standardised break magnitudes.
#'
#' @return A vector of risk values.
#'
#' @keywords internal
.fcst.risk.single.break <- function(eta, c, mu) {
  K <- length(mu)
  g <- (c - eta) / (1 - eta)
  drop(crossprod(mu)) * pmax(0, g)^2 + K / (1 - eta)
}


#' @rdname fcst.window
#' @order 3
#'
#' @param bp A break date. If `NULL` it is estimated by [fcst.window.LS].
#' @param n.grid Retained for compatibility. The risk is now evaluated at every
#' admissible integer window start, including full and post-break windows.
#'
#' @return `fcst.window.tradeoff` returns a list with the window start `start`,
#' the estimated standardised break magnitude `mu` and the risk profile.
#'
#' @details
#' `fcst.window.tradeoff` plugs the least squares estimates of the break date
#' and magnitude into the asymptotic risk and minimises it. This is the
#' finite-sample "trade-off" method of Pesaran and Timmermann (2007) in the
#' local asymptotic parameterisation of Hirano and Wright (2022).
#'
#' The estimated magnitude is
#' \eqn{\hat\mu = \sqrt{N} \, \Sigma_{xx}^{1/2}
#' (\hat\beta_{pre} - \hat\beta_{post}) / \hat\sigma}.
fcst.window.tradeoff <- function(
  y,
  x = NULL,
  trim = 0.15,
  bp = NULL,
  n.grid = 1000,
  time.index = NULL
) {
  time.index <- .fcst.time.index(time.index, NROW(y))
  .d <- .fcst.data(y, x)
  .fcst.require.regular(time.index[.d$rows], "fcst.window.tradeoff")
  y <- .d$y
  x <- .d$x
  N <- .d$N
  k <- .d$k
  .fcst.window.trim(trim)
  if (!is.numeric(n.grid) || length(n.grid) != 1L || !is.finite(n.grid) ||
      n.grid < 1 || n.grid != floor(n.grid)) {
    stop("ERROR! fcst.window.tradeoff: n.grid must be a positive integer")
  }

  if (is.null(bp)) {
    bp <- fcst.window.LS(y, x, trim)$bp
  }

  if (!is.numeric(bp) || length(bp) != 1L || !is.finite(bp) || bp != floor(bp) ||
      bp <= k || N - bp <= k) {
    stop("ERROR! fcst.window.tradeoff: both regimes need more than k observations")
  }

  m.pre <- OLS.reg(y[1:bp, , drop = FALSE], x[1:bp, , drop = FALSE])
  m.post <- OLS.reg(y[(bp + 1):N, , drop = FALSE], x[(bp + 1):N, , drop = FALSE])

  resid <- c(m.pre$residuals, m.post$residuals)
  sigma <- sqrt(mean(resid^2, na.rm = TRUE))
  if (!is.finite(sigma) || sigma <= 0) {
    stop("ERROR! fcst.window.tradeoff: positive residual variance required")
  }

  sigma.xx <- crossprod(x) / N
  ev <- eigen(sigma.xx, symmetric = TRUE)
  sqrt.xx <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), k) %*% t(ev$vectors)

  mu <- sqrt(N) * drop(sqrt.xx %*% (m.pre$coefficients - m.post$coefficients)) /
    sigma

  # eta is the fraction discarded, so a start at bp+1 has eta=bp/N.
  # Later starts cannot improve this risk, but are included up to the minimum
  # estimable window to make the complete admissible discrete domain explicit.
  starts <- seq_len(N - k)
  eta <- (starts - 1) / N
  risk <- .fcst.risk.single.break(eta, bp / N, mu)

  start <- starts[which.min(risk)]

  list(
    start = start,
    bp = bp,
    mu = mu,
    starts = starts,
    eta = eta,
    risk = risk
  )
}


#' @rdname fcst.window
#' @order 4
#'
#' @param omega The first observation of the pseudo out-of-sample evaluation
#' period, passed to [.fcst.cv.criterion].
#' @param restrict Whether to restrict the search to window start dates no
#' later than the first post-break observation (the `cv.pre` variant), as in
#' Pesaran and Timmermann's known-break candidate set.
#'
#' @return `fcst.window.cv` returns a list with the window start `start` and
#' the criterion profile `crit`.
#'
#' @details
#' `fcst.window.cv` reserves the tail of the sample and picks the window start
#' minimising the pseudo out-of-sample MSFE. In the simulations of Pesaran and
#' Timmermann (2007) this is one of the most reliable methods: when the break is
#' small, attempting to estimate its date is counterproductive.
fcst.window.cv <- function(
  y,
  x = NULL,
  trim = 0.15,
  omega = NULL,
  restrict = FALSE,
  bp = NULL,
  time.index = NULL
) {
  time.index <- .fcst.time.index(time.index, NROW(y))
  .d <- .fcst.data(y, x)
  .fcst.require.regular(time.index[.d$rows], "fcst.window.cv")
  y <- .d$y
  x <- .d$x
  N <- .d$N
  .fcst.window.trim(trim)
  if (!is.logical(restrict) || length(restrict) != 1L || is.na(restrict)) {
    stop("ERROR! fcst.window.cv: restrict must be TRUE or FALSE")
  }

  h.max <- round((1 - trim) * N)

  if (restrict) {
    if (is.null(bp)) {
      bp <- fcst.window.LS(y, x, trim)$bp
    }
    if (!is.numeric(bp) || length(bp) != 1L || !is.finite(bp) ||
        bp != floor(bp) || bp < 1 || bp >= N) stop("ERROR! fcst.window.cv: invalid break date")
    h.max <- min(h.max, bp + 1L)
  }

  crit <- .fcst.cv.criterion(y, x, h.max, omega)
  if (!any(is.finite(crit))) {
    stop("ERROR! fcst.window.cv: no admissible CV window")
  }

  list(
    start = which.min(crit),
    crit = crit,
    bp = bp
  )
}


#' @rdname fcst.window
#' @order 5
#'
#' @return `fcst.window.cvl` returns a list with the window start `start` and
#' the pseudo-posterior weights `posterior`.
#'
#' @details
#' `fcst.window.cvl` treats \eqn{\exp(-C(h) / 2\hat\sigma^2)} as a
#' pseudo-likelihood in the window start date and reports its posterior mean
#' under a flat prior:
#' \deqn{\hat h = \frac{\int h \, L(h) \, dh}{\int L(h) \, dh}.}
#'
#' The motivation is that the cross-validation criterion is noisy and the
#' underlying risk function is strongly asymmetric around its minimiser, so
#' taking the argmin discards information. Hirano and Wright (2022) find that
#' this dominates the plain cross-validation over most of their designs and
#' removes the pile-up at zero for exponential weighting.
fcst.window.cvl <- function(
  y,
  x = NULL,
  trim = 0.15,
  omega = NULL,
  restrict = FALSE,
  bp = NULL,
  time.index = NULL
) {
  time.index <- .fcst.time.index(time.index, NROW(y))
  .d <- .fcst.data(y, x)
  .fcst.require.regular(time.index[.d$rows], "fcst.window.cvl")
  y <- .d$y
  x <- .d$x
  N <- .d$N
  .fcst.window.trim(trim)
  if (!is.logical(restrict) || length(restrict) != 1L || is.na(restrict)) {
    stop("ERROR! fcst.window.cvl: restrict must be TRUE or FALSE")
  }

  h.max <- round((1 - trim) * N)
  crit <- .fcst.cv.criterion(y, x, h.max, omega)

  s.sq <- mean(OLS.reg(y, x)$residuals^2, na.rm = TRUE)

  if (restrict) {
    if (is.null(bp)) {
      bp <- fcst.window.LS(y, x, trim)$bp
    }
    if (!is.numeric(bp) || length(bp) != 1L || !is.finite(bp) ||
        bp != floor(bp) || bp < 1 || bp >= N) stop("ERROR! fcst.window.cvl: invalid break date")
    crit <- crit[seq_len(min(bp + 1L, length(crit)))]
  }
  if (!any(is.finite(crit))) {
    stop("ERROR! fcst.window.cvl: no admissible CV window")
  }

  # At zero scale the limiting posterior shares mass over exact minima.
  # Subtract before division to avoid Inf-Inf at very small positive scales.
  delta <- crit - min(crit[is.finite(crit)])
  post <- if (s.sq == 0) as.numeric(delta == 0) else exp(-0.5 * (delta / s.sq))
  post[!is.finite(post)] <- 0
  post <- post / sum(post)
  start <- max(1L, round(sum(post * seq_along(post))))
  if (!is.finite(crit[start])) {
    stop("ERROR! fcst.window.cvl: posterior mean selects an unidentified window")
  }

  list(
    start = start,
    posterior = post,
    scale = s.sq,
    degenerate.scale = s.sq == 0,
    crit = crit,
    bp = bp
  )
}


#' @title
#' Local linear estimator of the end-of-sample coefficient vector
#'
#' @description
#' Fits \eqn{y_t = x_t' \beta(1) + x_t' \beta^{(1)}(1) (t - T_f) / N +
#' \varepsilon_t} over the last `R0` observations and returns
#' \eqn{\tilde\beta(1)}.
#'
#' Inoue et al. (2017) use this pilot estimator in place of the unknown
#' \eqn{\beta(1)} when selecting the rolling window: its boundary bias vanishes
#' faster than that of the local constant (rolling) estimator, so the
#' substitution is asymptotically harmless.
#'
#' @param y A dependent variable.
#' @param x Explanatory variables.
#' @param R0 A pilot window size, at least 2*k+1.
#' @param time.index Original increasing integer time positions.
#' @param forecast.time Prediction endpoint in the same clock as `time.index`.
#' Defaults to the next period. In a direct h-step regression using response
#' dates as its clock, the temporal dispatcher supplies forecast origin + h.
#'
#' @return A vector of length `ncol(x)`.
#'
#' @keywords internal
.fcst.local.linear <- function(y, x, R0, time.index = NULL, forecast.time = NULL) {
  N <- nrow(y)
  k <- ncol(x)

  time.index <- .fcst.time.index(time.index, N)
  .fcst.require.regular(time.index, ".fcst.local.linear")
  if (is.null(forecast.time)) forecast.time <- max(time.index) + 1L
  if (!is.numeric(forecast.time) || length(forecast.time) != 1L ||
      !is.finite(forecast.time) || forecast.time != floor(forecast.time) ||
      forecast.time <= max(time.index)) {
    stop("ERROR! .fcst.local.linear: forecast.time must be an integer after the training sample")
  }
  if (!is.numeric(R0) || length(R0) != 1L || !is.finite(R0) ||
      R0 != floor(R0) || R0 < 2 * k + 1 || R0 > N) {
    stop("ERROR! .fcst.local.linear: R0 must be an integer in [2*k+1, N]")
  }
  idx <- (N - R0 + 1):N

  v <- (time.index[idx] - forecast.time) / (max(time.index) - min(time.index) + 1)
  xl <- cbind(x[idx, , drop = FALSE], x[idx, , drop = FALSE] * v)
  if (qr(xl)$rank < 2 * k) stop("ERROR! .fcst.local.linear: singular local-linear pilot design")

  .m <- OLS.reg(y[idx, , drop = FALSE], xl)
  .m$coefficients[seq_len(k)]
}


#' @rdname fcst.window
#' @order 6
#'
#' @param R0 A pilot window for the local linear estimator. If `NULL` it is set
#' by unrestricted Pesaran-Timmermann cross-validation after rejection of
#' parameter constancy. At least 2*k+1 observations are needed by the pilot.
#' @param R.range A two-element vector with the smallest and the largest
#' admissible window length. Defaults to
#' \eqn{[\max(1.5 N^{2/3}, 20), \min(c \, N^{2/3}, N)]} with `c.max` below.
#' @param c.max A multiplier of \eqn{N^{2/3}} bounding the window from above.
#' Inoue et al. (2017) report results for 4, 5 and 6.
#' @param cv.pretest Optional explicit critical value for the sup-Wald
#' pretest. NULL uses [get.cv.Wald] with the number of tested coefficients
#' and effective trimming; a supplied value overrides the calibrated cutoff.
#' @param level Nominal pretest size, default 5% as in Inoue et al. (2017).
#' @param vcov.type,bandwidth Settings of [fcst.test.Wald]. HAC is the default
#' to allow weakly dependent errors; the retained sample must remain regular.
#' @param pretest Whether to apply the article's stability pretest. FALSE
#' requests only the conditional pilot-based window rule.
#' @param omega First validation origin for pilot CV. NULL uses 75% of the
#' sample, following Pesaran and Timmermann (2007); the minimum window fraction
#' for this pilot is 10%. It does not condition on an estimated break date.
#' @param pilot.fun Internal zero-argument factory returning pilot CV details.
#' The temporal dispatcher supplies it to evaluate the requested horizon
#' using only responses observed at each validation origin.
#' @param forecast.time Prediction endpoint in the clock of `time.index`.
#' Defaults to the period after the original input sample, before removing
#' incomplete rows. For a direct h-step regression using response dates,
#' pass the actual forecast origin + h.
#'
#' @return `fcst.window.IJR` returns a list with the selected window length `R`,
#' the implied start date `start` and the criterion profile.
#'
#' @details
#' `fcst.window.IJR` selects the rolling window minimising the *conditional*
#' MSFE at the end of the sample,
#' \deqn{\hat R = \arg\min_R \left(\hat\beta_R(1) - \tilde\beta(1)\right)'
#' x_N x_N' \left(\hat\beta_R(1) - \tilde\beta(1)\right),}
#' with the pilot \eqn{\tilde\beta(1)} coming from [.fcst.local.linear]. The
#' optimal window is of order \eqn{N^{2/3}}, which is what bounds the search
#' range.
#'
#' Under the smoothness and dependence assumptions of the paper, the
#' conditional rule allows smoothly time-varying parameters, weakly dependent
#' errors, lagged dependent variables and direct forecast regressions. The
#' implementation follows Section 2: failure to reject stability returns the
#' full sample as the final window, without imposing the N^(2/3) search cap.
fcst.window.IJR <- function(
  y,
  x = NULL,
  x.new = NULL,
  trim = 0.15,
  R0 = NULL,
  R.range = NULL,
  c.max = 4,
  cv.pretest = NULL,
  level = 0.05,
  vcov.type = "HAC",
  bandwidth = NULL,
  time.index = NULL,
  omega = NULL,
  pretest = TRUE,
  pilot.fun = NULL,
  forecast.time = NULL
) {
  time.index <- .fcst.time.index(time.index, NROW(y))
  if (is.null(forecast.time)) forecast.time <- max(time.index) + 1L
  if (!is.numeric(forecast.time) || length(forecast.time) != 1L ||
      !is.finite(forecast.time) || forecast.time != floor(forecast.time) ||
      forecast.time <= max(time.index)) {
    stop("ERROR! fcst.window.IJR: forecast.time must be an integer after the input sample")
  }
  .d <- .fcst.data(y, x, x.new)
  y <- .d$y
  x <- .d$x
  x.new <- .d$x.new
  N <- .d$N
  k <- .d$k
  time.index <- time.index[.d$rows]
  .fcst.require.regular(time.index, "fcst.window.IJR")
  .fcst.window.trim(trim)
  .fcst.test.level(level)
  if (!is.logical(pretest) || length(pretest) != 1L || is.na(pretest)) {
    stop("ERROR! fcst.window.IJR: pretest must be TRUE or FALSE")
  }
  if (!is.null(pilot.fun) && !is.function(pilot.fun)) {
    stop("ERROR! fcst.window.IJR: pilot.fun must be a function")
  }
  if (!is.null(cv.pretest) && (!is.numeric(cv.pretest) || length(cv.pretest) != 1L ||
      !is.finite(cv.pretest) || cv.pretest < 0)) {
    stop("ERROR! fcst.window.IJR: cv.pretest must be finite and nonnegative")
  }
  if (!is.numeric(c.max) || length(c.max) != 1L || !is.finite(c.max) || c.max <= 0) {
    stop("ERROR! fcst.window.IJR: c.max must be finite and positive")
  }
  if (!is.null(R0) && (!is.numeric(R0) || length(R0) != 1L || !is.finite(R0) ||
      R0 != floor(R0) || R0 < 2 * k + 1 || R0 > N)) {
    stop("ERROR! fcst.window.IJR: R0 must be an integer in [2*k+1, N]")
  }
  test <- NULL
  if (pretest) {
    test <- fcst.test.Wald(y, x, trim = trim, vcov.type = vcov.type,
                          bandwidth = bandwidth, time.index = time.index, level = level)
    cutoff <- if (is.null(cv.pretest)) unname(test$critical.value["sup"]) else cv.pretest
    reject <- unname(test$statistics["sup"]) > cutoff
    test$critical.value["sup"] <- cutoff
    test$reject["sup"] <- reject
    if (!reject) {
      return(list(R = N, start = 1L, R0 = N, beta.pilot = OLS.reg(y, x)$coefficients,
                  R.all = N, crit = 0, pretest = test, full.only = TRUE,
                  pilot.cv = NULL, pilot.adjusted = FALSE, forecast.time = forecast.time))
    }
  }

  pilot.cv <- NULL
  pilot.adjusted <- FALSE
  if (is.null(R0)) {
    if (is.null(omega)) omega <- max(k + 1L, round(0.75 * N))
    # The temporal dispatcher supplies a forecast-origin-aware CV factory.
    # It is called only after rejection; direct h>1 must not validate rows
    # whose training responses would not yet have been observed.
    pilot.cv <- if (is.null(pilot.fun)) {
      fcst.window.cv(y, x, trim = 0.1, omega = omega, restrict = FALSE)
    } else pilot.fun()
    if (!is.list(pilot.cv) || length(pilot.cv$start) != 1L ||
        !is.finite(pilot.cv$start) || pilot.cv$start != floor(pilot.cv$start) ||
        pilot.cv$start < 1 || pilot.cv$start > N) {
      stop("ERROR! fcst.window.IJR: pilot CV returned an invalid start")
    }
    R0 <- N - pilot.cv$start + 1L
    pilot.adjusted <- R0 < 2 * k + 1L
    R0 <- max(R0, 2 * k + 1L)
  }

  beta.ll <- .fcst.local.linear(y, x, R0, time.index, forecast.time)

  if (is.null(R.range)) {
    R.range <- c(
      round(max(1.5 * N^(2 / 3), 20)),
      round(min(c.max * N^(2 / 3), N))
    )
  }
  if (!is.numeric(R.range) || length(R.range) != 2L || any(!is.finite(R.range)) ||
      any(R.range != floor(R.range)) || R.range[1] < 1 || R.range[1] > R.range[2]) {
    stop("ERROR! fcst.window.IJR: R.range must contain two ordered positive integers")
  }
  R.range[1] <- max(R.range[1], k + 1)
  R.range[2] <- min(R.range[2], N)

  if (R.range[1] > R.range[2]) {
    stop("ERROR! fcst.window.IJR: empty window search range")
  }

  R.all <- seq.int(R.range[1], R.range[2])
  crit <- sapply(R.all, function(R) {
    idx <- (N - R + 1):N
    if (qr(x[idx, , drop = FALSE])$rank < k) return(Inf)
    .m <- OLS.reg(y[idx, , drop = FALSE], x[idx, , drop = FALSE])
    drop(x.new %*% (.m$coefficients - beta.ll))^2
  })
  if (!any(is.finite(crit))) stop("ERROR! fcst.window.IJR: no identified candidate window")

  R <- R.all[which.min(crit)]

  list(
    R = R,
    start = N - R + 1,
    R0 = R0,
    forecast.time = forecast.time,
    beta.pilot = beta.ll,
    R.all = R.all,
    crit = crit,
    pretest = test,
    full.only = FALSE,
    pilot.cv = pilot.cv,
    pilot.adjusted = pilot.adjusted
  )
}


#' @rdname fcst.window
#' @order 7
#'
#' @param method One of `ls`, `bp`, `tradeoff`, `cv`, `cv.pre`, `cvl`,
#' `cvl.pre`, `ijr`.
#' @param ... Extra arguments passed to the underlying routine.
#'
#' @return An object of class `bt_fcstWindow`: a list with the window start
#' `start`, the window length `R`, the `method` used and the full output of the
#' underlying routine in `details`.
#'
#' @details
#' `fcst.window` is the common entry point. It reproduces, in a single call,
#' the window rules studied by Hirano and Wright (2022), with the IJR pretest
#' and pilot following the original Inoue et al. (2017) article.
#'
#' @export
fcst.window <- function(
  y,
  x = NULL,
  method = c("cv", "ls", "bp", "tradeoff", "cv.pre", "cvl", "cvl.pre", "ijr"),
  ...
) {
  method <- match.arg(method)

  .d <- .fcst.data(y, x)
  N <- .d$N

  res <- switch(
    method,
    ls = fcst.window.LS(y, x, ...),
    bp = fcst.window.BP(y, x, ...),
    tradeoff = fcst.window.tradeoff(y, x, ...),
    cv = fcst.window.cv(y, x, restrict = FALSE, ...),
    cv.pre = fcst.window.cv(y, x, restrict = TRUE, ...),
    cvl = fcst.window.cvl(y, x, restrict = FALSE, ...),
    cvl.pre = fcst.window.cvl(y, x, restrict = TRUE, ...),
    ijr = fcst.window.IJR(y, x, ...)
  )

  result <- list(
    start = res$start,
    R = N - res$start + 1,
    N = N,
    method = method,
    details = res
  )

  class(result) <- "bt_fcstWindow"
  result
}


# Shared validation of the trimming domain for window estimators.
.fcst.window.trim <- function(trim) {
  if (!is.numeric(trim) || length(trim) != 1L || !is.finite(trim) ||
      trim <= 0 || trim >= 0.5) {
    stop("ERROR! fcst.window: trim must be in (0, 0.5)")
  }
  invisible(trim)
}

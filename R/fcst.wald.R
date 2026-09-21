#' @title
#' Sup-Wald and exponential-Wald structural change tests
#'
#' @description
#' Tests a common unknown break in selected regression coefficients. Each
#' candidate date is fitted with [OLS.reg] and the module's covariance routines.
#' The statistic is in chi-square units, not divided by the number of tested
#' coefficients. Approximate p-values use Hansen's (1997) response surfaces.
#'
#' @param y,x An aligned regression; NULL x supplies a constant.
#' @param break.vars Original design columns allowed to change. NULL tests all.
#' @param trim,trim.exp Symmetric trimming fractions for sup and exp tests.
#' Dates strictly inside (trim*N, (1-trim)*N) are searched. The effective
#' endpoints are reported and used in the p-value approximation.
#' @param vcov.type "const", "HC0", "HC3" or "HAC". HC0 is the default.
#' @param bandwidth,time.index HAC bandwidth and original integer positions.
#' The estimation block must be contiguous; internal gaps are unsupported.
#' @param level Nominal rejection probability.
#'
#' @details
#' exp = log(mean(exp(W(date)/2))) is evaluated by log-sum-exp. The restricted
#' coefficients remain constant throughout the sample. All admissible dates
#' must be identified; singular dates are not silently removed from the test's
#' search region. bp maximizes Wald; bp.LS minimizes the unrestricted SSR.
#' They need not coincide with robust covariance. The sup and exp p-values
#' require the source asymptotic assumptions and are not exact finite-sample
#' probabilities under arbitrary nonstationarity or post-selection.
#'
#' @return A bt_fcstWald with named statistics, p.value, critical.value,
#' reject, estimated dates, effective trimming, and the complete profile.
#' `p.value.status` marks p-values outside the monotone response-surface
#' range. These are `NA`; decisions still compare the statistic to an
#' available critical value at the requested level.
#'
#' @references
#' Andrews (1993), Econometrica 61, 821-856, and Andrews-Ploberger (1994),
#' Econometrica 62, 1383-1414. Hansen (1997), "Approximate Asymptotic
#' P Values for Structural-Change Tests", JBES 15, 60-67.
#' https://users.ssc.wisc.edu/~behansen/progs/jbes_97.html.
#'
#' @export
fcst.test.Wald <- function(y, x = NULL, break.vars = NULL, trim = .15,
                           trim.exp = trim, vcov.type = "HC0", bandwidth = NULL,
                           time.index = NULL, level = .05) {
  time.index <- .fcst.time.index(time.index, NROW(y))
  d <- .fcst.data(y, x)
  .fcst.inference.data(d)
  cols <- .fcst.break.columns(break.vars, d$k)
  .fcst.wald.arguments(length(cols), trim)
  .fcst.test.level(level)
  time.index <- time.index[d$rows]
  .fcst.require.regular(time.index, "fcst.test.Wald")
  gs <- .fcst.wald.grid(d$N, d$k, trim, length(cols))
  ge <- .fcst.wald.grid(d$N, d$k, trim.exp, length(cols))
  dates <- sort(unique(c(gs, ge)))
  wald <- ssr <- numeric(length(dates))
  for (i in seq_along(dates)) {
    split <- .fcst.split.design(d$x, dates[i], cols)
    xx <- split$xs
    if (qr(xx)$rank < ncol(xx)) stop("ERROR! fcst.test.Wald: singular design inside search region")
    fit <- OLS.reg(d$y, xx)
    V <- .fcst.vcov(xx, fit$residuals, vcov.type, bandwidth, time.index)
    R <- matrix(0, length(cols), ncol(xx))
    R[cbind(seq_along(cols), split$idx.1)] <- 1
    R[cbind(seq_along(cols), split$idx.2)] <- -1
    difference <- R %*% fit$coefficients
    VD <- R %*% V %*% t(R)
    if (rcond(VD) < 1e-12 || any(!is.finite(VD)) ||
        min(eigen(VD, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
      stop("ERROR! fcst.test.Wald: singular coefficient covariance inside search region")
    }
    wald[i] <- drop(crossprod(difference, solve(VD, difference)))
    ssr[i] <- sum(fit$residuals^2)
  }
  sw <- wald[match(gs, dates)]
  ew <- wald[match(ge, dates)] / 2
  statistics <- c(sup = max(sw), exp = max(ew) + log(mean(exp(ew - max(ew)))))
  effective <- c(sup = .fcst.wald.trim(range(gs) / d$N),
                  exp = .fcst.wald.trim(range(ge) / d$N))
  p <- cv <- setNames(numeric(2), names(statistics))
  for (type in names(statistics)) {
    p[type] <- get.pv.Wald(statistics[type], length(cols), effective[type], type)
    cv[type] <- get.cv.Wald(length(cols), effective[type], level, type)
  }
  structure(list(statistics = statistics, p.value = p, critical.value = cv,
                  p.value.status = setNames(ifelse(is.na(p),
                    "outside monotone response-surface range", "Hansen approximation"), names(p)),
                  reject = statistics > cv, bp = gs[which.max(sw)],
                  bp.LS = gs[which.min(ssr[match(gs, dates)])],
                  profile = data.frame(bp = dates, Wald = wald, SSR = ssr,
                                       sup = dates %in% gs, exp = dates %in% ge),
                  trim = c(sup = trim, exp = trim.exp), effective.trim = effective,
                  level = level, df = length(cols), break.vars = cols,
                  vcov.type = vcov.type, bandwidth = bandwidth,
                  N = d$N, rows = d$rows, time.index = time.index), class = "bt_fcstWald")
}


#' @title
#' Hansen p-values and critical values for structural change tests
#' @param statistic A nonnegative sup-Wald or log-mean-exp(Wald/2) statistic.
#' @param df Number of tested coefficients, 1:40.
#' @param trim Symmetric trimming in `[.01,.5]`. For asymmetric endpoints
#' (a,b), use 1/(1+sqrt(b*(1-a)/(a*(1-b)))).
#' @param type "sup" or "exp".
#' @param level Upper-tail probability for a critical value.
#' @details Numerical response-surface tables are verified against both the
#' local GAUSS source and Hansen's author archive. Critical values invert the
#' same p-value function, avoiding typos in the local qa_crit/ap_crit tables.
#' At trim=.5 there is one date: sup is chi-square and exp is chi-square/2.
#' The latter factor corrects the degenerate-date branch in the old code.
#' Beyond the maximum of a relevant concave quadratic response surface,
#' the approximation would reverse the tail. P-values there are `NA`, not
#' a capped probability. Critical values are inverted only within the
#' monotone region and are `NA` if the requested tail cannot be bracketed.
#' @return A numeric p-value vector or scalar approximate critical value;
#' `NA` denotes a request outside the usable response-surface range.
#' @importFrom stats pchisq uniroot
#' @export
get.pv.Wald <- function(statistic, df = 1L, trim = .15, type = c("sup", "exp")) {
  type <- match.arg(type)
  .fcst.wald.arguments(df, trim)
  if (!is.numeric(statistic) || !length(statistic) || anyNA(statistic) || any(statistic < 0)) {
    stop("ERROR! get.pv.Wald: statistic must be nonnegative")
  }
  surface <- .fcst.wald.surface(df, trim, type)
  beta <- surface$beta
  single <- function(value) {
    if (is.infinite(value)) return(0)
    if (value == 0) return(1)
    limit <- pchisq(value * if (type == "exp") 2 else 1, df, lower.tail = FALSE)
    if (trim == .5) return(limit)
    if (value > surface$upper) return(NA_real_)
    m <- ncol(beta) - 1L
    polynomial <- beta[, m]
    for (j in seq.int(m - 1L, 1L)) polynomial <- polynomial * value + beta[, j]
    pp <- pchisq(pmax(polynomial, 0), beta[, m + 1L], lower.tail = FALSE)
    sum(surface$weights * pp) + surface$chi.weight * limit
  }
  vapply(statistic, single, 0.0, USE.NAMES = FALSE)
}


#' @rdname get.pv.Wald
#' @export
get.cv.Wald <- function(df = 1L, trim = .15, level = .05, type = c("sup", "exp")) {
  type <- match.arg(type)
  .fcst.wald.arguments(df, trim)
  .fcst.test.level(level)
  f <- function(q) get.pv.Wald(q, df, trim, type) - level
  bound <- .fcst.wald.surface(df, trim, type)$upper
  upper <- min(max(10, 2 * df), bound)
  while (is.finite(f(upper)) && f(upper) > 0 && upper < min(bound, 1e6)) {
    upper <- min(upper * 2, bound, 1e6)
  }
  if (!is.finite(f(upper)) || f(upper) > 0) return(NA_real_)
  uniroot(f, c(0, upper), tol = 1e-9)$root
}


# Select only rows contributing to the published linear interpolation.
# The upper bound excludes the declining side of a quadratic response surface;
# it is a numerical-domain bound, not an additional tail approximation.
.fcst.wald.surface <- function(df, trim, type) {
  beta <- if (type == "sup") .fcst.hansen.sup[[df]] else .fcst.hansen.exp[[df]]
  chi.weight <- 0
  if (trim >= .49) {
    rows <- 1L
    weights <- (.5 - trim) * 100
    chi.weight <- (trim - .49) * 100
  } else if (trim <= .01) {
    rows <- 25L
    weights <- 1
  } else {
    at <- (.51 - trim) * 50
    lower <- min(24L, floor(at))
    rows <- c(lower, lower + 1L)
    weights <- c(lower + 1 - at, at - lower)
  }
  active <- weights > 1e-12
  beta <- beta[rows[active], , drop = FALSE]
  weights <- weights[active]
  upper <- Inf
  if (ncol(beta) == 4L && nrow(beta)) {
    concave <- beta[, 3] < 0
    if (any(concave)) upper <- min(-beta[concave, 2] / (2 * beta[concave, 3]))
  }
  list(beta = beta, weights = weights, chi.weight = chi.weight, upper = upper)
}


.fcst.wald.arguments <- function(df, trim) {
  if (!is.numeric(df) || length(df) != 1L || !is.finite(df) ||
      df != floor(df) || df < 1 || df > 40) stop("ERROR! get.pv.Wald: df must be an integer in 1:40")
  if (!is.numeric(trim) || length(trim) != 1L || !is.finite(trim) || trim < .01 || trim > .5) {
    stop("ERROR! get.pv.Wald: trim must be in [.01,.5]")
  }
}


.fcst.wald.trim <- function(endpoints) {
  1 / (1 + sqrt(endpoints[2] * (1 - endpoints[1]) / (endpoints[1] * (1 - endpoints[2]))))
}


.fcst.wald.grid <- function(N, k, trim, n.break = k) {
  .fcst.wald.arguments(1, trim)
  # Stable coefficients are pooled: each side needs room for the changing
  # columns, not a separate estimate of every column. The expanded design
  # has k+n.break columns; its actual rank is checked at every candidate.
  if (N <= k + n.break) stop("ERROR! fcst.test.Wald: no residual degrees of freedom")
  first <- max(floor(trim * N) + 1L, n.break + 1L)
  last <- min(ceiling((1 - trim) * N) - 1L, N - n.break - 1L)
  if (first > last) stop("ERROR! fcst.test.Wald: trimming leaves no identified search dates")
  seq.int(first, last)
}


.fcst.test.level <- function(level) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    stop("ERROR! fcst.test: level must be in (0,1)")
  }
}


.fcst.break.columns <- function(cols, k, empty = FALSE) {
  if (is.null(cols)) return(seq_len(k))
  if (!is.numeric(cols) || (!empty && !length(cols)) || any(!is.finite(cols)) ||
      any(cols != floor(cols) | cols < 1 | cols > k) || anyDuplicated(cols)) {
    stop("ERROR! fcst.test: break.vars must index distinct design columns")
  }
  sort(as.integer(cols))
}


.fcst.inference.data <- function(d) {
  if (!is.numeric(d$y) || !is.numeric(d$x) || any(!is.finite(d$y)) ||
      any(!is.finite(d$x)) || any(!is.finite(d$x.new)) || nrow(d$x.new) != 1L ||
      qr(d$x)$rank < d$k) stop("ERROR! fcst.test: finite full-rank regression and one forecast row required")
}

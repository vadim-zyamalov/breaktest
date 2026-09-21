#' @title
#' Weighted GLS with post-break leave-one-out cross-validation
#'
#' @description
#' The feasible WGLS forecast of Lee, Parsaeian and Ullah (2022), equations
#' (11)-(13). The Aitchison-Aitken kernel discounts pre-break observations by
#' `gamma`; accounting for the regime variances gives the effective relative
#' weight \eqn{\gamma^* = \gamma / q^2}, where
#' \eqn{q = \sigma_{(1)} / \sigma_{(2)}} is a standard deviation ratio.
#'
#' @param y A dependent variable. Incomplete rows of `y` and `x` are removed
#' jointly. The retained observations must form a regular uninterrupted sample.
#' @param x Explanatory variables. If `NULL` a single constant is used.
#' @param x.new A (1 x k) vector of regressors for the forecast. Defaults to
#' the last complete row of `x`.
#' @param bp The last pre-break observation in the complete-case sample.
#' If `NULL` the single break is estimated by [fcst.window.LS]. Both regimes
#' must have full column rank and more observations than regressors.
#' @param trim A trimming fraction in (0, 0.5), used when `bp` is estimated.
#' @param time.index Original increasing integer observation positions.
#' @param gamma.grid Candidate kernel parameters in `[0, 1]`. By default the
#' grid is `seq(0, 1, by = 0.01)`, as in `ac/empirfunc.m`. Zero is included
#' as the continuous post-break-only limit of the paper's kernel.
#' @param gamma.star.grid Optional finite, nonnegative effective weights.
#' When supplied, implements the direct search in Remark 3 without estimating
#' `q`. Values above one are allowed. Do not also supply `gamma.grid`.
#'
#' @details
#' By default, regime OLS residuals give \eqn{\hat\sigma_j^2 = RSS_j / n_j}
#' and \eqn{\hat q^2 = \hat\sigma_1^2 / \hat\sigma_2^2}. These pilot estimates
#' and the break date are held fixed throughout cross-validation. Each
#' post-break row is omitted in turn from the weighted regression, while all
#' pre-break rows remain. The criterion returned is the mean squared LOO
#' error. Dividing it by the fixed post-break variance yields equation (13)
#' and leaves its minimiser unchanged. Direct search needs no variance scale.
#'
#' For fixed weights the LOO error is \eqn{e_s / (1-h_s)}, where
#' \eqn{h_s = w_s x_s'(X'WX)^{-1}x_s}. Thus only one [WLS.reg] fit per
#' candidate is needed. Candidates with singular or numerically unstable
#' fits, including any post-break deletion with leverage near one, receive
#' `Inf`. An error is raised if no candidate is valid. Grids are sorted and
#' deduplicated; exact ties select the smaller effective weight.
#'
#' This low-level routine takes an already aligned regression. Multi-step
#' direct and fixed-coefficient iterated forecasts use the same source LOO
#' through [fcst.forecast]; rolling-refit CV is not supported. In the
#' two-step version `gamma = 1` gives full-sample GLS; it gives OLS only when
#' the estimated regime variances coincide. `gamma.star = 1` always gives OLS.
#' The finite search grid in direct mode is chosen by the caller: Remark 3
#' does not prescribe an upper bound, and `gamma.star` need not be below one.
#'
#' @return An object of class `bt_fcstWGLS`, a list of:
#' * `forecast`, `coefficients`, `weights`: the prediction, fitted coefficients
#'   and observation weights normalised to sum to one,
#' * `bp`, `N`, `k`: the break date and complete-case sample dimensions,
#' * `gamma`, `gamma.star`: selected kernel and effective weights (`gamma` is
#'   `NA` in direct mode),
#' * `q`, `sigma`: pilot standard deviation ratio and named regime standard
#'   deviations (`NULL` in direct mode),
#' * `cv`: a data frame of `gamma`, `gamma.star` and the unscaled `criterion`,
#' * `cv.min`, `loo.errors`: the selected criterion and its post-break errors,
#' * `fit`: the [WLS.reg] fit on all observations at the selected weight.
#' The adaptive `fit` retains residual diagnostics, but its covariance,
#' standard errors and t-statistics are unavailable: fixed-weight inference
#' does not account for selection of the break date or weight.
#'
#' @references
#' Lee, Tae-Hwy, Shahnaz Parsaeian, and Aman Ullah.
#' "Forecasting under Structural Breaks Using Improved Weighted Estimation."
#' Oxford Bulletin of Economics and Statistics, 2022.
#'
#' @export
fcst.WGLS.cv <- function(
  y,
  x = NULL,
  x.new = NULL,
  bp = NULL,
  trim = 0.15,
  gamma.grid = seq(0, 1, by = 0.01),
  gamma.star.grid = NULL,
  time.index = NULL
) {
  setup <- .fcst.wgls.setup(y, x, x.new, bp, trim,
                             if (missing(gamma.grid)) NULL else gamma.grid,
                             gamma.star.grid, time.index)
  y <- setup$y
  x <- setup$x
  x.new <- setup$x.new
  N <- setup$N
  k <- setup$k
  bp <- setup$bp
  post <- setup$post
  grid <- setup$grid
  gamma <- setup$gamma
  gamma.star <- setup$gamma.star
  q <- setup$q
  sigma <- setup$sigma

  criterion <- rep(Inf, length(grid))
  best <- NULL
  best.value <- Inf
  for (i in seq_along(grid)) {
    w <- weights.AA.kernel(N, bp, gamma.star[i])
    candidate <- .fcst.wgls.loo(y, x, w, post)
    if (!is.null(candidate)) {
      criterion[i] <- mean(candidate$loo.errors^2)
      if (is.finite(criterion[i]) && criterion[i] < best.value) {
        best <- candidate
        best.value <- criterion[i]
      }
    }
  }
  if (is.null(best)) {
    stop("ERROR! fcst.WGLS.cv: no weight gives nonsingular post-break LOO fits")
  }
  selected <- which.min(criterion)

  result <- list(
    forecast = .fcst.predict(best$fit, x.new),
    coefficients = best$fit$coefficients,
    weights = best$fit$weights,
    bp = bp,
    N = N,
    k = k,
    gamma = gamma[selected],
    gamma.star = gamma.star[selected],
    q = q,
    sigma = sigma,
    cv = data.frame(gamma = gamma, gamma.star = gamma.star,
                    criterion = criterion),
    cv.min = best.value,
    loo.errors = best$loo.errors,
    fit = best$fit
  )

  class(result) <- "bt_fcstWGLS"
  .fcst.suppress.inference(result)
}


# Fixed-weight PRESS identity, evaluated only at the post-break rows.
# A singular full or deleted fit invalidates this candidate, not the grid.
.fcst.wgls.loo <- function(y, x, w, post) {
  fit <- tryCatch({
    model <- WLS.reg(y, x, w)
    xwx.inv <- solve(crossprod(x, x * w))
    list(model = model, xwx.inv = xwx.inv)
  }, error = function(e) NULL)
  if (is.null(fit)) {
    return(NULL)
  }

  x.post <- x[post, , drop = FALSE]
  h <- w[post] * rowSums((x.post %*% fit$xwx.inv) * x.post)
  if (any(!is.finite(h)) || any(1 - h <= sqrt(.Machine$double.eps))) {
    return(NULL)
  }
  loo.errors <- fit$model$residuals[post] / (1 - h)
  if (any(!is.finite(loo.errors))) {
    return(NULL)
  }
  list(fit = fit$model, loo.errors = loo.errors)
}


# Pilot estimates and weight grid for the paper's post-break LOO criterion.
.fcst.wgls.setup <- function(y, x, x.new, bp, trim,
                              gamma.grid = NULL, gamma.star.grid = NULL,
                              time.index = NULL) {
  direct <- !is.null(gamma.star.grid)
  if (direct && !is.null(gamma.grid)) {
    stop("ERROR! fcst.WGLS.cv: supply only one of gamma.grid and gamma.star.grid")
  }
  if (!direct && is.null(gamma.grid)) gamma.grid <- seq(0, 1, by = 0.01)
  grid <- if (direct) gamma.star.grid else gamma.grid
  if (!is.numeric(grid) || !is.null(dim(grid)) || length(grid) == 0 ||
      any(!is.finite(grid)) || any(grid < 0) || (!direct && any(grid > 1))) {
    stop("ERROR! fcst.WGLS.cv: invalid weight grid")
  }
  grid <- sort(unique(grid))

  time.index <- .fcst.time.index(time.index, NROW(y))
  .d <- .fcst.data(y, x, x.new)
  .fcst.require.regular(time.index[.d$rows], "fcst.WGLS.cv")
  y <- .d$y
  x <- .d$x
  x.new <- .d$x.new
  N <- .d$N
  k <- .d$k
  if (!is.numeric(y) || !is.numeric(x) || !is.numeric(x.new) ||
      any(!is.finite(y)) || any(!is.finite(x)) || any(!is.finite(x.new)) ||
      nrow(x.new) != 1 || k == 0) {
    stop("ERROR! fcst.WGLS.cv: finite numeric data and one row of x.new required")
  }

  if (is.null(bp)) {
    if (!is.numeric(trim) || length(trim) != 1 || !is.finite(trim) ||
        trim <= 0 || trim >= 0.5) {
      stop("ERROR! fcst.WGLS.cv: trim should be in (0, 0.5)")
    }
    bp <- fcst.window.LS(y, x, trim)$bp
  }
  if (!is.numeric(bp) || length(bp) != 1 || !is.finite(bp) ||
      bp != floor(bp) || bp <= k || N - bp <= k) {
    stop("ERROR! fcst.WGLS.cv: bp must leave more than k observations per regime")
  }
  pre <- seq_len(bp)
  post <- seq.int(bp + 1, N)
  x.pre <- x[pre, , drop = FALSE]
  x.post <- x[post, , drop = FALSE]
  if (qr(x.pre)$rank < k || qr(x.post)$rank < k) {
    stop("ERROR! fcst.WGLS.cv: each regime must have full column rank")
  }

  sigma <- q <- NULL
  if (direct) {
    gamma <- rep(NA_real_, length(grid))
    gamma.star <- grid
  } else {
    # OLS.reg from breaktest supplies the pilot residuals. Use RSS/n, as in
    # the MATLAB pilot estimates, but keep variances distinct from SDs.
    pre.fit <- OLS.reg(y[pre, , drop = FALSE], x.pre)
    post.fit <- OLS.reg(y[post, , drop = FALSE], x.post)
    sigma.sq <- c(pre = mean(pre.fit$residuals^2),
                  post = mean(post.fit$residuals^2))
    if (any(!is.finite(sigma.sq)) || any(sigma.sq <= 0)) {
      stop("ERROR! fcst.WGLS.cv: two-step WGLS requires positive regime variances")
    }
    sigma <- sqrt(sigma.sq)
    q.sq <- unname(sigma.sq["pre"] / sigma.sq["post"])
    if (!is.finite(q.sq) || q.sq <= 0) {
      stop("ERROR! fcst.WGLS.cv: the estimated variance ratio is not finite and positive")
    }
    q <- sqrt(q.sq)
    gamma <- grid
    gamma.star <- gamma / q.sq
    if (any(!is.finite(gamma.star))) {
      stop("ERROR! fcst.WGLS.cv: effective weights overflow; rescale the data")
    }
  }

  list(y = y, x = x, x.new = x.new, N = N, k = k, bp = bp, post = post,
       grid = grid, gamma = gamma, gamma.star = gamma.star, q = q, sigma = sigma)
}

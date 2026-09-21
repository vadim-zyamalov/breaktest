#' @title
#' Coerce and validate data for forecasting routines
#'
#' @description
#' Common entry point for all `fcst.*` functions. Converts `y` and `x` to
#' matrices, checks conformability and drops incomplete rows.
#'
#' @param y A dependent variable.
#' @param x Explanatory variables. If `NULL` a single constant is used.
#' @param x.new A (1 x k) vector of regressor values used to form the forecast.
#' If `NULL` the last row of `x` is reused.
#'
#' @return A list of:
#' * `y`: cleaned (N x 1) matrix,
#' * `x`: cleaned (N x k) matrix,
#' * `x.new`: (1 x k) matrix,
#' * `N`, `k`: dimensions.
#'
#' @importFrom Rfast rowAll
#'
#' @keywords internal
.fcst.data <- function(y, x = NULL, x.new = NULL) {
  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }
  if (!is.numeric(y) || ncol(y) != 1 || any(is.infinite(y))) {
    stop("ERROR! .fcst.data: finite numeric data required; y must be univariate without infinite values")
  }

  if (is.null(x)) {
    x <- matrix(.const(nrow(y)), ncol = 1)
    colnames(x) <- "const"
  }
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  if (nrow(x) != nrow(y)) {
    stop("ERROR! .fcst.data: y and x have different number of rows")
  }
  if (!is.numeric(x) || !ncol(x) || any(is.infinite(x))) {
    stop("ERROR! .fcst.data: finite numeric data required; x must conform to y without infinite values")
  }

  rows <- rowAll(!is.na(y)) & rowAll(!is.na(x))
  y <- y[rows, , drop = FALSE]
  x <- x[rows, , drop = FALSE]

  N <- nrow(y)
  k <- ncol(x)

  if (N <= k + 1) {
    stop("ERROR! .fcst.data: not enough observations")
  }

  if (is.null(x.new)) {
    x.new <- x[N, , drop = FALSE]
  } else {
    if (!is.matrix(x.new)) {
      x.new <- matrix(x.new, nrow = 1)
    }
    if (!is.numeric(x.new) || nrow(x.new) != 1L || ncol(x.new) != k ||
        any(!is.finite(x.new))) {
      stop("ERROR! .fcst.data: x.new must have one row of finite numeric data conformable with x")
    }
  }

  list(
    y = y,
    x = x,
    x.new = x.new,
    N = N,
    k = k,
    rows = which(rows)
  )
}


# A calendar-aware covariance is not a missing-observation theorem. Statistical
# procedures use a contiguous complete block; leading lag losses are allowed.
.fcst.require.regular <- function(time.index, fun) {
  if (length(time.index) > 1L && any(diff(time.index) != 1L)) {
    stop("ERROR! ", fun, ": consecutive complete observations required; internal gaps are unsupported")
  }
  invisible(time.index)
}


#' @title
#' Heteroskedasticity and autocorrelation consistent variance
#'
#' @description
#' Long-run variance (or covariance matrix) of a series used by the
#' Diebold-Mariano, Clark-West and Fluctuation tests.
#'
#' This function replaces three nearly identical routines found in the
#' replication codes (`HAC_est.m`, `nw.m` and `nwest.m` of Rossi's package and
#' the inline estimator of `dmtstat.m` by Hirano and Wright). Kernel weights are
#' taken from the package-wide `.lr.weight` so that all HAC calculations in
#' `breaktest` share a single implementation.
#'
#' @param y A series (or a matrix of series) of interest.
#' @param bandwidth Positive integer scale B. Compact kernels include lags
#' through B-1; Quadratic uses all observed lags with scale B. If `NULL`, the
#' heuristic \eqn{\lfloor N^{1/4} \rfloor} is used. It is not an automatic
#' optimal bandwidth selector.
#' @param kernel A kernel name passed to the package helper `.lr.weight`.
#' @param demean Whether the series should be demeaned first.
#' @param time.index Integer original-period positions. Gaps are preserved
#' when forming lag products, rather than compressed into adjacent rows.
#'
#' @return A scalar for a univariate `y`, a (k x k) matrix otherwise.
#'
#' @references
#' Newey, Whitney K., and Kenneth D. West.
#' “A Simple, Positive Semi-Definite, Heteroskedasticity and Autocorrelation
#' Consistent Covariance Matrix.”
#' Econometrica 55, no. 3 (1987): 703–8.
#' https://doi.org/10.2307/1913610.
#'
#' @keywords internal
HAC.variance <- function(
  y,
  bandwidth = NULL,
  kernel = "Bartlett",
  demean = TRUE,
  time.index = NULL
) {
  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }
  if (!is.numeric(y) || !ncol(y) || !is.logical(demean) ||
      length(demean) != 1L || is.na(demean)) {
    stop("ERROR! HAC.variance: numeric data and a scalar logical demean required")
  }
  time.index <- .fcst.time.index(time.index, nrow(y))
  keep <- rowAll(is.finite(y))
  y <- y[keep, , drop = FALSE]
  time.index <- time.index[keep]

  N <- nrow(y)
  if (N < 2) {
    stop("ERROR! HAC.variance: not enough observations")
  }

  if (is.null(bandwidth)) {
    bandwidth <- floor(N^(1 / 4))
  }
  if (!is.numeric(bandwidth) || length(bandwidth) != 1 ||
      !is.finite(bandwidth) || bandwidth < 1 || bandwidth != floor(bandwidth)) {
    stop("ERROR! HAC.variance: bandwidth must be a positive integer")
  }

  if (demean) {
    y <- sweep(y, 2, colMeans(y), "-")
  }

  wgtF <- .lr.weight(kernel)

  res <- crossprod(y) / N
  span <- max(time.index) - min(time.index)
  max.lag <- if (identical(kernel, "Quadratic")) span else min(bandwidth - 1L, span)
  kernel.scale <- if (identical(kernel, "Quadratic")) bandwidth else bandwidth - 1L
  if (max.lag > 0) {
    for (j in seq_len(max.lag)) {
      previous <- match(time.index - j, time.index)
      current <- which(!is.na(previous))
      if (!length(current)) next
      gamma <- crossprod(
        y[current, , drop = FALSE],
        y[previous[current], , drop = FALSE]
      ) / N
      res <- res + (gamma + t(gamma)) * wgtF(j, kernel.scale)
    }
  }

  if (ncol(res) == 1) drop(res) else res
}


.fcst.time.index <- function(time.index, N) {
  if (is.null(time.index)) return(seq_len(N))
  if (!is.numeric(time.index) || length(time.index) != N ||
      any(!is.finite(time.index)) || any(time.index != floor(time.index)) ||
      any(diff(time.index) <= 0)) {
    stop("ERROR! HAC.variance: time.index must contain increasing integer positions")
  }
  time.index
}


#' @title
#' Linear interpolation of a tabulated function on a finer grid
#'
#' @description
#' Auxiliary routine used to evaluate tabulated critical values and equal
#' performance boundaries at an arbitrary break fraction. Replaces
#' `extend_boundary.m` of Boot and Pick (2019), which performed the same
#' operation through an explicit interpolation matrix.
#'
#' @param values A vector of tabulated values.
#' @param grid A vector of arguments the values are tabulated at.
#' @param at A vector of points the values are needed at.
#'
#' @return A vector of interpolated values. Points outside of `grid` are
#' assigned the nearest boundary value.
#'
#' @keywords internal
.grid.interp <- function(values, grid, at) {
  if (length(values) != length(grid)) {
    stop("ERROR! .grid.interp: values and grid have different length")
  }

  keep <- is.finite(values) & is.finite(grid)
  if (sum(keep) == 0) {
    return(rep(NA_real_, length(at)))
  }
  if (sum(keep) == 1) {
    return(rep(values[keep], length(at)))
  }

  values <- values[keep]
  grid <- grid[keep]

  at <- pmin(pmax(at, min(grid)), max(grid))
  approx(x = grid, y = values, xout = at, method = "linear")$y
}


#' @title
#' Standardise a matrix column-wise keeping the moments
#'
#' @description
#' Analogue of MATLAB's `zscore` used in the empirical part of Boot and Pick
#' (2019). Forecasts are produced for the standardised series and then
#' transformed back, so the moments have to be returned as well.
#'
#' @param x A matrix to be standardised.
#'
#' @return A list of the standardised matrix `z`, the vector of means `mu` and
#' the vector of standard deviations `sigma`.
#'
#' @importFrom stats sd
#'
#' @keywords internal
.fcst.zscore <- function(x) {
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }

  mu <- colMeans(x, na.rm = TRUE)
  sigma <- apply(x, 2, sd, na.rm = TRUE)
  sigma[sigma == 0 | is.na(sigma)] <- 1

  list(
    z = sweep(sweep(x, 2, mu, "-"), 2, sigma, "/"),
    mu = mu,
    sigma = sigma
  )
}


#' @title
#' Split-sample design matrix for a break at a given date
#'
#' @param x A matrix of regressors.
#' @param bp A break point, the last observation of the first regime.
#' @param break.vars Indices of the columns subject to the break.
#' If `NULL` all the columns are allowed to break (pure structural change).
#'
#' @return A list of:
#' * `xs`: the expanded design matrix,
#' * `idx.1`, `idx.2`: positions of the pre- and post-break coefficients
#'   of the breaking regressors inside the expanded coefficient vector,
#' * `map.2`: positions producing the post-break coefficient vector in terms of
#'   the original regressors.
#'
#' @keywords internal
.fcst.split.design <- function(x, bp, break.vars = NULL) {
  N <- nrow(x)
  k <- ncol(x)

  if (is.null(break.vars)) {
    break.vars <- seq_len(k)
  }
  stable.vars <- setdiff(seq_len(k), break.vars)

  # .du(bp, N) is the package-wide break-in-constant dummy, 1 after the break
  d2 <- drop(.du(bp, N))
  d1 <- 1 - d2

  xs <- cbind(
    x[, break.vars, drop = FALSE] * d1,
    x[, break.vars, drop = FALSE] * d2,
    if (length(stable.vars) > 0) x[, stable.vars, drop = FALSE] else NULL
  )

  nb <- length(break.vars)

  map.2 <- numeric(k)
  map.2[break.vars] <- nb + seq_len(nb)
  if (length(stable.vars) > 0) {
    map.2[stable.vars] <- 2 * nb + seq_along(stable.vars)
  }

  list(
    xs = xs,
    idx.1 = seq_len(nb),
    idx.2 = nb + seq_len(nb),
    map.2 = map.2,
    break.vars = break.vars,
    stable.vars = stable.vars
  )
}


#' @title
#' Covariance matrix of OLS coefficients with optional HC correction
#'
#' @param x A design matrix.
#' @param resid A vector of residuals.
#' @param vcov.type One of `const` (homoskedastic), `HC0`, `HC3` or `HAC`.
#' @param bandwidth,time.index Arguments for [HAC.variance] in HAC mode.
#' The empirical part of Boot and Pick (2019) uses `HC3`.
#'
#' @return An estimated covariance matrix of the OLS estimator.
#'
#' @keywords internal
.fcst.vcov <- function(x, resid, vcov.type = "HC3", bandwidth = NULL,
                       time.index = NULL) {
  if (!vcov.type %in% c("const", "HC0", "HC3", "HAC")) {
    stop("ERROR! .fcst.vcov: unknown vcov.type")
  }

  N <- nrow(x)
  k <- ncol(x)

  xx.inv <- spdinv(crossprod(x))
  e2 <- as.numeric(resid)^2

  if (vcov.type == "HAC") {
    meat <- N * as.matrix(HAC.variance(x * as.numeric(resid), bandwidth,
                                       time.index = time.index))
    return(xx.inv %*% meat %*% xx.inv)
  }

  if (vcov.type == "const") {
    return(sum(e2) / (N - k) * xx.inv)
  }

  if (vcov.type == "HC3") {
    lev <- rowSums((x %*% xx.inv) * x)
    lev <- pmin(lev, 1 - 1e-8)
    e2 <- e2 / (1 - lev)^2
  }

  xx.inv %*% crossprod(x, x * e2) %*% xx.inv
}

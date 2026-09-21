#' @title
#' Weighted least squares with arbitrary observation weights
#'
#' @description
#' The workhorse behind every weighting-based forecast in this module:
#' \deqn{\hat\beta(w) = \left(\sum_t w_t x_t x_t'\right)^{-1}
#' \sum_t w_t x_t y_t.}
#'
#' Weights are *not* required to be positive. Pesaran et al. (2013) show that
#' under multiple breaks the MSFE-optimal weights can legitimately be negative,
#' so no restriction is imposed here.
#'
#' @param y A dependent variable.
#' @param x Explanatory variables.
#' @param w A vector of observation weights. If `NULL` equal weights are used,
#' which reproduces [OLS.reg].
#'
#' @return An object of class `bt_wls` inheriting from `bt_ols` with the usual
#' fields plus the `weights` used.
#' For nonnegative fixed weights, the residual scale uses the exact weighted
#' residual degrees of freedom under a common error variance. This working
#' covariance does not account for estimated weights, breaks or model selection.
#' With signed weights or no residual degrees of freedom, inference fields are
#' `NA`; the coefficient estimate and point forecast remain available.
#'
#' @references
#' Pesaran, M. Hashem, Andreas Pick, and Mikhail Pranovich.
#' “Optimal Forecasts in the Presence of Structural Breaks.”
#' Journal of Econometrics 177, no. 2 (2013): 134–52.
#'
#' @importFrom Rfast rowAll
#' @importFrom Rfast spdinv
#' @importFrom stats AIC BIC nobs
#'
#' @keywords internal
WLS.reg <- function(y, x, w = NULL) {
  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }

  N <- nrow(y)
  if (!is.numeric(y) || ncol(y) != 1L || !is.numeric(x) ||
      nrow(x) != N || !ncol(x) || any(is.infinite(y)) || any(is.infinite(x))) {
    stop("ERROR! WLS.reg: numeric univariate y and a conformable finite design required")
  }

  if (is.null(w)) {
    w <- rep(1 / N, N)
  }
  if (!is.numeric(w) || !is.null(dim(w)) || length(w) != N || any(is.infinite(w))) {
    stop("ERROR! WLS.reg: finite numeric weights conformable with y required")
  }

  rows <- rowAll(!is.na(y)) & rowAll(!is.na(x)) & !is.na(w)

  .y <- y[rows, , drop = FALSE]
  .x <- x[rows, , drop = FALSE]
  .w <- w[rows]
  if (!length(.w) || !any(.w != 0)) {
    stop("ERROR! WLS.reg: at least one nonzero complete observation weight required")
  }
  .w <- .w / max(abs(.w))

  xwx <- crossprod(.x, .x * .w)
  xwy <- crossprod(.x, .y * .w)

  signed <- any(.w < 0)
  if (signed) {
    cf <- drop(qr.solve(xwx, xwy))
    xwx.inv <- qr.solve(xwx)
  } else {
    weighted <- .x * sqrt(.w)
    decomposition <- qr(weighted)
    if (decomposition$rank < ncol(.x)) {
      stop("ERROR! WLS.reg: weighted design is rank deficient")
    }
    cf <- drop(qr.coef(decomposition, .y * sqrt(.w)))
    inverse <- chol2inv(qr.R(decomposition))
    xwx.inv <- inverse[order(decomposition$pivot), order(decomposition$pivot), drop = FALSE]
  }

  r <- drop(.y - .x %*% cf)
  meat <- crossprod(.x, .x * .w^2)
  df <- if (signed) NA_real_ else sum(.w) - sum(diag(xwx.inv %*% meat))
  s.sq <- NA_real_
  vcov <- matrix(NA_real_, ncol(.x), ncol(.x))
  inference <- "unavailable for signed weights"
  if (!signed) {
    inference <- "unavailable without residual degrees of freedom"
    if (is.finite(df) && df > sqrt(.Machine$double.eps) * sum(.w)) {
      s.sq <- sum(.w * r^2) / df
      vcov <- s.sq * xwx.inv %*% meat %*% xwx.inv
      vcov <- (vcov + t(vcov)) / 2
      inference <- "common error variance, conditional on fixed nonnegative weights and design"
    }
  }
  se.cf <- sqrt(pmax(diag(vcov), 0))
  t.stats <- ifelse(is.finite(se.cf) & se.cf > 0, cf / se.cf, NA_real_)

  resid <- rep(NA, N)
  resid[rows] <- r

  fitted <- rep(NA, N)
  fitted[rows] <- drop(.x %*% cf)

  result <- list(
    coefficients = cf,
    s.sq = s.sq,
    se.coefs = se.cf,
    t.stats = t.stats,
    vcov = vcov,
    residual.df = df,
    inference = inference,
    residuals = resid,
    fitted.values = fitted,
    weights = w,
    endog = y,
    exog = x
  )

  class(result) <- c("bt_wls", "bt_ols")
  result
}


# Selected weights/designs do not inherit fixed-weight WLS inference. Preserve
# residuals and their working scale for diagnostics, but not inferential fields.
.fcst.suppress.inference <- function(object) {
  if (inherits(object, "bt_wls")) {
    object$vcov[] <- NA_real_
    object$se.coefs[] <- NA_real_
    object$t.stats[] <- NA_real_
    object$inference <- "unavailable after data-dependent selection of weights, windows or regressors"
  } else if (is.list(object) && !is.data.frame(object)) {
    for (j in seq_along(object)) object[j] <- list(.fcst.suppress.inference(object[[j]]))
  }
  object
}


#' @title
#' Point forecast from a fitted regression object
#'
#' @param obj An object returned by [OLS.reg], [WLS.reg] or any other routine
#' carrying a `coefficients` field.
#' @param x.new A (1 x k) vector of regressor values.
#'
#' @return A scalar point forecast.
#'
#' @keywords internal
.fcst.predict <- function(obj, x.new) {
  if (!is.matrix(x.new)) {
    x.new <- matrix(x.new, nrow = 1)
  }
  drop(x.new %*% obj$coefficients)
}


#' @title
#' Sum of squared residuals of a two-regime regression
#'
#' @description
#' Concentrated SSR of the split-sample regression used by the least squares
#' break date estimator and by the quasi-likelihood ratio statistic of Koo and
#' Seo (2015). Builds the partially restricted split design and calls [OLS.reg].
#' Full-break profiles use `.fcst.ssr.matrix` to reuse the package recursion.
#'
#' @param y A dependent variable.
#' @param x Explanatory variables.
#' @param bp A break point, the last observation of the first regime.
#' @param break.vars Indices of the columns subject to the break.
#'
#' @return The SSR value.
#'
#' @keywords internal
.fcst.ssr.split <- function(y, x, bp, break.vars = NULL) {
  .s <- .fcst.split.design(x, bp, break.vars)
  if (qr(.s$xs)$rank < ncol(.s$xs)) return(Inf)
  .m <- OLS.reg(y, .s$xs)
  sum(.m$residuals^2, na.rm = TRUE)
}


# Reuse the package recursion, extending only an unidentified starting block.
# Singular shorter segments stay at Inf; they are not valid candidate regimes.
.fcst.ssr.recursive <- function(y, x, beg, end, width = 2) {
  y <- as.matrix(y)
  x <- as.matrix(x)
  result <- rep(Inf, nrow(y))
  first <- max(beg + width - 1L, beg + ncol(x))
  if (first > end) return(result)
  while (first <= end && qr(x[seq.int(beg, first), , drop = FALSE])$rank < ncol(x)) {
    first <- first + 1L
  }
  if (first > end) return(result)
  value <- tryCatch(SSR.recursive(y, x, beg, end, first - beg + 1L),
                    error = function(e) NULL)
  if (!is.null(value) && all(is.finite(value[seq.int(first, end)]))) return(value)
  # A QR fallback also handles numerically singular normal equations without
  # changing breaktest. Only used when its recursive estimator cannot run.
  for (last in seq.int(first, end)) {
    rows <- seq.int(beg, last)
    model <- qr(x[rows, , drop = FALSE])
    if (model$rank == ncol(x)) result[last] <- sum(qr.resid(model, y[rows, , drop = FALSE])^2)
  }
  result
}


.fcst.ssr.matrix <- function(y, x, width = 2) {
  y <- as.matrix(y)
  x <- as.matrix(x)
  N <- nrow(y)
  if (!is.numeric(width) || length(width) != 1L || !is.finite(width) ||
      width != floor(width) || width < 1 || width > N) {
    stop("ERROR! .fcst.ssr.matrix: invalid minimum segment width")
  }
  value <- tryCatch(SSR.matrix(y, x, width), error = function(e) NULL)
  if (!is.null(value) && !anyNA(value)) return(value)
  result <- matrix(Inf, N, N)
  for (first in seq_len(N - width + 1L)) {
    result[first, ] <- .fcst.ssr.recursive(y, x, first, N, width)
  }
  result
}

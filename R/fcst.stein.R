#' @title
#' Stein-like combination under the published iid GLS model
#'
#' @description
#' Combine the full-sample estimate and the estimate after a structural
#' break, using alpha = min(1, tau/H) as in Lee, Parsaeian and Ullah.
#' The implementation is restricted to a single coefficient break and
#' uncorrelated errors with constant variance within each regime.
#'
#' @param y,x,x.new An aligned regression as in [fcst.forecast].
#' @param bp Last pre-break observation in the complete-case sample. NULL
#' estimates a single break with [fcst.window.LS]. The model must have only
#' one actual break; supplying the last of several breaks is not supported.
#' @param tau Nonnegative shrinkage parameter. Default max(k-2,0).
#' @param trim Trimming for estimation of an unknown break date.
#' @param gls Must be `TRUE`: use inverse regime-variance weights in the
#' full-sample fit. The former OLS extension is not implemented.
#' @param vcov.type Must be `"iid"`. The former HAC shrinkage extension
#' has been removed because the published risk rule does not justify it.
#' @param bandwidth,kernel Unsupported legacy HAC settings; must be `NULL`.
#' @param time.index Original integer-period positions. The estimation block
#' must be contiguous after trimming leading or trailing missing rows.
#'
#' @details
#' Write B_F = (X'aX)^(-1), B_P = (X_P'X_P)^(-1), with scalar pilot weights
#' a_t. Rows of the linear maps are `L_F[t,]` = a_t*x_t'*B_F and
#' `L_P[t,]` = 1(t>bp)*x_t'*B_P. The error covariance is diagonal with
#' regime variances estimated from regime OLS residuals. Applying
#' these linear maps gives coefficient covariance blocks V_F,
#' V_P and C=Cov(beta_F,beta_P). Then
#' \deqn{V_D=V_P+V_F-C-C',\qquad H=(\hat\beta_P-\hat\beta_F)'V_D^{-1}
#' (\hat\beta_P-\hat\beta_F).}
#' No further factor N belongs in H. With GLS the difference covariance
#' reduces to V_D=V_P-V_F.
#'
#' The source paper assumes uncorrelated errors within and across regimes.
#' The `risk.eligible` flag checks only dimension and the range of tau in
#' the source result. Its model, estimation and quadratic-loss conditions
#' W=V_D^(-1) must also hold. It does not assert lower MSFE for an arbitrary
#' forecast direction, or verify the assumptions from observed data.
#' tau=0 returns the post-break fit when the method is identified. A singular
#' difference covariance produces an error, without a fallback strategy.
#'
#' @return A bt_fcstStein with forecast, coefficients, alpha, H, tau,
#' beta.full, beta.post, forecast.full, forecast.post, regime sigma,
#' vcov (joint/full/post/cross/difference), covariance.rank, `risk.eligible`
#' and the conditional interpretation `risk.conditions`. No unconditional
#' `dominance` result is returned.
#'
#' @references
#' Lee, T.-H., S. Parsaeian and A. Ullah (2022). "Optimal Forecast under
#' Structural Breaks." Journal of Applied Econometrics 37, 965-987.
#' Equations (4)-(6), Theorem 1 and Corollary 3.1.
#' https://doi.org/10.1002/jae.2908.
#'
#' @export
fcst.stein <- function(
  y, x = NULL, x.new = NULL, bp = NULL, tau = NULL, trim = 0.15,
  gls = TRUE, vcov.type = "iid", bandwidth = NULL,
  kernel = NULL, time.index = NULL
) {
  if (!identical(vcov.type, "iid")) {
    stop("ERROR! fcst.stein: only the published iid covariance model is supported")
  }
  if (!isTRUE(gls)) stop("ERROR! fcst.stein: only gls=TRUE is supported")
  if (!is.null(bandwidth) || !is.null(kernel)) {
    stop("ERROR! fcst.stein: bandwidth and kernel are unsupported HAC settings")
  }
  time.index <- .fcst.time.index(time.index, NROW(y))
  .d <- .fcst.data(y, x, x.new)
  y <- .d$y
  x <- .d$x
  x.new <- .d$x.new
  N <- .d$N
  k <- .d$k
  time.index <- time.index[.d$rows]
  .fcst.require.regular(time.index, "fcst.stein")
  if (!is.numeric(y) || !is.numeric(x) || any(!is.finite(y)) ||
      any(!is.finite(x)) || any(!is.finite(x.new)) || nrow(x.new) != 1L) {
    stop("ERROR! fcst.stein: finite numeric data and one forecast row required")
  }
  if (is.null(bp)) bp <- fcst.window.LS(y, x, trim)$bp
  if (!is.numeric(bp) || length(bp) != 1L || !is.finite(bp) ||
      bp != floor(bp) || bp <= k || N - bp <= k) {
    stop("ERROR! fcst.stein: each regime must contain more than k observations")
  }
  if (is.null(tau)) tau <- max(k - 2, 0)
  if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) || tau < 0) {
    stop("ERROR! fcst.stein: tau must be a finite nonnegative scalar")
  }
  pre <- seq_len(bp)
  post <- seq.int(bp + 1L, N)
  if (qr(x[pre, , drop = FALSE])$rank < k || qr(x[post, , drop = FALSE])$rank < k) {
    stop("ERROR! fcst.stein: both regime designs must have full rank")
  }
  m1 <- OLS.reg(y[pre, , drop = FALSE], x[pre, , drop = FALSE])
  m2 <- OLS.reg(y[post, , drop = FALSE], x[post, , drop = FALSE])
  variances <- c(pre = mean(m1$residuals^2), post = mean(m2$residuals^2))
  if (any(!is.finite(variances)) || any(variances <= 0)) {
    stop("ERROR! fcst.stein: positive residual variance required in each regime")
  }
  sigma.sq <- c(rep(variances[1], length(pre)), rep(variances[2], length(post)))
  a <- 1 / sigma.sq
  full <- WLS.reg(y, x, a / sum(a))
  L.full <- (x * a) %*% solve(crossprod(x, x * a))
  L.post <- matrix(0, N, k)
  L.post[post, ] <- x[post, , drop = FALSE] %*% solve(crossprod(x[post, , drop = FALSE]))
  maps <- cbind(L.full, L.post)
  joint <- crossprod(maps * sqrt(sigma.sq))
  f <- seq_len(k)
  p <- k + f
  v.full <- joint[f, f, drop = FALSE]
  v.post <- joint[p, p, drop = FALSE]
  cross <- joint[f, p, drop = FALSE]
  dv <- v.full + v.post - cross - t(cross)
  dv <- (dv + t(dv)) / 2
  ev <- eigen(dv, symmetric = TRUE)
  tolerance <- max(abs(ev$values)) * 1e-10
  if (any(ev$values < -tolerance)) {
    stop("ERROR! fcst.stein: difference covariance is indefinite")
  }
  positive <- ev$values > tolerance
  rank <- sum(positive)
  if (rank < k) stop("ERROR! fcst.stein: singular difference covariance; method is unidentified")
  b.full <- as.numeric(full$coefficients)
  b.post <- as.numeric(m2$coefficients)
  d <- b.post - b.full
  H <- sum(drop(crossprod(ev$vectors, d))^2 / ev$values)
  alpha <- if (tau == 0) 0 else if (H == 0) 1 else min(1, tau / H)
  cf <- alpha * b.full + (1 - alpha) * b.post
  risk.eligible <- k > 2 && tau > 0 && tau <= 2 * (k - 2)
  risk.conditions <- paste("dimension and tau eligibility only; requires the source single-break",
                           "iid GLS assumptions and quadratic loss W=V_D^(-1);",
                           "does not guarantee MSFE improvement for arbitrary x.new")
  structure(list(coefficients = cf, forecast = drop(x.new %*% cf),
                  alpha = alpha, H = H, tau = tau, risk.eligible = risk.eligible,
                  risk.conditions = risk.conditions,
                  bp = bp, beta.full = b.full, beta.post = b.post,
                  forecast.full = drop(x.new %*% b.full), forecast.post = drop(x.new %*% b.post),
                  sigma = sqrt(variances), vcov.type = vcov.type, gls = gls,
                  time.index = time.index,
                  covariance.rank = rank, diagnostic = risk.conditions,
                  vcov = list(full = v.full, post = v.post, cross = cross,
                              difference = dv, joint = joint)), class = "bt_fcstStein")
}

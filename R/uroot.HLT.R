#' @title
#' Unit root testing procedure under a single structural break.
#'
#' @param y A time series of interest.
#' @param const Whether a constant should be included.
#' @param trim The trimming parameter to find the lower and upper bounds of
#' possible break dates.
#'
#' @return The value of test statistic.
#'
#' @references
#' Harvey, David I., Stephen J. Leybourne, and A. M. Robert Taylor.
#' “Unit Root Testing under a Local Break in Trend.”
#' Journal of Econometrics 167, no. 1 (2012): 140–67.
#'
#' @importFrom Rfast spdinv
#' @export
uroot.HLT <- function(y, const = FALSE, trim = 0.15) {
  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }

  N <- nrow(y)
  m_xi <- ifelse(const, 1.052, 0.853)
  dy <- .msub(.diffn(y), -1)

  bp.min <- trunc(trim * N)
  bp.max <- trunc((1 - trim) * N)

  t0 <- -Inf
  t1 <- -Inf

  KPSS_y <- NA
  KPSS_dy <- NA

  for (bp in bp.min:bp.max) {
    du <- .du(bp, N)

    x <- cbind(
      .const(N),
      .trend(N),
      if (const) du else NULL,
      .dt(bp, N)
    )

    .model <- OLS.reg(y, x)
    r <- .model$residuals
    lrv_y <- .lr.var.bartlett(r)
    invxx <- spdinv(t(x) %*% x)

    t0_stats <- abs(
      .model$coefficients[ncol(x)] /
        sqrt(lrv_y * invxx[ncol(x), ncol(x)])
    )

    if (t0_stats > t0) {
      t0 <- t0_stats
      KPSS_y <- .kpss.statistic(r, lrv_y)
    }

    x <- cbind(
      .const(N),
      if (const) .diffn(du) else NULL,
      du
    )[2:N, ]

    .model <- OLS.reg(dy, x)
    r <- .model$residuals
    lrv_dy <- .lr.var.bartlett(r)
    invxx <- spdinv(t(x) %*% x)

    t1_stats <- abs(
      .model$coefficients[ncol(x)] /
        sqrt(lrv_dy * invxx[ncol(x), ncol(x)])
    )

    if (t1_stats > t1) {
      t1 <- t1_stats
      KPSS_dy <- .kpss.statistic(r, lrv_dy)
    }
  }

  lam_KPSS <- exp(-(500 * KPSS_y * KPSS_dy)^2)

  lam_KPSS * t0 + m_xi * (1 - lam_KPSS) * t1
}

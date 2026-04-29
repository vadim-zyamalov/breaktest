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
#' @export
ur.KPSS.HLT <- function(y, const = FALSE, trim = 0.15) {
  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }

  N <- nrow(y)
  m.ksi <- ifelse(const, 1.052, 0.853)
  dy <- diff(y)

  bp.min <- trunc(trim * N)
  bp.max <- trunc((1 - trim) * N)

  t0 <- -Inf
  t1 <- -Inf

  var.y <- NA
  var.dy <- NA

  for (bp in bp.min:bp.max) {
    du <- .du(bp, N)

    x <- cbind(
      .const(N),
      .trend(N),
      if (const) du else NULL,
      .dt(bp, N)
    )

    .model <- OLS.reg(y, x)

    .y.lr.var <- .lr.var.bartlett(.model$residuals)
    .xx.inv <- qr.solve(t(x) %*% x)

    .t0 <- abs(
      .model$coefficients[ncol(x)] /
        sqrt(.y.lr.var * .xx.inv[ncol(x), ncol(x)])
    )

    x <- cbind(
      .const(N),
      if (const) .diffn(du) else NULL,
      du
    )[-1, ]

    .model <- OLS.reg(dy, x)

    .dy.lr.var <- .lr.var.bartlett(.model$residuals)
    .xx.inv <- qr.solve(t(x) %*% x)

    .t1 <- abs(
      .model$coefficients[ncol(x)] /
        sqrt(.dy.lr.var * .xx.inv[ncol(x), ncol(x)])
    )

    if (.t0 > t0) {
      t0 <- .t0
      var.y <- .y.lr.var
    }
    if (.t1 > t1) {
      t1 <- .t1
      var.dy <- .dy.lr.var
    }
  }

  .kpss.y <- .kpss.statistic(.model$residuals, var.y)
  .kpss.dy <- .kpss.statistic(.model$residuals, var.dy)

  .kpss.lmb <- exp(-((500 * .kpss.y * .kpss.dy)^2))

  .kpss.lmb * t0 + m.ksi * (1 - .kpss.lmb) * t1
}

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
kpss_hlt <- function(y,
                     const = FALSE,
                     trim = 0.15) {
  if (!is.matrix(y)) y <- as.matrix(y)

  n_obs <- nrow(y)

  if (!const) {
    m_ksi <- 0.853
  } else {
    m_ksi <- 1.052
  }

  .dy <- diff(y)

  bp_min <- trunc(trim * n_obs)
  bp_max <- trunc((1 - trim) * n_obs)

  t0 <- -Inf
  t1 <- -Inf

  var_y <- NA
  var_dy <- NA

  for (bp in bp_min:bp_max) {
    du <- c(rep(0, bp), rep(1, n_obs - bp))
    dt <- du * (1:n_obs - bp)

    x <- cbind(
      rep(1, n_obs),
      1:n_obs,
      if (const) du else NULL,
      dt
    )

    .model <- .estimate_ols(y, x)

    .var_y_lr <- .variance_lr_bartlett(.model$residuals)
    .xx_inv <- qr.solve(t(x) %*% x)

    .t0 <- abs(.model$beta[ncol(x)] /
      sqrt(.var_y_lr * .xx_inv[ncol(x), ncol(x)]))

    x <- cbind(
      rep(1, n_obs - 1),
      if (const) diff(du) else NULL,
      du[2:n_obs]
    )

    .model <- .estimate_ols(.dy, x)

    .var_dy_lr <- .variance_lr_bartlett(.model$residuals)
    .xx_inv <- qr.solve(t(x) %*% x)

    .t1 <- abs(.model$beta[ncol(x)] /
      sqrt(.var_dy_lr * .xx_inv[ncol(x), ncol(x)]))

    if (.t0 > t0) {
      t0 <- .t0
      var_y <- .var_y_lr
    }
    if (.t1 > t1) {
      t1 <- .t1
      var_dy <- .var_dy_lr
    }
  }

  .kpss_y <- .kpss_stat(.model$residuals, var_y)
  .kpss_dy <- .kpss_stat(.model$residuals, var_dy)

  .kpss_lambda <- exp(-((500 * .kpss_y * .kpss_dy)^2))
  .kpss_lambda_t <- .kpss_lambda * t0 + m_ksi * (1 - .kpss_lambda) * t1

  .kpss_lambda_t
}

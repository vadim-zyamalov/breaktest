#' @title
#' Pre-calculate matrix of recursive SSR values.
#'
#' @param y A dependent variable.
#' @param x Explanatory variables.
#' @param width Minimum spacing between the breaks.
#'
#' @return A matrix of recursive SSR values.
#'
#' @keywords internal
ssr_matrix <- function(y,
                       x,
                       width = 2) {
  if (!is.matrix(y)) y <- as.matrix(y)
  if (!is.matrix(x)) x <- as.matrix(x)

  n_obs <- nrow(y)

  result <- matrix(data = Inf, nrow = n_obs, ncol = n_obs)

  for (i in 1:(n_obs - width + 1)) {
    result[i, 1:n_obs] <- ssr_recursive(
      y, x, i, n_obs, width
    )
  }

  result
}


#' @title
#' Calculate SSR recursively
#'
#' @param y A dependent variable.
#' @param x Explanatory variables.
#' @param beg,end The start and the end of SSR calculating period.
#' @param width Minimum spacing between the breaks.
#'
#' @return The vector of calculated recursive SSR.
#'
#' @references
#' Brown, R. L., J. Durbin, and J. M. Evans.
#' “Techniques for Testing the Constancy of Regression Relationships over Time.”
#' Journal of the Royal Statistical Society.
#' Series B (Methodological) 37, no. 2 (1975): 149–92.
#'
#' @keywords internal
ssr_recursive <- function(y,
                          x,
                          beg,
                          end,
                          width = 2) {
  if (!is.matrix(y)) y <- as.matrix(y)
  if (!is.matrix(x)) x <- as.matrix(x)

  if (beg < 1) beg <- 1
  if (nrow(y) < end) end <- nrow(y)

  n_obs <- nrow(y)

  result <- matrix(data = Inf, nrow = n_obs, ncol = 1)

  y0 <- y[beg:(beg + width - 1), , drop = FALSE]
  x0 <- x[beg:(beg + width - 1), , drop = FALSE]

  xx_inv0 <- qr.solve(t(x0) %*% x0)
  .model <- .estimate_ols(y0, x0)
  beta0 <- .model$beta
  resid0 <- .model$residuals
  rm(.model)

  result[beg + width - 1, 1] <- drop(t(resid0) %*% resid0)

  for (step in (beg + width):end) {
    if (step > end) break

    .x_loop <- x[step, , drop = FALSE]

    .resid_loop <- drop(y[step, , drop = FALSE] - .x_loop %*% beta0)

    denom <- 1 + .x_loop %*% xx_inv0 %*% t(.x_loop)
    denom <- drop(denom)

    xx_inv1 <- xx_inv0 -
      (xx_inv0 %*% t(.x_loop) %*% .x_loop %*% xx_inv0) / denom

    beta1 <- beta0 + xx_inv1 %*% t(.x_loop) * .resid_loop

    result[step, 1] <- result[step - 1, 1] + .resid_loop^2 / denom

    xx_inv0 <- xx_inv1

    beta0 <- beta1
  }

  result
}

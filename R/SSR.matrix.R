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
SSR.matrix <- function(y,
                       x,
                       width = 2) {
  if (!is.matrix(y)) y <- as.matrix(y)
  if (!is.matrix(x)) x <- as.matrix(x)

  n.obs <- nrow(y)

  result <- matrix(data = Inf, nrow = n.obs, ncol = n.obs)

  for (i in 1:(n.obs - width + 1)) {
    result[i, 1:n.obs] <- SSR.recursive(
      y, x, i, n.obs, width
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
SSR.recursive <- function(y,
                          x,
                          beg,
                          end,
                          width = 2) {
  if (!is.matrix(y)) y <- as.matrix(y)
  if (!is.matrix(x)) x <- as.matrix(x)

  n.obs <- nrow(y)

  beg <- max(beg, 1)
  end <- min(end, n.obs)

  result <- matrix(data = Inf, nrow = n.obs, ncol = 1)

  y0 <- y[beg:(beg + width - 1), , drop = FALSE]
  x0 <- x[beg:(beg + width - 1), , drop = FALSE]

  xx.inv0 <- qr.solve(t(x0) %*% x0)
  .model <- .OLS(y0, x0)
  beta0 <- .model$beta
  resid0 <- .model$residuals
  rm(.model)

  result[beg + width - 1, 1] <- drop(t(resid0) %*% resid0)

  for (step in (beg + width):end) {
    if (step > end) break

    .x.loop <- x[step, , drop = FALSE]

    .resid_loop <- drop(y[step, , drop = FALSE] - .x.loop %*% beta0)

    denom <- 1 + .x.loop %*% xx.inv0 %*% t(.x.loop)
    denom <- drop(denom)

    xx.inv1 <- xx.inv0 -
      (xx.inv0 %*% t(.x.loop) %*% .x.loop %*% xx.inv0) / denom

    beta1 <- beta0 + xx.inv1 %*% t(.x.loop) * .resid_loop

    result[step, 1] <- result[step - 1, 1] + .resid_loop^2 / denom

    xx.inv0 <- xx.inv1

    beta0 <- beta1
  }

  result
}

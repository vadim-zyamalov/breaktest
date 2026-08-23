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
SSR.matrix <- function(y, x, width = 2) {
  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }

  N <- nrow(y)

  result <- matrix(data = Inf, nrow = N, ncol = N)

  for (i in 1:(N - width + 1)) {
    result[i, 1:N] <- SSR.recursive(
      y,
      x,
      i,
      N,
      width
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
SSR.recursive <- function(y, x, beg, end, width = 2) {
  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }

  N <- nrow(y)

  beg <- max(beg, 1)
  end <- min(end, N)

  vecSSR <- rep(Inf, N)

  y0 <- .msub(y, beg:(beg + width - 1))
  x0 <- .msub(x, beg:(beg + width - 1))

  xx.inv <- solve(t(x0) %*% x0)
  .model <- OLS.reg(y0, x0)
  beta <- .model$coefficients
  residl <- .model$residuals
  vecSSR[beg + width - 1] <- sum(residl^2)

  for (step in (beg + width):end) {
    if (step > end) {
      break
    }

    z <- .msub(x, step)
    v <- drop(.msub(y, step) - z %*% beta)
    f <- drop(1 + z %*% xx.inv %*% t(z))

    vecSSR[step] <- vecSSR[step - 1] + v^2 / f

    beta <- beta + xx.inv %*% t(z) * v / f
    xx.inv <- xx.inv -
      (xx.inv %*% t(z) %*% z %*% xx.inv) / f
  }

  vecSSR
}

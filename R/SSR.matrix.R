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

  N <- nrow(y)

  result <- matrix(data = Inf, nrow = N, ncol = N)

  for (i in 1:(N - width + 1)) {
    result[i, 1:N] <- SSR.recursive(
      y, x, i, N, width
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

  N <- nrow(y)

  beg <- max(beg, 1)
  end <- min(end, N)

  result <- matrix(data = Inf, nrow = N, ncol = 1)

  y0 <- y[beg:(beg + width - 1), , drop = FALSE]
  x0 <- x[beg:(beg + width - 1), , drop = FALSE]

  xx.inv <- qr.solve(t(x0) %*% x0)
  .model <- OLS.reg(y0, x0)
  beta <- .model$coefficients
  residl <- .model$residuals
  rm(.model)

  result[beg + width - 1, 1] <- drop(t(residl) %*% residl)

  for (step in (beg + width):end) {
    if (step > end) break

    xl <- x[step, , drop = FALSE]
    residl <- drop(y[step, , drop = FALSE] - xl %*% beta)

    denom <- drop(1 + xl %*% xx.inv %*% t(xl))

    result[step, 1] <- result[step - 1, 1] + residl^2 / denom

    beta <- beta + xx.inv %*% t(xl) * residl
    xx.inv <- xx.inv -
      (xx.inv %*% t(xl) %*% xl %*% xx.inv) / denom
  }

  result
}

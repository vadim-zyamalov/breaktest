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

    return(result)
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

    if (beg < 1) beg <- 1
    if (nrow(y) < end) end <- nrow(y)

    n.obs <- nrow(y)

    result <- matrix(data = Inf, nrow = n.obs, ncol = 1)

    y.0 <- y[beg:(beg + width - 1), , drop = FALSE]
    x.0 <- x[beg:(beg + width - 1), , drop = FALSE]

    inv.XX.0 <- qr.solve(t(x.0) %*% x.0)
    tmp.OLS <- .estimate_ols(y.0, x.0)
    beta.0 <- tmp.OLS$beta
    resid.0 <- tmp.OLS$residuals
    rm(tmp.OLS)

    result[beg + width - 1, 1] <- drop(t(resid.0) %*% resid.0)

    for (step in (beg + width):end) {
        if (step > end) break

        step.x <- x[step, , drop = FALSE]

        step.resid <- y[step, , drop = FALSE] - step.x %*% beta.0
        step.resid <- drop(step.resid)

        denom <- 1 + step.x %*% inv.XX.0 %*% t(step.x)
        denom <- drop(denom)

        inv.XX.1 <- inv.XX.0 -
            (inv.XX.0 %*% t(step.x) %*% step.x %*% inv.XX.0) / denom

        beta.1 <- beta.0 + inv.XX.1 %*% t(step.x) * step.resid

        result[step, 1] <- result[step - 1, 1] + step.resid^2 / denom

        inv.XX.0 <- inv.XX.1

        beta.0 <- beta.1
    }

    return(result)
}

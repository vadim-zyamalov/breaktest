#' @title
#' Information criterions
#'
#' @details
#' Calculating the value of the following informational criterions:
#' * Akaike,
#' * Schwarz (Bayesian),
#' * Hannan-Quinn,
#' * Liu et al.
#'
#' @param resids Input residuals needed for estimating the values of
#' information criterions.
#' @param extra Number of extra parameters needed for estimating the punishment
#' term.
#' @param modification Whether the unit-root test modificaton is needed.
#' See Ng and Perron (2001) for further information.
#' @param alpha The coefficient \eqn{\alpha} of \eqn{y_{t-1}} in ADF model.
#' Needed only for criterion modification purposes.
#' @param y The vector of \eqn{y_{t-1}} in ADF model.
#' Needed only for criterion modification purposes.
#'
#' @return
#' The list of information criterions values.
#'
#' @references
#' Ng, Serena, and Pierre Perron. “Lag Length Selection and the Construction of
#' Unit Root Tests with Good Size and Power.”
#' Econometrica 69, no. 6 (2001): 1519–54.
#' https://doi.org/10.1111/1468-0262.00256.
#'
#' @export
info.criterion <- function(resids,
                           extra,
                           modification = FALSE,
                           alpha = 0,
                           y = NULL) {
    if (!is.matrix(resids)) resids <- as.matrix(resids)

    n.obs <- nrow(resids)

    if (modification) {
        s2 <- drop(t(resids) %*% resids) / n.obs
        tau <- (alpha ^ 2) * drop(t(y) %*% y) / s2
    } else {
        tau <- 0
    }

    log.RSS <- log(drop(t(resids) %*% resids) / n.obs)

    return(
        list(
            aic = log.RSS + 2 * (tau + extra) / n.obs,
            bic = log.RSS + (tau + extra) * log(n.obs) / n.obs,
            hq = log.RSS + 2 * (tau + extra) * log(log(n.obs)) / n.obs,
            lwz = log.RSS + 0.299 * (tau + extra) * (log(n.obs))^2.1
        )
    )
}

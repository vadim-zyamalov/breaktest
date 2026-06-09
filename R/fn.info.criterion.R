#' @export
HQIC <- function(obj, ...) UseMethod("HQIC")

#' @export
LWZ <- function(obj, ...) UseMethod("LWZ")


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
#' @param obj An object of one of the class `bt_ols`.
#' It should be pointed that `AR.reg` and `DOLS.reg` return object that belong
#' to subclasses of `bt_ols`, so they also can be used.
#' @param criterion One of the following string values (case insensitive):
#' * `"AIC"`,
#' * `"BIC"`,
#' * `"HQIC"`,
#' * `"LWZ"`.
#' @param alpha The coefficient \eqn{\alpha} of \eqn{y_{t-1}} in ADF model.
#' Needed only for criterion modification purposes,
#' see Ng and Perron (2001) for further information.
#' If NULL then no modification is applied.
#' @param y The vector of \eqn{y_{t-1}} in ADF model.
#' Needed only for criterion modification purposes.
#'
#' @return
#' The value of the desired information criterion.
#'
#' @references
#' Ng, Serena, and Pierre Perron. “Lag Length Selection and the Construction of
#' Unit Root Tests with Good Size and Power.”
#' Econometrica 69, no. 6 (2001): 1519–54.
#' https://doi.org/10.1111/1468-0262.00256.
#'
#' @export
info.criterions <- function(obj, criterion, alpha = NULL, y = NULL) {
  switch(
    toupper(criterion),
    AIC = AIC,
    BIC = BIC,
    HQIC = HQIC,
    LWZ = LWZ,
    stop("IC: unknown criterion'", criterion, "'")
  )(obj, alpha = alpha, y = y)
}

#' @rdname info.criterions
#' @exportS3Method
AIC.bt_ols <- function(obj, k = 2, alpha = NULL, y = NULL) {
  r <- na.omit(obj$residuals)
  N <- length(r)
  Nx <- length(obj$coefficients)

  s2 <- sum(r^2) / N

  tau <- if (!is.null(alpha) && !is.null(y)) {
    alpha^2 * sum(y^2) / s2
  } else {
    0
  }

  log(s2) + k * (Nx + tau) / N
}

#' @rdname info.criterions
#' @exportS3Method
BIC.bt_ols <- function(obj, alpha = NULL, y = NULL) {
  AIC(obj, k = log(nobs(obj)), alpha = alpha, y = y)
}

#' @rdname info.criterions
#' @exportS3Method
HQIC.bt_ols <- function(obj, alpha = NULL, y = NULL) {
  AIC(obj, k = 2 * log(log(nobs(obj))), alpha = alpha, y = y)
}

#' @rdname info.criterions
#' @exportS3Method
LWZ.bt_ols <- function(obj, alpha = NULL, y = NULL) {
  AIC(obj, k = 0.299 * (log(nobs(obj)))^2.1 * nobs(obj), alpha = alpha, y = y)
}

#' @exportS3Method
nobs.bt_ols <- function(obj) {
  length(na.omit(obj$residuals))
}

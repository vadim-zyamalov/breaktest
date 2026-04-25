#' @title
#' Custom OLS with extra information
#'
#' @description
#' Getting OLS estimates of betas, residuals, forecasted values and t-values.
#'
#' @param y A dependent variable.
#' @param x Explanatory variables.
#'
#' @return A list of:
#' * `beta`: estimates of coefficients,
#' * `resid`: estimated residuals,
#' * `predict`: forecasted values,
#' * `t.beta`: \eqn{t}-statistics for `beta`.
#'
#' @importFrom stats .lm.fit
#'
#' @keywords internal
OLS.reg <- function(y, x) {
  if (!is.matrix(y)) y <- as.matrix(y)
  if (!is.matrix(x)) x <- as.matrix(x)

  N <- length(y)

  yrows <- apply(y, 1, function(r) any(is.na(r)))
  xrows <- apply(x, 1, function(r) any(is.na(r)))
  rows <- !yrows & !xrows

  y <- y[rows, , drop = FALSE]
  x <- x[rows, , drop = FALSE]

  .model <- lm.fit(x, y)
  r <- .model$residuals
  cf <- .model$coefficients
  s.sq <- sum(r^2) / (nrow(x) - ncol(x))
  se.cf <- sqrt(diag(s.sq * qr.solve(t(x) %*% x)))
  t.beta <- cf / se.cf

  resid <- rep(NA, N)
  resid[rows] <- r

  result <- list(
    coefficients = as.matrix(cf),
    se.coefs = se.cf,
    t.stats = t.beta,
    residuals = resid,
    fitted.values = .model$fitted.values,
    endog = y,
    exog = x
  )

  class(result) <- "bt_ols"
  result
}


#' @title
#' Estimating DOLS regression for multiple known break points
#'
#' @param y A time series of interest.
#' @param x A matrix of explanatory stochastic regressors.
#' @param model A scalar or vector of break types:
#' * 1: for the break in const.
#' * 2: for the break in trend.
#' * 3: for the break in const and trend.
#' @param break.point An array of moments of structural breaks.
#' @param const,trend Whether a constant or trend are to be included.
#' @param k.lags,k.leads A number of lags and leads in DOLS regression.
#'
#' @return A list of:
#' * Estimates of coefficients,
#' * Estimates of residuals,
#' * A set of informational criterions values,
#' * \eqn{t}-statistics for the estimates of coefficients.
#'
#' @keywords internal
DOLS.many <- function(
  y,
  x,
  const = FALSE,
  trend = FALSE,
  break.type,
  break.point,
  break.coint = FALSE,
  n.lags,
  n.leads,
  ...
) {
  if (!is.matrix(y)) y <- as.matrix(y)
  if (is.null(x)) {
    stop("ERROR! DOLS.multiple: explanatory variables needed for DOLS")
  }
  if (!is.matrix(x)) x <- as.matrix(x)

  .vars_dols <- DOLS.mlt.regressors(
    y,
    x,
    const = FALSE,
    trend = FALSE,
    break.type,
    break.point,
    break.coint = FALSE,
    n.lags,
    n.leads
  )

  result <- OLS.reg(.vars_dols$yreg, .vars_dols$xreg)
  result$break.type <- break.type
  result$break.point <- break.point
  result$break.coint <- break.coint
  result$criterions <- info.criterions(result$residuals, ncol(.vars_dols$xreg))
  result$lags <- n.lags
  result$leads <- n.leads

  class(result) <- "bt_dols"
  result
}


#' @title
#' Custom GLS with extra information
#'
#' @description
#' Getting GLS estimates of betas, residuals, forecasted values and t-values.
#'
#' @param y A dependent variable.
#' @param x Explanatory variables.
#' @param c A coefficient for \eqn{\rho} calculation.
#'
#' @return A list of:
#' * `beta`: estimates of coefficients,
#' * `resid`: estimated residuals,
#' * `predict`: forecasted values,
#' * `t.beta`: \eqn{t}-statistics for `beta`.
#'
#' @keywords internal
GLS.reg <- function(y, z, c) {
  if (!is.matrix(y)) y <- as.matrix(y)
  if (!is.matrix(z)) z <- as.matrix(z)

  N <- nrow(y)
  Nc <- ncol(y)

  rho <- 1 + c / N

  y_hat <- y - rho * .lagn(y, 1)
  y_hat[1, ] <- y[1, ]

  z_hat <- z - rho * .lagn(z, 1)
  z_hat[1, ] <- z[1, ]

  betas <- NULL
  resids <- NULL
  fitted <- NULL
  t.betas <- NULL

  for (i in 1:Nc) {
    .model <- OLS.reg(y_hat[, i, drop = FALSE], z_hat)
    betas <- cbind(betas, .model$coefficients)
    t.betas <- cbind(t.betas, .model$t.stats)
    fitted <- cbind(fitted, z %*% .model$coefficients)
    resids <- cbind(resids, y[, i, drop = FALSE] - fitted)
  }

  result <- list(
    coefficients  = betas,
    t.stats       = drop(t.betas),
    residuals     = resids,
    fitted.values = fitted
  )
  class(result) <- "bt_gls"
  result
}


#' @title
#' Custom AR with extra information
#'
#' @param y A dependent variable.
#' @param x Explanatory variables.
#' @param max.lag The maximum number of lags.
#' @param criterion A criterion for lag number estimation.
#'
#' @return A list of:
#' * `beta`: estimates of coefficients,
#' * `residuals`: estimated residuals,
#' * `predict`: forecasted values,
#' * `t.beta`: \eqn{t}-statistics for `beta`,
#' * `lag`: estimated number of lags.
#'
#' @keywords internal
AR.reg <- function(
  y,
  x,
  max.lag,
  criterion = "aic"
) {
  if (!is.null(criterion)) {
    if (!criterion %in% c("bic", "aic", "lwz", "hq")) {
      stop("ERROR! .AR: Unknown criterion")
    }
  }

  if (!is.matrix(y)) y <- as.matrix(y)
  if (!is.null(x) && !is.matrix(x)) x <- as.matrix(x)

  if (!is.null(x)) {
    Nx <- ncol(x)
    .rhs <- x
  } else {
    Nx <- 0
    .rhs <- NULL
  }

  .rhs <- cbind(
    .rhs,
    apply(as.array(1:max.lag), 1, function(l) .lagn(y, l))
  )

  yrows <- apply(y, 1, function(r) any(is.na(r)))
  xrows <- apply(.rhs, 1, function(r) any(is.na(r)))
  rows <- !yrows & !xrows

  .lhs <- y[rows, , drop = FALSE]
  .rhs <- .rhs[rows, , drop = FALSE]

  if (is.null(criterion)) {
    .lag <- max.lag
    result <- OLS.reg(.lhs, .rhs[, 1:(Nx + .lag), drop = FALSE])
  } else {
    .lag <- 0

    result <- NULL
    .ic <- Inf

    for (l in 0:max.lag) {
      .model <- OLS.reg(.lhs, .rhs[, 1:(Nx + l), drop = FALSE])
      .model_ic <- info.criterions(.model$residuals, l)[[criterion]]

      if (.model_ic < .ic) {
        .ic <- .model_ic
        .lag <- l
        result <- .model
      }
    }
  }

  result$lag <- .lag
  result$criterion <- .ic
  result$criterion.name <- criterion

  class(result) <- "bt_ar"
  result
}


#' @title
#' Nadaraya–Watson kernel regression.
#'
#' @param y A dependent variable.
#' @param x Explanatory variables.
#' @param h A bandwidth parameter.
#' @param kernel Needed kernel, currently only `unif` and `gauss`:
#' * `unif`: \eqn{K(x) = \left\{\begin{array}{ll}
#' 1 & \frac{|x - x_i|}{h} \leq 1 \\
#' 0 & \textrm{otherwize}
#' \end{array}\right.}
#' * `gauss`: \eqn{\Phi(\frac{x - x_i}{h})}
#'
#' @return A list of arguments as well as the estimated coefficient vector and
#' residuals.
#'
#' @references
#' Harvey, David I., S. Leybourne, Stephen J., and Yang Zu.
#' “Nonparametric Estimation of the Variance Function
#' in an Explosive Autoregression Model.”
#' School of Economics. University of Nottingham, 2022.
#'
#' @keywords internal
NW.reg <- function(
  y,
  x,
  h,
  kernel = "unif"
) {
  if (!kernel %in% c("unif", "gauss")) {
    stop("ERROR! NW.reg: unknown kernel")
  }

  N <- length(y)

  rho <- numeric(N)

  for (k in 1:N) {
    .w <- NW.kernel(k, (1:N) / N, h, kernel)
    rho[k] <- sum(x * .w * y) / sum(x * .w * x)
  }

  list(
    my      = y,
    mx      = x,
    h       = h,
    kernel  = kernel,
    rr1.est = rho,
    u.hat   = y - rho * x
  )
}


#' @title
#' Nadaraya–Watson kernel volatility estimation
#'
#' @param e A series of interest.
#' @param h A bandwidth parameter.
#' @param kernel Needed kernel, currently only `unif` and `gauss`:
#' * `unif`: \eqn{K(x) = \left\{\begin{array}{ll}
#' 1 & \frac{|x - x_i|}{h} \leq 1 \\
#' 0 & \textrm{otherwize}
#' \end{array}\right.}
#' * `gauss`: \eqn{\Phi(\frac{x - x_i}{h})}
#'
#' @return A list of arguments as well as the estimated omega and s.e.
#'
#' @references
#' Cavaliere, Giuseppe, Peter C. B. Phillips, Stephan Smeekes,
#' and A. M. Robert Taylor.
#' “Lag Length Selection for Unit Root Tests in the Presence
#' of Nonstationary Volatility.”
#' Econometric Reviews 34, no. 4 (April 21, 2015): 512–36.
#' https://doi.org/10.1080/07474938.2013.808065.
#'
#' Harvey, David I., S. Leybourne, Stephen J., and Yang Zu.
#' “Nonparametric Estimation of the Variance Function
#' in an Explosive Autoregression Model.”
#' School of Economics. University of Nottingham, 2022.
#'
#' @keywords internal
NW.variance <- function(
  e,
  h,
  kernel = "unif"
) {
  if (!kernel %in% c("unif", "gauss")) {
    stop("ERROR! NW.variance: unknown kernel")
  }

  .e <- na.omit(e)
  N <- length(.e)
  omega2 <- numeric(N)

  for (k in 1:N) {
    .w <- NW.kernel(k, (1:N) / N, h, kernel)
    omega2[k] <- sum(.w * e^2) / sum(.w)
  }

  dN <- length(e) - length(.e)
  omega2 <- c(rep(NA, dN), omega2)

  list(
    me       = e,
    h        = h,
    kernel   = kernel,
    omega.sq = omega2,
    se       = sqrt(omega2)
  )
}


#' @title
#' LOO-CV for h in Nadaraya–Watson kernel regression.
#'
#' @param y A dependent variable.
#' @param x An explanatory variable.
#' @param kernel Needed kernel, currently only `unif` and `gauss`:
#' * `unif`: \eqn{K(x) = \left\{\begin{array}{ll}
#' 1 & \frac{|x - x_i|}{h} \leq 1 \\
#' 0 & \textrm{otherwize}
#' \end{array}\right.}
#' * `gauss`: \eqn{\Phi(\frac{x - x_i}{h})}
#'
#' @references
#' Harvey, David I., S. Leybourne, Stephen J., and Yang Zu.
#' “Nonparametric Estimation of the Variance Function
#' in an Explosive Autoregression Model.”
#' School of Economics. University of Nottingham, 2022.
#'
#' @return A list of arguments as well as the estimated bandwidth `h`.
#'
#' @keywords internal
NW.bandwidth <- function(y, x, kernel = "unif") {
  if (!kernel %in% c("unif", "gauss")) {
    stop("ERROR! NW.bandwidth: unknown kernel")
  }

  .y <- na.omit(y)
  .x <- x[!is.na(y)]
  N <- length(.y)

  h_candidates <- seq(N^(-0.5), N^(-0.3), by = 0.01)
  rss <- Inf
  h <- NULL

  for (.h in h_candidates) {
    rho <- numeric(N)
    for (k in 1:N) {
      .w <- NW.kernel(k, (1:N) / N, .h, kernel)
      .w[k] <- 0
      rho[k] <- sum(.x * .w * .y) / sum(.x * .w * .x)
    }

    .rss <- sum((y - rho * x)^2)
    if (.rss < rss) {
      rss <- .rss
      h <- .h
    }
  }

  list(
    my     = y,
    mx     = x,
    kernel = kernel,
    h      = h
  )
}


#' @title
#' Returning a kernel value for a specific point
#'
#' @param i An index of the base point.
#' @param x A series for kernel calculations.
#' @param h A bandwidth parameter.
#' @param kernel Needed kernel, currently only `unif` and `gauss`:
#' * `unif`: \eqn{K(x) = \left\{\begin{array}{ll}
#' 1 & \frac{|x - x_i|}{h} \leq 1 \\
#' 0 & \textrm{otherwize}
#' \end{array}\right.}
#' * `gauss`: \eqn{\Phi(\frac{x - x_i}{h})}
#'
#' @importFrom stats pnorm
#'
#' @keywords internal
NW.kernel <- function(i,
                      x,
                      h,
                      kernel = "unif") {
  switch(kernel,
    unif  = ifelse((abs((x - x[i]) / h) <= 1), 1, 0),
    gauss = pnorm((x - x[i]) / h)
  )
}

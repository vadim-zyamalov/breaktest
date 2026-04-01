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
.OLS <- function(y, x) {
  if (!is.matrix(y)) y <- as.matrix(y)
  if (!is.matrix(x)) x <- as.matrix(x)

  .model <- lm(y ~ x - 1)
  s.sq <- drop(t(.model$residuals) %*% .model$residuals) /
    (nrow(x) - ncol(x))
  t.beta <- .model$coefficients / sqrt(diag(s.sq * qr.solve(t(x) %*% x)))

  list(
    beta = as.matrix(.model$coefficients),
    residuals = .model$residuals,
    predict = .model$fitted.values,
    t.beta = t.beta
  )
}


#' @title
#' Estimating DOLS regression for a single known break point
#'
#' @param y A time series of interest.
#' @param x A matrix of explanatory stochastic regressors.
#' @param model See Carrion-i-Silvestre and Sansó (2006)
#' * 1: for model An,
#' * 2: for model A,
#' * 3: for model B,
#' * 4: for model C,
#' * 5: for model D,
#' * 6: for model E.
#' @param bp A position of the break point.
#' @param n.lags,n.leads A number of lags and leads in DOLS regression.
#'
#' @return A list of:
#' * Estimates of coefficients,
#' * Estimates of residuals,
#' * A value of BIC,
#' * \eqn{t}-statistics for the estimates of coefficients.
#'
#' @references
#' Carrion-i-Silvestre, Josep Lluís, and Andreu Sansó.
#' “Testing the Null of Cointegration with Structural Breaks.”
#' Oxford Bulletin of Economics and Statistics 68, no. 5 (October 2006): 623–46.
#' https://doi.org/10.1111/j.1468-0084.2006.00180.x.
#'
#' @keywords internal
.DOLS.single <- function(y,
                         x,
                         model,
                         bp,
                         n.lags,
                         n.leads) {
  if (is.null(x)) {
    stop("ERROR! DOLS.single: explanatory variables needed for DOLS")
  }
  if (!is.matrix(y)) y <- as.matrix(y)
  if (!is.matrix(x)) x <- as.matrix(x)

  N <- nrow(y)
  rows <- (n.lags + 2):(N - n.leads)
  dmodel <- c(1, 2, 3, 4, 1, 4)

  d.x <- .diffn(x)

  Ld.x <- if (n.lags != 0) {
    apply(as.array(1:n.lags), 1, function(l) .lagn(d.x, l))
  } else {
    NULL
  }

  Fd.x <- if (n.leads != 0) {
    apply(as.array(1:n.leads), 1, function(l) .lagn(d.x, -l))
  } else {
    NULL
  }

  deter <- if (model != 0) {
    trend.kpss.single(dmodel[model], N, bp)
  } else {
    NULL
  }

  xdu <- if (model %in% c(5, 6)) {
    sweep(x, 1, deter[, 2, drop = FALSE], `*`)
  } else {
    NULL
  }

  xreg <- cbind(deter, x, xdu, d.x, Ld.x, Fd.x)[rows, , drop = FALSE]
  y <- y[rows, 1, drop = FALSE]

  .res_ols <- .OLS(y, xreg)

  bic <- log(drop(t(.res_ols$residuals) %*% .res_ols$residuals) / nrow(xreg)) +
    ncol(xreg) * log(nrow(xreg)) / nrow(xreg)

  list(
    beta      = .res_ols$beta,
    residuals = .res_ols$residuals,
    bic       = bic,
    t.beta    = .res_ols$t.beta
  )
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
.DOLS.multiple <- function(
  y,
  x,
  model,
  bp,
  const = FALSE,
  trend = FALSE,
  n.lags,
  n.leads
) {
  if (!is.matrix(y)) y <- as.matrix(y)
  if (is.null(x)) {
    stop("ERROR! DOLS.multiple: explanatory variables needed for DOLS")
  }
  if (!is.matrix(x)) x <- as.matrix(x)

  .vars_dols <- variables.dols.multiple(
    y, x,
    model, bp,
    const, trend,
    n.lags, n.leads
  )

  .model <- .OLS(.vars_dols$yreg, .vars_dols$xreg)

  criterions <- .ic.values(.model$residuals, ncol(.vars_dols$xreg))

  list(
    beta       = .model$beta,
    residuals  = .model$residuals,
    criterions = criterions,
    t.beta     = .model$t.beta
  )
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
.GLS <- function(y, z, c) {
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
    .model <- .OLS(y_hat[, i, drop = FALSE], z_hat)
    betas <- cbind(betas, .model$beta)
    t.betas <- cbind(t.betas, .model$t.beta)
    fitted <- cbind(fitted, z %*% .model$beta)
    resids <- cbind(resids, y[, i, drop = FALSE] - fitted)
  }

  list(
    beta = betas,
    residuals = resids,
    predict = fitted,
    t.beta = drop(t.betas)
  )
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
.AR <- function(
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

  N <- nrow(y)
  rows <- (1 + max.lag):N

  .lhs <- y[rows, , drop = FALSE]

  if (!is.null(x)) {
    Nx <- ncol(x)
    .rhs <- x
  } else {
    Nx <- 0
    .rhs <- NULL
  }

  .rhs <- cbind(
    .rhs,
    apply(as.array(1:max.lag), 1, function(x) .lagn(y, l))
  )[rows, , drop = FALSE]

  if (is.null(criterion)) {
    .lag <- max.lag
    .model <- .OLS(.lhs, .rhs[, 1:(Nx + .lag), drop = FALSE])
    .beta <- .model$beta
    .resid <- .model$residuals
    .predict <- .model$predict
    .t_beta <- .model$t.beta
  } else {
    .lag <- 0

    .model <- NULL
    .beta <- NULL
    .resid <- NULL
    .predict <- NULL
    .t_beta <- NULL
    .ic <- Inf

    for (l in 0:max.lag) {
      .model <- .OLS(.lhs, .rhs[, 1:(Nx + l), drop = FALSE])
      .model_ic <- .ic.values(.model$residuals, l)[[criterion]]

      if (.model_ic < .ic) {
        .ic <- .model_ic
        .beta <- .model$beta
        .resid <- .model$residuals
        .predict <- .model$predict
        .t_beta <- .model$t.beta
        .lag <- l
      }
    }
  }

  list(
    beta      = .beta,
    residuals = .resid,
    predict   = .predict,
    t.beta    = .t_beta,
    lag       = .lag,
    criterion = .ic
  )
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
.NW.reg <- function(
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
    .w <- .NW.kernel(k, (1:N) / N, h, kernel)
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
.NW.variance <- function(
  e,
  h,
  kernel = "unif"
) {
  if (!kernel %in% c("unif", "gauss")) {
    stop("ERROR! NW.variance: unknown kernel")
  }

  N <- length(e)

  omega2 <- numeric(N)

  for (k in 1:N) {
    .w <- .NW.kernel(k, (1:N) / N, h, kernel)
    omega2[k] <- sum(.w * e^2) / sum(.w)
  }

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
.NW.bandwidth <- function(y, x, kernel = "unif") {
  if (!kernel %in% c("unif", "gauss")) {
    stop("ERROR! NW.bandwidth: unknown kernel")
  }

  N <- length(y)

  h_candidates <- seq(N^(-0.5), N^(-0.3), by = 0.01)
  rss <- Inf
  h <- NULL

  for (.h in h_candidates) {
    rho <- numeric(N)
    for (k in 1:N) {
      .w <- .NW.kernel(k, (1:N) / N, .h, kernel)
      .w[k] <- 0
      rho[k] <- sum(x * .w * y) / sum(x * .w * x)
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
.NW.kernel <- function(i,
                       x,
                       h,
                       kernel = "unif") {
  switch(kernel,
    unif  = if (abs((x - x[i]) / h) <= 1) 1 else 0,
    gauss = pnorm((x - x[i]) / h)
  )
}

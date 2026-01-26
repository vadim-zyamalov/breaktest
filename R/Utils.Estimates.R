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
.estimate_ols <- function(y, x) {
  if (!is.matrix(y)) y <- as.matrix(y)
  if (!is.matrix(x)) x <- as.matrix(x)

  .model <- .lm.fit(x, y)
  s_2 <- drop(t(.model$residuals) %*% .model$residuals) /
    (nrow(x) - ncol(x))
  t.beta <- .model$coefficients / sqrt(diag(s_2 * qr.solve(t(x) %*% x)))

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
#' @param break.point A position of the break point.
#' @param k.lags,k.leads A number of lags and leads in DOLS regression.
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
.estimate_dols_single <- function(y,
                                  x,
                                  model,
                                  tb,
                                  k_lags,
                                  k_leads) {
  if (is.null(x)) {
    stop("ERROR! Explanatory variables needed for DOLS")
  }
  if (!is.matrix(y)) y <- as.matrix(y)
  if (!is.matrix(x)) x <- as.matrix(x)

  n_obs <- nrow(y)

  .diff_x <- diff(x)
  .diff_x_lag <- .diff_x
  .diff_x_lead <- .diff_x

  for (i in 1:k_lags) {
    .diff_x_lag <- cbind(
      .diff_x_lag,
      lagn(.diff_x, i)
    )
  }

  for (i in 1:k_leads) {
    .diff_x_lead <- cbind(
      .diff_x_lead,
      lagn(.diff_x, -i)
    )
  }

  if (k_lags != 0 && k_leads != 0) {
    lags <- .diff_x_lag
    leads <- .diff_x_lead[, (ncol(x) + 1):(ncol(.diff_x_lead)), drop = FALSE]
    lags_leads <- cbind(lags, leads)
    lags_leads <-
      lags_leads[(k_lags + 1):(n_obs - 1 - k_leads), , drop = FALSE]
  } else if (k_lags != 0 && k_leads == 0) {
    lags <- .diff_x_lag
    lags_leads <- lags[(k_lags + 1):(n_obs - 1), , drop = FALSE]
  } else if (k_lags == 0 && k_leads != 0) {
    lags <- .diff_x_lag
    leads <- .diff_x_lead[, (ncol(x) + 1):(ncol(.diff_x_lead)), drop = FALSE]
    lags_leads <- cbind(lags, leads)
    lags_leads <- lags_leads[1:(n_obs - 1 - k_leads), , drop = FALSE]
  } else if (k_lags == 0 && k_leads == 0) {
    lags_leads <- .diff_x_lag
  }

  if (model == 0) {
    xreg <- cbind(
      x[(k_lags + 2):(n_obs - k_leads), , drop = FALSE],
      lags_leads
    )
  } else if (model >= 1 && model <= 4) {
    deter <- trend_kpss_single(model, n_obs, tb)
    xreg <- cbind(
      deter[(k_lags + 2):(n_obs - k_leads), , drop = FALSE],
      x[(k_lags + 2):(n_obs - k_leads), , drop = FALSE],
      lags_leads
    )
  } else if (model == 5) {
    deter <- trend_kpss_single(1, n_obs, tb)
    xdu <- sweep(x, 1, deter[, 2, drop = FALSE], `*`)
    xreg <- cbind(
      deter[(k_lags + 2):(n_obs - k_leads), , drop = FALSE],
      x[(k_lags + 2):(n_obs - k_leads), , drop = FALSE],
      xdu[(k_lags + 2):(n_obs - k_leads), , drop = FALSE],
      lags_leads
    )
  } else if (model == 6) {
    deter <- trend_kpss_single(4, n_obs, tb)
    xdu <- sweep(x, 1, deter[, 2, drop = FALSE], `*`)
    xreg <- cbind(
      deter[(k_lags + 2):(n_obs - k_leads), , drop = FALSE],
      x[(k_lags + 2):(n_obs - k_leads), , drop = FALSE],
      xdu[(k_lags + 2):(n_obs - k_leads), , drop = FALSE],
      lags_leads
    )
  }

  .res_ols <- .estimate_ols(
    y[(k_lags + 2):(n_obs - k_leads), 1, drop = FALSE],
    xreg
  )

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
.estimate_dols_multiple <- function(
  y,
  x,
  model,
  break_point,
  const = FALSE,
  trend = FALSE,
  n_lags,
  n_leads
) {
  if (!is.matrix(y)) y <- as.matrix(y)
  if (is.null(x)) {
    stop("ERROR! Explanatory variables needed for DOLS")
  }
  if (!is.matrix(x)) x <- as.matrix(x)

  .vars_dols <- variables_dols_multiple(
    y, x,
    model, break_point,
    const, trend,
    n_lags, n_leads
  )

  .model <- .estimate_ols(.vars_dols$yreg, .vars_dols$xreg)

  criterions <- .info_criterion(.model$residuals, ncol(.vars_dols$xreg))

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
.estimate_gls <- function(y, z, c) {
  if (!is.matrix(y)) y <- as.matrix(y)
  if (!is.matrix(z)) z <- as.matrix(z)

  n_obs <- nrow(y)
  n_var <- ncol(y)

  rho <- 1 + c / n_obs

  y_hat <- y - rho * lagn(y, 1)
  y_hat[1, ] <- y[1, ]

  z_hat <- z - rho * lagn(z, 1)
  z_hat[1, ] <- z[1, ]

  betas <- NULL
  resids <- NULL
  fitted <- NULL
  t.betas <- NULL

  for (i in 1:n_var) {
    .model <- .estimate_ols(y_hat[, i, drop = FALSE], z_hat)
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
.estimate_ar <- function(
  y,
  x,
  max_lag,
  criterion = "aic"
) {
  if (!is.null(criterion)) {
    if (!criterion %in% c("bic", "aic", "lwz", "hq")) {
      stop("WARNING! Unknown criterion, none is used")
    }
  }

  if (!is.matrix(y)) y <- as.matrix(y)
  n_obs <- nrow(y)
  .lhs <- y[(1 + max_lag):n_obs, , drop = FALSE]

  if (!is.null(x)) {
    if (!is.null(x) && !is.matrix(x)) x <- as.matrix(x)
    k <- ncol(x)
    .rhs <- x[(1 + max_lag):n_obs, , drop = FALSE]
  } else {
    k <- 0
    .rhs <- NULL
  }

  for (l in 1:max_lag) {
    if (l <= max_lag) {
      .rhs <- cbind(
        .rhs,
        lagn(y, l)[(1 + max_lag):n_obs, , drop = FALSE]
      )
    }
  }

  if (is.null(criterion)) {
    .lag <- max_lag
    .model <- .estimate_ols(.lhs, .rhs[, 1:(k + .lag), drop = FALSE])
    .beta <- .model$beta
    .resid <- .model$residuals
    .predict <- .model$predict
    .t_beta <- .model$t.beta
  } else {
    .lag <- 0

    if (!is.null(x)) {
      .model <- .estimate_ols(.lhs, .rhs[, 1:k, drop = FALSE])
      .beta <- .model$beta
      .resid <- .model$residuals
      .predict <- .model$predict
      .t_beta <- .model$t.beta
      .ic <- log(drop(t(.resid) %*% .resid) / (n_obs - max_lag))
    } else {
      .ic <- Inf
    }

    for (l in 1:max_lag) {
      if (l <= max_lag) {
        .model <- .estimate_ols(.lhs, .rhs[, 1:(k + l), drop = FALSE])
        .model_ic <- .info_criterion(.model$residuals, l)[[criterion]]

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
.estimate_nw <- function(
  y,
  x,
  h,
  kernel = "unif"
) {
  if (!kernel %in% c("unif", "gauss")) {
    stop("WARNING! Unknown kernel, unif is used instead")
  }

  n_obs <- length(y)

  rho <- rep(0, n_obs)
  for (k in 1:n_obs) {
    .w <- .kernel_nw(k, (1:n_obs) / n_obs, h, kernel)
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
.volatility_nw <- function(
  e,
  h,
  kernel = "unif"
) {
  if (!kernel %in% c("unif", "gauss")) {
    stop("WARNING! Unknown kernel, unif is used instead")
  }

  n_obs <- length(e)

  omega2 <- rep(0, n_obs)

  for (k in 1:n_obs) {
    .w <- .kernel_nw(k, (1:n_obs) / n_obs, h, kernel)
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
.bandwidth_nw <- function(y, x, kernel = "unif") {
  if (!kernel %in% c("unif", "gauss")) {
    stop("WARNING! Unknown kernel, unif is used instead")
  }

  n_obs <- length(y)

  h_candidates <- seq(n_obs^(-0.5), n_obs^(-0.3), by = 0.01)
  rss <- Inf

  for (.h in h_candidates) {
    rho <- rep(0, n_obs)
    for (k in 1:n_obs) {
      .w <- .kernel_nw(k, (1:n_obs) / n_obs, .h, kernel)
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
.kernel_nw <- function(i,
                       x,
                       h,
                       kernel = "unif") {
  switch(kernel,
    unif  = if (abs((x - x[i]) / h) <= 1) 1 else 0,
    gauss = pnorm((x - x[i]) / h)
  )
}

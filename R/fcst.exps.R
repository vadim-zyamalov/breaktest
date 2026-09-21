#' @title
#' Data-driven exponential discounting
#'
#' @description
#' Select the exponential weight using rolling-origin CV, Laplace CV, or
#' the local-level maximum-likelihood and trade-off rules of Hirano-Wright.
#' Regression fits reuse [OLS.reg] and [WLS.reg].
#'
#' @param y,x,x.new An aligned one-step regression, as in [fcst.forecast].
#' @param method "cv", "cvl", "ml" or "tradeoff".
#' @param gamma.grid Candidate finite-sample discount factors in `[0,1]`.
#' @param eta.grid Alternative local-parameter grid: gamma = 1 - eta/T,
#' with T the original history length. Supply only one grid. The default
#' is 301 equally spaced eta values from 0 to min(30, 0.99*T).
#' @param omega First inner forecast origin. Defaults to the last 10 percent
#' of the history, as in [fcst.window.cv].
#' @param prior For CVL, "midpoint" implements the uniform prior on the
#' continuous exponential kernel's median location in `[0.5,0.95]`, as in
#' Section 3.2 of the paper, integrated over grid cells. "uniform" assigns
#' equal mass to equally spaced eta candidates, as in mcprog.m. Irregular
#' grids are rejected for this prior. A truncated grid truncates the prior.
#' @param sigma.sq A fixed estimate of the equation-error variance for CVL.
#' By default use the unweighted direct-regression residual mean square.
#' It must estimate the same variance as in the source criterion; it is not
#' an arbitrary temperature parameter. Zero uses the limiting posterior on
#' criterion minima. This is not a posterior from
#' a fully specified likelihood for dependent forecast errors.
#'
#' @details
#' CVL averages eta and refits at the resulting gamma. It does not average
#' candidate forecasts. CV scores all candidates on the same usable origins;
#' a candidate failing on any such origin is inadmissible. Exact CV ties
#' prefer the largest gamma. Adaptive methods require a regular retained
#' sample without internal missing periods.
#'
#' ML and tradeoff require a complete, equally spaced, constant-only model.
#' They fit the exact zero-mean Gaussian MA(1) likelihood of diff(y), with
#' MA coefficient -gamma in `[-1,0]`, the invertible local-level branch.
#' The profile need not be unimodal: all local minima detected on a dense
#' grid, supplemented near the unit-MA boundary, are refined and compared
#' with both boundaries. The reported optimisation diagnostics describe
#' this deterministic numerical search; it is not a proof of uniqueness.
#' Unlike unconstrained MATLAB
#' optimisation, a positive MA coefficient cannot create negative weights.
#' The fitted variances are sigma.error^2 = gamma*sigma.MA^2 and
#' sigma.level^2 = (1-gamma)^2*sigma.MA^2. The local magnitude for (A.6) is
#' mu = T*sigma.level/sigma.error; T*(1-gamma) in MATLAB is its local
#' approximation. The trade-off rule minimises (A.6) on eta in `[0,T]`.
#'
#' For direct multiple horizons use
#' [fcst.forecast] with method="exps" and gamma.method set to one of these
#' methods. ML/tradeoff support h=1 or iterated constant-level forecasts;
#' direct h>1 and nonconstant regressions require CV/CVL. Iterated CV/CVL for
#' h>1 and nested predictor selection are not part of the source procedure.
#'
#' @return A bt_fcstExpS with forecast, coefficients, normalised weights,
#' gamma, eta, method, fit, and the CV profile/posterior or ML diagnostics.
#'
#' @references
#' Hirano, K. and J. H. Wright (2022). "Analyzing Cross-Validation for
#' Forecasting with Structural Instability." Journal of Econometrics 226,
#' 139-154. Sections 2.3, 3.2, 4 and equation (A.6).
#' https://doi.org/10.1016/j.jeconom.2020.10.009.
#'
#' @importFrom stats optimize
#' @export
fcst.ExpS <- function(y, x = NULL, x.new = NULL,
                      method = c("cv", "cvl", "ml", "tradeoff"),
                      gamma.grid = NULL, eta.grid = NULL, omega = NULL,
                      prior = c("midpoint", "uniform"), sigma.sq = NULL) {
  method <- match.arg(method)
  prior <- match.arg(prior)
  fcst.forecast(y, x, x.new, method = "exps", gamma.method = method,
                 gamma.grid = gamma.grid, eta.grid = eta.grid, omega = omega,
                 prior = prior, sigma.sq = sigma.sq)$details
}


.fcst.exps.grid <- function(N, gamma.grid = NULL, eta.grid = NULL) {
  if (!is.null(gamma.grid) && !is.null(eta.grid)) {
    stop("ERROR! fcst.ExpS: supply gamma.grid or eta.grid, not both")
  }
  if (!is.null(gamma.grid) && (!is.numeric(gamma.grid) ||
      !is.null(dim(gamma.grid)) || !length(gamma.grid) ||
      any(!is.finite(gamma.grid)) || any(gamma.grid < 0 | gamma.grid > 1))) {
    stop("ERROR! fcst.ExpS: gamma.grid must be a numeric vector in [0,1]")
  }
  if (is.null(eta.grid)) {
    eta.grid <- if (is.null(gamma.grid)) seq(0, min(30, 0.99 * N), length.out = 301L)
      else (1 - gamma.grid) * N
  }
  if (!is.numeric(eta.grid) || !is.null(dim(eta.grid)) || !length(eta.grid) ||
      any(!is.finite(eta.grid)) || any(eta.grid < 0 | eta.grid > N)) {
    stop("ERROR! fcst.ExpS: eta must be in [0,T] and gamma in [0,1]")
  }
  eta <- sort(unique(eta.grid))
  data.frame(eta = eta, gamma = 1 - eta / N)
}


.fcst.exps.prior <- function(eta, prior) {
  prior <- match.arg(prior, c("midpoint", "uniform"))
  n <- length(eta)
  if (prior == "uniform") {
    if (n > 2L && any(abs(diff(eta) - mean(diff(eta))) >
                         sqrt(.Machine$double.eps) * max(1, max(abs(eta))))) {
      stop("ERROR! fcst.ExpS: uniform prior requires an equally spaced eta grid")
    }
    return(rep(1 / n, n))
  }
  midpoint <- function(z) {
    out <- numeric(length(z))
    small <- z < 0.01
    out[small] <- 0.5 + z[small] / 8 - z[small]^3 / 192 + z[small]^5 / 2880
    out[!small] <- 1 + (log1p(exp(-z[!small])) - log(2)) / z[!small]
    out
  }
  if (n == 1L) {
    if (midpoint(eta) > 0.95) stop("ERROR! fcst.ExpS: grid has no mass under the midpoint prior")
    return(1)
  }
  edges <- c(eta[1], (eta[-1] + eta[-n]) / 2, eta[n])
  mass <- diff(pmin(1, pmax(0, (midpoint(edges) - 0.5) / 0.45)))
  if (!any(mass > 0)) stop("ERROR! fcst.ExpS: grid has no mass under the midpoint prior")
  mass / sum(mass)
}


# Equation (A.6), divided by K*sigma^2. A Taylor expansion and negative
# exponentials avoid cancellation at zero and overflow for large eta.
.fcst.exps.risk <- function(eta, mu) {
  small <- eta < 0.05
  bias <- variance <- numeric(length(eta))
  z <- eta[small]
  bias[small] <- 1 / 3 - z / 12 + z^2 / 180 + z^3 / 720 -
    z^4 / 5040 - z^5 / 30240 + z^6 / 151200
  variance[small] <- 1 + z^2 / 12 - z^4 / 720 + z^6 / 30240
  z <- eta[!small]
  a <- -expm1(-z)
  bias[!small] <- (1 - 4 * exp(-z) + (2 * z + 3) * exp(-2 * z)) / (2 * z * a^2)
  variance[!small] <- z / (2 * tanh(z / 2))
  mu^2 * bias + variance
}


# Exact MA(1) profile via tridiagonal LDL', O(T) per gamma. A vector of
# gamma values is evaluated together, using O(length(gamma)) storage. The
# determinant and quadratic form use the same factorisation even at gamma=1.
.fcst.ma1.profile <- function(gamma, dy) {
  n <- length(dy)
  diagonal <- 1 + gamma^2
  transformed <- rep(dy[1], length(gamma))
  quadratic <- transformed^2 / diagonal
  log.det <- log(diagonal)
  for (i in seq.int(2L, n)) {
    transformed <- dy[i] + gamma / diagonal * transformed
    diagonal <- 1 + gamma^2 - gamma^2 / diagonal
    quadratic <- quadratic + transformed^2 / diagonal
    log.det <- log.det + log(diagonal)
  }
  sigma.sq <- quadratic / n
  list(objective = n * log(sigma.sq) + log.det, sigma.sq = sigma.sq)
}


.fcst.exps.ml <- function(y, method) {
  dy <- diff(as.numeric(y))
  N <- length(y)
  if (all(dy == 0)) {
    return(list(gamma = 1, eta = 0, gamma.ml = 1, mu = 0, sigma.error.sq = 0,
                 sigma.level.sq = 0, sigma.MA.sq = 0, risk = NULL,
                 optimization = NULL, boundary = "constant series"))
  }
  likelihood <- function(g) .fcst.ma1.profile(g, dy)$objective
  # A single Brent search can miss a better local maximum of the MA(1)
  # likelihood. Cover the entire invertible branch and refine every basin,
  # including the O(1/T) neighbourhood of gamma=1 relevant to the paper.
  boundary.grid <- 10^seq(-12, 0, length.out = 129L)
  grid <- sort(unique(c(seq(0, 1, length.out = 1025L),
                        boundary.grid, 1 - boundary.grid)))
  values <- .fcst.ma1.profile(grid, dy)$objective
  if (any(!is.finite(values))) {
    stop("ERROR! fcst.ExpS: nonfinite MA profile; rescale the response")
  }
  middle <- seq.int(2L, length(grid) - 1L)
  minima <- middle[values[middle] <= values[middle - 1L] &
                    values[middle] <= values[middle + 1L]]
  refined <- vapply(minima, function(i) {
    optimize(likelihood, grid[c(i - 1L, i + 1L)], tol = 1e-10)$minimum
  }, 0.0)
  candidates <- unique(c(1, 0, grid[which.min(values)], refined))
  objectives <- vapply(candidates, likelihood, 0.0)
  gamma.ml <- candidates[which.min(objectives)]
  optimization <- list(method = "dense grid and refinement of every detected basin",
                        grid.size = length(grid), candidates = candidates,
                        objective = objectives, n.basins = length(minima))
  scale <- .fcst.ma1.profile(gamma.ml, dy)$sigma.sq
  sigma.error <- gamma.ml * scale
  sigma.level <- (1 - gamma.ml)^2 * scale
  mu <- if (gamma.ml == 0) Inf else N * (1 - gamma.ml) / sqrt(gamma.ml)
  eta <- (1 - gamma.ml) * N
  risk <- NULL
  if (method == "tradeoff" && is.finite(mu)) {
    risk.function <- function(z) .fcst.exps.risk(z, mu)
    opt <- optimize(risk.function, c(0, N), tol = 1e-8)
    candidates <- c(0, opt$minimum, N)
    eta <- candidates[which.min(vapply(candidates, risk.function, 0.0))]
    risk <- risk.function(eta)
  }
  list(gamma = 1 - eta / N, eta = eta, gamma.ml = gamma.ml, mu = mu,
       sigma.error.sq = sigma.error, sigma.level.sq = sigma.level,
       sigma.MA.sq = scale, risk = risk, optimization = optimization,
       boundary = if (gamma.ml == 1) "no level innovations" else
         if (gamma.ml == 0) "no observation noise" else "interior")
}


.fcst.exps.horizon <- function(history, data, h, type, x.future, options, row) {
  method <- match.arg(options$gamma.method, c("cv", "cvl", "ml", "tradeoff"))
  if (!is.null(history$selection)) {
    stop("ERROR! fcst.ExpS: nested predictor selection is not supported by the source procedure")
  }
  .fcst.require.regular(data$map$target, "fcst.ExpS")
  if (method %in% c("cv", "cvl") && type == "iterated" && h > 1L) {
    stop("ERROR! fcst.ExpS: iterated CV/CVL for h > 1 is not supported; use direct forecasts")
  }
  if (method != "cvl" && !is.null(options$sigma.sq)) {
    stop("ERROR! fcst.ExpS: sigma.sq is only used by CVL")
  }
  if (!is.null(options[["gamma"]])) stop("ERROR! fcst.ExpS: gamma is only used with gamma.method='fixed'")
  prior <- if (is.null(options$prior)) "midpoint" else match.arg(options$prior, c("midpoint", "uniform"))
  cv <- profile <- posterior <- ml <- sigma.sq <- NULL
  if (method %in% c("ml", "tradeoff")) {
    if (history$k != 1L || any(!is.finite(history$z)) ||
        any(history$z[, 1] != history$z[1, 1]) ||
        data$k != 1L || data$x[1, 1] == 0 || row[1, 1] != data$x[1, 1] ||
        any(data$x[, 1] != data$x[1, 1]) || any(diff(data$map$target) != 1L) ||
        data$N != history$N || (type == "direct" && h > 1L)) {
      stop("ERROR! fcst.ExpS: ml/tradeoff require a complete constant-only one-step model (iteration allowed)")
    }
    if (!is.null(options$gamma.grid) || !is.null(options$eta.grid) || !is.null(options$omega) ||
        !is.null(options$sigma.sq)) stop("ERROR! fcst.ExpS: grids/omega/sigma.sq apply to CV/CVL only")
    ml <- .fcst.exps.ml(data$y, method)
    gamma <- ml$gamma
    eta <- ml$eta
  } else {
    profile <- .fcst.exps.grid(history$N, options$gamma.grid, options$eta.grid)
    cv <- .fcst.cv.folds(history, h, type, options$omega, x.future)
    errors <- matrix(NA_real_, length(cv$folds), nrow(profile))
    for (i in seq_along(cv$folds)) {
      fold <- cv$folds[[i]]
      .fcst.require.regular(fold$data$map$target, "fcst.ExpS")
      for (j in seq_len(nrow(profile))) {
        errors[i, j] <- tryCatch({
          w <- weights.exponential(fold$data$N, profile$gamma[j], fold$data$map$target)
          fit <- WLS.reg(fold$data$y, fold$data$x, w)
          fold$actual - .fcst.fold.predict(fit$coefficients, fold, history, h, type)
        }, error = function(e) NA_real_)
      }
    }
    criterion <- .fcst.cv.score(errors)
    profile$criterion <- criterion
    if (method == "cvl") {
      sigma.sq <- options$sigma.sq
      if (is.null(sigma.sq)) {
        sigma.sq <- mean(OLS.reg(data$y, data$x)$residuals^2)
      }
      if (!is.numeric(sigma.sq) || length(sigma.sq) != 1L ||
          !is.finite(sigma.sq) || sigma.sq < 0) {
        stop("ERROR! fcst.ExpS: CVL requires finite nonnegative sigma.sq")
      }
      mass <- .fcst.exps.prior(profile$eta, prior)
      valid <- is.finite(criterion) & mass > 0
      if (!any(valid)) stop("ERROR! fcst.ExpS: no valid candidate with positive prior mass")
      delta <- criterion - min(criterion[valid])
      if (sigma.sq == 0) {
        posterior <- ifelse(valid & delta == 0, mass, 0)
      } else {
        log.mass <- log(mass) - 0.5 * (delta / sigma.sq)
        log.mass[!valid] <- -Inf
        posterior <- exp(log.mass - max(log.mass))
      }
      posterior <- posterior / sum(posterior)
      eta <- sum(profile$eta * posterior)
      profile$prior <- mass
      profile$posterior <- posterior
    } else {
      eta <- profile$eta[which.min(criterion)]
    }
    gamma <- 1 - eta / history$N
    cv$errors <- errors
  }
  weights <- weights.exponential(data$N, gamma, data$map$target)
  fit <- WLS.reg(data$y, data$x, weights)
  details <- list(forecast = .fcst.predict(fit, row), coefficients = fit$coefficients,
                  weights = weights, gamma = gamma, eta = eta, method = method,
                  N = data$N, eta.scale = history$N, fit = fit, cv = profile,
                  prior = if (method == "cvl") prior else NULL, sigma.sq = sigma.sq,
                  posterior = posterior, ml = ml)
  if (!is.null(cv)) {
    details$cv.origins <- cv$origins
    details$cv.targets <- cv$targets
    details$cv.errors <- cv$errors
    details$cv.excluded <- cv$excluded
    details$cv.selection <- lapply(cv$folds, function(f) f$data$selection)
  }
  class(details) <- "bt_fcstExpS"
  .fcst.suppress.inference(structure(
    list(method = "exps", coefficients = fit$coefficients,
         weights = weights, N = data$N, bp = NULL, details = details), class = "bt_fcst"))
}

#' @title
#' A set of residual based tests for cointegration
#'
#' @details
#' The code provided is the original GAUSS code by Perron and Rodríguez
#' ported to R.
#'
#' @param y,x Variables of interest. `x` can be a matrix of several variables.
#' @param deter A value equal to
#' * 1: quasi-demeaned y and x,
#' * 2: quasi-detrended y and x,
#' * 3: quasi-demeaned y and quasi-detrended x.
#' @param min.lag A minimum number of lags to be used in the tests.
#'
#' @return A list of:
#' * 7x1-matrix of test statistics values,
#' * estimated number of lags.
#'
#' @references
#' Perron, Pierre, and Gabriel Rodríguez.
#' “Residuals‐based Tests for Cointegration with Generalized Least‐squares
#' Detrended Data.”
#' The Econometrics Journal 19, no. 1 (February 1, 2016): 84–111.
#' https://doi.org/10.1111/ectj.12056.
#'
#' @export
coint.PR <- function(y,
                     x,
                     deter,
                     min.lag = 0) {
  if (!is.matrix(y)) y <- as.matrix(y)
  if (!is.matrix(x)) x <- as.matrix(x)

  N <- nrow(y)
  Nc <- ncol(x)

  opt.cbar <- .cbar_PR[[min(Nc, 5)]][deter]

  max.lag <- round(4 * (N / 100)^(1 / 4))

  zy <- if (deter == 1 || deter == 3) {
    .const(N)
  } else if (deter == 2) {
    cbind(.const(N), .trend(N))
  }

  zx <- if (deter == 1) {
    .const(N)
  } else if (deter == 2 || deter == 3) {
    cbind(.const(N), .trend(N))
  }

  y.d <- GLS.reg(y, zy, opt.cbar)$residuals
  x.d <- GLS.reg(x, zx, opt.cbar)$residuals

  model <- OLS.reg(y.d, x.d)
  ud.hat <- cbind(model$residuals)

  statistics.PR(ud.hat, min.lag, max.lag, opt.cbar, deter)
}


#' @title
#' Internal procedure for calculating test statistics from
#' Perron-Rodriguez (2016)
#'
#' @param ud A vector of residuals for testing.
#' @param min.lag,max.lag Minimum and maximum lag number.
#' @param c.bar A `c` parameter used for GLS detrending purposes.
#' @param deter A value equal to
#' * 1: quasi-demeaned y and x,
#' * 2: quasi-detrended y and x,
#' * 3: quasi-demeaned y and quasi-detrended x.
#'
#' @return A list of:
#' * 7x1-matrix of test statistics values,
#' * estimated number of lags.
#'
#' @keywords internal
statistics.PR <- function(ud,
                          min.lag,
                          max.lag,
                          c.bar,
                          deter) {
  N <- nrow(ud)

  if (ncol(ud) > 1) {
    stop("ERROR! resid.tests.PR: too many columns")
  }

  gls.tests <- matrix(0, 7, 1)

  d.ud <- .diffn(ud)
  sum.ud.sq <- sum(.lagn(ud, 1)^2, na.rm = TRUE)

  model.1 <- OLS.reg(ud, .lagn(ud, 1))

  rho.hat <- as.matrix(model.1$coefficients)
  omega <- na.omit(model.1$residuals)
  s2.ud <- sum(omega^2) / (nrow(omega) - 1)
  t.rho <- (rho.hat - 1) / sqrt(s2.ud / sum.ud.sq)

  fin.bic <- Inf
  fin.lag <- min.lag
  lag.bic <- min.lag
  for (lag.bic in min.lag:max.lag) {
    tmp.reg <- .lagn(ud, 1, na = 0)

    if (lag.bic > 0) {
      tmp.reg <- cbind(
        tmp.reg,
        apply(as.array(1:lag.bic), 1, function(l) .lagn(d.ud, l, na = 0))
      )
    }

    # tmp.reg <-
    #   tmp.reg[(lag.bic + 2):N, , drop = FALSE]

    model.2 <- OLS.reg(d.ud, tmp.reg)

    eta <- na.omit(model.2$residuals)
    s2.eta <- sum(eta^2) / (nrow(eta) - ncol(tmp.reg))
    xtx.inv <- solve(t(tmp.reg) %*% tmp.reg)

    sumb <- if (lag.bic == 0) {
      0
    } else {
      sum(model.2$coefficients[2:(lag.bic + 1)])
    }

    s2.adj <- s2.eta / ((1 - sumb)^2)

    cur.bic <- log(c(t(eta) %*% eta) / (N - max.lag)) +
      log(N - max.lag) * lag.bic / (N - max.lag)
    if (cur.bic < fin.bic) {
      fin.bic <- cur.bic
      fin.lag <- lag.bic

      gls.tests[1, 1] <- (ud[N, 1]^2 / N - s2.adj) /
        (2 * sum.ud.sq / N^2)
      gls.tests[2, 1] <- sqrt(2 * sum.ud.sq / (N^2 * s2.adj))
      gls.tests[3, 1] <- gls.tests[1, 1] * gls.tests[2, 1]
      gls.tests[4, 1] <- model.2$coefficients[1] / sqrt(s2.eta * xtx.inv[1, 1])
      gls.tests[5, 1] <- (N - 1) * (rho.hat - 1) -
        (s2.adj - s2.ud) / (2 * sum.ud.sq / N^2)
      gls.tests[6, 1] <- sqrt(s2.ud / s2.adj) * t.rho -
        (s2.adj - s2.ud) / sqrt(4 * s2.adj * sum.ud.sq / N^2)
    }

    lag.bic <- lag.bic + 1
  }

  if (deter == 1 || deter == 3) {
    gls.tests[7, 1] <- (c.bar^2 * sum.ud.sq / N^2 - c.bar * ud[N, 1]^2 / N) / s2.adj
  } else if (deter == 2) {
    gls.tests[7, 1] <- (c.bar^2 * sum.ud.sq / N^2 + (1 - c.bar) * ud[N, 1]^2 / N) / s2.adj
  } else {
    stop("ERROR! Unknown `det.comp` value")
  }

  rownames(gls.tests) <- c(
    "MZ(rho)",
    "MSB",
    "MZ(t.rho)",
    "ADF",
    "Z(rho)",
    "Z(t.rho)",
    if (deter == 1 || deter == 3) {
      "MP(T, demeaned)"
    } else if (deter == 2) {
      "MP(T, detrended)"
    }
  )

  result <- list(
    gls.tests = gls.tests,
    lag = fin.lag
  )
  class(result) <- "bt_cointPR"

  result
}

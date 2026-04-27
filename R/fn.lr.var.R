LR.variance.single <- function(
  u,
  kernel = "bartlett",
  limit.selector = "kpss-q",
  bw.limit = NULL,
  recolor = FALSE,
  max.lag = 0,
  criterion = "bic"
) {
  N <- length(u)

  if (recolor) {
    minIC <- log(sum(u^2) / (N - max.lag))

    arModel <- AR.reg(u, NULL, max.lag, criterion)
    arIC <- info.criterions(arModel$residuals, arModel$lag)[[
      criterion
    ]]

    if (minIC < arIC) {
      # Sul, Phillips and Choi (2003)
      return((sum(u^2) / N) * min(1, N * 0.15))
    } else {
      arCoefs <- arModel$coefficients
      u <- arModel$residuals
      N <- length(u)
    }
  }

  u <- u - mean(u)

  rho <- sum(u[1:(N - 1)] * u[2:N]) / sum(u[2:N]^2)

  lmtF <- .lr.lag.limit(limit.selector, .alpha.single, N)
  wgtF <- .lr.weight(kernel)

  bw <- lmtF(rho, N)
  if (!is.null(bw.limit)) {
    bw <- min(bw, lmtF(bw.limit, N))
  }
  bw <- trunc(bw)

  lrv <- sum(u^2) / N
  for (i in 1:bw) {
    lrv <- lrv + 2 * sum(u[1:(N - i)] * u[(1 + i):N]) * wgtF(i, bw) / N
  }

  if (recolor) {
    lrvR <- lrv / (1 - sum(arCoefs))^2
    return(min(lrvR, N * 0.15 * lrv))
  }

  lrv
}


#' @rdname LR.variance.single
#' @order 2
.lr.var.bartlett <- function(y) {
  LR.variance.single(
    y,
    kernel = "Bartlett",
    limit.selector = "kpss-q"
  )
}

#' @rdname LR.variance.single
#' @order 3
.lr.var.quad <- function(y) {
  LR.variance.single(
    y,
    kernel = "Quadratic",
    limit.selector = "Quadratic"
  )
}

#' @rdname LR.variance.single
#' @order 4
.lr.var.kurozumi <- function(y) {
  LR.variance.single(
    y,
    kernel = "Bartlett",
    limit.selector = "Bartlett",
    bw.limit = 0.8
  )
}

#' @rdname LR.variance.single
#' @order 5
.lr.var.spc <- function(
  y,
  max.lag = 0,
  kernel = "Bartlett",
  criterion = "bic"
) {
  LR.variance.single(
    y,
    max.lag = max.lag,
    kernel = kernel,
    criterion = criterion,
    recolor = TRUE
  )
}


.lr.lag.limit <- function(selector, alpha, N) {
  switch(
    selector,
    "Bartlett" = function(r, N) {
      1.1447 * (N * alpha(r)$q1)^(1 / 3)
    },
    "Parzen" = function(r, N) {
      2.6614 * (N * alpha(r)$q2)^(1 / 5)
    },
    "Tuckey-Hanning" = function(r, N) {
      1.7462 * (N * alpha(r)$q2)^(1 / 5)
    },
    "Quadratic" = function(r, N) {
      1.3221 * (N * alpha(r)$q2)^(1 / 5)
    },
    "kpss-q" = function(r, N) {
      4 * (N / 100)^(1 / 4)
    },
    "kpss-m" = function(r, N) {
      12 * (N / 100)^(1 / 4)
    },
    "Full" = function(r, N) {
      N - 1
    },
    stop("Unknown selector!")
  )
}

.lr.weight <- function(kernel) {
  switch(
    kernel,
    "truncated" = function(i, l) {
      if (abs(i / (l + 1)) <= 1) {
        return(1)
      }
      0
    },
    "Bartlett" = function(i, l) {
      x <- i / (l + 1)
      if (abs(x) < 1) {
        return(1 - abs(x))
      }
      0
    },
    "Parzen" = function(i, l) {
      x <- i / (l + 1)
      if (abs(x) <= 0.5) {
        return(1 - 6 * x^2 + 6 * abs(x)^3)
      } else if (abs(x) <= 1) {
        return(2 * (1 - abs(x))^3)
      }
      0
    },
    "Tukey-Hanning" = function(i, l) {
      x <- i / (l + 1)
      if (abs(x) <= 1) {
        return((1 + cos(pi * x)) / 2)
      }
      0
    },
    "Quadratic" = function(i, l) {
      delta <- 6 * pi * i / (5 * l)
      (3 / delta^2) * (sin(delta) / delta - cos(delta))
    },
    stop("Unknown kernel!")
  )
}

.alpha.single <- function(r) {
  list(
    q1 = 4 * r^2 / (1 + r)^2 / (1 - r)^2,
    q2 = 4 * r^2 / (1 - r)^4
  )
}

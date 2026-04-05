#' @title
#' Test of the cointegration rank with a possible break in trend
#'
#' @description
#' This procedure is aimed on the problem of testing for the cointegration rank
#' of a vector autoregressive process in the case where a trend break may
#' potentially be present in the data.
#'
#' The test is based on estimating the quasi log likelihood for two situations,
#' with break, and without it. The one with the smallest value is considered to
#' be the result.
#'
#' @details
#' The code provided is the original GAUSS code by Harris et al.
#' ported to R.
#'
#' @param y A matrix of \eqn{n} VAR variables.
#' @param r The cointegration rank tested against the alternative of \eqn{n}.
#' @param max.lag The maximum number of lags.
#' @param trim The trimming parameter to find the lower and upper bounds of
#' possible break dates.
#'
#' @return A list of:
#' * the indicator of the rejection of null.
#' * the estimated break point.
#' * the estimated lag number.
#'
#' @references
#' Harris, David, Stephen J. Leybourne, and A. M. Robert Taylor.
#' “Tests of the Co-Integration Rank in VAR Models in the Presence
#' of a Possible Break in Trend at an Unknown Point.”
#' Journal of Econometrics, Innovations in Multiple Time Series Analysis,
#' 192, no. 2 (June 1, 2016): 451–67.
#' https://doi.org/10.1016/j.jeconom.2016.02.010.
#'
#' @export
VECM.test <- function(y,
                      r,
                      max.lag,
                      trim = 0.15) {
  N <- nrow(y)
  Nc <- ncol(y)

  if (r >= Nc) {
    stop("ERROR: VECM.test: invalid values r or y")
  }

  first.break <- trunc(trim * N)
  last.break <- trunc((1 - trim) * N)
  breaks.list <- as.matrix(first.break:last.break)

  loglp <- NULL
  loglbp <- NULL

  for (i in 1:max.lag) {
    loglp <- cbind(loglp, VECM.logl(y, i))
    loglbp <- cbind(loglbp, VECM.break.logl(y, i, breaks.list))
  }

  bhatrp.index <- t(apply(loglbp, 2, which.min))

  rc <- 0:(max.lag - 1) * (Nc + 1) + (r + 1)
  nc <- 0:(max.lag - 1) * (Nc + 1) + (Nc + 1)

  p <- which.min(
    -2 * t(loglp[1, nc]) + log(N) * (Nc^2 * (1:max.lag))
  )
  pbr <- which.min(
    -2 * diag(loglbp[bhatrp.index[1, rc], nc]) +
      log(N) * (Nc^2 * (1:max.lag))
  )

  brphat.index <- bhatrp.index[rc[pbr]]

  tr0.phat <- ifelse(
    2 * (loglp[nc[p]] - loglp[rc[p]]) > .cval_VECM$trend[Nc - r],
    1, 0
  )
  tr1.VECM.phat <- ifelse(
    2 * (loglbp[brphat.index, nc[pbr]] - loglbp[brphat.index, rc[pbr]]) >
      .cval_VECM$br[Nc - r, round(10 * breaks.list[brphat.index] / N) - 1],
    1, 0
  )

  SC1r <- -2 * loglbp[brphat.index, rc[pbr]] +
    (Nc + r + 2 + (Nc^2) * pbr) * log(N)
  SC0r <- -2 * loglp[rc[p]] + (Nc^2) * p * log(N)

  result <- list()

  if (SC1r < SC0r) {
    result$tr1.SCVECM.phat <- tr1.VECM.phat
    result$b.SCVECM <- breaks.list[brphat.index]
    result$p.SCVECM <- pbr
  } else {
    result$tr1.SCVECM.phat <- tr0.phat
    result$b.SCVECM <- 0
    result$p.SCVECM <- p
  }

  result
}


#' @title
#' Quasi log likelihood for VECM without and with breaks
#' @order 1
#'
#' @param y A matrix of \eqn{n} VAR variables.
#' @param p Number of lags.
#' @param breaks.list Vector of possible break dates.
#'
#' @return The vector or matrix with the values of quasi log likelihood
#' for all possible values of cointegration rank and different break dates,
#' if the function with breaks is called.
#'
#' @keywords internal
VECM.logl <- function(y,
                      p) {
  if (!is.matrix(y)) y <- as.matrix(y)

  N <- nrow(y)
  d.y <- .diffn(y)

  z0 <- d.y[-seq_len(p), , drop = FALSE]
  z1 <- cbind(.lagn(y, 1), .trend(N))[-seq_len(p), , drop = FALSE]

  Xp <- as.matrix(rep(1, N))
  if (p > 0) {
    Xp <- cbind(
      Xp,
      apply(as.array(1:p), 1, function(l) .lagn(d.y, l))
    )
  }
  Xp <- Xp[-seq_len(p), , drop = FALSE]

  r0 <- z0 - Xp %*% solve(t(Xp) %*% Xp) %*% t(Xp) %*% z0
  r1 <- z1 - Xp %*% solve(t(Xp) %*% Xp) %*% t(Xp) %*% z1
  Li <- solve(t(chol(t(r1) %*% r1)))

  lam <- cbind(
    c(
      0,
      rev(sort(eigen(
        Li %*% t(r1) %*% r0 %*%
          solve.qr(t(r0) %*% r0) %*%
          t(r0) %*% r1 %*% t(Li)
      )$values))[2:nrow(Li)]
    )
  )
  logL <- -(N - p) / 2 * log(det((t(r0) %*% r0) / (N - p))) +
    cumsum(log(1 - lam))

  rownames(logL) <- NULL
  colnames(logL) <- NULL
  t(logL)
}


#' @rdname VECM.logl
#' @order 2
#' @keywords internal
VECM.break.logl <- function(y,
                            p,
                            breaks.list) {
  if (!is.matrix(y)) y <- as.matrix(y)
  N <- nrow(y)

  d.y <- .diffn(y)

  z0 <- d.y[-seq_len(p), , drop = FALSE]

  logL <- NULL

  for (bc in seq_len(nrow(breaks.list))) {
    b <- breaks.list[bc, 1]

    D <- NULL
    if (p > 0) {
      for (j in 1:p) D <- cbind(D, ifelse(1:N == (b + j), 1, 0))
    }

    E1 <- ifelse(1:N <= b, 1, 0)
    E2 <- ifelse(1:N > b, 1, 0)
    tE <- apply(cbind(E1, E2), 2, cumsum)

    z1 <- cbind(.lagn(y, 1), tE)[-seq_len(p), , drop = FALSE]

    Xp <- cbind(E1, E2, D)
    if (p > 0) {
      for (j in 1:p) Xp <- cbind(Xp, .lagn(d.y, j))
    }
    Xp <- Xp[-seq_len(p), , drop = FALSE]

    r0 <- z0 - Xp %*% solve(t(Xp) %*% Xp) %*% t(Xp) %*% z0
    r1 <- z1 - Xp %*% solve(t(Xp) %*% Xp) %*% t(Xp) %*% z1
    Li <- solve(t(chol(t(r1) %*% r1)))
    lam <- cbind(
      c(
        0,
        rev(sort(eigen(
          Li %*% t(r1) %*% r0 %*%
            Rfast::spdinv(t(r0) %*% r0) %*%
            t(r0) %*% r1 %*% t(Li)
        )$values))[2:nrow(Li)]
      )
    )
    logLb <- -(N - p) / 2 * log(det((t(r0) %*% r0) / (N - p))) +
      cumsum(log(1 - lam))

    logL <- rbind(logL, logLb)
  }

  rownames(logL) <- NULL
  colnames(logL) <- NULL

  logL
}

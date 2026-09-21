#' @rdname uroot.SADF
#' @order 3
#'
#' @details
#' Tests with time transformation are the modified versions of
#' the ordinary `SADF` and `GSADF` tests using Nadaraya-Watson residuals and
#' reindexing procedure by Cavaliere and Taylor (2008).
#'
#' @param omega.est Whether the variance of Nadaraya-Watson residuals should be
#' used.
#' @param truncated Whether the truncation of Nadaraya-Watson residuals is
#' needed.
#' @param is.reindex Whether the Cavaliere and Taylor (2008) time transformation
#' is needed.
#' @param ksi.input The value of the truncation parameter. Can be either `auto`
#' or the explicit numerical value. In the former case the numeric value is
#' estimated.
#' @param hc The scaling parameter for Nadaraya-Watson bandwidth.
#' @param pc The scaling parameter for the estimated truncation parameter value.
#' @param add.p.value Whether the p-value is to be returned. This argument is
#' needed to suppress the calculation of p-values during the precalculation of
#' tables needed for the p-values estimating.
#'
#' @return [uroot.STADF] returns an object of type `bt_SADF` and subclass `bt_STADF`.
#' It's a list of:
#' \item{y}{}
#' \item{N}{Number of observations,}
#' \item{trim}{}
#' \item{const}{}
#' \item{omega.est}{}
#' \item{truncated}{}
#' \item{is.reindex}{}
#' \item{new.index}{the vector of new indices,}
#' \item{ksi.input}{}
#' \item{hc}{}
#' \item{h.est}{}
#' \item{u.hat}{}
#' \item{pc}{}
#' \item{w.sq}{}
#' \item{t.values}{vector of \eqn{t}-values,}
#' \item{*.value}{the value of the corresponding test statistic,}
#' \item{u.hat.truncated}{truncated residuals if truncation was asked for,}
#' \item{ksi`, `sigma}{estimated values of the truncation parameter and resulting s.e. if `ksi.input` equals `auto`,}
#' \item{eta.hat}{the values of reindexing function if reindexing was asked for,}
#' \item{\eqn{p}-value}{if it was asked for.}
#'
#' @references
#' Cavaliere, Giuseppe, and A. M. Robert Taylor.
#' “Time-Transformed Unit Root Tests for Models with Non-Stationary Volatility.”
#' Journal of Time Series Analysis 29, no. 2 (March 2008): 300–330.
#' https://doi.org/10.1111/j.1467-9892.2007.00557.x.
#'
#' @importFrom stats sd
#'
#' @export
uroot.STADF <- function(
  y,
  trim = 0.01 + 1.8 / sqrt(length(y)),
  const = FALSE,
  omega.est = TRUE,
  truncated = TRUE,
  is.reindex = TRUE,
  ksi.input = "auto",
  hc = 1,
  pc = 1,
  add.p.value = TRUE
) {
  n.obs <- length(y)

  ## Part 4.1. NW estimation.
  ## Estimate kernel regression either on the basis of CV or for a fixed h.
  y.0 <- y - y[1]
  my <- diff(y.0)
  mx <- y.0[1:(n.obs - 1)]
  nw.model.cv <- NW.bandwidth(my, mx)
  nw.model <- NW.reg(my, mx, h = nw.model.cv$h)

  h.est <- hc * nw.model.cv$h
  u.hat <- nw.model$u.hat

  ## Truncating the residuals.
  if (truncated == TRUE) {
    if (ksi.input == "auto") {
      ## Calculate sigma.
      sigma <- 0
      bd <- round(0.1 * (n.obs - 1))
      for (s in bd:(n.obs - 1)) {
        sigma1 <- sd(u.hat[(s - bd + 1):s])
        if (sigma1 > sigma) {
          sigma <- sigma1
        }
      }
      ksi <- pc * sigma * (n.obs - 1)^(1 / 7)
    } else {
      ksi <- ksi.input
    }
    u.hat.truncated <- ifelse(abs(u.hat) < ksi, u.hat, 0)
    u.hat.star <- u.hat.truncated
  } else {
    u.hat.star <- u.hat
  }

  ## w.sq - the average of squares residues.
  if (omega.est == TRUE) {
    w.sq <- mean(u.hat.star^2)
  } else {
    w.sq <- 1
  }

  ## Part 4.2. Reindex.
  if (is.reindex == TRUE) {
    tmp.reindex <- reindex.CT(u.hat.star)
    eta.hat <- tmp.reindex$eta.hat
    new.index <- tmp.reindex$new.index
    rm(tmp.reindex)
  } else {
    new.index <- c(0:(n.obs - 1))
  }
  y.tt <- y[new.index + 1]

  ## Part 4.3. STADF test.
  t.values <- c()
  m <- 1

  for (j in (floor(trim * n.obs)):n.obs) {
    ## If we consider a model with a constant,
    ## we subtract the moving average.
    if (const) {
      y.tt.norm <- y.tt - mean(y.tt[1:j])
    } else {
      y.tt.norm <- y.tt - y.tt[1]
    }

    t.values[m] <- (y.tt.norm[j]^2 - y.tt.norm[1]^2 - w.sq * (j - 1)) /
      (w.sq^0.5 * 2 * sum(y.tt.norm[1:(j - 1)]^2)^0.5)
    m <- m + 1
  }

  ## Take the maximum of the calculated t-statistics.
  STADF.value <- max(t.values)

  if (add.p.value) {
    cr.values <- if (const) {
      .cval_SADF_with_const
    } else {
      .cval_SADF_without_const
    }

    p.value <- get.p.values.SADF(STADF.value, n.obs, cr.values)
  }

  result <- list(
    y = y,
    N = n.obs,
    trim = trim,
    const = const,
    omega.est = omega.est,
    truncated = truncated,
    is.reindex = is.reindex,
    new.index = new.index,
    ksi.input = ksi.input,
    hc = hc,
    h.est = h.est,
    u.hat = u.hat,
    pc = pc,
    w.sq = w.sq,
    t.values = t.values,
    STADF.value = STADF.value
  )

  class(result) <- c("bt_STADF", "bt_SADF")

  if (truncated) {
    result$u.hat.truncated <- u.hat.truncated
  }

  if (ksi.input == "auto") {
    result$ksi <- ksi
    result$sigma <- sigma
  }

  if (is.reindex) {
    result$eta.hat <- eta.hat
  }

  if (add.p.value) {
    result$p.value <- p.value
  }

  result
}


#' @rdname uroot.SADF
#' @order 4
#'
#' @return [uroot.GSTADF] returns an object of type `bt_SADF` and subclass `bt_GSTADF`.
#'
#' @importFrom stats sd
#'
#' @export
uroot.GSTADF <- function(
  y,
  trim = 0.01 + 1.8 / sqrt(length(y)),
  const = FALSE,
  omega.est = TRUE,
  truncated = TRUE,
  is.reindex = TRUE,
  ksi.input = "auto",
  hc = 1,
  pc = 1,
  add.p.value = TRUE
) {
  n.obs <- length(y)

  ## Part 4.1. NW estimation.
  ## Estimate kernel regression either on the basis of CV or for a fixed h.
  y.0 <- y - y[1]
  my <- diff(y.0)
  mx <- y.0[1:(n.obs - 1)]
  nw.model.cv <- NW.bandwidth(my, mx)
  nw.model <- NW.reg(my, mx, h = nw.model.cv$h)

  h.est <- hc * nw.model.cv$h
  u.hat <- nw.model$u.hat

  ## Truncating the residuals.
  if (truncated == TRUE) {
    if (ksi.input == "auto") {
      ## Calculate sigma.
      sigma <- 0
      bd <- round(0.1 * (n.obs - 1))
      for (s in bd:(n.obs - 1)) {
        sigma1 <- sd(u.hat[(s - bd + 1):s])
        if (sigma1 > sigma) {
          sigma <- sigma1
        }
      }
      ksi <- pc * sigma * (n.obs - 1)^(1 / 7)
    } else {
      ksi <- ksi.input
    }

    u.hat.truncated <- ifelse(abs(u.hat) < ksi, u.hat, 0)
    u.hat.star <- u.hat.truncated
  } else {
    u.hat.star <- u.hat
  }

  ## w.sq - the average of squares residues.
  if (omega.est == TRUE) {
    w.sq <- mean(u.hat.star^2)
  } else {
    w.sq <- 1
  }

  ## Part 4.2. Reindex.
  if (is.reindex == TRUE) {
    tmp.reindex <- reindex.CT(u.hat.star)
    eta.hat <- tmp.reindex$eta.hat
    new.index <- tmp.reindex$new.index
    rm(tmp.reindex)
  } else {
    new.index <- c(0:(n.obs - 1))
  }
  y.tt <- y[new.index + 1]

  ## Part 4.3. STADF test.
  t.values <- c()
  m <- 1

  for (i in 1:(n.obs - floor(trim * n.obs) + 1)) {
    for (j in (i + floor(trim * n.obs) - 1):n.obs) {
      ## If we consider a model with a constant,
      ## we subtract the moving average.
      if (const) {
        y.tt.norm <- y.tt - mean(y.tt[i:j])
      } else {
        y.tt.norm <- y.tt - y.tt[1]
      }

      t.values[m] <- (y.tt.norm[j]^2 - y.tt.norm[i]^2 - w.sq * (j - i)) /
        (w.sq^0.5 * 2 * sum(y.tt.norm[i:(j - 1)]^2)^0.5)
      m <- m + 1
    }
  }

  ## Take the maximum of the calculated t-statistics.
  GSTADF.value <- max(t.values)

  if (add.p.value) {
    if (const == TRUE) {
      cr.values <- .cval_GSADF_with_const
    } else {
      cr.values <- .cval_GSADF_without_const
    }

    p.value <- get.p.values.SADF(GSTADF.value, n.obs, cr.values)
  }

  result <- list(
    y = y,
    N = n.obs,
    trim = trim,
    const = const,
    omega.est = omega.est,
    truncated = truncated,
    is.reindex = is.reindex,
    new.index = new.index,
    ksi.input = ksi.input,
    hc = hc,
    h.est = h.est,
    u.hat = u.hat,
    pc = pc,
    w.sq = w.sq,
    t.values = t.values,
    GSTADF.value = GSTADF.value
  )

  class(result) <- c("bt_GSTADF", "bt_SADF")

  if (truncated) {
    result$u.hat.truncated <- u.hat.truncated
  }

  if (ksi.input == "auto") {
    result$ksi <- ksi
    result$sigma <- sigma
  }

  if (is.reindex) {
    result$eta.hat <- eta.hat
  }

  if (add.p.value) {
    result$p.value <- p.value
  }

  result
}


#' @title
#' A function that makes reindexing
#'
#' @description
#' The function is aimed to calculate the sequence of indices providing a new
#' "time transformed" time series as in Cavaliere and Taylor (2008).
#'
#' @param u The residuals series for reindexing.
#'
#' @references
#' Cavaliere, Giuseppe, and A. M. Robert Taylor.
#' “Time-Transformed Unit Root Tests for Models with Non-Stationary Volatility.”
#' Journal of Time Series Analysis 29, no. 2 (March 2008): 300–330.
#' https://doi.org/10.1111/j.1467-9892.2007.00557.x.
#'
#' Kurozumi, Eiji, Anton Skrobotov, and Alexey Tsarev.
#' “Time-Transformed Test for Bubbles under Non-Stationary Volatility.”
#' Journal of Financial Econometrics, April 23, 2022.
#' https://doi.org/10.1093/jjfinec/nbac004.
#'
#' @keywords internal
reindex.CT <- function(u) {
  N <- length(u)
  u.sq <- as.numeric(u^2)

  s <- (0:N) / N

  eta.hat <- rep(0, (N + 1))
  for (i in 2:(N + 1)) {
    sT <- floor(s[i] * N)
    eta.hat[i] <- (sum(u.sq[1:sT]) + (s[i] * N - sT) * u.sq[sT + 1]) /
      sum(u.sq)
  }
  eta.hat[N + 1] <- 1

  eta.hat.inv <- rep(0, (N + 1))
  for (i in 2:(N + 1)) {
    k <- length(eta.hat[eta.hat < s[i]])
    s0 <- (k - 1) / N
    eta.hat.inv[i] <- s0 +
      (s[i] - eta.hat[k]) / (eta.hat[k + 1] - eta.hat[k]) / N
  }
  eta.hat.inv[N + 1] <- 1

  list(
    eta.hat = eta.hat,
    eta.hat.inv = eta.hat.inv,
    new.index = floor(eta.hat.inv * N)
  )
}

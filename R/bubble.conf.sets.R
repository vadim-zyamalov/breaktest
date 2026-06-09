#' Confidence set for the emergence, collapse, and restore date of a bubble
#' @order 1
#'
#' @param y      Numeric vector (levels).
#' @param full_N Full sample size.
#' @param phi_a  Explosive AR coefficient.
#' @param s2     Residual variance.
#' @param trim   Trimming parameter.
#'
#' @details
#' Both the "estimated" and "true" parameter variants are handled via the `phi_a`, `phi_b`, `s2` arguments,
#' which are to be taken from the results of
#' [segments.AR1] or [segments.NBCN] functions.
#'
#' @return Named list of binary vectors.
#'
#' @references
#' Kurozumi, Eiji, Anton Skrobotov, and Alexey Tsarev. 2022.
#' "Time-Transformed Test for Bubbles under Non-Stationary Volatility".
#' Journal of Financial Econometrics, advance online publication, 23.
#' https://doi.org/10.1093/jjfinec/nbac004.
#'
#' Kurozumi, Eiji, and Anton Skrobotov. 2026.
#' "Confidence Sets for the Emergence, Collapse, and Recovery Dates of a Bubble".
#' arXiv:2511.16172.
#' Preprint, arXiv. https://doi.org/10.48550/arXiv.2511.16172.
#'
#' @export
cs.bubble.emerge <- function(
  y,
  full_N,
  phi_a,
  s2,
  trim = 0.1
) {
  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }

  N <- length(y) - 1
  Teps <- max(1L, floor(N * trim))
  rho_a <- phi_a - 1

  cset_EMa12 <- integer(N)
  cset_EMb12 <- integer(N)
  cset_LRa12 <- integer(N)
  cset_LRb12 <- integer(N)
  cset_EMa21 <- integer(N)
  cset_EMb21 <- integer(N)
  cset_LRa21 <- integer(N)

  # --- T1 < T2 ---
  for (T1 in 2:(N + 1 - Teps)) {
    T2_range <- (T1 + Teps):(N + 1)
    tstat <- numeric(length(T2_range))
    statLRa <- numeric(length(T2_range))
    statLRb <- numeric(length(T2_range))

    for (k in seq_along(T2_range)) {
      T2 <- T2_range[k]
      rv <- .nbcn.t.stat(y, T1, T2)

      tstat[k] <- rv$t.stat / sqrt(s2)
      statLRa[k] <- (2 * rv$S_y1dy - rho_a * rv$S_yL2) / s2
      statLRb[k] <- (y[T2]^2 - rho_a * rv$S_yL2) / s2
    }

    ind_minEM <- which.min(tstat)
    T2_EM <- ind_minEM + T1 + Teps - 1
    ind_minLRa <- which.min(statLRa)
    T2_LR <- ind_minLRa + T1 + Teps - 1
    ind_minLRb <- which.min(statLRb)
    T2_LRb <- ind_minLRb + T1 + Teps - 1

    EMb_stat <- min(tstat) /
      sqrt(abs(rho_a) * phi_a^(2 * (T2_EM - T1)) * full_N / 2)
    EMa_stat <- sum(tstat) /
      sqrt(abs(phi_a^(2 * (N - T1 + 1)) * full_N / 2 / rho_a))
    LRa_stat <- min(statLRa) / (full_N * phi_a^(2 * (T2_LR - T1)) / 2)
    LRb_stat <- min(statLRb) / (full_N * phi_a^(2 * (T2_LRb - T1)) / 2)

    cval_EM12 <- get.cv.emergence(T1, T2_EM, full_N, (T1 - 1) / full_N)
    cval_LRa <- get.cv.emergence(T1, T2_LR, full_N)
    cval_LRb <- get.cv.emergence(T1, T2_LRb, full_N)

    cset_EMa12[T1 - 1] <- as.integer(EMa_stat > cval_EM12)
    cset_EMb12[T1 - 1] <- as.integer(EMb_stat > cval_EM12)
    cset_LRa12[T1 - 1] <- as.integer(LRa_stat > cval_LRa)
    cset_LRb12[T1 - 1] <- as.integer(LRb_stat > cval_LRb)
  }

  # --- T2 < T1 ---
  for (T1 in (Teps + 2):(N + 1)) {
    T2_range <- 2:(T1 - Teps)
    tstat <- numeric(length(T2_range))
    statLR <- numeric(length(T2_range))

    for (k in seq_along(T2_range)) {
      T2 <- T2_range[k]
      rv <- .nbcn.t.stat(y, T1, T2)
      tstat[k] <- rv$t.stat / sqrt(s2)
      statLR[k] <- (2 * rv$S_y1dy - rho_a * rv$S_yL2) / s2
    }

    EMa_stat <- sum(tstat) / N
    EMb_stat <- max(tstat)
    LR_stat <- max(statLR) / (N^2 * rho_a)

    cval_EM21 <- get.cv.emergence(T1, T1 - Teps, full_N, (T1 - 1) / N)
    cval_LR21 <- get.cv.emergence(T1, T1 - Teps, N)

    cset_EMa21[T1 - 1] <- as.integer(EMa_stat < cval_EM21["EMa21"])
    cset_EMb21[T1 - 1] <- as.integer(EMb_stat < cval_EM21["EMb21"])
    cset_LRa21[T1 - 1] <- as.integer(LR_stat < cval_LR21)
  }

  list(
    cset_EMa = pmin(cset_EMa12, cset_EMa21),
    cset_EMb = pmin(cset_EMb12, cset_EMb21),
    cset_LRa = pmin(cset_LRa12, cset_LRa21),
    cset_LE = pmin(cset_LRb12, cset_EMa21),
    cset_EMa12 = cset_EMa12,
    cset_EMb12 = cset_EMb12,
    cset_LRa12 = cset_LRa12,
    cset_LRb12 = cset_LRb12,
    cset_LE12 = cset_LRb12,
    cset_EMa21 = cset_EMa21,
    cset_EMb21 = cset_EMb21,
    cset_LRa21 = cset_LRa21,
    cset_LE21 = cset_EMa21
  )
}


#' @rdname cs.bubble.emerge
#' @order 2
#' @param lambda_e Scalar, fraction for critical value scaling.
#' @param phi_a    AR coefficient before collapse (explosive phase).
#' @param phi_b    AR coefficient after collapse (recovery phase).
#' @export
cs.bubble.collapse <- function(
  y,
  lambda_e,
  full_N,
  phi_a,
  phi_b,
  s2,
  trim = 0.1
) {
  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }

  N <- length(y) - 1
  Teps <- max(1L, floor(N * trim))

  rho_a <- phi_a - 1
  rho_b <- 1 - phi_b

  # T1 < T2
  cval_EM_lt <- get.cv.collapse(lambda_e, tail = "right")
  cval_LR_lt <- get.cv.collapse(lambda_e, tail = "right", LR = TRUE)

  # T2 < T1
  cval_EM_gt <- get.cv.collapse(lambda_e, tail = "left")
  cval_LR_gt <- get.cv.collapse(lambda_e, tail = "left", LR = TRUE)

  cset_EMa12 <- integer(N)
  cset_EMb12 <- integer(N)
  cset_LRa12 <- integer(N)

  cset_EMa21 <- integer(N)
  cset_EMb21 <- integer(N)
  cset_LRa21 <- integer(N)

  # --- T1 < T2 ---
  for (T1 in 2:(N + 1 - Teps)) {
    T2_range <- (T1 + Teps):(N + 1)
    tstat <- numeric(length(T2_range))
    statLR <- numeric(length(T2_range))

    for (k in seq_along(T2_range)) {
      rv <- .nbcn.t.stat(y, T1, T2_range[k])
      tstat[k] <- rv$t.stat / sqrt(s2)
      statLR[k] <- (2 * rv$S_y1dy + (2 - phi_a - phi_b) * rv$S_yL2) / s2
    }

    scale_EM <- sqrt(full_N * phi_a^(2 * (T1 - 1)) * abs(rho_b) / 2)
    scale_LR <- full_N * (phi_a - phi_b) * phi_a^(2 * (T1 - 1)) / (2 * rho_b)

    cset_EMa12[T1 - 1] <- as.integer(mean(tstat) / scale_EM < cval_EM_lt)
    cset_EMb12[T1 - 1] <- as.integer(max(tstat) / scale_EM < cval_EM_lt)
    cset_LRa12[T1 - 1] <- as.integer(max(statLR) / scale_LR < cval_LR_lt)
  }

  # --- T2 < T1 ---
  for (T1 in (Teps + 2):(N + 1)) {
    T2_range <- 2:(T1 - Teps)
    tstat <- numeric(length(T2_range))
    statLR <- numeric(length(T2_range))

    for (k in seq_along(T2_range)) {
      rv <- .nbcn.t.stat(y, T1, T2_range[k])
      tstat[k] <- rv$t.stat / sqrt(s2)
      statLR[k] <- (2 * rv$S_y1dy + (2 - phi_a - phi_b) * rv$S_yL2) / s2
    }

    scale_EM <- sqrt(full_N * phi_a^(2 * (T1 - 1)) * abs(rho_a) / 2)
    scale_LR <- full_N * (phi_a - phi_b) * phi_a^(2 * (T1 - 1)) / (2 * rho_a)

    idx <- T1 - 1
    cset_EMa21[idx] <- as.integer(mean(tstat) / scale_EM > cval_EM_gt)
    cset_EMb21[idx] <- as.integer(min(tstat) / scale_EM > cval_EM_gt)
    cset_LRa21[idx] <- as.integer(min(statLR) / scale_LR > cval_LR_gt)
  }

  list(
    cset_EMa = pmin(cset_EMa12, cset_EMa21),
    cset_EMb = pmin(cset_EMb12, cset_EMb21),
    cset_LRa = pmin(cset_LRa12, cset_LRa21),
    cset_EMa12 = cset_EMa12,
    cset_EMb12 = cset_EMb12,
    cset_LRa12 = cset_LRa12,
    cset_EMa21 = cset_EMa21,
    cset_EMb21 = cset_EMb21,
    cset_LRa21 = cset_LRa21
  )
}

#' @rdname cs.bubble.emerge
#' @order 3
#' @param lambda_c  Scalar lambda_c = Tc / T_yall.
#' @export
cs.bubble.recovery <- function(
  y,
  full_N,
  lambda_e,
  lambda_c,
  phi_a,
  phi_b,
  s2,
  trim = 0.1
) {
  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }

  N <- length(y) - 1
  Teps <- max(1L, floor(N * trim))
  rho_b <- 1 - phi_b

  exp_ac <- floor(full_N * (lambda_c - lambda_e))

  cset_LRa12 <- integer(N)
  cset_EMa12s <- integer(N)
  cset_EMb12s <- integer(N)
  cset_LRa21 <- integer(N)
  cset_LRb21 <- integer(N)
  cset_EMa21s <- integer(N)
  cset_EMb21s <- integer(N)

  # --- T1 < T2 ---
  for (T1 in 2:(N + 1 - Teps)) {
    T2_range <- (T1 + Teps):(N + 1)
    tstat <- numeric(length(T2_range))
    statLR <- numeric(length(T2_range))

    for (k in seq_along(T2_range)) {
      rv <- .nbcn.t.stat(y, T1, T2_range[k])
      tstat[k] <- rv$t.stat / sqrt(s2)
      statLR[k] <- (2 * rv$S_y1dy + rho_b * rv$S_yL2) / s2
    }

    ind_minLR <- which.min(statLR)
    T2_LR <- ind_minLR + T1 + Teps - 1
    scale_LRa <- full_N *
      (T2_LR - T1) *
      rho_b *
      phi_a^(2 * exp_ac) *
      phi_b^(2 * (T1 - 1))

    LRa_stat <- min(statLR) / scale_LRa
    EMas_stat <- sum(tstat) / N
    EMbs_stat <- min(tstat)

    cv_EM <- get.cv.recovery(lambda_e, (T1 - 1) / N, tail = "right")
    cv_LR <- get.cv.recovery(lambda_e, (T1 - 1) / N, tail = "right", LR = TRUE)

    cset_LRa12[T1 - 1] <- as.integer(LRa_stat > cv_LR["LRas"])
    cset_EMa12s[T1 - 1] <- as.integer(EMas_stat > cv_EM["EMas"])
    cset_EMb12s[T1 - 1] <- as.integer(EMbs_stat > cv_EM["EMbs"])
  }

  # --- T2 < T1 ---
  for (T1 in (2 * Teps + 1):(N + 1)) {
    T2_range <- (Teps + 1):(T1 - Teps)
    tstat <- numeric(length(T2_range))
    statLR <- numeric(length(T2_range))
    statLRb <- numeric(length(T2_range))

    for (k in seq_along(T2_range)) {
      T2 <- T2_range[k]
      rv <- .nbcn.t.stat(y, T1, T2)
      tstat[k] <- rv$t.stat / sqrt(s2)
      statLR[k] <- (2 * rv$S_y1dy + rho_b * rv$S_yL2) / s2
      statLRb[k] <- (-y[T2 - 1]^2 - rv$S_dy2 + rho_b * rv$S_yL2) / s2
    }

    ind_maxLR <- which.max(statLR)
    T2_LR <- ind_maxLR + Teps
    ind_maxLRb <- which.max(statLRb)
    T2_LRb <- ind_maxLRb + Teps
    ind_maxEM <- which.max(tstat)
    T2_EM <- ind_maxEM + Teps

    scale_LRa <- phi_a^(2 * exp_ac) * phi_b^(2 * (T2_LR - 1)) * full_N / 2
    scale_LRb <- phi_a^(2 * exp_ac) * phi_b^(2 * (T2_LRb - 1)) * full_N / 2
    scale_EMa <- phi_a^(exp_ac) * phi_b^(Teps) * sqrt(full_N / abs(rho_b) / 2)
    scale_EMb <- phi_a^(exp_ac) *
      phi_b^(T2_EM - 1) *
      sqrt(full_N * abs(rho_b) / 2)

    LRa_stat <- max(statLR) / scale_LRa
    LRb_stat <- max(statLRb) / scale_LRb
    EMas_stat <- sum(tstat) / scale_EMa
    EMbs_stat <- max(tstat) / scale_EMb

    cv_EM <- get.cv.recovery(lambda_e, lambda1 = 0, tail = "left")
    cv_LR <- get.cv.recovery(lambda_e, lambda1 = 0, tail = "left", LR = TRUE)

    cset_EMa21s[T1 - 1] <- as.integer(EMas_stat < cv_EM)
    cset_EMb21s[T1 - 1] <- as.integer(EMbs_stat < cv_EM)
    cset_LRa21[T1 - 1] <- as.integer(LRa_stat < cv_LR)
    cset_LRb21[T1 - 1] <- as.integer(LRb_stat < cv_LR)
  }

  list(
    cset_EMas = pmin(cset_EMa12s, cset_EMa21s),
    cset_EMbs = pmin(cset_EMb12s, cset_EMb21s),
    cset_LRa = pmin(cset_LRa12, cset_LRa21),
    cset_LE = pmin(cset_EMb12s, cset_LRb21),
    cset_EMa12s = cset_EMa12s,
    cset_EMb12s = cset_EMb12s,
    cset_LRa12 = cset_LRa12,
    cset_LE12 = cset_EMb12s,
    cset_EMa21s = cset_EMa21s,
    cset_EMb21s = cset_EMb21s,
    cset_LRa21 = cset_LRa21,
    cset_LRb21 = cset_LRb21,
    cset_LE21 = cset_LRb21
  )
}

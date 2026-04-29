#' @title
#' A wrapping function around [MDF.1br].
#'
#' @details
#' The code provided is the original Ox code by Skrobotov (2018)
#' ported to R.
#'
#' @param y A time series of interest.
#' @param const,trend Whether the constant term and trend should be included.
#' @param season Whether the seasonal adjustment is needed.
#' @param trim Trimming value for a possible break date bounds.
#'
#' @export
robust.tests.single <- function(
  y,
  const = FALSE,
  trend = FALSE,
  season = FALSE,
  trim = 0.15,
  criterion = "bic"
) {
  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }

  N <- nrow(y)

  if (season) {
    SEAS <- cbind(
      .const(N),
      seasonal.dummies(N)
    )
    y <- OLS.reg(y, SEAS)$residuals
  }

  result <- MDF.1br(
    y = y,
    const = const,
    trend = trend,
    trim = trim,
    criterion = criterion
  )
  result$season <- season
  result
}


#' @title
#' MDF procedure for a single unknown break.
#'
#' @details
#' The code provided is the original Ox code by Skrobotov (2018)
#' ported to R.
#'
#' @param y A time series of interest.
#' @param const Whether the constant term should be included.
#' @param trend Whether the trend term should be included.
#' @param trim Trimming value for a possible break date bounds.
#'
#' @return A list of sublists each containing
#' * The value of statistic: \eqn{MDF-GLS}, \eqn{MDF-OLS},
#' * The asymptotic critical values.
#' \eqn{UR} values are included as well.
MDF.1br <- function(
  y,
  const = FALSE,
  trend = FALSE,
  trim = 0.15,
  criterion = "bic"
) {
  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }
  model <- ifelse(const && trend, "ct", "t")

  result <- list(
    const = const,
    trend = trend
  )

  ## Magic numbers ##
  cv.DF.GLS.m <- .cval_MDF_single[[model]]$cv.DF.GLS.m
  cv.DF.OLS.m <- .cval_MDF_single[[model]]$cv.DF.OLS.m
  cv.DF.GLS.t <- .cval_MDF_single[[model]]$cv.DF.GLS.t
  cv.MDF.GLS <- .cval_MDF_single[[model]]$cv.MDF.GLS
  cv.DF.OLS.t <- .cval_MDF_single[[model]]$cv.DF.OLS.t
  cv.MDF.GLS.lib <- .cval_MDF_single[[model]]$cv.MDF.GLS.lib
  cv.MDF.OLS <- .cval_MDF_single[[model]]$cv.MDF.OLS
  cv.MDF.OLS.lib <- .cval_MDF_single[[model]]$cv.MDF.OLS.lib
  cv.HLT <- .cval_MDF_single[[model]]$cv.HLT
  sap.ur <- .cval_MDF_single[[model]]$sap.ur
  sap.ur1 <- .cval_MDF_single[[model]]$sap.ur1
  sap.ur2 <- .cval_MDF_single[[model]]$sap.ur2
  sap.ur3 <- .cval_MDF_single[[model]]$sap.ur3
  sap.ur4 <- .cval_MDF_single[[model]]$sap.ur4
  sap.ur5 <- .cval_MDF_single[[model]]$sap.ur5
  sap.cv.ur.k0 <- .cval_MDF_single[[model]]$sap.cv.ur.k0
  sap.cv.A.k0 <- .cval_MDF_single[[model]]$sap.cv.A.k0
  sap.cv.ur.k30 <- .cval_MDF_single[[model]]$sap.cv.ur.k30

  ## Start ##
  N <- nrow(y)

  max.lag <- trunc(12 * (N / 100)^(1 / 4))

  first.break <- trunc(trim * N)
  last.break <- trunc((1 - trim) * N)

  tb <- segments.GLS(
    y,
    const,
    trend,
    1,
    first.break,
    last.break,
    trim
  )
  tb <- drop(tb)
  result$break.time <- tb

  x <- cbind(
    .const(N),
    .trend(N),
    if (const) .du(tb, N) else NULL,
    if (trend) .dt(tb, N) else NULL
  )

  ## OLS/GLS Part ##
  ## Mean case
  r_OLS_m <- OLS.reg(y, x[, 1, drop = FALSE])$residuals
  k_m <- max(
    1,
    ADF.test(
      r_OLS_m,
      const = FALSE,
      trend = FALSE,
      max.lag = max.lag,
      criterion = "aic",
      modified.criterion = TRUE
    )$lag
  )

  DF_OLS_m <- ADF.test(
    r_OLS_m,
    const = FALSE,
    trend = FALSE,
    max.lag = k_m,
    criterion = NULL
  )$t.alpha

  r_GLS_m <- GLS.reg(y, x[, 1, drop = FALSE], -7)$residuals
  DF_GLS_m <- ADF.test(
    r_GLS_m,
    const = FALSE,
    trend = FALSE,
    max.lag = k_m,
    criterion = NULL
  )$t.alpha

  ## Trend case
  r_OLS_t <- OLS.reg(y, x[, 1:2])$residuals
  k_t <- max(
    1,
    ADF.test(
      r_OLS_t,
      const = FALSE,
      trend = FALSE,
      max.lag = max.lag,
      criterion = "aic",
      modified.criterion = TRUE
    )$lag
  )

  DF_OLS_t <- ADF.test(
    r_OLS_t,
    const = FALSE,
    trend = FALSE,
    max.lag = k_t,
    criterion = NULL
  )$t.alpha

  r_GLS_t <- GLS.reg(y, x[, 1:2], -13.5)$residuals
  DF_GLS_t <- ADF.test(
    r_GLS_t,
    const = FALSE,
    trend = FALSE,
    max.lag = k_t,
    criterion = NULL
  )$t.alpha

  ## MDF ##
  ## One break ##
  MDF_OLS <- Inf
  MDF_GLS <- Inf
  for (tb1 in first.break:last.break) {
    z <- cbind(
      .const(N),
      .trend(N),
      if (const) .du(tb1, N) else NULL,
      if (trend) .dt(tb1, N) else NULL
    )

    r_OLS <- OLS.reg(y, z)$residuals
    k_tb <- max(
      1,
      ADF.test(
        r_OLS,
        const = FALSE,
        trend = FALSE,
        max.lag = max.lag,
        criterion = "aic",
        modified.criterion = TRUE
      )$lag
    )

    DF1 <- ADF.test(
      r_OLS,
      const = FALSE,
      trend = FALSE,
      max.lag = k_tb,
      criterion = NULL
    )
    DF1_tb <- N * DF1$alpha / (1 - sum(DF1$coefficients) + DF1$alpha)

    r_GLS <- GLS.reg(y, z, -17.6)$residuals
    DF2_tb <- ADF.test(
      r_GLS,
      const = FALSE,
      trend = FALSE,
      max.lag = k_tb,
      criterion = NULL
    )$t.alpha

    MDF_OLS <- min(MDF_OLS, DF1_tb)
    MDF_GLS <- min(MDF_GLS, DF2_tb)
  }

  t_HLT <- ur.KPSS.HLT(y, const, trim)

  tmp.PY <- PY.statistic(y, const, trend, criterion, trim, max.lag)
  t_PY <- tmp.PY$statistic
  cv_PY <- tmp.PY$critical.value
  rm(tmp.PY)

  ## init value test ##
  t_alpha <- as.numeric(OLS.reg(y, x)$t.stats[1])
  t_alpha_id <- as.numeric(abs(t_alpha) > 1.96)

  ## trend test ##
  x <- cbind(.const(N), .trend(N))
  LRV_y <- sum(r_OLS_t^2) / N
  t_y <- as.numeric(OLS.reg(y, x)$t.stats[2])

  dy <- y[2:N] - (1 - 30 / N) * y[1:(N - 1)]
  dx <- x[2:N, ] - (1 - 30 / N) * x[1:(N - 1), ]
  t_dy <- as.numeric(OLS.reg(dy, dx)$t.stats[2])

  KPSS_t <- .kpss.statistic(r_OLS_t, LRV_y)
  lam <- exp(-0.00025 * ((DF_GLS_t / KPSS_t)^2))
  t_lam <- (1 - lam) * t_y + lam * t_dy

  ## UR-HLT
  t_lambda_id <- as.numeric(t_HLT > cv.HLT)

  ur_id_sa <- as.numeric(
    (DF_GLS_t < (sap.ur * sap.cv.ur.k0 * cv.DF.GLS.t)) ||
      (MDF_GLS < (sap.ur * sap.cv.ur.k0 * cv.MDF.GLS)) ||
      (DF_OLS_t < (sap.ur * sap.cv.ur.k0 * cv.DF.OLS.t)) ||
      (MDF_OLS < (sap.ur * sap.cv.ur.k0 * cv.MDF.OLS))
  )
  ur1_id_sa <- as.numeric(
    (DF_OLS_t < (sap.ur1 * sap.cv.ur.k0 * cv.DF.OLS.t)) ||
      (MDF_OLS < (sap.ur1 * sap.cv.ur.k0 * cv.MDF.OLS))
  )
  ur2_id_sa <- as.numeric(
    (MDF_GLS < (sap.ur2 * sap.cv.ur.k30 * cv.MDF.GLS.lib)) ||
      (MDF_OLS < (sap.ur2 * sap.cv.ur.k30 * cv.MDF.OLS.lib))
  )
  ur3.id.sa <- as.numeric(MDF_OLS < sap.cv.ur.k0 * cv.MDF.OLS.lib)
  UR_HLT <-
    (1 - t_lambda_id) *
    (1 - t_alpha_id) *
    ur_id_sa +
    (1 - t_lambda_id) * t_alpha_id * ur1_id_sa +
    t_lambda_id * (1 - t_alpha_id) * ur2_id_sa +
    t_lambda_id * t_alpha_id * ur3.id.sa

  ## UR-PY
  t_lambda_id <- as.numeric(t_PY > cv_PY[2])
  ur_id_sa <- as.numeric(
    (DF_GLS_t < (sap.ur * sap.cv.ur.k0 * cv.DF.GLS.t)) ||
      (MDF_GLS < (sap.ur * sap.cv.ur.k0 * cv.MDF.GLS)) ||
      (DF_OLS_t < (sap.ur * sap.cv.ur.k0 * cv.DF.OLS.t)) ||
      (MDF_OLS < (sap.ur * sap.cv.ur.k0 * cv.MDF.OLS))
  )
  ur1_id_sa <- as.numeric(
    (DF_OLS_t < (sap.ur1 * sap.cv.ur.k0 * cv.DF.OLS.t)) ||
      (MDF_OLS < (sap.ur1 * sap.cv.ur.k0 * cv.MDF.OLS))
  )
  ur2_id_sa <- as.numeric(
    (MDF_GLS < (sap.ur2 * sap.cv.ur.k30 * cv.MDF.GLS.lib)) ||
      (MDF_OLS < (sap.ur2 * sap.cv.ur.k30 * cv.MDF.OLS.lib))
  )
  ur3_id_sa <- as.numeric(MDF_OLS < sap.cv.ur.k30 * cv.MDF.OLS.lib)
  UR_PY <-
    (1 - t_lambda_id) *
    (1 - t_alpha_id) *
    ur_id_sa +
    (1 - t_lambda_id) * t_alpha_id * ur1_id_sa +
    t_lambda_id * (1 - t_alpha_id) * ur2_id_sa +
    t_lambda_id * t_alpha_id * ur3_id_sa

  ## A-HLT
  t_lambda_id <- as.numeric(t_HLT > cv.HLT)
  ur2_id_sa <- as.numeric(
    (MDF_GLS < (sap.ur2 * sap.cv.ur.k30 * cv.MDF.GLS.lib)) ||
      (MDF_OLS < (sap.ur2 * sap.cv.ur.k30 * cv.MDF.OLS.lib))
  )
  ur2a_id_sa <- as.numeric(MDF_OLS < sap.cv.ur.k30 * cv.MDF.OLS.lib)
  ur3_id_sa <- as.numeric(
    (MDF_GLS < (sap.ur3 * sap.cv.A.k0 * cv.MDF.GLS)) ||
      (MDF_OLS < (sap.ur3 * sap.cv.A.k0 * cv.MDF.OLS))
  )
  ur3a_id_sa <- as.numeric(MDF_OLS < sap.cv.A.k0 * cv.MDF.OLS)

  A_HLT <-
    (1 - t_lambda_id) *
    (1 - t_alpha_id) *
    ur3_id_sa +
    (1 - t_lambda_id) * t_alpha_id * ur3a_id_sa +
    t_lambda_id * (1 - t_alpha_id) * ur2_id_sa +
    t_lambda_id * t_alpha_id * ur2a_id_sa

  ## A-PY
  t_lambda_id <- as.numeric(t_PY > cv_PY[2])
  ur2_id_sa <- as.numeric(
    (MDF_GLS < (sap.ur2 * sap.cv.ur.k30 * cv.MDF.GLS.lib)) ||
      (MDF_OLS < (sap.ur2 * sap.cv.ur.k30 * cv.MDF.OLS.lib))
  )
  ur2a_id_sa <- as.numeric(MDF_OLS < sap.cv.ur.k30 * cv.MDF.OLS.lib)
  ur3_id_sa <- as.numeric(
    (MDF_GLS < (sap.ur3 * sap.cv.A.k0 * cv.MDF.GLS)) ||
      (MDF_OLS < (sap.ur3 * sap.cv.A.k0 * cv.MDF.OLS))
  )
  ur3_id_sa <- as.numeric(MDF_OLS < sap.cv.A.k0 * cv.MDF.OLS)

  A_PY <- (1 - t_lambda_id) *
    (1 - t_alpha_id) *
    ur3_id_sa +
    (1 - t_lambda_id) * t_alpha_id * ur3_id_sa +
    t_lambda_id * (1 - t_alpha_id) * ur2_id_sa +
    t_lambda_id * t_alpha_id * ur2a_id_sa

  ## lambda-HLT
  t_lam_t_id <- abs(t_lam) > 1.96
  t_lambda_id <- as.numeric(t_HLT > cv.HLT)

  ur_id_sa <- as.numeric(
    (DF_GLS_t < (sap.ur * sap.cv.ur.k0 * cv.DF.GLS.t)) ||
      (MDF_GLS < (sap.ur * sap.cv.ur.k0 * cv.MDF.GLS)) ||
      (DF_OLS_t < (sap.ur * sap.cv.ur.k0 * cv.DF.OLS.t)) ||
      (MDF_OLS < (sap.ur * sap.cv.ur.k0 * cv.MDF.OLS))
  )
  ur1_id_sa <- as.numeric(
    (DF_OLS_t < (sap.ur1 * sap.cv.ur.k0 * cv.DF.OLS.t)) ||
      (MDF_OLS < (sap.ur1 * sap.cv.ur.k0 * cv.MDF.OLS))
  )
  ur2_id_sa <- as.numeric(
    (MDF_GLS < (sap.ur2 * sap.cv.ur.k30 * cv.MDF.GLS.lib)) ||
      (MDF_OLS < (sap.ur2 * sap.cv.ur.k30 * cv.MDF.OLS.lib))
  )
  ur4_id_sa <- as.numeric(
    (DF_OLS_m < (sap.ur4 * sap.cv.ur.k0 * cv.DF.OLS.m)) ||
      (DF_OLS_t < (sap.ur4 * sap.cv.ur.k0 * cv.DF.OLS.t)) ||
      (MDF_OLS < (sap.ur4 * sap.cv.ur.k0 * cv.MDF.OLS))
  )
  ur5_id_sa <- as.numeric(
    (DF_GLS_m < (sap.ur5 * sap.cv.ur.k0 * cv.DF.GLS.m)) ||
      (DF_OLS_m < (sap.ur5 * sap.cv.ur.k0 * cv.DF.OLS.m)) ||
      (DF_GLS_t < (sap.ur5 * sap.cv.ur.k0 * cv.DF.GLS.t)) ||
      (MDF_GLS < (sap.ur5 * sap.cv.ur.k0 * cv.MDF.GLS)) ||
      (DF_OLS_t < (sap.ur5 * sap.cv.ur.k0 * cv.DF.OLS.t)) ||
      (MDF_OLS < (sap.ur5 * sap.cv.ur.k0 * cv.MDF.OLS))
  )

  t_lam_UR_HLT <- (1 - t_lambda_id) *
    ((1 - t_alpha_id) * ur_id_sa + t_alpha_id * ur1_id_sa) +
    t_lambda_id *
      (t_lam_t_id *
        (1 - t_alpha_id) *
        ur_id_sa +
        t_lam_t_id * t_alpha_id * ur1_id_sa +
        (1 - t_lam_t_id) * (1 - t_alpha_id) * ur5_id_sa +
        (1 - t_lam_t_id) * t_alpha_id * ur4_id_sa)

  ## lambda-PY
  t_lambda_id <- as.numeric(t_PY > cv_PY[2])
  t_lam_UR_PY <- (1 - t_lambda_id) *
    ((1 - t_alpha_id) * ur_id_sa + t_alpha_id * ur1_id_sa) +
    t_lambda_id *
      (t_lam_t_id *
        (1 - t_alpha_id) *
        ur_id_sa +
        t_lam_t_id * t_alpha_id * ur1_id_sa +
        (1 - t_lam_t_id) * (1 - t_alpha_id) * ur5_id_sa +
        (1 - t_lam_t_id) * t_alpha_id * ur4_id_sa)

  ## results
  result$HLT <- list(
    statistic = drop(t_HLT),
    cv = cv.HLT
  )
  result$PY <- list(
    statistic = t_PY,
    cv = cv_PY[2]
  )
  result$DF.GLS <- list(
    statistic = DF_GLS_t,
    cv = cv.DF.GLS.t
  )
  result$DF.OLS <- list(
    statistic = DF_OLS_t,
    cv = cv.DF.OLS.t
  )
  result$MDF.GLS <- list(
    statistic = MDF_GLS,
    cv = cv.MDF.GLS
  )
  result$MDF.OLS <- list(
    statistic = MDF_OLS,
    cv = cv.MDF.OLS
  )
  result$A.HLT <- A_HLT
  result$A.PY <- A_PY
  result$UR.HLT <- UR_HLT
  result$UR.PY <- UR_PY
  result$URR.HLT <- t_lam_UR_HLT
  result$URR.PY <- t_lam_UR_PY

  class(result) <- "bt_mdfHLT"

  result
}


#' @title
#' Perron-Yabu (2009) statistic for break at unknown date.
#'
#' @details
#' The code provided is the original Ox code by Skrobotov (2018)
#' ported to R.
#'
#' @param y A time series of interest.
#' @param const,trend Allowing the break in constant or trend.
#' @param criterion Needed information criterion: aic, bic, hq or lwz.
#' @param trim A trimming parameter to determine the lower and upper bounds for
#' a possible break point.
#' @param max.lag The maximum possible lag in the model.
#'
#' @return A list of the estimated Wald statistic as well as its c.v.
#'
#' @references
#' Perron, Pierre, and Tomoyoshi Yabu.
#' “Testing for Shifts in Trend With an Integrated or
#' Stationary Noise Component.”
#' Journal of Business & Economic Statistics 27, no. 3 (July 2009): 369–96.
#' https://doi.org/10.1198/jbes.2009.07268.
#'
#' @export
PY.statistic <- function(
  y,
  const = FALSE,
  trend = FALSE,
  criterion = "aic",
  trim = 0.15,
  max.lag = trunc(12 * (length(y) / 100)^(1 / 4))
) {
  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }
  if (!trim %in% c(0.01, 0.05, 0.10, 0.15, 0.25)) {
    stop("ERROR! PY.single: Illegal trim value")
  }

  model <- paste0(
    ifelse(const, "c", ""),
    ifelse(trend, "t", "")
  )
  if (model == "") {
    stop("ERROR! PY.single: Illegal model")
  }

  VR <- .cval_PY_single[[model]]$VR
  v.t <- .cval_PY_single[[model]]$v.t
  c.v <- .cval_PY_single[[model]]$c.v

  N <- nrow(y)
  first.break <- max(trunc(trim * N), max.lag + 2)
  last.break <- trunc((1 - trim) * N)

  vect1 <- NULL

  for (tb in first.break:last.break) {
    lambda <- (tb - 1) / N

    reg <- cbind(
      .const(N),
      .trend(N),
      if (const) .du(tb, N) else NULL,
      if (trend) .dt(tb, N) else NULL
    )

    # Estimation of alpha
    khat <- max(1, AR.reg(y, reg, max.lag, criterion)$lag)
    u <- OLS.reg(y, reg)$residuals
    du <- .diffn(u, na = 0)
    depu <- u[khat:N]
    regu <- .lagn(u, 1, na = 0)
    for (i in seq_len(khat - 1)) {
      regu <- cbind(regu, .lagn(du, i, na = 0))
    }
    regu <- regu[khat:N, , drop = FALSE]

    tmp.OLS <- OLS.reg(depu, regu)
    b <- tmp.OLS$coefficients
    ehat <- tmp.OLS$residuals

    VCV <- solve(t(regu) %*% regu) * sum(ehat^2) / length(ehat)

    ahat <- b[1]
    vahat <- VCV[1, 1]
    tau1 <- (ahat - 1) / sqrt(vahat)

    # Upper Biased Estimator
    tau05 <- v.t[ceiling(lambda * 10)]

    IP <- trunc((khat + 1) / 2)
    r <- ncol(reg)
    k <- 10

    c1 <- sqrt((1 + r) * N)
    c2 <- ((1 + r) * N - tau05^2 * (IP + N)) / (tau05 * (tau05 + k) * (IP + N))

    if (tau1 > tau05) {
      ctau <- -tau1
    } else if (tau1 <= tau05 && tau1 > -k) {
      ctau <- IP * tau1 / N - (r + 1) / (tau1 + c2 * (tau1 + k))
    } else if (tau1 <= -k && tau1 > -c1) {
      ctau <- IP * tau1 / N - (r + 1) / tau1
    } else if (tau1 <= -c1) {
      ctau <- 0
    }

    amus <- ahat + ctau * sqrt(vahat)
    if (amus >= 1) {
      amus <- 1
    } else if (amus <= -1) {
      amus <- -0.99
    }

    CR <- sqrt(N) * abs(amus - 1)
    if (CR <= 1) {
      amus <- 1
    }

    gdep <- rbind(
      y[1, , drop = FALSE],
      y[2:N, , drop = FALSE] -
        amus * y[1:(N - 1), , drop = FALSE]
    )
    greg <- rbind(
      reg[1, , drop = FALSE],
      reg[2:N, , drop = FALSE] -
        amus * reg[1:(N - 1), , drop = FALSE]
    )

    tmp.OLS <- OLS.reg(gdep, greg)
    b <- tmp.OLS$coefficients
    v <- tmp.OLS$residuals

    if (khat == 1) {
      h0 <- sum(v^2) / length(v)
    } else {
      if (amus == 1) {
        regv <- NULL
        for (ki in seq_len(khat - 1)) {
          regv <- cbind(regv, .lagn(v, ki, na = 0))
        }

        depv <- v[(khat - 1):N]
        regv <- regv[(khat - 1):N, , drop = FALSE]

        tmp.OLS <- OLS.reg(depv, regv)
        beta <- tmp.OLS$coefficients
        e <- tmp.OLS$residuals

        if (const && !trend) {
          vbeta <- matrix(0, nrow = khat - 1, ncol = 3)
          for (ki in seq_len(khat - 1)) {
            DUki <- .du(tb - ki, N)
            regki <- cbind(.const(N), .trend(N), DUki)
            gdepki <- rbind(
              y[1, , drop = FALSE],
              y[2:N, , drop = FALSE] - amus * y[1:(N - 1), , drop = FALSE]
            )
            gregki <- rbind(
              reg[1, ],
              regki[2:N, ] - amus * regki[1:(N - 1), ]
            )
            vbeta[ki, ] <- drop(OLS.reg(gdepki, gregki)$coefficients)
          }
          b[3] <- b[3] - drop(t(vbeta[, 3]) %*% beta)
          h0 <- sum(e^2) / (N - khat)
        }

        if (!const && trend) {
          h0 <- (sum(e^2) / (N - khat)) /
            ((1 - sum(beta))^2)
        }

        if (const && trend) {
          vbeta <- matrix(0, nrow = khat - 1, ncol = 4)
          for (ki in seq_len(khat - 1)) {
            regki <- cbind(
              .const(N),
              .trend(N),
              .du(tb - ki, N),
              .du(tb - ki, N) * (.trend(N) - tb)
            )
            gdepki <- rbind(
              y[1, , drop = FALSE],
              y[2:N, , drop = FALSE] - amus * y[1:(N - 1), , drop = FALSE]
            )
            gregki <- rbind(
              reg[1, ],
              regki[2:N, ] - amus * regki[1:(N - 1), ]
            )
            vbeta[ki, ] <- drop(OLS.reg(gdepki, gregki)$coefficients)
          }

          sige <- sum(e^2) / (N - khat)

          h0 <- sige / ((1 - sum(beta))^2)
          b[3] <- (sqrt(h0) / sqrt(sige)) *
            (b[3] - drop(t(vbeta[, 3]) %*% beta))
        }
      }

      if (abs(amus) < 1) {
        h0 <- .lr.var.quad(v)
      }
    }

    VCV <- h0 * solve(t(greg) %*% greg)
    vect1 <- c(
      vect1,
      t(VR %*% b) %*% solve(VR %*% VCV %*% t(VR)) %*% (VR %*% b)
    )
  }

  wald <- log(sum(exp(vect1 / 2)) / N)

  if (trim == 0.01) {
    cv <- c.v[1, ]
  }
  if (trim == 0.05) {
    cv <- c.v[2, ]
  }
  if (trim == 0.10) {
    cv <- c.v[3, ]
  }
  if (trim == 0.15) {
    cv <- c.v[4, ]
  }
  if (trim == 0.25) {
    cv <- c.v[5, ]
  }

  result <- list(
    statistic = wald,
    critical.value = cv
  )
  class(result) <- "bt_PY"

  result
}

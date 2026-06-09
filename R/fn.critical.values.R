#' @title
#' Critical values for break point confidence intervals
#'
#' @description
#' Auxiliary function returning pre-calculated critical values for
#' [coint.conf.sets]
#'
#' @details
#' The function is not intended to be used directly so it's not exported.
#'
#' @param lambda Relative break point position.
#' @param trend Whether thern is to be included.
#' @param conf.level Confidense level.
#' @param n.zb Number of variables with breaks.
#' @param n.zf Number of variables without breaks.
#'
#' @keywords internal
get.cv.coint.conf.sets <- function(lambda, trend, conf.level, n.zb, n.zf) {
  if (conf.level == 0.9) {
    if (n.zf == 0) {
      values_table <- .cval_break_date_cset[[1]][[1]]
    } else if (n.zb %in% 0:3) {
      values_table <- .cval_break_date_cset[[1]][[2 + n.zb]]
    } else {
      stop("ERROR: Invalid value of variable p_zb or p_zf")
    }
  } else if (conf.level == 0.9) {
    if (n.zf == 0) {
      values_table <- .cval_break_date_cset[[2]][[1]]
    } else if (n.zb %in% 0:3) {
      values_table <- .cval_break_date_cset[[2]][[2 + n.zb]]
    } else {
      stop("ERROR: Invalid value of variable p_zb or p_zf")
    }
  } else {
    stop("ERROR: Invalid value of variable conf.level")
  }

  lambda.d <- abs(lambda - 0.5)

  if (n.zf == 0) {
    if (!trend) {
      coef_sup <- values_table$sup_all[, n.zb]
      coef_avg <- values_table$avg_all[, n.zb]
      coef_exp <- values_table$exp_all[, n.zb]
    } else if (trend) {
      coef_sup <- values_table$sup_all[, (4 + n.zb)]
      coef_avg <- values_table$avg_all[, (4 + n.zb)]
      coef_exp <- values_table$exp_all[, (4 + n.zb)]
    } else {
      stop("ERROR: Invalid value of variable trend")
    }
  } else if (n.zb %in% 0:3) {
    if (!trend) {
      coef_sup <- values_table$sup_all[, n.zf]
      coef_avg <- values_table$avg_all[, n.zf]
      coef_exp <- values_table$exp_all[, n.zf]
    } else if (trend) {
      coef_sup <- values_table$sup_all[, (4 - n.zb + n.zf)]
      coef_avg <- values_table$avg_all[, (4 - n.zb + n.zf)]
      coef_exp <- values_table$exp_all[, (4 - n.zb + n.zf)]
    }
  } else {
    stop("ERROR: Invalid value of variable p_zb or p_zf")
  }

  cval_sup <- coef_sup[1] +
    coef_sup[2] / (lambda.d + 1) +
    coef_sup[3] * lambda.d +
    coef_sup[4] * lambda.d^2 +
    coef_sup[5] * lambda.d^3
  cval_avg <- coef_avg[1] +
    coef_avg[2] / (lambda.d + 1) +
    coef_avg[3] * lambda.d +
    coef_avg[4] * lambda.d^2 +
    coef_avg[5] * lambda.d^3
  cval_exp <- coef_exp[1] +
    coef_exp[2] / (lambda.d + 1) +
    coef_exp[3] * lambda.d +
    coef_exp[4] * lambda.d^2 +
    coef_exp[5] * lambda.d^3

  list(
    cval_sup = cval_sup,
    cval_avg = cval_avg,
    cval_exp = cval_exp
  )
}

#' @title
#' Critical values for SADF-type tests
#'
#' @description
#' Interpolating p-value for intermediate observation numbers for SADF-type
#' tests.
#'
#' @details
#' The function is not intended to be used directly so it's not exported.
#'
#' @param statistic The statistic value.
#' @param n.obs The number of observations.
#' @param cr.values The set of precalculated tables.
#'
#' @keywords internal
get.p.values.SADF <- function(statistic, n.obs, cr.values) {
  N.table.obs <- as.numeric(names(cr.values))

  if (n.obs < min(N.table.obs)) {
    warning("Too little number of observations, using data for T = 30")
    i.0 <- min(N.table.obs)
  } else {
    i.0 <- max(N.table.obs[N.table.obs <= n.obs])
  }
  p.0 <- sum(cr.values[[as.character(i.0)]] > statistic) /
    length(cr.values[[as.character(i.0)]])

  if (n.obs > max(N.table.obs)) {
    i.1 <- max(N.table.obs)
  } else {
    i.1 <- min(N.table.obs[N.table.obs >= n.obs])
  }
  p.1 <- sum(cr.values[[as.character(i.1)]] > statistic) /
    length(cr.values[[as.character(i.1)]])

  if (i.0 != i.1) {
    p.value <- p.0 + (p.1 - p.0) * (n.obs - i.0) / (i.1 - i.0)
  } else {
    p.value <- p.0
  }

  p.value
}


get.cv.collapse <- function(lambda_e, tail = c("right", "left"), LR = FALSE) {
  tail <- match.arg(tail)
  pwr <- if (LR) 1 else 0.5
  base <- (lambda_e * qchisq(0.05, df = 1))^pwr
  if (tail == "right") -base else base
}


get.cv.emergence <- function(T1, T2, N, lambda_e = NULL) {
  if (!is.null(lambda_e)) {
    if (T1 < T2) {
      sqrt(lambda_e * qchisq(0.05, df = 1))
    } else {
      c(
        EMa21 = .nbcn.cval(.cval_NBCN$EMa21_e, lambda_e),
        EMb21 = .nbcn.cval(.cval_NBCN$EMb21_e, lambda_e)
      )
    }
  } else {
    if (T1 < T2) {
      (T1 - 1) / N * qchisq(0.05, df = 1)
    } else {
      .nbcn.cval(.cval_NBCN$LRb21_e, (T1 - 1) / N)
    }
  }
}


get.cv.recovery <- function(
  lambda_e,
  lambda1 = NULL,
  tail = c("right", "left"),
  LR = FALSE
) {
  tail <- match.arg(tail)
  pwr <- if (LR) 1 else 0.5
  if (tail == "right") {
    if (!LR) {
      c(
        EMa = .nbcn.cval(.cval_NBCN$EMa12_r, lambda1),
        EMas = .nbcn.cval(.cval_NBCN$EMa12_rs, lambda1),
        EMb = .nbcn.cval(.cval_NBCN$EMb12_r, lambda1),
        EMbs = .nbcn.cval(.cval_NBCN$EMb12_rs, lambda1)
      )
    } else {
      c(
        LRa = lambda_e * qchisq(0.05, df = 1),
        LRas = .nbcn.cval(.cval_NBCN$LRa12_rs, lambda1)
      )
    }
  } else {
    -(lambda_e * qchisq(0.05, df = 1))^pwr
  }
}

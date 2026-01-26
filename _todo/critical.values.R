#' @title
#' Critical values for KPSS test with 1 break
#'
#' @description
#' Auxiliary function returning pre-calculated critical values for
#' [KPSS.1.break]
#'
#' @details
#' The function is not intended to be used directly so it's not exported.
#'
#' @param model A scalar equal to
#' * 1: for model An,
#' * 2: for model A,
#' * 3: for model B,
#' * 4: for model C,
#' * 5: for model D,
#' * 6: for model E.
#' @param break.point Position of the break point.
#' @param n.obs Number of observations.
#' @param k Number of RHS variables.
#'
#' @keywords internal
#' @export
get.cv.KPSS.1p <- function(model,
                           break.point,
                           n.obs,
                           k) {
    lambda <- round(break.point / n.obs, 1)
    i <- trunc(lambda * 10)
    if (i == 0) i <- 1
    if (i == 10) i <- 9

    m_vc <- .cval_kpss_1p[[model]][[k]]

    return(m_vc[, i, drop = FALSE])
}

#' @title
#' Critical values for KPSS test with 2 breaks
#'
#' @description
#' Auxiliary function returning pre-calculated critical values for
#' [KPSS.2.breaks]
#'
#' @details
#' The function is not intended to be used directly so it's not exported.
#'
#' @param model A scalar equal to
#' * 1: for the AA (without trend) model,
#' * 2: for the AA (with trend) model,
#' * 3: for the BB model,
#' * 4: for the CC model,
#' * 5: for the AC-CA model.
#' @param break.point Position of the break point.
#' @param n.obs Number of observations.
#'
#' @keywords internal
#' @export
get.cv.KPSS.2p <- function(model,
                           break.point,
                           n.obs) {
    lambda1 <- round(break.point[1] / n.obs, 1)
    i <- trunc(lambda1 * 10)

    lambda2 <- round(break.point[2] / n.obs, 1)
    j <- trunc(lambda2 * 10)
    if (1 <= model && model <= 4) {
        j <- j - 1
    }

    v_vc <- .cval_kpss_2p[[model]]

    return(
        c(
            v_vc[1][i, j],
            v_vc[2][i, j],
            v_vc[3][i, j],
            v_vc[4][i, j]
        )
    )
}

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
#' @param p_zb Number of variables with breaks.
#' @param p_zf Number of variables without breaks.
#'
#' @keywords internal
#' @export
get.cv.coint.conf.sets <- function(lambda,
                                   trend,
                                   conf.level,
                                   p_zb,
                                   p_zf) {
    if (conf.level == 0.9) {
        if (p_zf == 0) {
            values_table <- .cval_break_date_cset[[1]][[1]]
        } else if (p_zb %in% 0:3) {
            values_table <- .cval_break_date_cset[[1]][[2 + p_zb]]
        } else {
            stop("ERROR: Invalid value of variable p_zb or p_zf")
        }
    } else if (conf.level == 0.9) {
        if (p_zf == 0) {
            values_table <- .cval_break_date_cset[[2]][[1]]
        } else if (p_zb %in% 0:3) {
            values_table <- .cval_break_date_cset[[2]][[2 + p_zb]]
        } else {
            stop("ERROR: Invalid value of variable p_zb or p_zf")
        }
    } else {
        stop("ERROR: Invalid value of variable conf.level")
    }

    lambda.d <- abs(lambda - 0.5)

    if (p_zf == 0) {
        if (!trend) {
            coef_sup <- values_table$sup_all[, p_zb]
            coef_avg <- values_table$avg_all[, p_zb]
            coef_exp <- values_table$exp_all[, p_zb]
        } else if (trend) {
            coef_sup <- values_table$sup_all[, (4 + p_zb)]
            coef_avg <- values_table$avg_all[, (4 + p_zb)]
            coef_exp <- values_table$exp_all[, (4 + p_zb)]
        } else {
            stop("ERROR: Invalid value of variable trend")
        }
    } else if (p_zb %in% 0:3) {
        if (!trend) {
            coef_sup <- values_table$sup_all[, p_zf]
            coef_avg <- values_table$avg_all[, p_zf]
            coef_exp <- values_table$exp_all[, p_zf]
        } else if (trend) {
            coef_sup <- values_table$sup_all[, (4 - p_zb + p_zf)]
            coef_avg <- values_table$avg_all[, (4 - p_zb + p_zf)]
            coef_exp <- values_table$exp_all[, (4 - p_zb + p_zf)]
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

    return(
        list(
            cval_sup = cval_sup,
            cval_avg = cval_avg,
            cval_exp = cval_exp
        )
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
#' @export
get.p.values.SADF <- function(statistic,
                              n.obs,
                              cr.values) {
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
    return(p.value)
}

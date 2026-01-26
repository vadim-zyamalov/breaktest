#' @title
#' Custom functions for printing results in a nice way.
#'
#' @param x An jbject containing results.
#' @param ... Any additional arguments for [print] function.
#'
#' @keywords internal
#' @export
print.sadf <- function(x, ...) {
    test_name <- NULL

    if ("SADF.value" %in% names(x)) {
        test_statistic <- x$SADF.value
        test_name <- "SADF"
    } else if ("STADF.value" %in% names(x)) {
        test_statistic <- x$STADF.value
        test_name <- "STADF"
    } else if ("GSADF.value" %in% names(x)) {
        test_statistic <- x$GSADF.value
        test_name <- "GSADF"
    } else if ("GSTADF.value" %in% names(x)) {
        test_statistic <- x$GSTADF.value
        test_name <- "GSTADF"
    }

    if (!is.null(test_name)) {
        cat(test_name, "test statistic:", test_statistic, "\n")
    }

    if ("supBZ.value" %in% names(x)) {
        cat("supBZ statistic:", x$supBZ.value, "\n")
    }

    if ("p.value" %in% names(x)) {
        cat("P-value:", x$p.value, "\n")
    }

    if ("is.explosive" %in% names(x)) {
        cat(
            "Current process is",
            ifelse(x$is.explosive == 1, "", "not"),
            "explosive\n"
        )
    }
}


#' @rdname print.sadf
#' @importFrom stringr str_split
#' @keywords internal
#' @export
print.mdfHLT <- function(x, ...) {
    if (x$const && !x$trend) {
        cat("Model 0: Structural change in intercept\n")
        cat("Y{t}=a0+a1*DU+b0*t+e{t}\n",
            "    where DU=1(t>TB)\n")
    } else if (!x$const && x$trend) {
        cat("Model 1: Structural change in slope\n")
        cat("Y{t}=a0+b0*t+b1*DT+e{t}\n",
            "    where DT=1(t>TB)*(t-TB)\n")
    } else if (x$const && x$trend) {
        cat("Model 2: Structural change in both intercept and slope\n")
        cat("Y{t}=a0+a1*DU+b0*t+b1*DT+e{t}\n",
            "    where DU=1(t>TB) and DT=1(t>TB)*(t-TB)\n")
    }
    cat("\n")

    cat("The break date is estimated in", x$break.time, "\n")
    cat("\n")

    cat("Robust tests for the break with uncertaint over errors",
        "(integrated or stationary):\n\n",
        "\tstat\tc.v.\n")
    for (v in c("HLT", "PY")) {
        cat(sprintf("%-7s\t%.4f\t%.4f\n", v, x[[v]]$stat, x[[v]]$cv))
        if (x[[v]]$stat > x[[v]]$cv)
            cat("reject\n")
        else
            cat("fails to reject\n")
        cat("\n")
    }
    cat("\n")

    cat("Unit root tests:\n\n",
        "\tstat\tc.v.\n")
    for (v in c("DF.GLS", "DF.OLS", "MDF.GLS", "MDF.OLS", "MDF.t")) {
        cat(sprintf("%-7s\t%.4f\t%.4f\n", v, x[[v]]$stat, x[[v]]$cv))
        if (x[[v]]$stat < x[[v]]$cv)
            cat("reject\n")
        else
            cat("fails to reject\n")
        cat("\n")
    }
    cat("\n")

    cat("Testing strategies:\n\n")
    for (v in c("A.HLT", "A.PY", "UR.HLT", "UR.PY")) {
        tmp_str <- stringr::str_split(v, "\\.")
        cat(tmp_str[[1]][1], "*(t_", tmp_str[[1]][2], ", s_alpha): ", sep = "")
        if (x[[v]] == 1)
            cat("reject\n")
        else
            cat("fails to reject\n")
        cat("\n")
    }
}


#' @rdname print.sadf
#' @importFrom stringr str_split
#' @keywords internal
#' @export
print.mdfHLTN <- function(x, ...) {
    cat("\t\tstat\tc.v.\n\n")

    for (v in c(
        "MDF.GLS.1",
        "MDF.GLS.2",
        if (x$breaks == 3) "MDF.GLS.3" else NULL,
        "MDF.OLS.1",
        "MDF.OLS.2",
        if (x$breaks == 3) "MDF.OLS.3" else NULL
    )) {
        cat(sprintf("%-9s:\t%.4f\t%.4f\n", v, x[[v]]$stat, x[[v]]$cv))
        if (x[[v]]$stat < x[[v]]$cv)
            cat("reject\n")
        else
            cat("fails to reject\n")
        cat("\n")
    }
    cat("\n")

    cat(sprintf("UR^%d(s.alpha): ", x$breaks))
    if (x$UR1 == 1)
        cat("reject\n")
    else
        cat("fails to reject\n")

    cat(sprintf("UR^%d(s.alpha, %d): ", x$breaks, x$breaks.star))
    if (x$UR == 1)
        cat("reject\n")
    else
        cat("fails to reject\n")
    cat("\n")
}


#' @rdname print.sadf
#' @keywords internal
#' @export
print.mdfCHLT <- function(x, ...) {
    cat("\t\tstat\tc.v.\t wild c.v.\n\n")
    for (v in c("MZa", "MSB", "MZt", "ADF")) {
        cat(
            sprintf(
                "%s stat:\t%.4f\t%.4f\t%.4f\n",
                v,
                x[[v]]$stat,
                x[[v]]$cv,
                x[[v]]$cv.bootstrap
            )
        )
    }
}


#' @rdname print.sadf
#' @keywords internal
#' @export
print.cointGH <- function(x, ...) {
    cat("Gregory-Hansen tests\n")
    cat("\t\tstat\tc.v.\tasymptotic\n\n")
    for (v in c("Za", "Zt", "ADF")) {
        cat(
            sprintf(
                "%s stat:\t%.4f\t%.4f\t%.4f\n",
                v,
                x[[v]]$stat,
                x[[v]]$cv,
                x[[v]]$asy.cv
            )
        )
    }
}


#' @rdname print.sadf
#' @keywords internal
#' @export
print.robustUR <- function(x, ...) {
    cat(sprintf("Estimated break moment: %d\n", x$break.time))
    cat("\t\tstat\tc.v.\t\n\n")
    for (v in c("HLT", "PY", "DF.GLS", "DF.OLS",
                "MDF.GLS", "MDF.OLS", "MDF.t")) {
        cat(
            sprintf(
                "%s stat:\t%.4f\t%.4f\n",
                v,
                x[[v]]$stat,
                x[[v]]$cv
            )
        )
    }
}


#' @rdname print.sadf
#' @keywords internal
#' @export
print.robustURN <- function(x, ...) {
    cat(sprintf("Estimated break moment: %d\n", x$breaks.star))
    cat("\t\tstat\tc.v.\n\n")
    for (v in c("MDF.GLS.1", "MDF.OLS.1", "MDF.GLS.2", "MDF.OLS.2")) {
        cat(
            sprintf(
                "%s stat:\t%.4f\t%.4f\n",
                v,
                x[[v]]$stat,
                x[[v]]$cv
            )
        )
    }
}

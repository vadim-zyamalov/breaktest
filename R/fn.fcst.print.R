#' @title
#' Printing the results of the forecasting routines
#' @name fcst.print
#'
#' @param x An object of the corresponding class.
#' @param ... Unused, kept for compatibility with the generic.
#'
#' @keywords internal
NULL


#' @rdname fcst.print
#' @export
print.bt_fcstWald <- function(x, ...) {
  cat("\nSup/exp-Wald structural change tests\n")
  cat("Tested coefficients:", x$df, "  covariance:", x$vcov.type, "\n")
  print(data.frame(statistic = x$statistics, p.value = x$p.value,
                   critical.value = x$critical.value, reject = x$reject))
  if (anyNA(x$p.value)) cat("Tail approximation:", paste(unique(x$p.value.status), collapse = "; "), "\n")
  cat("Wald date:", x$bp, "  least-squares date:", x$bp.LS, "\n\n")
  invisible(x)
}


#' @rdname fcst.print
#' @export
print.bt_fcstFGLS <- function(x, ...) {
  cat("\nAltansukh-Osborn break inference and forecast\n")
  cat("Method:", x$inference, "  stepwise:", x$stepwise, "\n")
  cat("Coefficient break:", x$bp, "  variance break:", x$variance.bp, "\n")
  cat("Changing columns:", if (length(x$break.vars)) paste(x$break.vars, collapse = ", ") else "none", "\n")
  if (isTRUE(x$average)) {
    cat("Confidence-set averaging:", if (isTRUE(x$averaged)) length(x$break.set) else 0L,
        "dates;", x$conf.method, "confidence", x$conf.level, "\n")
    if (!isTRUE(x$averaged)) cat("  Full-sample forecast: coefficient stability was not rejected.\n")
    if (isTRUE(x$averaged) && x$stepwise) cat("  Changing columns above refer to the estimated-date fit; see component.break.vars.\n")
  }
  cat("Forecast:", format(x$forecast, digits = 6), "\n")
  cat(x$inference.note, "\n\n")
  invisible(x)
}


#' @rdname fcst.print
#' @export
print.bt_fcstWindow <- function(x, ...) {
  cat("\nEstimation window selection\n")
  cat("===========================\n")
  cat("Method             :", x$method, "\n")
  cat("Sample size        :", x$N, "\n")
  cat("Window starts at   :", x$start, "\n")
  cat("Window length      :", x$R, "\n")

  if (!is.null(x$details$bp)) {
    cat("Break date used    :", x$details$bp, "\n")
  }
  if (!is.null(x$details$F.stat)) {
    cat("Wald-scaled F      :", format(x$details$F.stat, digits = 4), "\n")
  }
  if (!is.null(x$details$n.breaks)) {
    cat("Selected breaks    :", x$details$n.breaks, "\n")
  }
  cat("\n")

  invisible(x)
}


#' @rdname fcst.print
#' @export
print.bt_fcstBP <- function(x, ...) {
  .yn <- function(z) if (is.na(z)) "not calibrated"
    else if (isTRUE(z)) "reject" else "do not reject"

  cat("\nBoot-Pick test of equal forecast accuracy\n")
  cat("=========================================\n")
  cat("H0: ", switch(
    x$target,
    postbreak = "post-break and full-sample forecasts are equally accurate",
    combined = "combined and full-sample forecasts are equally accurate"
  ), "\n\n", sep = "")

  cat("Estimated break    :", x$bp, " (tau = ",
      format(x$tau, digits = 3), ")\n", sep = "")
  cat("Boundary zeta^1/2  :", format(x$zeta.sqrt, digits = 4), "\n")
  cat("Nominal size       :", x$level, "\n")
  cat("Covariance         :", x$vcov.type, "\n\n")
  if (!is.null(x$inference)) cat("Calibration        :", x$inference, "\n\n")

  cat(sprintf(
    "  %-4s %10s %10s   %s\n",
    "stat", "value", "cv", "decision"
  ))
  cat(sprintf(
    "  %-4s %10.4f %10.4f   %s\n",
    "W", x$W, x$cv.W, .yn(x$reject.W)
  ))
  cat(sprintf(
    "  %-4s %10.4f %10.4f   %s\n",
    "S", x$S, x$cv.S, .yn(x$reject.S)
  ))

  cat("\nForecasts\n")
  cat("  full sample      :", format(x$forecast["full"], digits = 6), "\n")
  cat("  post-break       :", format(x$forecast["postbreak"], digits = 6), "\n")
  cat("  combined         :", format(x$forecast["combined"], digits = 6), "\n")
  cat("  weight on full   :", format(x$weight.full, digits = 4), "\n\n")

  invisible(x)
}


#' @rdname fcst.print
#' @export
print.bt_fcstConfset <- function(x, ...) {
  cat("\nEo-Morley confidence set for a break date\n")
  cat("========================================\n")
  cat("Confidence level   :", x$level, "\n")
  cat("Estimated break    :", x$bp, "\n")
  cat("Error variance     :", x$variance, "\n")
  cat("Scale covariance   :", x$vcov.type, "\n")
  cat("LR critical value  :", format(x$critical.value, digits = 6), "\n")
  cat("Accepted dates     :", length(x$set), "of", nrow(x$profile), "\n")
  cat("Contiguous components (last pre-break dates):\n")
  print(x$intervals, row.names = FALSE)
  cat("\n")
  invisible(x)
}


#' @rdname fcst.print
#' @export
print.bt_fcstWGLS <- function(x, ...) {
  cat("\nWGLS with post-break leave-one-out CV\n")
  cat("=====================================\n")
  cat("Sample size        :", x$N, "\n")
  cat("Break date         :", x$bp, "\n")
  if (!is.null(x$q)) {
    cat("Sigma pre / post   :", format(x$q, digits = 4), "\n")
    cat("Kernel gamma       :", format(x$gamma, digits = 4), "\n")
  } else {
    cat("Weight selection   : direct gamma.star search\n")
  }
  cat("Effective gamma*   :", format(x$gamma.star, digits = 4), "\n")
  cat("CV mean sq. error  :", format(x$cv.min, digits = 6), "\n")
  cat("Forecast           :", format(x$forecast, digits = 6), "\n\n")

  invisible(x)
}


#' @rdname fcst.print
#' @export
print.bt_fcstStein <- function(x, ...) {
  cat("\nStein-like combined estimator\n")
  cat("=============================\n")
  cat("Break date         :", x$bp, "\n")
  cat("Hausman statistic  :", format(x$H, digits = 4), "\n")
  cat("Shrinkage tau      :", x$tau, "\n")
  cat("Weight on full     :", format(x$alpha, digits = 4), "\n")
  if (!is.null(x$vcov.type)) cat("Covariance         :", x$vcov.type, "\n")
  cat("Risk eligibility   :", x$risk.eligible, "\n")
  cat("Risk conditions    :", x$risk.conditions, "\n")
  cat("  Eligibility does not establish dominance for these data.\n")
  cat("Sigma pre / post   :",
      format(x$sigma["pre"], digits = 4), "/",
      format(x$sigma["post"], digits = 4), "\n\n")

  cat("Forecasts\n")
  cat("  full sample      :", format(x$forecast.full, digits = 6), "\n")
  cat("  post-break       :", format(x$forecast.post, digits = 6), "\n")
  cat("  combined         :", format(x$forecast, digits = 6), "\n\n")

  invisible(x)
}


#' @rdname fcst.print
#' @export
print.bt_fcstCombine <- function(x, ...) {
  cat("\nCombined forecast\n")
  cat("=================\n")
  cat("Number of models   :", length(x$forecasts), "\n")

  if (!is.null(x$R.all)) {
    cat("Window lengths     :", min(x$R.all), "...", max(x$R.all), "\n")
  }
  if (!is.null(x$set)) {
    dates <- if (is.null(x$break.set)) x$set else x$break.set
    cat("Break date set     :", min(dates), "...", max(dates),
        " (", length(dates), " dates)\n", sep = "")
  }

  cat("Forecast           :", format(x$forecast, digits = 6), "\n")
  cat("Spread of forecasts:",
      format(min(x$forecasts), digits = 4), "...",
      format(max(x$forecasts), digits = 4), "\n\n")

  invisible(x)
}


#' @rdname fcst.print
#' @export
print.bt_fcst <- function(x, ...) {
  cat("\nForecast under structural breaks\n")
  cat("================================\n")
  cat("Method             :", x$method, "\n")
  cat("Sample size        :", x$N, "\n")
  if (!is.null(x$h)) cat("Horizon / type     :", x$h, "/", x$type, "\n")

  if (!is.null(x$bp)) {
    cat("Break date         :", x$bp, "\n")
  }
  if (!is.null(x$weights)) {
    eff <- 1 / sum(x$weights^2)
    cat("Effective n. obs   :", format(eff, digits = 4), "\n")
  }

  cat("Forecast           :", format(x$forecast, digits = 6), "\n\n")

  invisible(x)
}


#' @rdname fcst.print
#' @export
print.bt_fcstDM <- function(x, ...) {
  cat("\nEqual predictive accuracy test\n")
  cat("==============================\n")
  cat("Statistic          :", format(x$statistic, digits = 4), "\n")
  cat("Two-sided p-value  :", format(x$p.value, digits = 4), "\n")
  cat("Mean loss diff.    :", format(x$mean.diff, digits = 4), "\n")
  cat("Observations       :", x$n.obs, "\n\n")

  invisible(x)
}


#' @rdname fcst.print
#' @export
print.bt_fcstCW <- function(x, ...) {
  cat("\nClark-West adjusted predictive accuracy test\n")
  cat("===========================================\n")
  cat("Alternative        : larger model improves forecast accuracy\n")
  cat("Statistic          :", format(x$statistic, digits = 4), "\n")
  cat("One-sided p-value  :", format(x$p.value, digits = 4), "\n")
  cat("Mean adjusted diff.:", format(x$mean.diff, digits = 4), "\n")
  cat("Observations       :", x$n.obs, "\n\n")
  invisible(x)
}


#' @rdname fcst.print
#' @export
print.bt_fcstFluct <- function(x, ...) {
  cat("\nFluctuation-type test\n")
  cat("=====================\n")
  cat("Rolling window     :", x$window, " (mu = ",
      format(x$mu, digits = 3), ")\n", sep = "")
  cat("Sup statistic      :", format(x$statistic, digits = 4), "\n")
  cat("Attained at        :", x$at, "\n")
  cat("Critical value     :", format(x$cv, digits = 4),
      " (level ", x$level, ")\n", sep = "")
  if (!is.null(x$inference)) cat("Interpretation     :", x$inference, "\n")
  cat("Decision           :",
      if (isTRUE(x$reject)) "reject" else "do not reject", "\n\n")

  invisible(x)
}


#' @rdname fcst.print
#' @export
print.bt_fcstCompare <- function(x, ...) {
  cat("\nPseudo out-of-sample comparison\n")
  cat("===============================\n")
  cat("Forecast origins   :", length(x$origins), "\n")
  if (!is.null(x$h)) cat("Horizon            :", x$h, "\n")
  if (!is.null(x$n.common)) cat("Common valid dates :", x$n.common, "\n")
  if (!is.null(x$failures)) cat("Failed fits        :", nrow(x$failures), "\n")
  cat("Benchmark          :", x$benchmark, "\n\n")

  tbl <- data.frame(
    RMSFE = round(x$RMSFE, 6),
    rel.RMSFE = round(x$rel.RMSFE, 4),
    DM.vs.benchmark = round(x$DM, 3)
  )
  if (!is.null(x$n.available)) tbl$n.available <- x$n.available
  print(tbl)
  cat("\nA relative RMSFE below one means the method beats the benchmark.\n")
  cat("A negative DM statistic favours the method over the benchmark.\n\n")

  invisible(x)
}


#' @rdname fcst.print
#' @export
print.bt_fcstMulti <- function(x, ...) {
  cat("\nMulti-step forecasts under structural breaks\n")
  cat("Method / type:", x$method, "/", x$type, "\n")
  print(data.frame(h = x$h, forecast = unname(x$forecasts)), row.names = FALSE)
  invisible(x)
}


#' @rdname fcst.print
#' @export
print.bt_fcstCompareMulti <- function(x, ...) {
  cat("\nComparison by forecast horizon\n")
  cat("Benchmark:", x$benchmark, "\n")
  print(x$summary, row.names = FALSE)
  cat("n.common is the same evaluation sample for every method within a horizon.\n")
  invisible(x)
}


#' @rdname fcst.print
#' @export
print.bt_fcstExpS <- function(x, ...) {
  cat("\nData-driven exponential discounting\n")
  cat("Method             :", x$method, "\n")
  cat("Gamma / eta        :", format(x$gamma, digits = 6), "/",
      format(x$eta, digits = 6), "\n")
  if (!is.null(x$prior)) cat("CVL prior          :", x$prior, "\n")
  if (!is.null(x$cv.origins)) cat("Validation origins :", length(x$cv.origins), "\n")
  if (!is.null(x$ml)) cat("ML boundary        :", x$ml$boundary, "\n")
  cat("Forecast           :", format(x$forecast, digits = 6), "\n")
  invisible(x)
}


#' @rdname fcst.print
#' @export
print.bt_fcstSelect <- function(x, ...) {
  cat("\nSingle-stage OCMT predictor selection\n")
  cat("Observations       :", x$N, "\n")
  cat("Selected / total   :", length(x$selected), "/", x$k, "\n")
  cat("Covariance         :", x$vcov.type, "\n")
  cat("Threshold          :", format(x$threshold, digits = 5), "\n")
  cat("Selected columns   :", paste(x$selected.names, collapse = ", "), "\n")
  invisible(x)
}


#' @rdname fcst.print
#' @export
print.bt_fcstDGP <- function(x, ...) {
  cat("\nHirano-Wright simulated series\n")
  cat("Design / N / rho   :", x$design, "/", x$N, "/", x$rho, "\n")
  cat("Horizons           :", paste(x$h, collapse = ", "), "\n")
  cat("Jump dates         :", paste(x$breaks, collapse = ", "), "\n")
  invisible(x)
}


#' @rdname fcst.print
#' @export
print.bt_fcstSimulation <- function(x, ...) {
  cat("\nForecast simulation experiment\n")
  cat("Replications       :", nrow(x$grid), "\n")
  cat("Seed               :", if (is.null(x$settings$seed)) "NULL" else x$settings$seed, "\n")
  cat("Unconditional risks require every replication to succeed.\n")
  cat("Failed fits        :", nrow(x$failures), "\n")
  if (!is.null(x$score.scope)) cat("Score scope        :", x$score.scope, "\n")
  print(x$summary[, c("N", "design", "rho", "h", "method", "n.common", "RMSFE",
                      "scaled.risk", "MCSE.scaled.risk")], row.names = FALSE)
  invisible(x)
}

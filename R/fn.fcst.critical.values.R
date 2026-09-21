#' @title
#' Tabulated quantities of the Boot and Pick (2020) forecast break test
#' @name fcst.BP.tables
#'
#' @description
#' Boot and Pick (2020) test the null of *equal forecast accuracy* between the
#' post-break (or the combined) forecast and the full-sample forecast, rather
#' than the null of no break. Under the MSFE loss the critical break magnitude
#' is non-zero, so the critical values differ substantially from the usual
#' sup-Wald ones of Andrews (1993).
#'
#' Two objects are needed at every break fraction \eqn{\tau_b}:
#' * `zeta.sqrt`: the standardised break magnitude \eqn{\zeta^{1/2}} at which
#'   the two forecasts achieve equal MSFE, used to form the `S` statistic,
#' * critical values of the `W` and `S` statistics.
#'
#' The tables below are Tables 1 and 2 of the published paper, obtained by the
#' authors from 50 000 simulated Brownian motions on a grid with trimming 0.15
#' and step 0.05. They are used by default; [fcst.BP.simulate] regenerates them
#' for a different trimming or a finer grid.
#'
#' @references
#' Boot, Tom, and Andreas Pick.
#' “Does Modeling a Structural Break Improve Forecast Accuracy?”
#' Journal of Econometrics 215 (2020): 35-59. Local author version: 2019.
#' https://doi.org/10.1016/j.jeconom.2019.07.007.
#'
#' Andrews, Donald W. K.
#' “Tests for Parameter Instability and Structural Change
#' with Unknown Change Point.”
#' Econometrica 61, no. 4 (1993): 821–56.
#'
#' @keywords internal
NULL


#' @rdname fcst.BP.tables
.bp.tau.grid <- seq(0.15, 0.85, by = 0.05)


#' @rdname fcst.BP.tables
#'
#' @details
#' `.bp.cv.postbreak` corresponds to Table 1: post-break versus full sample.
.bp.cv.postbreak <- list(
  zeta.sqrt = c(
    2.99, 2.73, 2.55, 2.41, 2.28, 2.17, 2.06, 1.95,
    1.84, 1.75, 1.64, 1.54, 1.43, 1.31, 1.18
  ),
  W = rbind(
    "0.01" = c(
      30.54, 28.84, 27.29, 26.07, 25.00, 24.06, 23.15, 22.29,
      21.46, 20.70, 19.89, 19.08, 18.22, 17.23, 15.82
    ),
    "0.05" = c(
      23.71, 22.29, 20.99, 19.95, 19.04, 18.24, 17.46, 16.74,
      16.03, 15.38, 14.71, 14.02, 13.30, 12.48, 11.37
    ),
    "0.10" = c(
      20.44, 19.16, 17.99, 17.05, 16.22, 15.49, 14.79, 14.13,
      13.49, 12.91, 12.30, 11.68, 11.04, 10.32, 9.36
    )
  ),
  S = rbind(
    "0.01" = c(
      2.76, 2.79, 2.81, 2.83, 2.85, 2.86, 2.86, 2.87,
      2.86, 2.86, 2.85, 2.83, 2.80, 2.74, 2.60
    ),
    "0.05" = c(
      2.12, 2.16, 2.18, 2.20, 2.21, 2.22, 2.23, 2.23,
      2.22, 2.22, 2.20, 2.18, 2.14, 2.08, 1.94
    ),
    "0.10" = c(
      1.78, 1.82, 1.84, 1.86, 1.87, 1.88, 1.89, 1.89,
      1.88, 1.87, 1.86, 1.83, 1.80, 1.73, 1.59
    )
  )
)


#' @rdname fcst.BP.tables
#'
#' @details
#' `.bp.cv.combined` corresponds to Table 2: forecast combination versus full
#' sample.
.bp.cv.combined <- list(
  zeta.sqrt = c(
    2.81, 2.54, 2.36, 2.21, 2.08, 1.97, 1.85, 1.75,
    1.64, 1.55, 1.45, 1.35, 1.25, 1.15, 1.03
  ),
  W = rbind(
    "0.01" = c(
      28.74, 27.08, 25.57, 24.35, 23.30, 22.40, 21.53, 20.74,
      19.95, 19.23, 18.53, 17.81, 17.03, 16.18, 15.02
    ),
    "0.05" = c(
      22.15, 20.78, 19.51, 18.48, 17.60, 16.84, 16.10, 15.43,
      14.76, 14.16, 13.57, 12.97, 12.34, 11.64, 10.74
    ),
    "0.10" = c(
      19.01, 17.78, 16.63, 15.71, 14.92, 14.22, 13.56, 12.95,
      12.34, 11.81, 11.28, 10.74, 10.19, 9.58, 8.82
    )
  ),
  S = rbind(
    "0.01" = c(
      2.82, 2.85, 2.87, 2.89, 2.90, 2.90, 2.91, 2.91,
      2.90, 2.89, 2.87, 2.85, 2.82, 2.76, 2.63
    ),
    "0.05" = c(
      2.18, 2.22, 2.24, 2.25, 2.26, 2.27, 2.27, 2.27,
      2.26, 2.25, 2.23, 2.20, 2.17, 2.11, 1.98
    ),
    "0.10" = c(
      1.85, 1.88, 1.90, 1.92, 1.93, 1.93, 1.93, 1.93,
      1.91, 1.90, 1.88, 1.86, 1.82, 1.76, 1.63
    )
  )
)


#' @title
#' Critical values and equal performance boundary of the Boot-Pick test
#'
#' @param tau A break fraction (or a vector of them) the values are needed at.
#' @param level A nominal size present in the table's probability rows.
#' Built-in tables support 0.01, 0.05 and 0.10.
#' @param target Which comparison is tested:
#' * `postbreak`: post-break forecast versus full sample (Table 1),
#' * `combined`: forecast combination versus full sample (Table 2).
#' @param table The full output of [fcst.BP.simulate], or a target-specific
#' list with `tau`, `zeta.sqrt`, `W`, `S` and `metadata`. Metadata must identify
#' `trim` and `target`; row names of W/S identify upper-tail probabilities.
#' @param trim Requested search trimming. NULL uses the table's metadata.
#'
#' @return A list of `zeta.sqrt`, `cv.W` and `cv.S`, each of the same length as
#' `tau`, plus `calibrated`, `inference` and `metadata`. Compatible tables
#' are linearly interpolated. Missing metadata, unsupported levels, a changed
#' trimming or dates outside the table return NA values and a diagnostic.
#'
#' @keywords internal
get.cv.BootPick <- function(
  tau,
  level = 0.05,
  target = c("postbreak", "combined"),
  table = NULL,
  trim = NULL
) {
  target <- match.arg(target)
  .fcst.test.level(level)
  if (!is.numeric(tau) || !length(tau) || any(!is.finite(tau)) ||
      any(tau <= 0 | tau >= 1)) stop("ERROR! get.cv.BootPick: tau must be in (0,1)")
  if (!is.null(trim) && (!is.numeric(trim) || length(trim) != 1L ||
      !is.finite(trim) || trim <= 0 || trim >= .5)) {
    stop("ERROR! get.cv.BootPick: trim must be in (0,.5)")
  }

  if (is.null(table)) {
    tbl <- switch(
      target,
      postbreak = .bp.cv.postbreak,
      combined = .bp.cv.combined
    )
    tbl$tau <- .bp.tau.grid
    tbl$metadata <- list(trim = .15, target = target, source = "Boot and Pick (2020), Tables 1-2",
                         levels = c(.01, .05, .10), boundary.method = "exact-process")
  } else {
    if (!is.list(table)) stop("ERROR! get.cv.BootPick: table must be a list")
    tbl <- if (!is.null(table$postbreak) || !is.null(table$combined)) table[[target]] else table
  }
  if (!is.list(tbl) || !all(c("tau", "zeta.sqrt", "W", "S") %in% names(tbl))) {
    stop("ERROR! get.cv.BootPick: target table requires tau, zeta.sqrt, W and S")
  }
  grid <- tbl$tau
  if (!is.numeric(grid) || length(grid) < 2L || any(!is.finite(grid)) ||
      any(grid <= 0 | grid >= 1) || is.unsorted(grid, strictly = TRUE) ||
      !is.numeric(tbl$zeta.sqrt) || length(tbl$zeta.sqrt) != length(grid) ||
      any(is.infinite(tbl$zeta.sqrt)) || any(tbl$zeta.sqrt < 0, na.rm = TRUE)) {
    stop("ERROR! get.cv.BootPick: invalid tau grid or equal-risk boundary")
  }
  probabilities <- NULL
  for (type in c("W", "S")) {
    values <- tbl[[type]]
    if (!is.matrix(values) || !is.numeric(values) || !nrow(values) ||
        ncol(values) != length(grid) || any(is.infinite(values)) ||
        is.null(rownames(values))) {
      stop("ERROR! get.cv.BootPick: W and S must be numeric matrices with probability row names")
    }
    rows <- suppressWarnings(as.numeric(rownames(values)))
    if (any(!is.finite(rows)) || any(rows <= 0 | rows >= 1) || anyDuplicated(rows) ||
        (type == "W" && any(values < 0, na.rm = TRUE))) {
      stop("ERROR! get.cv.BootPick: invalid probability rows or negative W critical values")
    }
    if (is.null(probabilities)) probabilities <- rows
    else if (!identical(rows, probabilities)) {
      stop("ERROR! get.cv.BootPick: W and S probability rows must agree")
    }
    ordered <- values[order(rows), , drop = FALSE]
    if (nrow(ordered) > 1L && any(apply(ordered, 2, diff) > 1e-10, na.rm = TRUE)) {
      stop("ERROR! get.cv.BootPick: critical values must decrease as tail probability increases")
    }
  }
  metadata <- tbl$metadata
  reason <- NULL
  if (!is.list(metadata) || !all(c("trim", "target") %in% names(metadata))) {
    reason <- "custom table lacks metadata identifying trim and target"
  } else if (!is.numeric(metadata$trim) || length(metadata$trim) != 1L ||
             !is.finite(metadata$trim) || metadata$trim <= 0 || metadata$trim >= .5 ||
             !is.character(metadata$target) || length(metadata$target) != 1L ||
             is.na(metadata$target) || !metadata$target %in% c("postbreak", "combined")) {
    stop("ERROR! get.cv.BootPick: invalid table metadata")
  } else if (metadata$target != target) {
    reason <- "table target differs from the requested comparison"
  } else if (!is.null(trim) && abs(metadata$trim - trim) > 1e-10) {
    reason <- "table trimming differs from the requested search region"
  } else if (grid[1] > metadata$trim + 1e-10 ||
             tail(grid, 1L) < 1 - metadata$trim - 1e-10) {
    reason <- "table grid does not cover its complete search region"
  }
  row <- which(abs(probabilities - level) <= 1e-12)
  if (is.null(reason) && !length(row)) reason <- "requested significance level is not tabulated"
  if (is.null(reason) && any(tau < min(grid) - 1e-10 | tau > max(grid) + 1e-10)) {
    reason <- "requested break fraction lies outside the table grid"
  }
  if (is.null(reason) && (anyNA(tbl$zeta.sqrt) || anyNA(tbl$W[row, ]) || anyNA(tbl$S[row, ]))) {
    reason <- "equal-risk boundary or requested critical values are incomplete"
  }
  if (!is.null(reason)) {
    return(list(zeta.sqrt = rep(NA_real_, length(tau)), cv.W = rep(NA_real_, length(tau)),
                cv.S = rep(NA_real_, length(tau)), calibrated = FALSE,
                inference = paste("not calibrated:", reason), metadata = metadata))
  }

  list(
    zeta.sqrt = .grid.interp(tbl$zeta.sqrt, grid, tau),
    cv.W = .grid.interp(tbl$W[row, ], grid, tau),
    cv.S = .grid.interp(tbl$S[row, ], grid, tau),
    calibrated = TRUE,
    inference = "asymptotic single-break equal-forecast-risk calibration",
    metadata = metadata
  )
}


#' @title
#' Critical values of the Giacomini-Rossi Fluctuation test
#'
#' @description
#' Two-sided critical values of \eqn{\sup_t |F_{t,m}|} tabulated by Giacomini
#' and Rossi (2010) as a function of \eqn{\mu = m / P}, the ratio of the rolling
#' window to the number of out-of-sample observations.
#'
#' In the replication package this table is duplicated: it appears both inside
#' `FluctuationL.m` and as a standalone `GiacominiRossiCV.m`. Here it is stored
#' once.
#'
#' @param mu A ratio of the rolling window size to the evaluation sample size.
#' @param level A nominal size, either 0.05 or 0.10.
#'
#' @return A scalar critical value.
#' @details Values between the tabulated window fractions are linearly
#' interpolated. Only mu in `[0.1,0.9]` is supported; extrapolation is refused.
#' Interpolation and the finite Monte Carlo tables give approximate, not
#' exact, asymptotic calibration.
#'
#' @references
#' Giacomini, Raffaella, and Barbara Rossi.
#' “Forecast Comparisons in Unstable Environments.”
#' Journal of Applied Econometrics 25, no. 4 (2010): 595–620.
#'
#' @keywords internal
get.cv.fluctuation <- function(mu, level = 0.05) {
  tbl <- rbind(
    "0.1" = c(3.393, 3.170),
    "0.2" = c(3.179, 2.948),
    "0.3" = c(3.012, 2.766),
    "0.4" = c(2.890, 2.626),
    "0.5" = c(2.779, 2.500),
    "0.6" = c(2.634, 2.356),
    "0.7" = c(2.560, 2.252),
    "0.8" = c(2.433, 2.130),
    "0.9" = c(2.248, 1.950)
  )
  colnames(tbl) <- c("0.05", "0.10")

  .fcst.evaluation.table.args(mu, level)
  .grid.interp(tbl[, if (level == .05) 1L else 2L], seq(.1, .9, .1), mu)
}


#' @title
#' Critical values of the Fluctuation rationality test
#'
#' @description
#' Critical values of \eqn{\sup_t W_{t,m}} for the rolling-window forecast
#' rationality test of Rossi and Sekhposyan (2016), tabulated for 5% and 10%
#' significance as a function of the number of restrictions and of
#' \eqn{\mu = m / P}. These are the survey/model-free values in Table A.1c,
#' Panels B and A, of the authors' November 2014 paper and appendix.
#' Table 1 Panel C contains the 5% values, despite the source MATLAB
#' comment describing them as 10%.
#'
#' @param n.reg A number of coefficients tested (rows of the table).
#' @param mu A ratio of the rolling window size to the evaluation sample size.
#' @param level A nominal size, either 0.05 or 0.10.
#'
#' @return A scalar critical value. Linear interpolation between tabulated
#' mu values in `[0.1,0.9]` gives approximate asymptotic calibration. The
#' original Monte Carlo table is retained, including its simulation noise;
#' the values are not silently replaced by monotone fits.
#'
#' @references
#' Rossi, Barbara, and Tatevik Sekhposyan.
#' “Forecast Rationality Tests in the Presence of Instabilities,
#' with Applications to Federal Reserve and Survey Forecasts.”
#' Journal of Applied Econometrics 31, no. 3 (2016): 507–32.
#' https://doi.org/10.1002/jae.2440.
#' https://crei.cat/wp-content/uploads/users/working-papers/Rossi_forecast_rat.pdf.
#'
#' @keywords internal
get.cv.rationality <- function(n.reg, mu, level = 0.05) {
  .fcst.evaluation.table.args(mu, level)
  if (!is.numeric(n.reg) || length(n.reg) != 1L || !is.finite(n.reg) ||
      n.reg != floor(n.reg) || n.reg < 1 || n.reg > 10) {
    stop("ERROR! get.cv.rationality: n.reg must be an integer between 1 and 10")
  }
  tbl <- if (level == .05) rbind(
    c(11.8290, 10.5637, 8.9252, 8.1468, 8.1409, 7.2803, 6.4978, 6.0837, 5.4695),
    c(14.9966, 13.0846, 12.8141, 10.9084, 11.1314, 9.9386, 9.1724, 9.0589, 7.8305),
    c(17.6768, 15.7548, 15.0608, 13.4383, 13.2113, 12.6018, 10.9597, 10.8426, 9.4727),
    c(19.8434, 17.6051, 17.0158, 16.3186, 15.1404, 14.7573, 13.5928, 13.1087, 10.8243),
    c(21.7091, 20.4659, 18.7186, 18.2152, 17.1092, 15.6317, 15.4842, 13.9418, 13.6335),
    c(24.2721, 22.4870, 20.9717, 20.2839, 20.2971, 17.8602, 16.5583, 15.4633, 14.4789),
    c(26.2869, 24.2644, 22.8543, 21.6818, 20.5974, 20.1200, 19.0697, 17.7064, 15.9126),
    c(28.3030, 25.7461, 24.3315, 23.4497, 22.4328, 21.1563, 20.3632, 19.1440, 18.1475),
    c(29.5489, 27.9249, 26.8101, 25.2662, 24.2510, 22.7821, 21.7109, 20.2745, 19.7147),
    c(31.7548, 29.4709, 27.5980, 27.0357, 25.3011, 25.3250, 23.4556, 22.6180, 21.6647)
  ) else rbind(
    c(10.0909, 8.8274, 7.7116, 6.9555, 6.4272, 5.8410, 4.9404, 4.8508, 4.0096),
    c(13.2456, 11.4773, 10.7955, 9.6482, 9.3648, 8.3442, 7.7478, 7.4669, 6.2243),
    c(15.9915, 14.2049, 13.3396, 11.6461, 11.4939, 10.4839, 9.3900, 8.9699, 7.9423),
    c(18.4447, 15.5897, 15.1254, 13.8661, 13.2415, 12.7312, 11.5331, 10.7335, 9.3509),
    c(19.9690, 18.2447, 16.7190, 15.7116, 15.0672, 14.1355, 13.1798, 12.1317, 11.4857),
    c(22.3413, 20.3183, 19.1924, 18.0867, 17.3395, 15.6658, 14.7640, 13.4408, 12.3823),
    c(24.3462, 22.3182, 21.0891, 20.0822, 18.4705, 17.7375, 16.8276, 15.6473, 13.9220),
    c(26.5930, 23.8825, 22.4070, 21.6763, 20.1078, 18.7631, 18.1456, 17.0475, 15.8832),
    c(27.8409, 26.0117, 24.5040, 23.3912, 21.4939, 20.6265, 19.2173, 18.0560, 17.2398),
    c(29.2933, 27.1103, 25.8483, 24.6502, 23.4647, 22.5659, 20.6670, 20.0081, 18.4020)
  )
  .grid.interp(tbl[n.reg, ], seq(.1, .9, .1), mu)
}


.fcst.evaluation.table.args <- function(mu, level) {
  if (!is.numeric(mu) || length(mu) != 1L || !is.finite(mu) || mu < .1 || mu > .9) {
    stop("ERROR! forecast evaluation: tabulated window fraction must be in [0.1,0.9]")
  }
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      !level %in% c(.05, .10)) {
    stop("ERROR! forecast evaluation: level must be 0.05 or 0.10")
  }
  invisible(NULL)
}

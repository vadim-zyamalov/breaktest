#' @title
#' One-covariate-at-a-time selection before forecast estimation
#'
#' @description
#' The single-stage OCMT procedure of Chudik, Pesaran and Sharifvaghefi:
#' regress y on each candidate and the always-kept conditioning variables,
#' without observation down-weighting. Retain candidates with absolute
#' t-ratios above the multiple-testing threshold.
#'
#' @param y A numeric univariate response.
#' @param x A numeric design, including any constant explicitly. The number
#' of candidates may exceed the number of observations.
#' @param always Conditioning columns, by index or unique name. NULL keeps
#' nonzero constant columns. integer(0) specifies no conditioning variables.
#' @param p Nominal multiple-testing size in (0, 1).
#' @param delta Exponent in qnorm(1 - p / (2 * n.candidates^delta)).
#' Must be at least 1 for const and strictly greater than 1 for robust
#' covariance. An omitted delta defaults to 1 for const and 1.1 otherwise;
#' delta.source records which default, or the user's explicit value, was used.
#' @param vcov.type "const" uses RSS/N, as in the paper. "HC0", "HC3" and
#' "HAC" use the module's sandwich covariance. Robust standard errors with
#' delta > 1 are explicitly proposed in the source's Section 5.2, footnote 3.
#' @param bandwidth,time.index HAC settings as in [HAC.variance]. Original
#' integer time positions are retained after complete-case filtering.
#'
#' @details
#' This is the single-stage procedure in Section 3, equations (4)-(5), of
#' the parameter-instability paper, which excludes hidden signals. It is not
#' the multi-stage OCMT for hidden signals in Chudik et al. (2018), nor the
#' grouped Ibragimov-Muller test mentioned in comments in ac/application_*.m.
#' The paper's normal critical value replaces MATLAB's t quantile with 10000
#' degrees of freedom. HAC robustness does not itself establish the paper's
#' variable-selection guarantees under an arbitrary dependence process.
#' The source's signal separation, conditioning-set, exogeneity, tail and
#' dependence assumptions remain necessary. The robust recommendation is
#' not a general high-dimensional HAC theorem. Selection is unweighted;
#' subsequent down-weighted forecast estimation has direct support in
#' Section 6.2. Repeating OCMT inside other named adaptive procedures needs
#' separate support and is restricted by [fcst.forecast].
#'
#' All marginal regressions use the same complete cases. Collinearity with
#' the conditioning set is reported, not silently treated as a significant
#' coefficient. The selected joint design can still be rank deficient;
#' its rank is returned and the forecast interface refuses unidentified fits.
#' Retained observations must be consecutive on the original time grid.
#' Leading/trailing incomplete rows may be trimmed; internal gaps are refused.
#'
#' @return A bt_fcstSelect with selected column indices/names, conditioning
#' and candidate columns, the full table of t-ratios, threshold and reasons,
#' complete-case rows, and the rank of the selected joint design.
#'
#' @references
#' Chudik, A., M. H. Pesaran and M. Sharifvaghefi (2024).
#' "Variable selection in high dimensional linear regressions with
#' parameter instability." Journal of Econometrics 246, 105900.
#' Section 3, equations (4)-(5), and Section 5.2, footnote 3, of the July 2024
#' revision (Dallas Fed Globalization Institute Working Paper 394, version 3;
#' cover dated August 2024)
#' (earlier title: "Variable Selection and Forecasting in High Dimensional
#' Linear Regressions with Parameter Instability").
#' https://doi.org/10.1016/j.jeconom.2024.105900.
#' https://www.dallasfed.org/~/media/documents/research/international/wpapers/2020/0394r3.pdf.
#'
#' @importFrom stats qnorm
#' @export
fcst.select <- function(y, x, always = NULL, p = 0.1, delta = 1,
                        vcov.type = c("const", "HC0", "HC3", "HAC"),
                        bandwidth = NULL, time.index = NULL) {
  vcov.type <- match.arg(vcov.type)
  delta.source <- if (missing(delta)) {
    if (vcov.type == "const") "default.const" else "default.robust"
  } else "user"
  if (missing(delta) && vcov.type != "const") delta <- 1.1
  y <- as.matrix(y)
  x <- as.matrix(x)
  if (!is.numeric(y) || ncol(y) != 1L || !is.numeric(x) || !ncol(x) ||
      nrow(x) != nrow(y) || any(is.infinite(y)) || any(is.infinite(x))) {
    stop("ERROR! fcst.select: finite numeric y and a conformable design required")
  }
  time.index <- .fcst.time.index(time.index, nrow(y))
  rows <- which(is.finite(y[, 1]) & rowSums(!is.finite(x)) == 0L)
  y <- y[rows, , drop = FALSE]
  x <- x[rows, , drop = FALSE]
  time.index <- time.index[rows]
  .fcst.require.regular(time.index, "fcst.select")
  N <- nrow(x)
  if (N < 3L) stop("ERROR! fcst.select: not enough complete observations")
  k <- ncol(x)
  labels <- colnames(x)
  if (is.null(labels)) labels <- paste0("x", seq_len(k))
  if (anyNA(labels) || any(labels == "") || anyDuplicated(labels)) {
    stop("ERROR! fcst.select: design column names must be unique and nonempty")
  }
  if (!is.numeric(p) || length(p) != 1L || !is.finite(p) || p <= 0 || p >= 1 ||
      !is.numeric(delta) || length(delta) != 1L || !is.finite(delta) || delta < 1) {
    stop("ERROR! fcst.select: p must be in (0,1) and delta must be at least 1")
  }
  if (vcov.type != "const" && delta <= 1) {
    stop("ERROR! fcst.select: robust covariance requires delta greater than 1")
  }
  if (is.null(always)) {
    always <- which(vapply(seq_len(k), function(j) {
      x[1, j] != 0 && all(x[, j] == x[1, j])
    }, TRUE))
  } else if (is.character(always)) {
    always <- match(always, labels)
  }
  if (!is.numeric(always) || anyNA(always) || any(!is.finite(always)) ||
      any(always != floor(always)) || any(always < 1L | always > k) ||
      anyDuplicated(always)) stop("ERROR! fcst.select: invalid always columns")
  always <- sort(as.integer(always))
  if (N <= length(always) + 2L ||
      (length(always) && qr(x[, always, drop = FALSE])$rank < length(always))) {
    stop("ERROR! fcst.select: conditioning design is too large or rank deficient")
  }
  candidates <- setdiff(seq_len(k), always)
  nc <- length(candidates)
  # Upper-tail log probability avoids rounding 1 - p/(2*N^delta) to one.
  threshold <- if (nc) qnorm(log(p / 2) - delta * log(nc),
                             lower.tail = FALSE, log.p = TRUE) else Inf
  stats <- data.frame(column = candidates, variable = labels[candidates],
                      coefficient = rep(NA_real_, nc), se = rep(NA_real_, nc),
                      statistic = rep(NA_real_, nc), selected = rep(FALSE, nc),
                      reason = rep(NA_character_, nc))
  for (i in seq_along(candidates)) {
    cols <- c(always, candidates[i])
    xx <- x[, cols, drop = FALSE]
    if (qr(xx)$rank < ncol(xx)) {
      stats$reason[i] <- "collinear with conditioning variables"
      next
    }
    fit <- OLS.reg(y, xx)
    V <- .fcst.vcov(xx, fit$residuals, vcov.type, bandwidth, time.index)
    if (vcov.type == "const") V <- V * (N - ncol(xx)) / N
    j <- ncol(xx)
    variance <- V[j, j]
    if (!is.finite(variance) || variance < 0) {
      stats$reason[i] <- "invalid coefficient variance"
      next
    }
    b <- as.numeric(fit$coefficients[j])
    se <- sqrt(variance)
    statistic <- if (se > 0) b / se else if (b == 0) 0 else sign(b) * Inf
    stats$coefficient[i] <- b
    stats$se[i] <- se
    stats$statistic[i] <- statistic
    stats$selected[i] <- abs(statistic) > threshold
  }
  selected <- sort(c(always, candidates[stats$selected]))
  rank <- if (length(selected)) qr(x[, selected, drop = FALSE])$rank else 0L
  structure(list(selected = selected, selected.names = labels[selected],
                  always = always, candidates = candidates, statistics = stats,
                  threshold = threshold, p = p, delta = delta, vcov.type = vcov.type,
                  delta.source = delta.source, assumptions.verified = FALSE,
                  bandwidth = bandwidth, time.index = time.index, rows = rows,
                  N = N, k = k, rank = rank), class = "bt_fcstSelect")
}


.fcst.selection.options <- function(selection) {
  if (is.null(selection)) return(NULL)
  if (!is.list(selection)) stop("ERROR! fcst.forecast: selection must be NULL or a named list")
  .fcst.check.options(selection, c("always", "p", "delta", "vcov.type", "bandwidth"))
  selection
}


# Keep the original column map for recursion, whose lag positions must not
# change when a fold selects a different set of predictors.
.fcst.select.data <- function(data, selection, h = 1L, type = "direct") {
  data$columns <- seq_len(data$k)
  if (is.null(selection)) return(data)
  args <- selection
  if (is.null(args$vcov.type)) args$vcov.type <- if (h > 1L && type == "direct") "HAC" else "const"
  if (args$vcov.type == "HAC" && is.null(args$bandwidth)) {
    args$bandwidth <- max(if (type == "direct") h else 1L, floor(data$N^(1 / 4)))
  }
  selected <- do.call(fcst.select, c(list(y = data$y, x = data$x,
                                         time.index = data$map$target), args))
  if (!length(selected$selected)) {
    stop("ERROR! fcst.forecast: selection retained no regressors; supply an always-kept constant")
  }
  if (selected$rank < length(selected$selected) ||
      data$N <= length(selected$selected) + 1L) {
    stop("ERROR! fcst.forecast: selected joint design is too large or rank deficient")
  }
  data$columns <- selected$selected
  data$x <- data$x[, data$columns, drop = FALSE]
  data$k <- ncol(data$x)
  data$selection <- selected
  data
}


.fcst.expand.coef <- function(coef, columns, k) {
  result <- numeric(k)
  result[columns] <- as.numeric(coef)
  result
}

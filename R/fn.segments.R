#' @title
#' Procedure to minimize the SSR for 1 break point
#'
#' @param beg Start of the sample.
#' @param end End of the sample.
#' @param first.break First possible break point.
#' @param last.break Last possible break point.
#' @param N Total number of observations.
#' @param SSR.data The matrix of recursive SSR values.
#'
#' @return A list of:
#' * `SSR`: Optimal SSR value,
#' * `break.point`: The point of possible break.
#'
#' @references
#' Carrion-i-Silvestre, Josep Lluís, and Andreu Sansó.
#' “Testing the Null of Cointegration with Structural Breaks.”
#' Oxford Bulletin of Economics and Statistics 68, no. 5 (October 2006): 623–46.
#' https://doi.org/10.1111/j.1468-0084.2006.00180.x.
#'
#' @keywords internal
segments.CSS <- function(beg, end, bp.min, bp.max, N, SSR.data) {
  SSR <- matrix(data = Inf, nrow = N, ncol = 1)

  for (bp in bp.min:bp.max) {
    SSR[bp] <- SSR.data[beg, bp] + SSR.data[bp + 1, end]
  }

  list(
    SSR = min(SSR),
    break.point = which.min(SSR)
  )
}


#' @title
#' Find \eqn{m + 1} optimal partitions using sequential procedure
#'
#' @param y (Tx1)-vector of the dependent variable.
#' @param x (Txk)-vector of the explanatory stochastic regressors.
#' @param m Number of breaks.
#' @param width Minimum spacing between the breaks.
#' @param SSR.data Optional matrix of recursive SSR's.
#'
#' @details
#' The sequential procedure by Bai & Perron (2003) works as follows.
#' First we find the first break point by minimizing SSR on the fraction of initial sample.
#' The upper bound is calculated using `width` parameter ensuring that there will be
#' enough observations for all breaks.
#' Then we find optimal two-break segmentations for all possible upper bounds.
#' Repeating this procedure we get an optimal \eqn{m}-breaks segmentation.
#'
#' @return A list of:
#' * optimal SSR,
#' * the vector of break points.
#'
#' @references
#' Bai, Jushan, and Pierre Perron.
#' “Computation and Analysis of Multiple Structural Change Models.”
#' Journal of Applied Econometrics 18, no. 1 (2003): 1–22.
#' https://doi.org/10.1002/jae.659.
#'
#' @keywords internal
#' @export
segments.BP <- function(
  y,
  x,
  m = 1,
  h = 2,
  SSR.data = NULL
) {
  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }

  N <- nrow(y)

  if (is.null(SSR.data)) {
    SSR.data <- SSR.matrix(y, x, h)
  }

  # For one break use the procedure above
  if (m == 1) {
    return(
      segments.CSS(1, N, h, N - h, N, SSR.data)
    )
  }

  optdat <- matrix(
    data = NA,
    nrow = N,
    ncol = m
  )
  optSSR <- matrix(
    data = Inf,
    nrow = N,
    ncol = m
  )

  for (ib in 1:m) {
    if (ib == 1) {
      for (j1 in (2 * h):N) {
        .segments <- segments.CSS(1, j1, h, j1 - h, j1, SSR.data)
        optSSR[j1, 1] <- .segments$SSR
        optdat[j1, 1] <- .segments$break.point
      }
    } else if (ib == m) {
      dvec <- rep(Inf, N)
      for (jb in (ib * h):(N - h)) {
        dvec[jb] <- optSSR[jb, ib - 1] + SSR.data[jb + 1, N]
      }
      optSSR[N, ib] <- min(dvec)
      optdat[N, ib] <- which.min(dvec)
    } else {
      for (jlast in ((ib + 1) * h):N) {
        dvec <- rep(Inf, N)
        for (jb in (ib * h):(jlast - h)) {
          dvec[jb] <- optSSR[jb, ib - 1] + SSR.data[jb + 1, jlast]
        }
        optSSR[jlast, ib] <- min(dvec)
        optdat[jlast, ib] <- which.min(dvec)
      }
    }
  }

  datevec <- numeric(m)
  datevec[m] <- optdat[N, m]
  for (i in seq_len(m - 1)) {
    xx <- m - i
    datevec[xx] <- optdat[datevec[xx + 1], xx]
  }

  list(
    SSR = optSSR[N, m],
    break.point = datevec
  )
}


#' @title
#' Procedure to minimize the GLS-SSR for 1 break point
#'
#' @param y Variable of interest.
#' @param const Whether there is a break in the constant.
#' @param trend Whether there is a break in the trend.
#' @param breaks Number of breaks.
#' @param first.break First possible break point.
#' @param last.break Last possible break point.
#' @param trim Trim value to calculate `first.break` and `last.break`
#' if not provided.
#'
#' @return The point of possible break.
#'
#' @references
#' Skrobotov, Anton.
#' “On Trend Breaks and Initial Condition in Unit Root Testing.”
#' Journal of Time Series Econometrics 10, no. 1 (2018): 1–15.
#' https://doi.org/10.1515/jtse-2016-0014.
#'
#' @keywords internal
segments.GLS <- function(
  y,
  const = FALSE,
  trend = FALSE,
  breaks = 1,
  bp_min = NULL,
  bp_max = NULL,
  trim = 0.15
) {
  if (!is.matrix(y)) {
    y <- as.matrix(y)
  }

  if (breaks < 1) {
    stop("At least one break is needed!")
  }
  if (breaks > 3) {
    stop("More than three breaks are not supported at the moment!")
  }

  N <- nrow(y)

  if (is.null(bp_min)) {
    bp_min <- floor(trim * N)
  }
  if (is.null(bp_max)) {
    bp_max <- floor((1 - trim) * N)
  }
  width <- bp_min

  steps <- c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.975, 1)

  SSR <- Inf
  resBreaks <- rep(0, breaks)

  for (alpha in steps) {
    if (breaks == 1) {
      for (bp1 in bp_min:bp_max) {
        x <- cbind(
          .const(N),
          .trend(N),
          if (const) .du(bp1, N) else NULL,
          if (trend) .dt(bp1, N) else NULL
        )

        c_bar <- N * (alpha - 1)
        r_GLS <- GLS.reg(y, x, c_bar)$residuals
        loopSSR <- sum(r_GLS^2)

        if (loopSSR < SSR) {
          SSR <- loopSSR
          resBreaks <- c(bp1)
        }
      }
    } else if (breaks == 2) {
      for (bp1 in bp_min:(bp_max - width)) {
        for (bp2 in (bp1 + width):bp_max) {
          x <- cbind(
            .const(N),
            .trend(N),
            if (const) .du(bp1, N) else NULL,
            if (trend) .dt(bp1, N) else NULL,
            if (const) .du(bp2, N) else NULL,
            if (trend) .dt(bp2, N) else NULL
          )

          c_bar <- N * (alpha - 1)
          r_GLS <- GLS.reg(y, x, c_bar)$residuals
          loopSSR <- sum(r_GLS^2)

          if (loopSSR < SSR) {
            SSR <- loopSSR
            resBreaks <- c(bp1, bp2)
          }
        }
      }
    } else if (breaks == 3) {
      for (bp1 in bp_min:(bp_max - 2 * width)) {
        for (bp2 in (bp1 + width):(bp_max - width)) {
          for (bp3 in (bp2 + width):bp_max) {
            x <- cbind(
              .const(N),
              .trend(N),
              if (const) .du(bp1, N) else NULL,
              if (trend) .dt(bp1, N) else NULL,
              if (const) .du(bp2, N) else NULL,
              if (trend) .dt(bp2, N) else NULL,
              if (const) .du(bp3, N) else NULL,
              if (trend) .dt(bp3, N) else NULL
            )

            c_bar <- N * (alpha - 1)
            r_GLS <- GLS.reg(y, x, c_bar)$residuals
            loopSSR <- sum(r_GLS^2)

            if (loopSSR < SSR) {
              SSR <- loopSSR
              resBreaks <- c(bp1, bp2, bp3)
            }
          }
        }
      }
    }
  }

  resBreaks
}

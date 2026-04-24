#' @title
#' Procedure to minimize the SSR for 1 break point
#'
#' @param beg Start of the sample.
#' @param end End of the sample.
#' @param first.break First possible break point.
#' @param last.break Last possible break point.
#' @param len Total number of observations.
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
segments.CSS <- function(beg,
                         end,
                         bp.min,
                         bp.max,
                         len,
                         SSR.data) {
  .rss <- matrix(data = Inf, nrow = len, ncol = 1)

  for (bp in bp.min:bp.max) {
    .rss[bp] <- SSR.data[beg, bp] + SSR.data[bp + 1, end]
  }

  final.rss <- min(.rss[bp.min:bp.max])
  final.bp <- (bp.min - 1) + which.min(.rss[bp.min:bp.max])

  list(
    SSR         = final.rss,
    break.point = final.bp
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
  width = 2,
  SSR.data = NULL
) {
  if (!is.matrix(y)) y <- as.matrix(y)
  if (!is.matrix(x)) x <- as.matrix(x)

  N <- nrow(y)

  if (is.null(SSR.data)) {
    SSR.data <- SSR.matrix(y, x, width)
  }

  # For one break use the procedure above
  if (m == 1) {
    return(
      segments.CSS(1, N, width, N - width, N, SSR.data)
    )
  }

  cNvars <- N - (m + 1) * width + 1
  vSSR <- rep(Inf, cNvars)
  mBreaks <- matrix(
    data = 0,
    nrow = cNvars,
    ncol = m
  )

  for (step in 1:m) {
    loopSSR <- rep(Inf, cNvars)

    if (step == 1) {
      for (v in 1:cNvars) {
        upperBorder <- 2 * width + v - 1
        .segments <- segments.CSS(
          1,
          upperBorder,
          width,
          upperBorder - width,
          upperBorder, SSR.data
        )
        vSSR[v] <- .segments$SSR
        mBreaks[v, 1] <- .segments$break.point
      }
    } else if (step == m) {
      for (v in 1:cNvars) {
        loopSSR[v] <- vSSR[v] + SSR.data[step * width + v, N]
      }
      finalSSR <- min(loopSSR)
      finalIdx <- which.min(loopSSR)
      finalBreaks <- mBreaks[finalIdx, ]
      finalBreaks[m] <- step * width + finalIdx - 1
    } else {
      vNewSSR <- rep(Inf, cNvars)
      mNewBreaks <- matrix(
        data = 0,
        nrow = cNvars,
        ncol = m
      )

      # Looping through the possible upperBounds for step-breaks segmentation.
      for (upperBorder in ((step + 1) * width):(N - (m - step) * width)) {
        searchIdx <- upperBorder - (step + 1) * width + 1

        # For every v we calculate a new SSR value as the sum of step-1 breaks
        # segmentation with last break at v and SSR of the rest part till upperBound.
        for (v in 1:cNvars) {
          loopSSR[v] <- vSSR[v] + SSR.data[step * width + v, upperBorder]
        }

        # Look for the minimum loopSSR which corresponds to the optimal step-breaks
        # segmentation for current upper border.
        vNewSSR[searchIdx] <- min(loopSSR)
        minIdx <- which.min(loopSSR)
        mNewBreaks[searchIdx, 1:m] <- mBreaks[minIdx, ]
        mNewBreaks[searchIdx, step] <- step * width + minIdx - 1
      }

      # Update vSSR and mBreaks with optimal step-breaks segments.
      vSSR <- vNewSSR
      mBreaks <- mNewBreaks
    }
  }

  list(
    SSR         = finalSSR,
    break.point = finalBreaks
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
segments.GLS <- function(y,
                         const = FALSE,
                         trend = FALSE,
                         breaks = 1,
                         bp_min = NULL,
                         bp_max = NULL,
                         trim = 0.15) {
  if (!is.matrix(y)) y <- as.matrix(y)

  if (breaks < 1) {
    stop("At least one break is needed!")
  }
  if (breaks > 3) {
    stop("More than three breaks are not supported at the moment!")
  }

  N <- nrow(y)

  if (is.null(bp_min)) bp_min <- floor(trim * N) + 1
  if (is.null(bp_max)) bp_max <- floor((1 - trim) * N) + 1
  width <- bp_min - 1

  steps <- c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.975, 1)

  .rss <- Inf
  .bps <- rep(0, breaks)

  for (alpha in steps) {
    if (breaks == 1) {
      for (bp1 in bp_min:bp_max) {
        du1 <- as.numeric(trend > bp1)
        dt1 <- du1 * (trend - bp1)

        x <- cbind(
          .const(N),
          .trend(N),
          if (const) du1 else NULL,
          if (trend) dt1 else NULL
        )

        c_bar <- N * (alpha - 1)
        resids <- GLS.reg(y, x, c_bar)$residuals

        .rss_loop <- drop(t(resids) %*% resids)

        if (.rss_loop < .rss) {
          .rss <- .rss_loop
          .bps <- c(bp1)
        }
      }
    } else if (breaks == 2) {
      for (bp1 in bp_min:(bp_max - width)) {
        for (bp2 in (bp1 + width):bp_max) {
          du1 <- as.numeric(trend > bp1)
          dt1 <- du1 * (trend - bp1)
          du2 <- as.numeric(trend > bp2)
          dt2 <- du2 * (trend - bp2)

          x <- cbind(
            .const(N),
            .trend(N),
            if (const) du1 else NULL,
            if (trend) dt1 else NULL,
            if (const) du2 else NULL,
            if (trend) dt2 else NULL
          )

          c_bar <- N * (alpha - 1)
          resids <- GLS.reg(y, x, c_bar)$residuals

          .rss_loop <- drop(t(resids) %*% resids)

          if (.rss_loop < .rss) {
            .rss <- .rss_loop
            .bps <- c(bp1, bp2)
          }
        }
      }
    } else if (breaks == 3) {
      for (bp1 in bp_min:(bp_max - 2 * width)) {
        for (bp2 in (bp1 + width):(bp_max - width)) {
          for (bp3 in (bp2 + width):bp_max) {
            du1 <- as.numeric(trend > bp1)
            dt1 <- du1 * (trend - bp1)
            du2 <- as.numeric(trend > bp2)
            dt2 <- du2 * (trend - bp2)
            du3 <- as.numeric(trend > bp3)
            dt3 <- du3 * (trend - bp3)

            x <- cbind(
              .const(N),
              .trend(N),
              if (const) du1 else NULL,
              if (trend) dt1 else NULL,
              if (const) du2 else NULL,
              if (trend) dt2 else NULL,
              if (const) du3 else NULL,
              if (trend) dt3 else NULL
            )

            c_bar <- N * (alpha - 1)
            resids <- GLS.reg(y, x, c_bar)$residuals

            .rss_loop <- drop(t(resids) %*% resids)

            if (.rss_loop < .rss) {
              .rss <- .rss_loop
              .bps <- c(bp1, bp2, bp3)
            }
          }
        }
      }
    }
  }

  .bps
}

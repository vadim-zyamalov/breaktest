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
segments.ols.single <- function(beg,
                                end,
                                bp.min,
                                bp.max,
                                len,
                                SSR.data) {
  .rss <- matrix(data = Inf, nrow = len, ncol = 1)

  for (bp in bp.min:bp.max) {
    .rss[bp] <- SSR.data[beg, bp] + SSR.data[bp + 1, end]
  }

  .rss.min <- min(.rss[bp.min:bp.max])
  .bp <- (bp.min - 1) + which.min(.rss[bp.min:bp.max])

  list(
    SSR         = .rss.min,
    break.point = .bp
  )
}


#' @title
#' Procedure to minimize the SSR for 2 break points
#'
#' @param y A time series of interest.
#' @param model A scalar equal to
#' * 1: for the AA (without trend) model,
#' * 2: for the AA (with trend) model,
#' * 3: for the BB model,
#' * 4: for the CC model,
#' * 5: for the AC-CA model.
#'
#' @return A list of
#' * resid: (Tx1) vector of estimated OLS residuals,
#' * tb1: The first break point,
#' * tb2: The second break point.
#'
#' @references
#' Carrion-i-Silvestre, Josep Lluís, and Andreu Sansó.
#' “The KPSS Test with Two Structural Breaks.”
#' Spanish Economic Review 9, no. 2 (May 16, 2007): 105–27.
#' https://doi.org/10.1007/s10108-006-9017-8.
#'
#' @keywords internal
segments.ols.double <- function(y, model) {
  if (!is.matrix(y)) y <- as.matrix(y)

  n.obs <- nrow(y)

  .resids <- 0
  .rss <- Inf
  .bp1 <- 0
  .bp2 <- 0

  if (1 <= model && model <= 4) {
    for (bp1 in 2:(n.obs - 4)) {
      for (bp2 in (bp1 + 2):(n.obs - 2)) {
        z <- trend.kpss.double(model, n.obs, c(bp1, bp2))
        resids <- .lm.fit(z, y)$residuals
        ssr <- drop(t(resids) %*% resids)
        if (ssr < .rss) {
          .resids <- resids
          .rss <- ssr
          .bp1 <- bp1
          .bp2 <- bp2
        }
      }
    }
  } else if (5 <= model && model <= 7) {
    for (bp1 in 2:(n.obs - 4)) {
      for (bp2 in (bp1 + 2):(n.obs - 2)) {
        z <- trend.kpss.double(model, n.obs, c(bp1, bp2))
        resids <- .lm.fit(z, y)$residuals
        ssr <- drop(t(resids) %*% resids)
        if (ssr < .rss) {
          .resids <- resids
          .rss <- ssr
          .bp1 <- bp1
          .bp2 <- bp2
        }
      }
    }
    for (bp2 in 2:(n.obs - 4)) {
      for (bp1 in (bp2 + 2):(n.obs - 2)) {
        z <- trend.kpss.double(model, n.obs, c(bp1, bp2))
        resids <- .lm.fit(z, y)$residuals
        ssr <- drop(t(resids) %*% resids)
        if (ssr < .rss) {
          .resids <- resids
          .rss <- ssr
          .bp1 <- bp1
          .bp2 <- bp2
        }
      }
    }
  }

  list(
    residuals = .resids,
    tb1       = .bp1,
    tb2       = .bp2
  )
}


#' @title
#' Find \eqn{m + 1} optimal partitions
#'
#' @param y (Tx1)-vector of the dependent variable.
#' @param x (Txk)-vector of the explanatory stochastic regressors.
#' @param m Number of breaks.
#' @param width Minimum spacing between the breaks.
#' @param SSR.data Optional matrix of recursive SSR's.
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
segments.ols.mulitiple <- function(
  y,
  x,
  m = 1,
  width = 2,
  SSR.data = NULL
) {
  if (!is.matrix(y)) y <- as.matrix(y)
  if (!is.matrix(x)) x <- as.matrix(x)

  n.obs <- nrow(y)

  if (is.null(SSR.data)) {
    SSR.data <- SSR.matrix(y, x, width)
  }

  if (m == 1) {
    .segments <- segments.ols.single(
      1, n.obs,
      width, n.obs - width,
      n.obs, SSR.data
    )
    .rss_final <- .segments$SSR
    .bp_final <- .segments$break.point
  } else {
    n_variants <- n.obs - (m + 1) * width + 1
    .rss <- matrix(
      data = Inf,
      nrow = n_variants,
      ncol = 1
    )
    .bp <- matrix(
      data = 0,
      nrow = n_variants,
      ncol = m
    )
    for (step in 1:m) {
      .rss_loop <- matrix(
        data = Inf,
        nrow = n_variants,
        ncol = 1
      )
      if (step == 1) {
        for (v in 1:n_variants) {
          .last_step <- 2 * width + v - 1
          .segments <- segments.ols.single(
            1,
            .last_step,
            width,
            .last_step - width,
            .last_step, SSR.data
          )
          .rss[v, 1] <- .segments$SSR
          .bp[v, 1] <- .segments$break.point
        }
      } else if (step == m) {
        for (v in 1:n_variants) {
          .rss_loop[v, 1] <- .rss[v, 1] +
            SSR.data[step * width + v, n.obs]
        }
        .rss_final <- min(.rss_loop)
        .idx_final <- which.min(.rss_loop)
        .bp_final <- .bp[.idx_final, ]
        .bp_final[m] <- step * width + .idx_final - 1
      } else {
        .rss_new <- matrix(
          data = Inf,
          nrow = n_variants,
          ncol = 1
        )
        .bp_new <- matrix(
          data = 0,
          nrow = n_variants,
          ncol = m
        )
        for (.last_step in ((step + 1) * width):(n.obs - (m - step) * width)) {
          .v_new <- .last_step - (step + 1) * width + 1
          for (v in 1:n_variants) {
            .rss_loop[v, 1] <- .rss[v, 1] +
              SSR.data[step * width + v, .last_step]
          }
          .rss_new[.v_new, 1] <- min(.rss_loop)
          .idx_new <- which.min(.rss_loop)
          .bp_new[.v_new, 1:m] <- .bp[.idx_new, ]
          .bp_new[.v_new, step] <- step * width + .idx_new - 1
        }
        .rss <- .rss_new
        .bp <- .bp_new
      }
    }
  }

  list(
    SSR         = .rss_final,
    break.point = .bp_final
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
segments.gls <- function(y,
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

  n.obs <- nrow(y)
  const <- rep(1, n.obs)
  trend <- 1:n.obs

  if (is.null(bp_min)) bp_min <- floor(trim * n.obs) + 1
  if (is.null(bp_max)) bp_max <- floor((1 - trim) * n.obs) + 1
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
          const,
          trend,
          if (const) du1 else NULL,
          if (trend) dt1 else NULL
        )

        c_bar <- n.obs * (alpha - 1)
        resids <- .GLS(y, x, c_bar)$residuals

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
            const,
            trend,
            if (const) du1 else NULL,
            if (trend) dt1 else NULL,
            if (const) du2 else NULL,
            if (trend) dt2 else NULL
          )

          c_bar <- n.obs * (alpha - 1)
          resids <- .GLS(y, x, c_bar)$residuals

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
              const,
              trend,
              if (const) du1 else NULL,
              if (trend) dt1 else NULL,
              if (const) du2 else NULL,
              if (trend) dt2 else NULL,
              if (const) du3 else NULL,
              if (trend) dt3 else NULL
            )

            c_bar <- n.obs * (alpha - 1)
            resids <- .GLS(y, x, c_bar)$residuals

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

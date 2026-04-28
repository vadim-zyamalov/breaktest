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
robust.tests.single <- function(y,
                                const = FALSE,
                                trend = FALSE,
                                season = FALSE,
                                trim = 0.15) {
  if (!is.matrix(y)) y <- as.matrix(y)

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
    trim = trim
  )

  result <- append(result, list(season = season), 2)
  class(result) <- "robustUR"

  result
}


#' @title
#' A wrapping function around [KP] and [MDF.mlt].
#'
#' @details
#' The code provided is the original Ox code by Skrobotov (2018)
#' ported to R.
#'
#' @param y A series of interest.
#' @param const Whether the constant term should be included.
#' @param season Whether the seasonal adjustment is needed.
#' @param breaks Number of breaks.
#' @param trim Trimming value for a possible break date bounds.
#'
#' @export
robust.tests.multiple <- function(y,
                                  const = FALSE,
                                  season = FALSE,
                                  breaks = 2,
                                  trim = 0.15) {
  if (!is.matrix(y)) y <- as.matrix(y)

  ## Start ##
  N <- nrow(y)

  if (season) {
    SEAS <- cbind(
      .const(N),
      seasonal.dummies(N)
    )
    y <- OLS.reg(y, SEAS)$residuals
  }

  m.star <- KP(
    y = y,
    const = const, breaks = breaks,
    criterion = "aic", trim = trim
  )

  result <- MDF.mlt(
    y = y,
    const = const,
    breaks = breaks,
    breaks.star = m.star,
    trim = trim,
    ZA = FALSE
  )

  result <- append(result, list(season = season), 1)
  class(result) <- "robustURN"

  result
}


#' @title
#' Kejrival-Perron procedure of breaks number detection
#'
#' @details
#' The code provided is the original Ox code by Skrobotov (2018)
#' ported to R.
#'
#' @param y An input series of interest.
#' @param const Whether the break in constant is allowed.
#' @param breaks Number of breaks.
#' @param criterion Needed information criterion: aic, bic, hq or lwz.
#' @param trim A trimming value for a possible break date bounds.
#'
#' @return The estimated optimal break point.
#'
#' @references
#' Kejriwal, Mohitosh, and Pierre Perron.
#' “A Sequential Procedure to Determine the Number of Breaks in Trend
#' with an Integrated or Stationary Noise Component:
#' Determination of Number of Breaks in Trend.”
#' Journal of Time Series Analysis 31, no. 5 (September 2010): 305–28.
#' https://doi.org/10.1111/j.1467-9892.2010.00666.x.
#'
#' @export
KP <- function(y,
               const = FALSE,
               breaks = 1,
               criterion = "aic",
               trim = 0.15) {
  if (!is.matrix(y)) y <- as.matrix(y)

  N <- nrow(y)
  k.max <- trunc(12 * (N / 100)^(1 / 4))

  model <- as.numeric(const) + 1
  trim.pos <- which(c(0.01, 0.05, 0.1, 0.15, 0.25) == trim)

  res <- 0

  for (l in 0:(breaks - 1)) {
    test.stat <- KP.seq.statistic(y, const, l, criterion, trim, k.max)
    c.v <- .cval_KP[[model]][[trim.pos]]

    if (test.stat < c.v[2, l + 1]) {
      res <- l
      break
    }
  }

  res
}

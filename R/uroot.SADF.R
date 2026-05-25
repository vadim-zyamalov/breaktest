#' @title
#' Supremum ADF tests
#' @order 1
#'
#' @details
#' This set of unit root tests are developed for testing for bubbles
#' under time-varying non-stationary volatility.
#' Because the limiting distribution of the seminal Phillips et al. (2011) test
#' depends on the variance function and usually requires a bootstrap implementation
#' under heteroskedasticity, Kurozumi et al. (2022) construct the test
#' based on a deformation of the time domain.
#' The proposed test is asymptotically pivotal under the null hypothesis
#' and its limiting distribution coincides with that of the standard test
#' under homoskedasticity, so that the test does not require computationally
#' extensive methods for inference.
#'
#' `uroot.SADF` is a test statistic equal to the minimum value of [uroot.ADF] for
#' subsamples starting at \eqn{t = 1}.
#'
#' `GSADF.test` is a generalized version of `SADF.test`. Subsamples are allowed
#' to start at any point between 1 and \eqn{T (1 - trim)}.
#'
#' @param y A time series of interest.
#' @param trim A trimming parameter to determine the lower and upper bounds for
#' a possible break point.
#' @param const Whether the constant needs to be included.
#' @param add.p.value Whether the p-value is to be returned. This argument is
#' needed to suppress the calculation of p-values during the precalculation of
#' tables needed for the p-values estimating.
#'
#' @return [uroot.SADF] returns an object of type `bt_SADF`. It's a list of:
#' * `y`,
#' * `trim`,
#' * `const`,
#' * vector of \eqn{t}-values,
#' * the value of the corresponding test statistic,
#' * \eqn{p}-value if it was asked for.
#'
#' @references
#' Kurozumi, Eiji, Anton Skrobotov, and Alexey Tsarev.
#' “Time-Transformed Test for Bubbles under Non-Stationary Volatility.”
#' Journal of Financial Econometrics, April 23, 2022.
#' https://doi.org/10.1093/jjfinec/nbac004.
#'
#' @export
uroot.SADF <- function(
  y,
  trim = 0.01 + 1.8 / sqrt(length(y)),
  const = TRUE,
  add.p.value = TRUE,
  boot.p = FALSE,
  boot.iter = 999
) {
  N <- length(y)

  if (!const) {
    y <- y - y[1]
  }

  t.values <- c()
  m <- 1
  for (j in (floor(trim * N)):N) {
    model <- uroot.ADF(y[1:j], const = const)
    t.values[m] <- model$t.alpha
    m <- m + 1
  }

  ## Take the maximum of the calculated t-statistics.
  SADF.value <- max(t.values)

  result <- list(
    y = y,
    trim = trim,
    const = const,
    t.values = t.values,
    SADF.value = SADF.value
  )
  class(result) <- "bt_SADF"

  if (add.p.value) {
    cr.values <- ifelse(const, .cval_SADF_with_const, .cval_SADF_without_const)
    result$p.value <- get.p.values.SADF(SADF.value, N, cr.values)
  }

  if (boot.p) {
    result$boot.p.value <- bootstrap(result, boot.iter)
  }

  result
}


#' @rdname uroot.SADF
#' @order 2
#'
#' @return [uroot.GSADF] returns an object of class `bt_SADF` and subclass `bt_GSADF`.
#'
#' @export
uroot.GSADF <- function(
  y,
  trim = 0.01 + 1.8 / sqrt(length(y)),
  const = TRUE,
  add.p.value = TRUE,
  boot.p = FALSE,
  boot.iter = 999
) {
  N <- length(y)

  if (const == FALSE) {
    y <- y - y[1]
  }

  t.values <- c()
  m <- 1
  for (i in 1:(N - floor(trim * N) + 1)) {
    for (j in (i + floor(trim * N) - 1):N) {
      model <- uroot.ADF(y[i:j], const = const)
      t.values[m] <- model$t.alpha
      m <- m + 1
    }
  }

  ## Take the maximum of the calculated t-statistics.
  GSADF.value <- max(t.values)

  result <- list(
    y = y,
    trim = trim,
    const = const,
    t.values = t.values,
    GSADF.value = GSADF.value
  )
  class(result) <- c("bt_GSADF", "bt_SADF")

  if (add.p.value) {
    cr.values <- ifelse(
      const,
      .cval_GSADF_with_const,
      .cval_GSADF_without_const
    )
    result$p.value <- get.p.values.SADF(GSADF.value, N, cr.values)
  }

  if (boot.p) {
    result$boot.p.value <- bootstrap(result, boot.iter)
  }

  result
}

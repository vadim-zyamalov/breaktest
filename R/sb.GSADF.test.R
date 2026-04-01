#' @title
#' Sign-based SADF test
#'
#' @details
#' Refactored original code by Kurozumi et al.
#'
#' @param y A time series of interest.
#' @param trim A trimming parameter to determine the lower and upper bounds for
#' a possible break point.
#' @param const Whether the constant needs to be included.
#' @param alpha Needed level of significance.
#' @param iter Number of bootstrapping iterations.
#' @param urs Use union of rejections strategy if `TRUE`.
#' @param seed The seed parameter for the random number generator.
#'
#' @references
#' Harvey, David I., Stephen J. Leybourne, and Yang Zu.
#' “Sign-Based Unit Root Tests for Explosive Financial Bubbles
#' in the Presence of Deterministically Time-Varying Volatility.”
#' Econometric Theory 36, no. 1 (February 2020): 122–69.
#' https://doi.org/10.1017/S0266466619000057.
#'
#' Kurozumi, Eiji, Anton Skrobotov, and Alexey Tsarev.
#' “Time-Transformed Test for Bubbles under Non-Stationary Volatility.”
#' Journal of Financial Econometrics, April 23, 2022.
#' https://doi.org/10.1093/jjfinec/nbac004.
#'
#' @import doSNOW
#' @import foreach
#' @import parallel
#' @importFrom stats quantile
#' @importFrom stats rnorm
#' @importFrom stats sd
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#'
#' @export
sb.GSADF.test <- function(y,
                          trim = 0.01 + 1.8 / sqrt(length(y)),
                          const = TRUE,
                          alpha = 0.05,
                          iter = 999,
                          urs = TRUE,
                          seed = round(10^4 * sd(y))) {
  N <- length(y)

  ## Find supSBADF_value.
  supSBADF.model <- supSBADF.statistic(y, trim)

  ## Do parallel.
  cores <- detectCores()

  progress.bar <- txtProgressBar(max = iter, style = 3)
  progress <- function(n) setTxtProgressBar(progress.bar, n)

  cluster <- makeCluster(max(cores - 1, 1))
  clusterExport(cluster, c(
    "ADF.test",
    "GSADF.test",
    "supSBADF.statistic",
    ".cval_GSADF_without_const",
    ".cval_GSADF_with_const",
    ".diffn"
  ))
  registerDoSNOW(cluster)

  GSADF.supSBADF.bootstrap.values <- foreach( # nolint
    step = 1:iter,
    .combine = rbind,
    .options.snow = list(progress = progress)
  ) %dopar% {
    y.star <- cumsum(rnorm(N - 1) * .diffn(y))
    tmp.GSADF.value <- NA
    supSBADF.value <- NA
    if (urs) {
      gsadf.model <- GSADF.test(y.star, trim, const)
      tmp.GSADF.value <- gsadf.model$sadf.value
    }
    supSBADF.model <- supSBADF.statistic(y.star, trim)
    tmp.supSBADF.value <- supSBADF.model$supSBADF.value
    c(tmp.GSADF.value, tmp.supSBADF.value)
  }

  stopCluster(cluster)

  ## Get sadf_supSBADF_bootstrap_values
  supSBADF.bootstrap.values <- GSADF.supSBADF.bootstrap.values[, 2]

  ## Find critical value.
  supSBADF.cr.value <- as.numeric(quantile(
    supSBADF.bootstrap.values,
    1 - alpha
  ))

  ## A union of rejections strategy.
  if (urs) {
    ## Find sadf_value.
    gsadf.model <- GSADF.test(y, trim, const)
    t.values <- gsadf.model$t.values
    GSADF.value <- gsadf.model$GSADF.value

    ## Get sadf_supSBADF_bootstrap_values
    GSADF.bootstrap.values <- GSADF.supSBADF.bootstrap.values[, 1]

    ## Find critical value.
    GSADF.cr.value <- as.numeric(quantile(
      GSADF.bootstrap.values,
      1 - alpha
    ))

    ## Calculate U value.
    U.value <- max(
      GSADF.value,
      GSADF.cr.value / supSBADF.cr.value * supSBADF.value
    )

    ## Find U_bootstrap_values.
    U.bootstrap.values <- c()
    for (b in 1:iter) {
      U.bootstrap.values[b] <- max(
        GSADF.bootstrap.values[b],
        GSADF.cr.value /
          supSBADF.cr.value * supSBADF.bootstrap.values[b]
      )
    }

    ## Find critical value.
    U.cr.value <- as.numeric(quantile(U.bootstrap.values, 1 - alpha))

    p.value <- round(sum(U.bootstrap.values > U.value) / iter, 4)

    is.explosive <- ifelse(U.value > U.cr.value, 1, 0)
  } else {
    p.value <- round(sum(supSBADF.bootstrap.values > supSBADF.value) / iter, 4)

    is.explosive <- ifelse(supSBADF.value > supSBADF.cr.value, 1, 0)
  }

  result <- c(
    list(
      y = y,
      trim = trim,
      const = const,
      alpha = alpha,
      iter = iter,
      urs = urs,
      seed = seed,
      SBADF.values = supSBADF.model$SBADF.values,
      supSBADF.value = supSBADF.model$supSBADF_value,
      supSBADF.bootstrap.values = supSBADF.bootstrap.values,
      supSBADF.cr.value = supSBADF.cr.value,
      p.value = p.value,
      is.explosive = is.explosive
    ),
    if (urs) {
      list(
        t.values = t.values,
        GSADF.value = GSADF.value,
        GSADF.bootstrap.values = GSADF.bootstrap.values,
        GSADF.cr.value = GSADF.cr.value,
        U.value = U.value,
        U.bootstrap.values = U.bootstrap.values,
        U.cr.value = U.cr.value
      )
    } else {
      NULL
    }
  )

  class(result) <- "sadf"

  result
}


#' @title
#' Calculate superior sign-based SADF statistic.
#'
#' @param y The series of interest.
#' @param trim Trimming parameter to determine the lower and upper bounds.
#' @param generalized Whether to calculate generalized statistic value.
#'
#' @return A list of
#' * `y`,
#' * `trim`,
#' * `C.t`: the cumulative sum of "signs" (1 or -1) of the first difference of
#' `y`,
#' * `SBADF.values`: series of sign-based ADF statistics,
#' * `supSBADF.value`: the maximum of `SBADF.values`.
#'
#' @references
#' Harvey, David I., Stephen J. Leybourne, and Yang Zu.
#' “Sign-Based Unit Root Tests for Explosive Financial Bubbles in the Presence
#' of Deterministically Time-Varying Volatility.”
#' Econometric Theory 36, no. 1 (February 2020): 122–69.
#' https://doi.org/10.1017/S0266466619000057.
#'
#' @keywords internal
supSBADF.statistic <- function(y,
                               trim = 0.01 + 1.8 / sqrt(length(y)),
                               generalized = FALSE) {
  N <- length(y)

  ## Calculate C.t.
  C.t <- cumsum(sign(diff(y)))

  SBADF.values <- c()
  m <- 1

  if (!generalized) {
    for (j in (floor(trim * N)):N) {
      t.beta <- .OLS(.diffn(C.t, na = 0)[1:j], C.t[1:j])$t.beta
      SBADF.values[m] <- drop(t.beta)
      m <- m + 1
    }
  } else {
    for (i in 1:(N - floor(trim * N) + 1)) {
      for (j in (i + floor(trim * N) - 1):N) {
        t.beta <- .OLS(.diffn(C.t, na = 0)[i:j], C.t[i:j])$t.beta
        SBADF.values[m] <- drop(t.beta)
        m <- m + 1
      }
    }
  }

  supSBADF.value <- max(SBADF.values)

  list(
    y = y,
    trim = trim,
    C.t = C.t,
    SBADF.values = SBADF.values,
    supSBADF.value = supSBADF.value
  )
}

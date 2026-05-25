#' @rdname uroot.SADF
#' @order 5
#'
#' @details
#' Refactored original code by Kurozumi et al.
#'
#' @param alpha Needed level of significance.
#' @param iter Number of bootstrapping iterations.
#' @param urs Use union of rejections strategy if `TRUE`.
#' @param seed The seed parameter for the random number generator.
#'
#' @return [uroot.sb.GSADF] returns an object of class `bt_SADF` and subclass `bt_sbGSADF`. It's a list of:
#' * `y`,
#' * main parameters, such as `trim`, `const`, `alpha`, `iter`,
#  * `urs`,
#' * `seed`,
#' supremum SBADF-statistic values incl. bootstrapped ones,
#' * p-value,
#' * indicator of explosive process.
#'
#' If `urs=TRUE` then some additional values included:
#' * t-values,
#' * GSADF bootstrapped statistics,
#' * bootstrapped U-values.
#'
#' @references
#' Harvey, David I., Stephen J. Leybourne, and Yang Zu.
#' “Sign-Based Unit Root Tests for Explosive Financial Bubbles
#' in the Presence of Deterministically Time-Varying Volatility.”
#' Econometric Theory 36, no. 1 (February 2020): 122–69.
#' https://doi.org/10.1017/S0266466619000057.
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
uroot.sb.GSADF <- function(
  y,
  trim = 0.01 + 1.8 / sqrt(length(y)),
  const = TRUE,
  alpha = 0.05,
  iter = 999,
  urs = TRUE,
  seed = round(10^4 * sd(y))
) {
  N <- length(y)

  ## Find supSBADF_value.
  supSBADF.model <- supSBADF.statistic(y, trim)
  supSBADF.value <- supSBADF.model$supSBADF.value

  ## Do parallel.
  cores <- detectCores()

  progress.bar <- txtProgressBar(max = iter, style = 3)
  progress <- function(n) setTxtProgressBar(progress.bar, n)

  cluster <- makeCluster(max(cores - 1, 1))
  clusterExport(
    cluster,
    c(
      "ADF.test",
      "GSADF.test",
      # "supSBADF.statistic",
      ".cval_GSADF_without_const",
      ".cval_GSADF_with_const",
      ".diffn"
    )
  )
  registerDoSNOW(cluster)

  GSADF.bootstrap.values <- foreach(
    i = 1:iter,
    .combine = rbind,
    .options.snow = list(progress = progress)
  ) %dopar%
    {
      y.star <- cumsum(rnorm(N) * .diffn(y, na = 0))
      res <- c(NA, NA)
      if (urs) {
        gsadf.model <-
          res[1] <- uroot.GSADF(y.star, trim, const)$GSADF.value
      }
      res[2] <- supSBADF.statistic(y.star, trim)$statistic
      res
    }

  stopCluster(cluster)

  ## Get sadf_supSBADF_bootstrap_values
  supSBADF.bootstrap.values <- GSADF.bootstrap.values[, 2]

  ## Find critical value.
  supSBADF.cr.value <- as.numeric(quantile(
    supSBADF.bootstrap.values,
    1 - alpha
  ))

  ## A union of rejections strategy.
  if (urs) {
    ## Find sadf_value.
    gsadf.model <- uroot.GSADF(y, trim, const)
    t.values <- gsadf.model$t.values
    GSADF.value <- gsadf.model$GSADF.value

    ## Get sadf_supSBADF_bootstrap_values
    GSADF.bootstrap.values <- GSADF.bootstrap.values[, 1]

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
          supSBADF.cr.value *
          supSBADF.bootstrap.values[b]
      )
    }

    ## Find critical value.
    U.cr.value <- as.numeric(quantile(U.bootstrap.values, 1 - alpha))

    p.value <- sum(U.bootstrap.values > U.value) / iter

    is.explosive <- ifelse(U.value > U.cr.value, 1, 0)
  } else {
    p.value <- sum(supSBADF.bootstrap.values > supSBADF.value) / iter

    is.explosive <- ifelse(supSBADF.value > supSBADF.cr.value, 1, 0)
  }

  result <- list(
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
  )

  class(result) <- c("bt_sbGSADF", "bt_SADF")

  if (urs) {
    result <- c(
      result,
      list(
        t.values = t.values,
        GSADF.value = GSADF.value,
        GSADF.bootstrap.values = GSADF.bootstrap.values,
        GSADF.cr.value = GSADF.cr.value,
        U.value = U.value,
        U.bootstrap.values = U.bootstrap.values,
        U.cr.value = U.cr.value
      )
    )
  }

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
supSBADF.statistic <- function(
  y,
  trim = 0.01 + 1.8 / sqrt(length(y)),
  generalized = FALSE
) {
  N <- length(y)

  ## Calculate C.t.
  C.t <- cumsum(sign(.diffn(y, na = 0)))

  SBADF.values <- c()
  m <- 1

  if (!generalized) {
    for (j in (floor(trim * N)):N) {
      t.beta <- OLS.reg(.diffn(C.t, na = 0)[1:j], C.t[1:j])$t.stats
      SBADF.values[m] <- drop(t.beta)
      m <- m + 1
    }
  } else {
    for (i in 1:(N - floor(trim * N) + 1)) {
      for (j in (i + floor(trim * N) - 1):N) {
        t.beta <- OLS.reg(.diffn(C.t, na = 0)[i:j], C.t[i:j])$t.stats
        SBADF.values[m] <- drop(t.beta)
        m <- m + 1
      }
    }
  }

  list(
    SBADF.values = SBADF.values,
    statistic = max(SBADF.values)
  )
}

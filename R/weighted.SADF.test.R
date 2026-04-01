#' @title
#' Weighted supremum ADF test
#' @order 1
#'
#' @details
#' Refactored original code by Kurozumi et al.
#'
#' @param y A time series of interest.
#' @param trim The trimming parameter to find the lower and upper bounds of
#' possible break dates.
#' @param const Whether the constant needs to be included.
#' @param alpha A significance level of interest.
#' @param iter Nnumber of iterations.
#' @param urs Use `union of rejections` strategy.
#' @param seed A seed parameter for the random number generator.
#'
#' @return An object of type `sadf`. It's a list of:
#' * `y`,
#' * `trim`,
#' * `const`,
#' * `alpha`,
#' * `iter`,
#' * `urs`,
#' * `seed`,
#' * `sigma.sq`: the estimated variance,
#' * `BZ.values`: a series of BZ-statistic,
#' * `supBZ.value`: the maximum of `supBZ.values`,
#' * `supBZ.bootstsrap.values`: bootstrapped supremum BZ values,
#' * `supBZ.cr.value`: supremum BZ \eqn{\alpha} critical value,
#' * `p.value`,
#' * `is.explosive`: 1 if `supBZ.value` is greater than `supBZ.cr.value`.
#'
#' if `urs` is `TRUE` the following items are also included:
#' * vector of \eqn{t}-values,
#' * the value of the SADF test statistic,
#' * `SADF.bootstrap.values`: bootstrapped SADF values,
#' * `U.value`: union test statistic value,
#' * `U.bootstrap.values`: bootstrapped series of `U.value`,
#' * `U.cr.value`: critical value of `U.value`.
#'
#' @references
#' Harvey, David I., Stephen J. Leybourne, and Yang Zu.
#' “Testing Explosive Bubbles with Time-Varying Volatility.”
#' Econometric Reviews 38, no. 10 (November 26, 2019): 1131–51.
#' https://doi.org/10.1080/07474938.2018.1536099.
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
weighted.SADF.test <- function(y,
                               trim = 0.01 + 1.8 / sqrt(length(y)),
                               const = TRUE,
                               alpha = 0.05,
                               iter = 4 * 200,
                               urs = TRUE,
                               seed = round(10^4 * sd(y))) {
  N <- length(y)

  ## Find supBZ.value.
  supBZ.model <- supBZ.statistic(y, trim)
  sigma.sq <- supBZ.model$sigma.sq
  BZ.values <- supBZ.model$BZ.values
  supBZ.value <- supBZ.model$supBZ.value

  ## Do parallel.
  cores <- detectCores()

  progress.bar <- txtProgressBar(max = iter, style = 3)
  progress <- function(n) setTxtProgressBar(progress.bar, n)

  cluster <- makeCluster(max(cores - 1, 1))
  clusterExport(cluster, c(
    "ADF.test",
    "SADF.test",
    "supBZ.statistic",
    ".cval_SADF_without_const",
    ".cval_SADF_with_const",
    ".diffn"
  ))
  registerDoSNOW(cluster)

  SADF.supBZ.bootstrap.values <- foreach(
    i = 1:iter,
    .combine = rbind,
    .options.snow = list(progress = progress)
  ) %dopar% {
    y.star <- cumsum(rnorm(N - 1) * .diff(y, na = 0))
    tmp.SADF.value <- NA
    if (urs) {
      tmp.sadf.model <- SADF.test(y.star, trim, const)
      tmp.SADF.value <- tmp.sadf.model$SADF.value
    }
    tmp.supBZ.model <- supBZ.statistic(y.star, trim, sigma.sq)
    tmp.supBZ.value <- tmp.supBZ.model$supBZ.value
    c(tmp.SADF.value, tmp.supBZ.value)
  }

  stopCluster(cluster)

  ## Get sadf_supBZ_bootstrap.values.
  supBZ.bootstrap.values <- SADF.supBZ.bootstrap.values[, 2]

  ## Find critical value.
  supBZ.cr.value <- as.numeric(quantile(
    supBZ.bootstrap.values,
    1 - alpha
  ))

  ## A union of rejections strategy.
  if (urs == TRUE) {
    ## Find SADF.value.
    sadf.model <- SADF.test(y, trim, const)
    t.values <- sadf.model$t.values
    SADF.value <- sadf.model$SADF.value

    ## Get sadf_supBZ.bootstrap.values.
    SADF.bootstrap.values <- SADF.supBZ.bootstrap.values[, 1]

    ## Find critical value.
    SADF.cr.value <- as.numeric(quantile(SADF.bootstrap.values, 1 - alpha))

    ## Calculate U value.
    U.value <- max(
      SADF.value,
      SADF.cr.value / supBZ.cr.value * supBZ.value
    )

    ## Find U_bootstrap.values.
    U.bootstrap.values <- c()
    for (b in 1:iter) {
      U.bootstrap.values[b] <- max(
        SADF.bootstrap.values[b],
        SADF.cr.value /
          supBZ.cr.value *
          supBZ.bootstrap.values[b]
      )
    }

    ## Find critical value.
    U.cr.value <- as.numeric(quantile(U.bootstrap.values, 1 - alpha))

    p.value <- round(sum(U.bootstrap.values > U.value) / iter, 4)

    is.explosive <- ifelse(U.value > U.cr.value, 1, 0)
  } else {
    p.value <- round(sum(supBZ.bootstrap.values > supBZ.value) / iter, 4)

    is.explosive <- ifelse(supBZ.value > supBZ.cr.value, 1, 0)
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
      sigma.sq = sigma.sq,
      BZ.values = BZ.values,
      supBZ.value = supBZ.value,
      supBZ.bootstrap.values = supBZ.bootstrap.values,
      supBZ.cr.value = supBZ.cr.value,
      p.value = p.value,
      is.explosive = is.explosive
    ),
    if (urs) {
      list(
        t.values = t.values,
        SADF.value = SADF.value,
        SADF.bootstrap.values = SADF.bootstrap.values,
        SADF.cr.value = SADF.cr.value,
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


#' @rdname weighted.SADF.test
#' @order 2
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
weighted.GSADF.test <- function(y,
                                trim = 0.01 + 1.8 / sqrt(length(y)),
                                const = TRUE,
                                alpha = 0.05,
                                iter = 4 * 200,
                                urs = TRUE,
                                seed = round(10^4 * sd(y))) {
  N <- length(y)

  ## Find supBZ.value.
  supBZ.model <- supBZ.statistic(y, trim)
  sigma.sq <- supBZ.model$sigma.sq
  BZ.values <- supBZ.model$BZ.values
  supBZ.value <- supBZ.model$supBZ.value

  ## Do parallel.
  cores <- detectCores()

  progress.bar <- txtProgressBar(max = iter, style = 3)
  progress <- function(n) setTxtProgressBar(progress.bar, n)

  cluster <- makeCluster(max(cores - 1, 1))
  clusterExport(cluster, c(
    "ADF.test",
    "GSADF.test",
    "supBZ.statistic",
    ".cval_GSADF_without_const",
    ".cval_GSADF_with_const",
    ".diffn"
  ))
  registerDoSNOW(cluster)

  SADF.supBZ.bootstrap.values <- foreach(
    step = 1:iter,
    .combine = rbind,
    .options.snow = list(progress = progress)
  ) %dopar% {
    y.star <- cumsum(rnorm(N - 1) * .diffn(y, na = 0))
    tmp.GSADF.value <- NA
    if (urs) {
      gsadf.model <- GSADF.test(y.star, trim, const)
      tmp.GSADF.value <- gsadf.model$GSADF.value
    }
    supBZ.model <- supBZ.statistic(y.star, trim, sigma.sq)
    tmp.supBZ.value <- supBZ.model$supBZ.value
    c(tmp.GSADF.value, tmp.supBZ.value)
  }

  stopCluster(cluster)

  ## Get sadf_supBZ.bootstsrap.values.
  supBZ.bootstsrap.values <- SADF.supBZ.bootstrap.values[, 2]

  ## Find critical value.
  supBZ.cr.value <- as.numeric(quantile(
    supBZ.bootstsrap.values,
    1 - alpha
  ))

  ## A union of rejections strategy.
  if (urs) {
    ## Find sadf.value.
    gsadf.model <- GSADF.test(y, trim, const)
    t.values <- gsadf.model$t.values
    GSADF.value <- gsadf.model$GSADF.value

    ## Get sadf_supBZ.bootstsrap.values.
    GSADF.bootstsrap.values <- SADF.supBZ.bootstrap.values[, 1]

    ## Find critical value.
    GSADF.cr.value <- as.numeric(quantile(
      GSADF.bootstsrap.values,
      1 - alpha
    ))

    ## Calculate U value.
    U.value <- max(
      GSADF.value,
      GSADF.cr.value / supBZ.cr.value * supBZ.value
    )

    ## Find U.bootstsrap.values.
    U.bootstsrap.values <- c()
    for (b in 1:iter) {
      U.bootstsrap.values[b] <- max(
        GSADF.bootstsrap.values[b],
        GSADF.cr.value / supBZ.cr.value * supBZ.bootstsrap.values[b]
      )
    }

    ## Find critical value.
    U.cr.value <- as.numeric(quantile(U.bootstsrap.values, 1 - alpha))

    p.value <- round(sum(U.bootstsrap.values > U.value) / iter, 4)

    is.explosive <- ifelse(U.value > U.cr.value, 1, 0)
  } else {
    p.value <- round(sum(supBZ.bootstsrap.values > supBZ.value) / iter, 4)

    is.explosive <- ifelse(supBZ.value > supBZ.cr.value, 1, 0)
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
      sigma.sq = sigma.sq,
      BZ.values = BZ.values,
      supBZ.value = supBZ.value,
      supBZ.bootstsrap.values = supBZ.bootstsrap.values,
      supBZ.cr.value = supBZ.cr.value,
      p.value = p.value,
      is.explosive = is.explosive
    ),
    if (urs) {
      list(
        t.values = t.values,
        GSADF.value = GSADF.value,
        GSADF.bootstsrap.values = GSADF.bootstsrap.values,
        GSADF.cr.value = GSADF.cr.value,
        U.value = U.value,
        U.bootstsrap.values = U.bootstsrap.values,
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
#' Calculate supBZ statistic
#'
#' @param y A time series of interest.
#' @param trim The trimming parameter to find the lower and upper bounds of
#' possible break dates.
#' @param sigma.sq Local non-parametric estimates of variance. If `NULL` they
#' will be estimated via Nadaraya-Watson procedure.
#' @param generalized Whether to calculate generalized statistic value.
#'
#' @return A list of:
#' * `y`,
#' * `trim`,
#' * `sigma.sq`,
#' * `BZ.values`: a series of BZ-statistic,
#' * `supBZ.value`: the maximum of `supBZ.values`,
#' * `h.est`: the estimated value of bandwidth if `sigma.sq` is `NULL`.
#'
#' @references
#' Harvey, David I., Stephen J. Leybourne, and Yang Zu.
#' “Testing Explosive Bubbles with Time-Varying Volatility.”
#' Econometric Reviews 38, no. 10 (November 26, 2019): 1131–51.
#' https://doi.org/10.1080/07474938.2018.1536099.
#'
#' @keywords internal
supBZ.statistic <- function(y,
                            trim = 0.01 + 1.8 / sqrt(length(y)),
                            sigma.sq = NULL,
                            generalized = FALSE) {
  n.obs <- length(y)

  if (is.null(sigma.sq)) {
    ## NW estimation.
    my <- (diff(y))^2
    mx <- rep(1, n.obs - 1)
    nw.loocv.model <- .NW.bandwidth(my, mx, kernel = "gauss")
    h.est <- nw.loocv.model$h
    nw.model <- .NW.variance(
      my,
      kernel = "gauss",
      h = nw.loocv.model$h
    )
    sigma.sq <- nw.model$omega.sq
  }

  y <- y - y[1]
  d.y <- diff(y)
  l.y <- y[1:(n.obs - 1)]

  BZ.values <- c()
  m <- 1

  if (!generalized) {
    for (j in (floor(trim * n.obs)):n.obs) {
      BZ.values[m] <-
        sum(d.y[1:(j - 1)] * l.y[1:(j - 1)] / sigma.sq[1:(j - 1)]) /
          (sum(l.y[1:(j - 1)]^2 / sigma.sq[1:(j - 1)]))^0.5
      m <- m + 1
    }
  } else {
    for (i in 1:(n.obs - floor(trim * n.obs) + 1)) {
      for (j in (i + floor(trim * n.obs) - 1):n.obs) {
        BZ.values[m] <-
          sum(d.y[i:(j - 1)] * l.y[i:(j - 1)] / sigma.sq[i:(j - 1)]) /
            (sum(l.y[i:(j - 1)]^2 / sigma.sq[i:(j - 1)]))^0.5
        m <- m + 1
      }
    }
  }

  supBZ.value <- max(BZ.values)

  c(
    list(
      y = y,
      trim = trim,
      sigma.sq = sigma.sq,
      BZ.values = BZ.values,
      supBZ.value = supBZ.value
    ),
    if (exists("h.est")) {
      list(h.est = h.est)
    } else {
      NULL
    }
  )
}

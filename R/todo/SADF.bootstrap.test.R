#' @title
#' Supremum ADF tests with wild bootstrap.
#' @order 1
#'
#' @description
#' `SADF.bootstrap.test` is a wild bootstrapping procedure for estimating
#' critical and \eqn{p}-values for [SADF.test].
#'
#' `GSADF.bootstrap.test` is the same procedure but for [GSADF.test].
#'
#' @details
#' Refactored original code by Kurozumi et al.
#'
#' @param y A time series of interest.
#' @param trim A trimming parameter to determine the lower and upper bounds for
#' a possible break point.
#' @param const Whether the constant needs to be included.
#' @param alpha The significance level of interest.
#' @param iter The number of iterations.
#' @param seed The seed parameter for the random number generator.
#'
#' @return An object of type `sadf`. It's a list of:
#' * `y`,
#' * `trim`,
#' * `const`,
#' * `alpha`,
#' * `iter`,
#' * `seed`,
#' * vector of \eqn{t}-values,
#' * the value of the corresponding test statistic,
#' * series of bootstrapped test statistics,
#' * bootstrapped critical values,
#' * \eqn{p}-value.
#'
#' @references
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
SADF.bootstrap.test <- function(y,
                                trim = 0.01 + 1.8 / sqrt(length(y)),
                                const = TRUE,
                                alpha = 0.05,
                                iter = 999,
                                seed = round(10^4 * sd(y))) {
    n.obs <- length(y)

    ## Find SADF.value.
    model <- SADF.test(y, trim, const)
    t.values <- model$t.values
    SADF.value <- model$SADF.value

    ## Do parallel.
    cores <- detectCores()

    progress.bar <- txtProgressBar(max = iter, style = 3)
    progress <- function(n) setTxtProgressBar(progress.bar, n)

    cluster <- makeCluster(max(cores - 1, 1))
    clusterExport(cluster, c("SADF.test"))
    registerDoSNOW(cluster)

    SADF.bootstrap.values <- foreach(
        step = 1:iter,
        .combine = c,
        .options.snow = list(progress = progress)
    ) %dopar% {
        y.star <- cumsum(c(0, rnorm(n.obs - 1) * diff(y)))
        model <- SADF.test(y.star, trim, const)
        model$SADF.value
    }

    stopCluster(cluster)

    ## Find critical value.
    cr.value <- as.numeric(quantile(SADF.bootstrap.values, 1 - alpha))

    p.value <- round(sum(SADF.bootstrap.values > SADF.value) / iter, 4)

    result <- list(
        y = y,
        trim = trim,
        const = const,
        alpha = alpha,
        iter = iter,
        seed = seed,
        t.values = t.values,
        SADF.value = SADF.value,
        SADF.bootstrap.values = SADF.bootstrap.values,
        cr.value = cr.value,
        p.value = p.value
    )

    class(result) <- "sadf"

    return(result)
}

#' @rdname SADF.bootstrap.test
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
GSADF.bootstrap.test <- function(y,
                                 trim = 0.01 + 1.8 / sqrt(length(y)),
                                 const = TRUE,
                                 alpha = 0.05,
                                 iter = 4 * 200,
                                 seed = round(10^4 * sd(y))) {
    n.obs <- length(y)

    ## Find SADF.value.
    model <- GSADF.test(y, trim, const)
    t.values <- model$t.values
    GSADF.value <- model$GSADF.value

    ## Do parallel.
    cores <- detectCores()

    progress.bar <- txtProgressBar(max = iter, style = 3)
    progress <- function(n) setTxtProgressBar(progress.bar, n)

    cluster <- makeCluster(max(cores - 1, 1))
    clusterExport(cluster, c("GSADF.test"))
    registerDoSNOW(cluster)

    GSADF.bootstsrap.values <- foreach(
        step = 1:iter,
        .combine = c,
        .options.snow = list(progress = progress)
    ) %dopar% {
        y.star <- cumsum(c(0, rnorm(n.obs - 1) * diff(y)))
        model <- GSADF.test(y.star, trim, const)
        model$GSADF.value
    }
    stopCluster(cluster)

    ## Find critical value.
    cr.value <- as.numeric(quantile(GSADF.bootstsrap.values, 1 - alpha))

    p.value <- round(sum(GSADF.bootstsrap.values > GSADF.value) / iter, 4)

    result <- list(
        y = y,
        trim = trim,
        const = const,
        alpha = alpha,
        iter = iter,
        seed = seed,
        t.values = t.values,
        GSADF.value = GSADF.value,
        GSADF.bootstsrap.values = GSADF.bootstsrap.values,
        cr.value = cr.value,
        p.value = p.value
    )

    class(result) <- "sadf"

    return(result)
}

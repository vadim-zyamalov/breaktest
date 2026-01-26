#' @title
#' Andrews-Kim (2006) test
#'
#' @description
#' Test for structural break at the end of the sample.
#'
#' @details
#' See Andrews and Kim (2006) for the detailed description.
#'
#' @param eq Base model formula. At the moment all the variables included should
#' be defined explicitly, dynamic regressors (i.e. functions etc.) are not
#' supported.
#' @param m Post-break period length.
#' @param dataset Source of the data.
#'
#' @return A list of
#' * `m`,
#' * estimated values of P- and R-tests,
#' * sequences of auxiliary statistics \eqn{P_j} and \eqn{R_j},
#' * the corresponding p-values.
#'
#' @references
#' Andrews, D. W. K.
#' “End-of-Sample Instability Tests.”
#' Econometrica 71, no. 6 (2003): 1661–94.
#' https://doi.org/10.1111/1468-0262.00466.
#'
#' Andrews, Donald W. K, and Jae-Young Kim.
#' “Tests for Cointegration Breakdown Over a Short Time Period.”
#' Journal of Business & Economic Statistics 24, no. 4 (2006): 379–94.
#' https://doi.org/10.1198/073500106000000297.
#'
#' @importFrom stats lm
#' @importFrom stats predict
#'
#' @export
eos.break.test <- function(eq,
                           m,
                           dataset) {
    result <- list(m = m)

    dep.var <- eq[[2]]

    n.obs <- dim(dataset)[1] - m

    A <- matrix(NA, m, m)
    for (i in 1:m) for (j in 1:m) A[i, j] <- min(i, j)

    tmp.model <- lm(formula = eq, data = dataset)
    tmp.resid <- as.matrix(tmp.model$residuals[(n.obs + 1):(n.obs + m)])

    result$P <- drop(t(tmp.resid) %*% tmp.resid)
    result$R <- drop(t(tmp.resid) %*% A %*% tmp.resid)

    result$Pj <- seq(0, 0, length.out = n.obs - m + 1)
    result$Rj <- seq(0, 0, length.out = n.obs - m + 1)

    for (j in seq_len(n.obs - m + 1)) {
        low <- j
        high <- j + ceiling(m / 2) - 1
        tmp.data <-
            dataset[1:n.obs, , drop = FALSE][-(low:high), , drop = FALSE]

        tmp.model <- lm(formula = eq, data = tmp.data)

        tmp.resid <- dataset[[dep.var]][j:(j + m - 1)] -
            predict(tmp.model, newdata = dataset[j:(j + m - 1), , drop = FALSE])
        tmp.resid <- as.matrix(tmp.resid)

        result$Pj[j] <- drop(t(tmp.resid) %*% tmp.resid)
        result$Rj[j] <- drop(t(tmp.resid) %*% A %*% tmp.resid)
    }

    result$p.value <- (1 / (n.obs - m + 1)) * sum(I(result$P <= result$Pj))
    result$r.value <- (1 / (n.obs - m + 1)) * sum(I(result$R <= result$Rj))

    return(result)
}

cat('\014')

library(zeallot)
library(parallel)
library(doSNOW)
library(breaktest)

N_sim <- 999
N_obs <- c(50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 1000)

progress.bar <- txtProgressBar(max = N_sim, style = 3)
progress <- function(n) setTxtProgressBar(progress.bar, n)

cores <- detectCores()
cluster <- makeCluster(max(cores - 1, 1), type = "SOCK")
registerDoSNOW(cluster)

data <- list()
data$level <- list()
data$`level-trend` <- list()
data$regime <- list()
for (i in 1:6) {
    data$level[[i]] <- list(Za = NULL, Zt = NULL, ADF = NULL)
    data$`level-trend`[[i]] <- list(Za = NULL, Zt = NULL, ADF = NULL)
    data$regime[[i]] <- list(Za = NULL, Zt = NULL, ADF = NULL)
}

for (m in c("level", "level-trend", "regime")) {
    for (N in N_obs) {
        cat("\nCalculating values for", N, "observations (", m, ")\n")

        tmp <- foreach(
            i = 1:N_sim,
            .combine = rbind,
            .inorder = FALSE,
            .errorhandling = "remove",
            .packages = c("breaktest"),
            .options.snow = list(progress = progress)
        ) %dopar% {
            y <- cumsum(rnorm(N))

            for (m in 1:6) {
                y <- cbind(y, cumsum(rnorm(N)))
            }

            r1 <- coint.test.GH(
                y[, 1], y[, 2],
                shift = "level"
            )
            r2 <- coint.test.GH(
                y[, 1], y[, 2], y[, 3],
                shift = "level"
            )
            r3 <- coint.test.GH(
                y[, 1], y[, 2], y[, 3], y[, 4],
                shift = "level"
            )
            r4 <- coint.test.GH(
                y[, 1], y[, 2], y[, 3], y[, 4], y[, 5],
                shift = "level"
            )
            r5 <- coint.test.GH(
                y[, 1], y[, 2], y[, 3], y[, 4], y[, 5], y[, 6],
                shift = "level"
            )
            r6 <- coint.test.GH(
                y[, 1], y[, 2], y[, 3], y[, 4], y[, 5], y[, 6], y[, 7],
                shift = "level"
            )

            c(
                r1$Za, r2$Za, r3$Za, r4$Za, r5$Za, r6$Za,
                r1$Zt, r2$Zt, r3$Zt, r4$Zt, r5$Zt, r6$Zt,
                r1$ADF, r2$ADF, r3$ADF, r4$ADF, r5$ADF, r6$ADF
            )
        }

        for (i in 1:6) {
            data[[m]][[i]]$Za <- rbind(
                data[[m]][[i]]$Za,
                c(
                    quantile(tmp[, i], probs = c(0.01, 0.025, 0.05, 0.1, 0.975)),
                    1,
                    1 / N,
                    1 / N^2
                )
            )
            data[[m]][[i]]$Zt <- rbind(
                data[[m]][[i]]$Zt,
                c(
                    quantile(tmp[, 6 + i], probs = c(0.01, 0.025, 0.05, 0.1, 0.975)),
                    1,
                    1 / N,
                    1 / N^2
                )
            )
            data[[m]][[i]]$ADF <- rbind(
                data[[m]][[i]]$ADF,
                c(
                    quantile(tmp[, 12 + i], probs = c(0.01, 0.025, 0.05, 0.1, 0.975)),
                    1,
                    1 / N,
                    1 / N^2
                )
            )
        }
    }
}
stopCluster(cluster)

for (m in c("level", "level-trend", "regime")) {
    .cval_coint_hg[[m]] <- list(
        Za = list(
            b0 = matrix(0, nrow = 6, ncol = 5),
            b1 = matrix(0, nrow = 6, ncol = 5)
        ),
        Zt = list(
            b0 = matrix(0, nrow = 6, ncol = 5),
            b1 = matrix(0, nrow = 6, ncol = 5)
        ),
        ADF = list(
            b0 = matrix(0, nrow = 6, ncol = 5),
            b1 = matrix(0, nrow = 6, ncol = 5)
        )
    )

    for (i in 1:6) {
        for (j in 1:5) {
            c(b, ., ., .) %<-% OLS(
                                   data[[m]][[i]]$Za[, j, drop = FALSE],
                                   data[[m]][[i]]$Za[, 6:7, drop = FALSE]
                               )
            .cval_coint_hg[[m]]$Za$b0[i, j] <- b[1]
            .cval_coint_hg[[m]]$Za$b1[i, j] <- b[2]

            c(b, ., ., .) %<-% OLS(
                                   data[[m]][[i]]$Zt[, j, drop = FALSE],
                                   data[[m]][[i]]$Zt[, 6:7, drop = FALSE]
                               )
            .cval_coint_hg[[m]]$Zt$b0[i, j] <- b[1]
            .cval_coint_hg[[m]]$Zt$b1[i, j] <- b[2]

            c(b, ., ., .) %<-% OLS(
                                   data[[m]][[i]]$ADF[, j, drop = FALSE],
                                   data[[m]][[i]]$ADF[, 6:7, drop = FALSE]
                               )
            .cval_coint_hg[[m]]$ADF$b0[i, j] <- b[1]
            .cval_coint_hg[[m]]$ADF$b1[i, j] <- b[2]
        }
    }
}
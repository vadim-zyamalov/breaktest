library(parallel)
library(doSNOW)
library(breaktest)

N_sim <- 999
N_obs <- c(
  seq(from = 30, to = 99, by = 10),
  seq(from = 100, to = 199, by = 20),
  seq(from = 200, to = 399, by = 50),
  seq(from = 400, to = 799, by = 100),
  seq(from = 800, to = 1599, by = 200),
  seq(from = 1600, to = 3200, by = 400)
)

## 30-100 шаг 10 наблюдений
## 100-200 шаг 20 наблюдений
## 200-400 шаг 50 наблюдений
## 400-800 шаг 100 наблюдений
## 800-1600 шаг 200 наблюдений
## 1600-3200 шаг 400 (или сразу 800 бахнуть) наблюдений

progress.bar <- txtProgressBar(max = N_sim, style = 3)
progress <- function(n) setTxtProgressBar(progress.bar, n)

cores <- detectCores()
cluster <- makeCluster(max(cores - 1, 1), type = "SOCK")
registerDoSNOW(cluster)

.cval_SADF_without_const <- list()
.cval_SADF_with_const <- list()
.cval_GSADF_without_const <- list()
.cval_GSADF_with_const <- list()

for (step in N_obs) {
  cat("\nCalculating tables for", step, "observations\n")

  cat("without CONST\n")
  tmp <- NULL
  tmp <- foreach(
    i = 1:N_sim,
    .combine = c,
    .inorder = FALSE,
    .errorhandling = "remove",
    .packages = c("breaktest"),
    .options.snow = list(progress = progress)
  ) %dopar% {
    repeat {
      y <- rnorm(step)
      y <- cumsum(y)

      model <- STADF.test(y, const = FALSE, add.p.value = FALSE)

      if (!is.na(model$STADF.value)) {
        break
      }
    }
    model$STADF.value
  }
  tmp <- sort(tmp)
  names(tmp) <- NULL
  .cval_SADF_without_const[[as.character(step)]] <- tmp
  cat("\nLen:", length(tmp), "\n")

  cat("with CONST\n")
  tmp <- NULL
  tmp <- foreach(
    i = 1:N_sim,
    .combine = c,
    .inorder = FALSE,
    .errorhandling = "remove",
    .packages = c("breaktest"),
    .options.snow = list(progress = progress)
  ) %dopar% {
    repeat {
      y <- rnorm(step)
      y <- cumsum(y)

      model <- STADF.test(y, const = TRUE, add.p.value = FALSE)

      if (!is.na(model$STADF.value)) {
        break
      }
    }
    model$STADF.value
  }
  tmp <- sort(tmp)
  names(tmp) <- NULL
  .cval_SADF_with_const[[as.character(step)]] <- tmp
  cat("\nLen:", length(tmp), "\n")

  cat("without CONST\n")
  tmp <- NULL
  tmp <- foreach(
    i = 1:N_sim,
    .combine = c,
    .inorder = FALSE,
    .errorhandling = "remove",
    .packages = c("breaktest"),
    .options.snow = list(progress = progress)
  ) %dopar% {
    repeat {
      y <- rnorm(step)
      y <- cumsum(y)

      model <- GSTADF.test(y, const = FALSE, add.p.value = FALSE)

      if (!is.na(model$GSTADF.value)) {
        break
      }
    }
    model$GSTADF.value
  }
  tmp <- sort(tmp)
  names(tmp) <- NULL
  .cval_GSADF_without_const[[as.character(step)]] <- tmp
  cat("\nLen:", length(tmp), "\n")

  cat("with CONST\n")
  tmp <- NULL
  tmp <- foreach(
    i = 1:N_sim,
    .combine = c,
    .inorder = FALSE,
    .errorhandling = "remove",
    .packages = c("breaktest"),
    .options.snow = list(progress = progress)
  ) %dopar% {
    repeat {
      y <- rnorm(step)
      y <- cumsum(y)

      model <- GSTADF.test(y, const = TRUE, add.p.value = FALSE)

      if (!is.na(model$GSTADF.value)) {
        break
      }
    }
    model$GSTADF.value
  }
  tmp <- sort(tmp)
  names(tmp) <- NULL
  .cval_GSADF_with_const[[as.character(step)]] <- tmp
  cat("\nLen:", length(tmp), "\n")

  save(
    .cval_SADF_without_const,
    .cval_SADF_with_const,
    .cval_GSADF_without_const,
    .cval_GSADF_with_const,
    file = "res_tmp.rda"
  )
}

stopCluster(cluster)

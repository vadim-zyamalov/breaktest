hbm.sim <- function(N,
                    mu = 0,
                    tau,
                    sigma0,
                    sigma1,
                    x0 = 0,
                    vs = "shift",
                    tau.break,
                    delta.break) {
  e <- rnorm(N)

  if (vs == "shift") {
    sigma <- sigma0 + (sigma1 - sigma0) *
      ifelse((1:N) / N > tau, 1, 0)
  } else if (vs == "double_shift") {
    sigma <- sigma0 + (sigma1 - sigma0) *
      ifelse(((1:N) / N > 0.4) & ((1:N) / N <= 0.6), 1, 0)
  } else if (vs == "lst") {
    sigma <- sigma0 + (sigma1 - sigma0) *
      1 / (1 + exp(-50 * ((1:N) / N - 0.5)))
  } else if (vs == "trend") {
    sigma <- sigma0 + (sigma1 - sigma0) *
      ((1:N) / N)
  }
  eps <- sigma * e
  u <- c(x0)
  for (i in 2:(N + 1)) {
    if (i <= floor(tau.break[1] * N)) {
      u[i] <- u[i - 1] + eps[i - 1]
    } else if (i > floor(tau.break[1] * N) & i <= floor(tau.break[2] * N)) {
      u[i] <- (1 + delta.break[1]) * u[i - 1] + eps[i - 1]
    } else if (i > floor(tau.break[2] * N) & i <= floor(tau.break[3] * N)) {
      u[i] <- (1 - delta.break[2]) * u[i - 1] + eps[i - 1]
    } else {
      u[i] <- u[i - 1] + eps[i - 1]
    }
  }
  y <- mu + u
  if (y[floor(tau.break[2] * N)] < y[floor(tau.break[1] * N)]) {
    y <- (-1) * y
  }
  return(
    list(
      N = N,
      mu = mu,
      sigma0 = sigma0,
      sigma1 = sigma1,
      x0 = x0,
      tau.break = tau.break,
      delta.break = delta.break,
      tau = tau,
      sigma = sigma,
      eps = eps,
      u = u,
      y = y
    )
  )
}

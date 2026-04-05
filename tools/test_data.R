source("d:/git/breaktest/_origin/_other/BP2003-master/R/Utilities.R", echo = FALSE)
source("d:/git/breaktest/_origin/_other/BP2003-master/R/EstDate.R", echo = FALSE)

alpha_0 <- 5
alpha_1 <- 6
beta_0 <- 0.5
beta_1 <- 0.6
sd <- 0.5

steps <- 10000
N <- 1000
m <- 50

e <- rnorm(N + m + 10, mean = 0, sd = sd)

library(breaktest)

# 1 break

end_0 <- N + 10
end_1 <- N + m + 10

Y_a <- seq(0, 0, length.out = end_1)
Y_b <- seq(0, 0, length.out = end_1)
Y_ab <- seq(0, 0, length.out = end_1)
Y_u <- seq(0, 0, length.out = end_1)
X <- seq(0, 0, length.out = end_1)

for (t in 2:(end_0)) {
  Y_a[t] <- Y_a[t - 1] * beta_0 + alpha_0 + e[t]
  Y_b[t] <- Y_b[t - 1] * beta_0 + alpha_0 + e[t]
  Y_ab[t] <- Y_ab[t - 1] * beta_0 + alpha_0 + e[t]
  X[t] <- X[t - 1] * beta_0 + alpha_0 + e[t]
}
for (t in (end_0 + 1):(end_1)) {
  Y_a[t] <- Y_a[t - 1] * beta_0 + alpha_1 + e[t]
  Y_b[t] <- Y_b[t - 1] * beta_1 + alpha_0 + e[t]
  Y_ab[t] <- Y_ab[t - 1] * beta_1 + alpha_1 + e[t]
  X[t] <- X[t - 1] * beta_0 + alpha_0 + e[t]
}
data <- data.frame(
  C = rep(1, 1050), Trend = 1:1050,
  Y_a = Y_a[11:(end_1)], dY_a = Y_a[10:(end_1 - 1)],
  Y_b = Y_b[11:(end_1)], dY_b = Y_b[10:(end_1 - 1)],
  Y_ab = Y_ab[11:(end_1)], dY_ab = Y_ab[10:(end_1 - 1)],
  X = X[11:(end_1)], dX = X[10:(end_1 - 1)]
)

# 2 breaks

end_0 <- N / 2 + 10
end_1 <- N + 10
end_2 <- N + m + 10

Y_a <- seq(0, 0, length.out = end_1)
Y_b <- seq(0, 0, length.out = end_1)
Y_ab <- seq(0, 0, length.out = end_1)
X <- seq(0, 0, length.out = end_1)

for (t in 2:(end_0)) {
  Y_a[t] <- Y_a[t - 1] * beta_0 + alpha_0 + e[t]
  Y_b[t] <- Y_b[t - 1] * beta_0 + alpha_0 + e[t]
  Y_ab[t] <- Y_ab[t - 1] * beta_0 + alpha_0 + e[t]
  X[t] <- X[t - 1] * beta_0 + alpha_0 + e[t]
}
for (t in (end_0 + 1):(end_1)) {
  Y_a[t] <- Y_a[t - 1] * beta_0 + alpha_1 + e[t]
  Y_b[t] <- Y_b[t - 1] * beta_1 + alpha_0 + e[t]
  Y_ab[t] <- Y_ab[t - 1] * beta_1 + alpha_1 + e[t]
  X[t] <- X[t - 1] * beta_0 + alpha_0 + e[t]
}
for (t in (end_1 + 1):(end_2)) {
  Y_a[t] <- Y_a[t - 1] * beta_0 + alpha_0 + e[t]
  Y_b[t] <- Y_b[t - 1] * beta_0 + alpha_0 + e[t]
  Y_ab[t] <- Y_ab[t - 1] * beta_0 + alpha_0 + e[t]
  X[t] <- X[t - 1] * beta_0 + alpha_0 + e[t]
}
data2 <- data.frame(
  C = rep(1, 1050),
  Y_a = Y_a[11:(end_2)], dY_a = Y_a[10:(end_2 - 1)],
  Y_b = Y_b[11:(end_2)], dY_b = Y_b[10:(end_2 - 1)],
  Y_ab = Y_ab[11:(end_2)], dY_ab = Y_ab[10:(end_2 - 1)],
  X = X[11:(end_2)], dX = X[10:(end_2 - 1)]
)


#################
# Tests 1 break #
#################
SSR.1 <- breaktest:::SSR.matrix(data$Y_ab, cbind(data$C, data$dY_ab), width = 4)

res.1.1 <- dating(as.matrix(data$Y_ab), as.matrix(cbind(data$C, data$dY_ab)), 4, 1, 1, 1005)
res.1.2 <- breaktest:::segments.OLS.mlt(data$Y_ab, cbind(data$C, data$dY_ab), 1, width = 4, SSR.data = SSR.1)

res.1.1$datevec
res.1.1$glb

res.1.2$break.point
res.1.2$SSR

# Carrion-i-Silvestre and Sansó 2006 OBES
KPSS.1.2 <- KPSS.1br(data$Y_ab, data$dY_ab, 1, res.1.2$break.point, FALSE, 4)
KPSS.1.2b <- KPSS.mlt(data$Y_ab,
  data$dY_ab,
  model = 1,
  bp = res.1.2$break.point,
  const = TRUE,
  trend = FALSE,
  lags.init = 4,
  leads.init = 4,
  kernel = NULL,
  weakly.exog = FALSE
)
KPSS.1.2boot <- KPSS.mlt.bootstrap(data$Y_ab,
  data$dY_ab,
  model = 1,
  bp = res.1.2$break.point,
  const = TRUE,
  trend = FALSE,
  lags.init = 4,
  leads.init = 4,
  kernel = NULL,
  weakly.exog = FALSE,
  max.lag = 0
)

##################
# Tests 2 breaks #
##################

res.2.1 <- dating(as.matrix(data$Y_ab), as.matrix(cbind(data$C, data$dY_ab)), 4, 2, 1, 1005)
res.2.2 <- breaktest:::segments.OLS.mlt(data2$Y_ab, cbind(data2$C, data2$dY_ab), 2, width = 4, SSR.data = SSR.1)

res.2.1$datevec
res.2.1$glb

res.2.2$break.point
res.2.2$SSR

# Carrion-i-Silvestre and Sansó 2007
KPSS.2.2 <- KPSS.2br(as.matrix(data$Y_ab), 1, res.2.2$break.point, 4, NULL)
KPSS.2.2b <- KPSS.mlt(data2$Y_ab,
  data2$dY_ab,
  model = c(1, 1),
  bp = res.2.2$break.point,
  const = FALSE,
  trend = TRUE,
  lags.init = 4,
  leads.init = 4,
  kernel = NULL,
  weakly.exog = FALSE
)
KPSS.2.2boot <- KPSS.mlt.bootstrap(data$Y_ab,
  data$dY_ab,
  model = c(1, 1),
  bp = res.2.2$break.point,
  const = FALSE,
  trend = TRUE,
  lags.init = 4,
  leads.init = 4,
  kernel = NULL,
  weakly.exog = FALSE,
  max.lag = 0
)

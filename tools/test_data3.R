source("tools/hbm.sim.R")

N_sim <- 456
# tau_sim - Moment of volatility shift.
tau_sim <- 0.5
# sigma0, sigma1.
sigma0 <- 1
sigma1 <- 1
# delta_b1_sim - Bubble size.
delta_b1_sim <- 0.02

Y <- hbm.sim(
  N = N_sim, tau = tau_sim, sigma0 = sigma0, sigma1 = sigma1,
  tau.break = c(0.4, 0.6, 1), delta.break = c(delta_b1_sim, 0)
)$y
Z <- cumsum(rnorm(N_sim))

res.Y.sadf <- SADF.test(Y)
res.Y.stadf <- STADF.test(Y)
res.Y.sadf.b <- SADF.bootstrap.test(Y)
res.Y.gsadf <- GSADF.test(Y)
res.Y.gstadf <- GSTADF.test(Y)

res.Z.sadf <- SADF.test(Z)
res.Z.stadf <- STADF.test(Z)
res.Z.sadf.b <- SADF.bootstrap.test(Z)
res.Z.gsadf <- GSADF.test(Z)
res.Z.gstadf <- GSTADF.test(Z)

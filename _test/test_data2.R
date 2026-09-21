Y <- rnorm(100)
Y2 <- Y
Y2[51:100] <- Y2[51:100] + 2
X <- rnorm(100)
data <- data.frame(Y, Y2, X)

SSR.1 <- SSR.matrix(data$Y2, 1:100, width = 4)
res.1.2 <-
  segs.SSR.N.breaks(data$Y2, 1:100, 1, width = 4, SSR.data = SSR.1)
KPSS.1.2a <- KPSS.N.breaks(
  data$Y,
  NULL,
  model = 1,
  break.point = res.1.2$break.point,
  const = FALSE,
  trend = TRUE,
  lags.init = 4,
  leads.init = 4,
  corr.max = 0,
  kernel = NULL,
  weakly.exog = TRUE
)
KPSS.1.2a.boot <- KPSS.N.breaks.bootstrap(
  data$Y,
  NULL,
  model = 1,
  break.point = res.1.2$break.point,
  const = FALSE,
  trend = TRUE,
  lags.init = 4,
  leads.init = 4,
  corr.max = 0,
  kernel = NULL,
  weakly.exog = TRUE
)

KPSS.1.2b <- KPSS.N.breaks(
  data$Y,
  data$X,
  model = 1,
  break.point = res.1.2$break.point,
  const = FALSE,
  trend = TRUE,
  lags.init = 4,
  leads.init = 4,
  corr.max = 0,
  kernel = NULL,
  weakly.exog = FALSE
)
KPSS.1.2b.boot <- KPSS.N.breaks.bootstrap(
  data$Y,
  data$X,
  model = 1,
  break.point = res.1.2$break.point,
  const = FALSE,
  trend = TRUE,
  lags.init = 4,
  leads.init = 4,
  corr.max = 0,
  kernel = NULL,
  weakly.exog = FALSE
)

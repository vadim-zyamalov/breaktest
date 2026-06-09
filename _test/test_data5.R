library(breaktest)
data <- readxl::read_excel("tools/data_asp.xlsx")

wss1 <- uroot.ADF(rnorm(100), max.lag = 12, trend = TRUE, criterion = "aic", modified.criterion = TRUE, boot.iter = 1999, boot.p = TRUE)
urs1 <- uroot.ADF(data$gdp, max.lag = 12, trend = TRUE, criterion = "aic", modified.criterion = TRUE, boot.iter = 1999, boot.p = TRUE)

wss2 <- uroot.ADF(rnorm(100), max.lag = 12, trend = TRUE, criterion = "aic", modified.criterion = TRUE, recursive = TRUE, boot.p = TRUE, boot.iter = 1999)
urs2 <- uroot.ADF(data$gdp, max.lag = 12, trend = TRUE, criterion = "aic", modified.criterion = TRUE, recursive = TRUE, boot.p = TRUE, boot.iter = 1999)

wss3 <- coint.CSS(rnorm(100), const = TRUE, break.type = "c", break.point = 50, boot.p = TRUE, boot.type = "Rademacher")
urs3 <- coint.CSS(data$gdp, const = TRUE, break.type = "c", break.point = 50, boot.p = TRUE, boot.type = "Rademacher")

wss4 <- coint.CSS(rnorm(100), rnorm(100), const = TRUE, weakly.exog = FALSE, max.lags = 10, max.leads = 10, break.type = "c", break.point = 50, boot.p = TRUE)
urs4 <- coint.CSS(data$gdp, data$cons, const = TRUE, weakly.exog = FALSE, max.lags = 10, max.leads = 10, break.type = "c", break.point = 50, boot.p = TRUE)

wss5 <- uroot.robust(rnorm(100), TRUE, TRUE)
wss5
urs5 <- uroot.robust(data$gdp, TRUE, TRUE)
urs5

wss6 <- uroot.robust.mlt(rnorm(100), TRUE, FALSE)
wss6
urs6 <- uroot.robust.mlt(data$gdp, TRUE, FALSE)
urs6

wss7 <- coint.PR(rnorm(100), data$cons, 1)
wss7
urs7 <- coint.PR(data$gdp, data$cons, 1)
urs7

wss8 <- coint.conf.sets(rnorm(100), zb = rnorm(100))
wss8
urs8 <- coint.conf.sets(data$gdp, zb = data$cons)
urs8

wss9 <- uroot.CHLT(rnorm(100), boot.cv = TRUE)
wss9
urs9 <- uroot.CHLT(data$gdp, boot.cv = TRUE)
urs9

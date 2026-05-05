library(breaktest)
data <- readxl::read_excel("D:/data_asp.xlsx")

res1 <- uroot.ADF(data$gdp, max.lag = 12, trend = TRUE, criterion = "aic", modified.criterion = TRUE, boot.iter = 1999, boot.p = TRUE)
res2 <- uroot.ADF(data$gdp, max.lag = 12, trend = TRUE, criterion = "aic", modified.criterion = TRUE, recursive = TRUE, boot.p = TRUE, boot.iter = 1999)

res3 <- coint.CSS(data$gdp, data$cons, const = TRUE, break.type = "c", break.point = 50, boot.p = TRUE, kernel = "Bartlett", max.lag = 10)

res4 <- breaktest::PY.statistic(data$gdp, criterion = "bic", const = TRUE, trend = TRUE)
res5 <- breaktest::KP.seq.statistic(data$gdp, breaks = 0, criterion = "bic", const = TRUE)

res6 <- uroot.robust(data$gdp, TRUE, TRUE)
res6

res7 <- uroot.robust.mlt(data$gdp, TRUE, FALSE)
res7

res8 <- coint.PR(data$gdp, data$cons, 1)
res8

res9 <- coint.conf.sets(data$gdp, zb = data$cons)
res9

res10 <- uroot.CHLT(data$gdp, boot.cv = TRUE)
res10

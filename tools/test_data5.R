library(breaktest)
data <- readxl::read_excel("D:/data_asp.xlsx")

res1 <- ADF.test(data$gdp, max.lag = 12, trend = TRUE, criterion = "aic", modified.criterion = TRUE, boot.iter = 1999, boot.p = TRUE)
res2 <- ADF.test(data$gdp, max.lag = 12, trend = TRUE, criterion = "aic", modified.criterion = TRUE, recursive = TRUE, boot.p = TRUE, boot.iter = 1999)

res3 <- coint.CSS(data$gdp, data$cons, const = TRUE, break.type = "c", break.point = 50, boot.p = TRUE, kernel = "Bartlett", max.lag = 10)

res4 <- breaktest::PY.statistic(data$gdp, criterion = "bic", const = TRUE, trend = TRUE)
res5 <- breaktest::KP.seq.statistic(data$gdp, breaks = 0, criterion = "bic", const = TRUE)

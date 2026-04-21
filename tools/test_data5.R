library(breaktest)
data <- readxl::read_excel("D:/data_asp.xlsx")

res1 <- ADF.test(data$gdp, max.lag = 12, trend = TRUE, criterion = "aic", modified.criterion = TRUE, iter = 1999, bootstrap.p = TRUE)
res2 <- ADF.test(data$gdp, max.lag = 12, trend = TRUE, criterion = "aic", modified.criterion = TRUE, recursive = TRUE, bootstrap.p = TRUE, iter = 1999)

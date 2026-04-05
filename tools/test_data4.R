library(devtools)
library(breaktest)
library(openxlsx)

df <- read.xlsx("D:/cryptos.xlsx")
# Removing the date column.
j <- 1
df <- df[(1 + j * 30):(399 + j * 30), 2:dim(df)[2]]

yt <- df[1:100, 1]

SADF.test(yt)
GSADF.test(yt)

STADF.test(yt)
GSTADF.test(yt)

sb.GSADF.test(yt, iter = 5)

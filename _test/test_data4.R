library(breaktest)
library(openxlsx)

df <- read.xlsx("D:/data_asp.xlsx")
# Removing the date column.
j <- 1
yt <- df[1:100, 1]

SADF.test(yt)
GSADF.test(yt)

STADF.test(yt)
GSTADF.test(yt)

sb.GSADF.test(yt, iter = 5)

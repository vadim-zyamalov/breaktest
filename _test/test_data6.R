library(breaktest)
data <- readxl::read_excel("tools/data_asp.xlsx")
rdata <- list(y = rnorm(100),
              zb = rnorm(100),
              zf = rnorm(100))

wss8 <- coint.conf.sets(rdata$y, zb = rdata$zb)
wss8
urs8 <- coint.conf.sets(data$gdp, zb = data$cons)
urs8

wss8a <- breaktest:::cset_break_coint(
  rdata$y,
  FALSE,
  zb = rdata$zb,
  z_ld = -1,
  z_lg = -1,
  c_level = 0.9
)
wss8a

wss8o <- breaktest:::cset_break_optimal(
  rdata$y,
  zb = rdata$zb,
  c_level = 0.9
)
wss8o

urs8a <- breaktest:::cset_break_coint(
  data$gdp,
  FALSE,
  zb = data$cons,
  z_ld = -1,
  z_lg = -1,
  c_level = 0.9
)

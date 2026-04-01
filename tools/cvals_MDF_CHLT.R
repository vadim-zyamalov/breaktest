g.brk.MZ <- 6
g.brk.ADF <- 3

cv.MZa.lim <- -16.62
cv.MSB.lim <- 0.171
cv.MZt.lim <- -2.85
cv.ers.lim <- -2.85

MZa.cv.vals.lim <- c(
  -23.06, -23.49, -23.83, -23.95, -23.98,
  -23.95, -23.90, -23.77, -23.52, -23.25,
  -22.90, -22.46, -21.83, -21.16, -20.24
)
MSB.cv.vals.lim <- c(
  0.146, 0.144, 0.143, 0.143, 0.143,
  0.143, 0.143, 0.144, 0.145, 0.146,
  0.147, 0.148, 0.151, 0.153, 0.157
)
MZt.cv.vals.lim <- c(
  -3.37, -3.40, -3.42, -3.43, -3.43,
  -3.44, -3.44, -3.42, -3.41, -3.39,
  -3.37, -3.34, -3.29, -3.24, -3.17
)
ADF.cv.vals.lim <- MZt.cv.vals.lim

cbar.vals <- c(
  -17.6, -17.8, -18.2, -18.4, -18.6,
  -18.4, -18.4, -18.2, -18.0, -17.6,
  -17.4, -17.0, -16.6, -16.0, -15.2
)

tau.cbar.MZa.cv.lim <- cbind(
  seq(from = 0.15, by = 0.05, length.out = 15),
  cbar.vals,
  MZa.cv.vals.lim
)
tau.cbar.MSB.cv.lim <- cbind(
  seq(from = 0.15, by = 0.05, length.out = 15),
  cbar.vals,
  MSB.cv.vals.lim
)
tau.cbar.MZt.cv.lim <- cbind(
  seq(from = 0.15, by = 0.05, length.out = 15),
  cbar.vals,
  MZt.cv.vals.lim
)
tau.cbar.ADF.cv.lim <- cbind(
  seq(from = 0.15, by = 0.05, length.out = 15),
  cbar.vals,
  ADF.cv.vals.lim
)

.cval_MDF_CHLT <- list(
  g.brk.MZ = g.brk.MZ,
  g.brk.ADF = g.brk.ADF,
  cv.MZa.lim = cv.MZa.lim,
  cv.MSB.lim = cv.MSB.lim,
  cv.MZt.lim = cv.MZt.lim,
  cv.ers.lim = cv.ers.lim,
  MZa.cv.vals.lim = MZa.cv.vals.lim,
  MSB.cv.vals.lim = MSB.cv.vals.lim,
  MZt.cv.vals.lim = MZt.cv.vals.lim,
  ADF.cv.vals.lim = MZt.cv.vals.lim,
  cbar.vals = cbar.vals,
  tau.cbar.MZa.cv.lim = tau.cbar.MZa.cv.lim,
  tau.cbar.MSB.cv.lim = tau.cbar.MSB.cv.lim,
  tau.cbar.MZt.cv.lim = tau.cbar.MZt.cv.lim,
  tau.cbar.ADF.cv.lim = tau.cbar.ADF.cv.lim
)

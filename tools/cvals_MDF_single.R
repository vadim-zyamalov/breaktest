.commons <- list(
  cv.DF.GLS.m = -1.93078,
  cv.DF.OLS.m = -2.85706,
  cv.DF.GLS.t = -2.84317,
  cv.MDF.GLS = -3.84632,
  cv.DF.OLS.t = -3.39735,
  cv.MDF.GLS.lib = -3.69209,
  sap.ur = 1.1159,
  sap.ur1 = 1.0743,
  sap.ur2 = 1.0408,
  sap.ur3 = 1.0357,
  sap.ur4 = 1.1199,
  sap.ur5 = 1.1987,
  sap.cv.ur.k0 = 1.0171,
  sap.cv.A.k0 = 1.0059,
  sap.cv.ur.k30 = 1
)

.cval_MDF_single <- list(
  t = c(
    .commons,
    list(
      cv.MDF.OLS = -34.7163,
      cv.MDF.OLS.lib = -32.3989,
      cv.HLT = 2.563
    )
  ),
  ct = c(
    .commons,
    list(
      cv.MDF.OLS = -39.99,
      cv.MDF.OLS.lib = -32.8742,
      cv.HLT = 3.162
    )
  )
)

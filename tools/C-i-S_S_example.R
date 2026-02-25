library(breaktest)

x <- cumsum(rnorm(200))
y <- cumsum(rnorm(200))

#################
# Tests 1 break #
#################
SSR.1 <- breaktest:::ssr_matrix(y, x, width = 4)
res.1 <- breaktest:::segments_ols_mulitiple(y, x, 1, width = 4, rss_values = SSR.1)

res.1$break.point
res.1$SSR

# Carrion-i-Silvestre and Sansó 2006 OBES
KPSS.1u <- kpss_single_unknown(y, x, 1, TRUE, 4)
cv.1u <- cvalues_kpss_single(1, KPSS.1u["min(stat)", "tb"], length(y), 3)

KPSS.1 <- kpss_single(y, x, 1, res.1$break.point, TRUE, 4)
cv.1 <- cvalues_kpss_single(1, KPSS.1$break.point, length(KPSS.1$residuals), length(KPSS.1$beta))

#################
# Tests 2 break #
#################
KPSS.2b <- kpss_double_unknown(y, 4, 4, "bartlett")
cv.2b <- cvalues_kpss_double(4, c(KPSS.2b$tb1, KPSS.2b$tb2), 200)

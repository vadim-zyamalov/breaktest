source("tools/cvals_kpss_1p.R")
source("tools/cvals_kpss_2p.R")
source("tools/cvals_gsadf.R")
source("tools/cvals_kp.R")
source("tools/cvals_coint_gh.R")
source("tools/cvals_break_date_cset.R")
source("tools/cvals_PY_single.R")
source("tools/cvals_PY_sequential.R")

save(.cval_kpss_1p,
  .cval_kpss_2p,
  .cval_SADF_without_const,
  .cval_SADF_with_const,
  .cval_GSADF_without_const,
  .cval_GSADF_with_const,
  .cval_KP,
  .cval_coint_gh,
  .cval_break_date_cset,
  .cval_PY_single,
  .cval_PY_sequential,
  file = "sysdata.rda",
  compress = "xz"
)

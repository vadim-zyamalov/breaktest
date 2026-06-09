source("tools/cvals_kpss_1p.R")
source("tools/cvals_kpss_2p.R")
source("tools/cvals_gsadf.R")
source("tools/cvals_kp.R")
source("tools/cvals_coint_gh.R")
source("tools/cvals_break_date_cset.R")
source("tools/cvals_PY_single.R")
source("tools/cvals_PY_sequential.R")
source("tools/cvals_VECM.R")
source("tools/cvals_MDF_single.R")
source("tools/cvals_MDF_multiple.R")
source("tools/cvals_MDF_CHLT.R")
source("tools/cbar_PR.R")
source("tools/cvals_PR.R")
source("tools/cvals_NBCN.R")

save(
  .cval_kpss_1p,
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
  .cval_VECM,
  .cval_MDF_single,
  .cval_MDF_multiple,
  .cval_MDF_CHLT,
  .cbar_PR,
  .cval_PR,
  .cval_NBCN,
  file = "sysdata.rda",
  compress = "xz"
)

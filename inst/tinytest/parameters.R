parameters <- c(
  list(cv_dir = cs_dir_tst),
  yaml::read_yaml(file.path(cs_dir_ctl, "testParameters.yaml"))
)

parameters$nv_tfc_inl <- unlist(parameters$nv_tfc_inl)
parameters$nv_enr_rto <- unlist(parameters$nv_enr_rto)
parameters$nm_olc <- matrix(parameters$nm_olc, ncol = 2L)
parameters$df_cdt <- as.data.frame(parameters$df_cdt, stringsAsFactors = FALSE)

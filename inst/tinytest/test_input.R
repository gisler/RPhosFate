#### demoProject ####
expect_warning(
  {
    demoProject()
    cs_dir_tst <- demoProject()
  },
  info = 'existing "demoProject" returns warning'
)

#### clean-up ####
unlink(cs_dir_tst, recursive = TRUE)

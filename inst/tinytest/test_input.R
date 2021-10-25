#### demoProject ####
expect_warning(
  {
    cs_dir_tst <- demoProject()
    demoProject()
  },
  info = 'existing "demoProject" returns warning'
)

#### clean-up ####
unlink(cs_dir_tst, recursive = TRUE)

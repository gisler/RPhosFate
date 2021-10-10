control <- RPhosFate(
  cv_dir = system.file("tinytest", "testProject", package = "RPhosFate"),
  ls_ini = TRUE
)

cs_dir_tst <- demoProject()

x <- RPhosFate(
  cv_dir = cs_dir_tst,
  ls_ini = TRUE
)



unlink(cs_dir_tst, recursive = TRUE)

#### preparations ####
cs_wd <- getwd()

cs_dir_ctl <- system.file("tinytest", "testProject", package = "RPhosFate")
control <- RPhosFate(
  cv_dir = cs_dir_ctl,
  ls_ini = TRUE
)

source("parameters.R") # nolint

cs_dir_tst <- demoProject()
x <- RPhosFate(
  cv_dir = cs_dir_tst,
  ls_ini = TRUE
)

#### getLayer ####

#### getParameter ####

#### setParameter ####

#### clean-up ####
expect_identical(
  cs_wd,
  getwd(),
  info = "working directory is left untouched (getters and setters)"
)

unlink(cs_dir_tst, recursive = TRUE)

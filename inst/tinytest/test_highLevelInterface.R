#### preparations ####
cs_wd <- getwd()

cs_dir_ctl <- system.file("tinytest", "testProject", package = "RPhosFate")
control <- RPhosFate(
  cv_dir = cs_dir_ctl,
  ls_ini = TRUE
)

cs_dir_tst <- demoProject()
x <- RPhosFate(
  cv_dir = cs_dir_tst,
  ls_ini = TRUE
)

#### RPhosFate and catchment ####
parameters <- list(
  cv_dir = cs_dir_tst,
  ns_slp_min = 11.1,
  ns_slp_max = 99.9,
  ns_rhy_a = 0.11,
  ns_rhy_b = 0.99,
  ns_cha_rto = 0.9,
  ns_man_rip = 0.99,
  ns_man_cha = 0.11,
  ns_dep_ovl = 99.9e-4,
  ns_dep_cha = 11.1e-4,
  nv_enr_rto = c(PP = 11.0),
  nv_tfc_inl = c(SS = 0.9, PP = 0.9),
  iv_fDo = rev(c(32L, 16L, 8L, 64L, 0L, 4L, 128L, 1L, 2L)),
  nm_olc = matrix(c(4704255, 2795195), 1L),
  df_cdt = read.table(
    file.path(cs_dir_tst, "cdt.txt"),
    header = TRUE,
    stringsAsFactors = FALSE
  )
)

y <- do.call(RPhosFate, parameters)
z <- do.call(catchment, parameters)

expect_identical(
  y,
  z,
  info = '"RPhosFate" and "catchment" yield identical results (creation)'
)

expect_identical(
  getParameter(y),
  parameters[-1L],
  info = "parameters are created correctly"
)

file.copy(
  file.path(cs_dir_ctl, "testParameters.yaml"),
  file.path(cs_dir_tst, "parameters.yaml"),
  overwrite = TRUE
)

y <- RPhosFate(
  cv_dir = cs_dir_tst,
  ls_ini = TRUE
)
z <- catchment(
  cv_dir = cs_dir_tst,
  ls_ini = TRUE
)

expect_identical(
  y,
  z,
  info = '"RPhosFate" and "catchment" yield identical results (loading)'
)

expect_identical(
  getParameter(y),
  parameters[-1L],
  info = "parameters are loaded correctly"
)

#### firstRun ####
x <- firstRun(x, "SS")

layers <- c("inl", "LFa", "rhy", "rip", "SFa", "slp_cap", "ero")
for (layer in layers) {
  expect_true(
    raster::all.equal(
      getLayer(x      , layer),
      getLayer(control, layer)
    ),
    info = 'substance independent "firstRun" outputs are correct'
  )
}

layers <- c("xxr", "xxt", "xxt_cld", "xxt_ctf", "xxt_inp", "xxt_out")
substances <- slotNames(control@substances)
layers <- list(
  layer = c(rep("xxe", length(substances) - 1L), layers),
  substance = c(setdiff(substances, "SS"), rep("SS", length(layers)))
)
for (i in seq_along(layers$layer)) {
  expect_true(
    raster::all.equal(
      getLayer(x      , layers$layer[i], layers$substance[i]),
      getLayer(control, layers$layer[i], layers$substance[i])
    ),
    info = 'substance dependent "firstRun" outputs are correct'
  )
}

#### subsequentRun ####
layers <- c("xxr", "xxt", "xxt_cld", "xxt_ctf", "xxt_inp", "xxt_out")
for (emissiveSubstance in setdiff(slotNames(control@substances), "SS")) {
  x <- subsequentRun(x, emissiveSubstance)
  for (layer in layers) {
    expect_true(
      raster::all.equal(
        getLayer(x      , layer, emissiveSubstance),
        getLayer(control, layer, emissiveSubstance)
      ),
      info = sprintf(
        '%s "subsequentRun" outputs are correct',
        emissiveSubstance
      )
    )
  }
}

#### snapGauges ####
df_cdt <- getParameter(snapGauges(x), "df_cdt")

expect_identical(
  df_cdt$ID,
  c("G3", "G2", "G1"),
  info = "IDs of gauges are left untouched"
)

expect_identical(
  df_cdt$x,
  c(4704255, 4704195, 4704065),
  info = "x-coordinates are correct"
)

expect_identical(
  df_cdt$y,
  c(2795195, 2795375, 2795585),
  info = "y-coordinates are correct"
)

#### calibrationQuality ####

#### autoCalibrate ####

#### saveState ####
saveState(x)

expect_identical(
  readRDS(file.path(cs_dir_tst, "order.rds")),
  readRDS(file.path(cs_dir_ctl, "order.rds")),
  info = '"order.rds" is written correctly'
)

expect_identical(
  yaml::read_yaml(file.path(cs_dir_tst, "parameters.yaml"))[-1L],
  yaml::read_yaml(file.path(cs_dir_ctl, "parameters.yaml"))[-1L],
  info = '"parameters.yaml" is written correctly'
)

#### clean-up ####
expect_identical(
  cs_wd,
  getwd(),
  info = "working directory is left untouched (high level interface)"
)

unlink(cs_dir_tst, recursive = TRUE)

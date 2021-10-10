#### preparations ####
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

#### RPhosFate ####

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

#### calibrationQuality ####

#### autoCalibrate ####

#### saveState ####
saveState(x)

expect_identical(
  readRDS(file.path(cs_dir_tst, "order.rds")),
  readRDS(file.path(cs_dir_ctl, "order.rds"))
)

expect_identical(
  yaml::read_yaml(file.path(cs_dir_tst, "parameters.yaml"))[-1L],
  yaml::read_yaml(file.path(cs_dir_ctl, "parameters.yaml"))[-1L]
)

#### clean-up ####
unlink(cs_dir_tst, recursive = TRUE)

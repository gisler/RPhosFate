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

source("parameters.R", local = TRUE) # nolint

#### RPhosFate and catchment ####
y <- do.call(RPhosFate, parameters)
z <- do.call(catchment, parameters)

expect_identical(
  getParameter(y),
  parameters[-(1:2)],
  info = "parameters are created correctly"
)

expect_identical(
  y,
  z,
  info = '"RPhosFate" and "catchment" yield identical results (creation)'
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
  getParameter(y),
  parameters[-(1:2)],
  info = "parameters are loaded correctly"
)

expect_identical(
  y,
  z,
  info = '"RPhosFate" and "catchment" yield identical results (loading)'
)

y <- RPhosFate(
  cv_dir = cs_dir_tst,
  ls_ini = TRUE,
  ns_slp_min = 33.3
)
parameters$ns_slp_min <- 33.3

expect_identical(
  getParameter(y, "ns_slp_min"),
  33.3,
  info = "overriding saved parameter works"
)

expect_identical(
  getParameter(y),
  parameters[-(1:2)],
  info = "correct parameter is overridden"
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

#### snapGauges ####
x <- snapGauges(x)

df_cdt <- getParameter(x, "df_cdt")

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
for (substance in substances) {
  expect_identical(
    calibrationQuality(x, substance, sprintf("%s_load", substance)),
    readRDS(file.path(cs_dir_ctl, "calibrationQuality.rds"))[[substance]],
    info = "calibration quality is assessed correctly"
  )

  expect_stdout(
    calibrationQuality(x, substance, sprintf("%s_load", substance)),
    info = "calibration quality is printed"
  )
}

#### Monte Carlo simulation mode ####
x <- RPhosFate(
  cv_dir = c(cs_dir_tst, cs_dir_ctl),
  ls_ini = TRUE,
  is_MCi = 1L
)

expect_identical(
  basename(x@erosion@rl_CFa@file@name),
  "CFa1.tif",
  info = "Monte Carlo input data is detected"
)

x <- subsequentRun(x, "PP", erosion = TRUE, emission = TRUE)

layers <- c("ero1", "ppe1", "ppt1", "ppt_cld1")
for (layer in layers) {
  expect_true(
    file.exists(file.path(cs_dir_tst, "Result", sprintf("%s.tif", layer))),
    info = "Monte Carlo simulation mode outputs exist"
  )
}

#### clean-up ####
expect_identical(
  getwd(),
  cs_wd,
  info = "working directory is left untouched (high level interface)"
)

unlink(cs_dir_tst, recursive = TRUE)

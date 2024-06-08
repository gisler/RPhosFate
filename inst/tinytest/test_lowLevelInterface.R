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

#### erosionPrerequisites ####
x <- erosionPrerequisites(x)

layers <- c("LFa", "SFa", "slp_cap")
for (layer in layers) {
  expect_true(
    terra::all.equal(
      getLayer(x      , layer),
      getLayer(control, layer),
      maxcell = Inf
    ),
    info = '"erosionPrerequisites" outputs are correct'
  )
}

#### erosion ####
x <- erosion(x)

expect_true(
  terra::all.equal(
    getLayer(x      , "ero"),
    getLayer(control, "ero"),
    maxcell = Inf
  ),
  info = '"erosion" output is correct'
)

#### emission ####
for (emissiveSubstance in setdiff(slotNames(control@substances), "SS")) {
  x <- emission(x, emissiveSubstance)

  expect_true(
    terra::all.equal(
      getLayer(x      , "xxe", emissiveSubstance),
      getLayer(control, "xxe", emissiveSubstance),
      maxcell = Inf
    ),
    info = sprintf('%s "emission" output is correct', emissiveSubstance)
  )
}

#### transportPrerequisites ####
x <- transportPrerequisites(x)

layers <- c("inl", "rhy", "rip")
for (layer in layers) {
  expect_true(
    terra::all.equal(
      getLayer(x      , layer),
      getLayer(control, layer),
      maxcell = Inf
    ),
    info = '"transportPrerequisites" outputs are correct'
  )
}

#### transport ####
layers <- c("xxr", "xxt", "xxt_cld", "xxt_ctf", "xxt_inp", "xxt_out")
for (substance in slotNames(control@substances)) {
  x <- transport(x, substance)

  for (layer in layers) {
    expect_true(
      terra::all.equal(
        getLayer(x      , layer, substance),
        getLayer(control, layer, substance),
        maxcell = Inf
      ),
      info = sprintf('%s "transport" outputs are correct', substance)
    )
  }
}

#### clean-up ####
expect_identical(
  getwd(),
  cs_wd,
  info = "working directory is left untouched (low-level interface)"
)

unlink(cs_dir_tst, recursive = TRUE)

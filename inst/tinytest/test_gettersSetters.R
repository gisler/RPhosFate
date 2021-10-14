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

#### getLayer ####
layers <- list.files(cs_dir_ctl, "\\.img$", full.names = TRUE, recursive = TRUE)
for (layer in layers) {
  expect_true(
    raster::all.equal(
      raster::raster(layer),
      {
        layer <- sub("\\.img$", "", basename(layer))
        substance <- regmatches(layer, regexpr(sprintf(
          "(^%s)",
          paste(tolower(slotNames(control@substances)), collapse = "|^")
        ), layer))
        if (length(substance) == 0L) {
          getLayer(control, layer)
        } else {
          getLayer(control, sub(substance, "xx", layer), toupper(substance))
        }
      }
    ),
    info = '"getLayer" works correctly'
  )
}

#### setParameter and getParameter ####
x <- do.call(setParameter, c(list(x), parameters[-(1:2)]))

expect_identical(
  getParameter(x),
  parameters[-(1:2)],
  info = "setting and getting all parameters at once works correctly"
)

for (parameter in names(parameters[-(1:2)])) {
  expect_identical(
    getParameter(x, parameter),
    parameters[[parameter]],
    info = sprintf('getting parameter "%s" works correctly', parameter)
  )
}

#### clean-up ####
expect_identical(
  cs_wd,
  getwd(),
  info = "working directory is left untouched (getters and setters)"
)

unlink(cs_dir_tst, recursive = TRUE)

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

#### getLayer and [ extract operator ####
layers <- list.files(
  cs_dir_ctl,
  "^\\D+\\.tif$",
  full.names = TRUE,
  recursive = TRUE
)
for (layer in layers) {
  expect_true(
    raster::all.equal(
      raster::raster(layer),
      {
        layer <- sub("\\.tif$", "", basename(layer))
        substance <- regmatches(layer, regexpr(sprintf(
          "(^%s)",
          paste(tolower(slotNames(control@substances)), collapse = "|^")
        ), layer))
        if (length(substance) == 0L) {
          rl <- getLayer(control, layer)
        } else {
          rl <- getLayer(control, sub(substance, "xx", layer), toupper(substance))
        }

        rl
      }
    ),
    info = '"getLayer" works correctly'
  )

  expect_identical(
    if (length(substance) == 0L) {
      control[layer]@file@name
    } else {
      control[sub(substance, "xx", layer), toupper(substance)]@file@name
    },
    rl@file@name,
    info = '"[" and "getLayer" refer to the same raster layer'
  )

  expect_error(
    getLayer(x, "non-existent"),
    info = "non-existent layer returns error"
  )
}

#### setParameter and getParameter ####
x <- do.call(setParameter, c(list(x), rev(parameters[-(1:2)])))

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

x <- setParameter(x, nv_tfc_inl = c(Hg = 0.1, PP = 0.1))

expect_identical(
  getParameter(x, "nv_tfc_inl"),
  c(SS = 0.9, PP = 0.1, Hg = 0.1),
  info = "updating and adding substance parameter values works correctly"
)

#### clean-up ####
expect_identical(
  getwd(),
  cs_wd,
  info = "working directory is left untouched (getters and setters)"
)

unlink(cs_dir_tst, recursive = TRUE)

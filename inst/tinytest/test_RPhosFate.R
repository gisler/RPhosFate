control <- RPhosFate(
  cv_dir = system.file("tinytest", "testProject", package = "RPhosFate"),
  ls_ini = TRUE
)

cs_dir_tst <- demoProject()

x <- RPhosFate(
  cv_dir = cs_dir_tst,
  ls_ini = TRUE
)

x <- firstRun(x, "SS")

layers <- list(
  layer = c("inl", "LFa", "rhy", "rip", "SFa", "slp_cap", "ero"),
  substance = NULL
)
for (i in seq_along(layers$layer)) {
  expect_true(
    raster::all.equal(
      getLayer(x      , layers$layer[i], layers$substance[i]),
      getLayer(control, layers$layer[i], layers$substance[i])
    )
  )
}

layers <- list(
  layer = c("xxe", "xxr", "xxt", "xxt_cld", "xxt_ctf", "xxt_inp", "xxt_out"),
  substance = c("PP", rep("SS", 6L))
)
for (i in seq_along(layers$layer)) {
  expect_true(
    raster::all.equal(
      getLayer(x      , layers$layer[i], layers$substance[i]),
      getLayer(control, layers$layer[i], layers$substance[i])
    )
  )
}

x <- subsequentRun(x)

layers <- list(
  layer = c("xxr", "xxt", "xxt_cld", "xxt_ctf", "xxt_inp", "xxt_out"),
  substance = rep("PP", 6L)
)
for (i in seq_along(layers$layer)) {
  expect_true(
    raster::all.equal(
      getLayer(x      , layers$layer[i], layers$substance[i]),
      getLayer(control, layers$layer[i], layers$substance[i])
    )
  )
}

unlink(cs_dir_tst, recursive = TRUE)

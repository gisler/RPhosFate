#### demoProject ####
expect_warning(
  {
    cs_dir_tst <- demoProject()
    demoProject()
  },
  info = 'existing "demoProject" returns warning'
)

unlink(cs_dir_tst, recursive = TRUE)

#### DEMrelatedInput ####
if (.Platform$OS.type == "windows" && tinytest::at_home()) {
  cs_dir_ctl <- system.file("tinytest", "testProject", package = "RPhosFate")
  control <- RPhosFate(
    cv_dir = cs_dir_ctl,
    ls_ini = TRUE
  )

  cv_dir <- normalizePath(
    tempfile("cmt"),
    winslash = .Platform$file.sep,
    mustWork = FALSE
  )
  cs_dir_lrg <- system.file("tinytest", "largeData", package = "RPhosFate")

  nm_olc <- DEMrelatedInput(
    cv_dir = cv_dir,
    cs_dem = file.path(cs_dir_lrg, "dem_lrg.tif"),
    cs_cha = file.path(cs_dir_lrg, "cha_lrg.tif"),
    sp_msk = raster::shapefile(file.path(cs_dir_lrg, "msk.shp")),
    sp_olp = raster::shapefile(file.path(cs_dir_lrg, "olp.shp")),
    sp_sds = raster::shapefile(file.path(cs_dir_lrg, "sds.shp")),
    cs_rds = file.path(cs_dir_lrg, "rds_lrg.tif"),
    cs_wgs = file.path(cs_dir_lrg, "wgs_lrg.tif"),
    ls_tmp = TRUE
  )

  layers <- list.files(
    file.path(cv_dir, "Input"),
    pattern = "^\\D+\\.tif$",
    full.names = TRUE
  )
  for (layer in layers) {
    expect_true(
      raster::all.equal(
        raster::raster(layer),
        getLayer(control, sub("\\.tif$", "", basename(layer)))
      ),
      info = '"DEMrelatedInput" outputs are correct (standard use case)'
    )
  }

  unlink(cv_dir, recursive = TRUE)
}

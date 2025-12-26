devtools::clean_dll()
devtools::load_all()

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
  sp_msk = terra::vect(file.path(cs_dir_lrg, "msk.shp")),
  sp_olp = terra::vect(file.path(cs_dir_lrg, "olp.shp")),
  sp_sds = terra::vect(file.path(cs_dir_lrg, "sds.shp")),
  cs_rds = file.path(cs_dir_lrg, "rds_lrg.tif"),
  ls_tmp = TRUE
)

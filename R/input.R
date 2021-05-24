adjustExtent <- function(rl, ex) {
  extend(crop(rl, ex), ex)
}

#' @export
DEMrelatedInput <- function(
  cv_dir,
  cs_dem,
  cs_cha,
  sp_msk,
  sp_olp,
  sp_sds,
  cs_rds = NULL,
  cs_wgs = NULL,
  ns_brn = 50,
  is_adj = 1L,
  is_ths = 1L,
  ls_tmp = FALSE
) {
  if (!requireNamespace("whitebox", quietly = TRUE)) {
    stop(paste(
      'Package "whitebox" must be installed',
      'from "http://R-Forge.R-project.org" for this functionality.'
    ), call. = FALSE)
  }
  whitebox::wbt_init()
  if (Sys.which("mpiexec") == "" || Sys.which("AreaD8") == "") {
    stop(paste(
      '"TauDEM" must be installed and added to "PATH" environment variable',
      'for this functionality.'
    ), call. = FALSE)
  }
  qassert(cv_dir, "S+")
  qassert(cs_dem, "S1")
  qassert(cs_cha, "S1")
  assertClass(sp_msk, "SpatialPolygonsDataFrame")
  assertClass(sp_olp, "SpatialPointsDataFrame")
  assertClass(sp_sds, "SpatialPointsDataFrame")
  if (!is.null(cs_rds)) {
    qassert(cs_rds, "S1")
  }
  if (!is.null(cs_wgs)) {
    qassert(cs_wgs, "S1")
  }
  qassert(ns_brn, "N1[0,)")
  qassert(is_adj, "X1[0,)")
  qassert(is_ths, "X1[1,)")
  qassert(ls_tmp, "B1")

  cs_dir_old <- setwd(cv_dir[1L])
  on.exit(setwd(cs_dir_old))

  dir.create("Input", showWarnings = FALSE)
  dir.create("temp", showWarnings = FALSE)
  setwd("temp")

  # Extract oversized DEM by mask
  rl_dem_ovr <- raster(cs_dem)
  rl_dem_ovr <- adjustExtent(rl_dem_ovr, sp_msk)
  rl_dem_ovr <- mask(
    rl_dem_ovr,
    sp_msk,
    filename = "dem_ovr.tif",
    datatype = "FLT4S",
    options = "COMPRESSED=YES",
    overwrite = TRUE
  )

  # Burn streams (oversized DEM)
  rl_cha_map <- raster(cs_cha)
  rl_cha_map <- adjustExtent(rl_cha_map, sp_msk)

  rl_dem_bnt <- overlay(
    x = rl_dem_ovr,
    y = rl_cha_map,
    fun = function(x, y) {
      ifelse(is.na(y), x, x - ns_brn)
    }
  )

  for (i in seq_len(is_adj)) {
    rl_cha_map[adjacent(
      rl_cha_map,
      Which(!is.na(rl_cha_map), cells = TRUE),
      directions = 8,
      pairs = FALSE,
      include = TRUE
    )] <- 1L

    rl_dem_bnt <- overlay(
      x = rl_dem_bnt,
      y = rl_cha_map,
      fun = function(x, y) {
        ifelse(is.na(y), x, x - ns_brn)
      }
    )
  }

  writeRaster(
    rl_dem_bnt,
    filename = "dem_bnt.tif",
    datatype = "FLT4S",
    options = "COMPRESSED=YES",
    overwrite = TRUE
  )
  rl_dem_bnt <- raster("dem_bnt.tif")
  rm(rl_cha_map)

  # Breach depressions (oversized DEM)
  whitebox::wbt_breach_depressions(
    dem = file.path(normalizePath("."), "dem_bnt.tif"),
    output = file.path(normalizePath("."), "dem_brd.tif")
  )

  # Calculate D8 flow directions (oversized DEM)
  whitebox::wbt_d8_pointer(
    dem = file.path(normalizePath("."), "dem_brd.tif"),
    output = file.path(normalizePath("."), "dir_ovr.tif"),
    esri_pntr = TRUE
  )

  # Identify watershed
  shapefile(sp_olp, "olp.shp", overwrite = TRUE)
  whitebox::wbt_watershed(
    d8_pntr = file.path(normalizePath("."), "dir_ovr.tif"),
    pour_pts = file.path(normalizePath("."), "olp.shp"),
    output = file.path(normalizePath("."), "wsh_ovr.tif"),
    esri_pntr = TRUE
  )

  rl_wsh <- trim(
    raster("wsh_ovr.tif")
  )
  writeRaster(
    rl_wsh,
    filename = "wsh.tif",
    datatype = "INT1U",
    options = "COMPRESSED=YES",
    overwrite = TRUE
  )
  rl_wsh <- raster("wsh.tif")

  # Extract DEM by watershed
  rl_dem <- mask(
    crop(raster("dem_ovr.tif"), rl_wsh),
    rl_wsh,
    filename = "dem.tif",
    datatype = "FLT4S",
    options = "COMPRESSED=YES",
    overwrite = TRUE
  )

  # Extract flow directions by watershed
  rl_dir <- mask(
    crop(raster("dir_ovr.tif"), rl_wsh),
    rl_wsh,
    filename = "dir.tif",
    datatype = "INT4S",
    options = "COMPRESSED=YES",
    overwrite = TRUE
  )

  # Channel cells
  shapefile(sp_sds, "sds.shp", overwrite = TRUE)
  whitebox::wbt_trace_downslope_flowpaths(
    seed_pts = file.path(normalizePath("."), "sds.shp"),
    d8_pntr = file.path(normalizePath("."), "dir.tif"),
    output = file.path(normalizePath("."), "cha.tif"),
    esri_pntr = TRUE
  )

  rl_cha <- raster("cha.tif")
  rl_cha[!is.na(rl_cha)] <- 1L
  writeRaster(
    rl_cha,
    filename = "cha.tif",
    datatype = "INT1U",
    options = "COMPRESSED=YES",
    overwrite = TRUE
  )
  rl_cha <- raster("cha.tif")

  # Road cells
  if (!is.null(cs_rds)) {
    rl_rds <- raster(cs_rds)
    rl_rds <- adjustExtent(rl_rds, rl_wsh)
    rl_rds[!rl_rds %in% c(0L, 1L)] <- NA_integer_
    rl_rds <- mask(
      rl_rds,
      rl_wsh,
      filename = "rds.tif",
      datatype = "INT1U",
      options = "COMPRESSED=YES",
      overwrite = TRUE
    )
  } else {
    rl_rds <- rl_wsh
    rl_rds[] <- NA_integer_
    writeRaster(
      rl_rds,
      filename = "rds.tif",
      datatype = "INT1U",
      options = "COMPRESSED=YES",
      overwrite = TRUE
    )
    rl_rds <- raster("rds.tif")
  }

  # Flow accumulation
  rl_dir_tau <- subs(
    rl_dir,
    data.frame(
      by    = c(1L, 2L, 4L, 8L, 16L, 32L, 64L, 128L),
      which = c(1L, 8L, 7L, 6L,  5L,  4L,  3L,   2L)
    ),
    filename = "dir_tau.tif",
    datatype = "INT1U",
    options = "COMPRESSED=YES",
    overwrite = TRUE
  )

  system2(
    "mpiexec",
    sprintf(
      "-n %s AreaD8 -nc -p %s -ad8 %s",
      is_ths,
      shQuote(file.path(normalizePath("."), "dir_tau.tif")),
      shQuote(file.path(normalizePath("."), "acc.tif"))
    )
  )

  if (!is.null(cs_wgs)) {
    rl_wgs <- raster(cs_wgs)
    rl_wgs <- adjustExtent(rl_wgs, rl_wsh)
    rl_wgs <- mask(
      rl_wgs,
      rl_wsh,
      filename = "wgs.tif",
      datatype = "FLT4S",
      options = "COMPRESSED=YES",
      overwrite = TRUE
    )

    system2(
      "mpiexec",
      sprintf(
        "-n %s AreaD8 -nc -p %s -ad8 %s -wg %s",
        is_ths,
        shQuote(file.path(normalizePath("."), "dir_tau.tif")),
        shQuote(file.path(normalizePath("."), "acc_wtd.tif")),
        shQuote(file.path(normalizePath("."), "wgs.tif"))
      )
    )
  } else {
    file.copy("acc.tif", "acc_wtd.tif", overwrite = TRUE)
  }

  if (!is.null(cs_rds)) {
    rl_dir_rds <- overlay(
      x = rl_dir_tau,
      y = raster("cha.tif"),
      z = rl_rds,
      fun = function(x, y, z) {
        ifelse(is.na(y), ifelse(is.na(z), x, NA_integer_), x)
      },
      filename = "dir_tau_rds.tif",
      datatype = "INT1U",
      options = "COMPRESSED=YES",
      overwrite = TRUE
    )

    system2(
      "mpiexec",
      sprintf(
        "-n %s AreaD8 -nc -p %s -ad8 %s",
        is_ths,
        shQuote(file.path(normalizePath("."), "dir_tau_rds.tif")),
        shQuote(file.path(normalizePath("."), "acc_rds.tif"))
      )
    )

    rl_acc <- overlay(
      x = rl_cha,
      y = raster("acc_rds.tif"),
      z = raster("acc.tif"),
      fun = function(x, y, z) {
        ifelse(is.na(x), y, z)
      },
      filename = "acc.tif",
      datatype = "INT4S",
      options = "COMPRESSED=YES",
      overwrite = TRUE
    )

    if (!is.null(cs_wgs)) {
      system2(
        "mpiexec",
        sprintf(
          "-n %s AreaD8 -nc -p %s -ad8 %s -wg %s",
          is_ths,
          shQuote(file.path(normalizePath("."), "dir_tau_rds.tif")),
          shQuote(file.path(normalizePath("."), "acc_wtd_rds.tif")),
          shQuote(file.path(normalizePath("."), "wgs.tif"))
        )
      )
    } else {
      file.copy("acc_rds.tif", "acc_wtd_rds.tif", overwrite = TRUE)
    }

    rl_acc_wtd <- overlay(
      x = rl_cha,
      y = raster("acc_wtd_rds.tif"),
      z = raster("acc_wtd.tif"),
      fun = function(x, y, z) {
        ifelse(is.na(x), y, z)
      },
      filename = "acc_wtd.tif",
      datatype = "FLT4S",
      options = "COMPRESSED=YES",
      overwrite = TRUE
    )
  } else {
    rl_acc <- raster("acc.tif")
    rl_acc_wtd <- raster("acc_wtd.tif")
  }

  # Calculate D8 slope (oversized DEM)
  nm_slp_ovr <- D8slope(
    im_dir = as.matrix(raster("dir_ovr.tif")),
    nm_dem = as.matrix(rl_dem_ovr),
    im_fDo = matrix(
      c(32L, 16L, 8L, 64L, 0L, 4L, 128L, 1L, 2L),
      3L
    ),
    ns_fpl = xres(rl_dir),
    is_ths = as.integer(is_ths)
  )

  rl_slp <- mask(
    crop(raster(nm_slp_ovr, template = rl_dem_ovr), rl_wsh),
    rl_wsh,
    filename = "slp.tif",
    datatype = "FLT4S",
    options = "COMPRESSED=YES",
    overwrite = TRUE
  )
  rm(nm_slp_ovr)

  # Copy to input
  toInput <- list(
    acc     = rl_acc    ,
    acc_wtd = rl_acc_wtd,
    cha     = rl_cha    ,
    rds     = rl_rds    ,
    dem     = rl_dem    ,
    dir     = rl_dir    ,
    slp     = rl_slp    ,
    wsh     = rl_wsh
  )
  mapply(
    function(rl, filename) {
      writeRaster(
        rl,
        filename = filename,
        datatype = dataType(rl),
        options = "COMPRESSED=YES",
        overwrite = TRUE
      )
    },
    toInput, file.path("..", "Input", sprintf("%s.%s", names(toInput), "img"))
  )

  # Outlet coordinates
  nm_olc <- xyFromCell(
    rl_acc,
    Which(rl_acc == cellStats(rl_acc, max), cells = TRUE)
  )

  # Clean up
  if (!ls_tmp) {
    setwd("..")
    unlink("temp", recursive = TRUE)
  }

  nm_olc
}

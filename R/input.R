dataLicense <- readLines("inst/tinytest/LICENSE.md")[-1]
#' Demonstration project
#'
#' @description
#' Copies a demonstration project to an existing or a temporary directory.
#'
#' The demonstration project data are a derivative of the
#' `r dataLicense`
#'
#' @param cs_dir An optional character string specifying an existing directory.
#'
#' @return A character string containing the demonstration project root
#'   directory.
#'
#' @seealso [`RPhosFate`], [`catchment`]
#'
#' @examples
#' \dontrun{
#' demoProject()
#' }
#'
#' @export
demoProject <- function(cs_dir = tempdir(TRUE)) {
  assertDirectoryExists(cs_dir, access = "w")

  demoRoot <- file.path(cs_dir, "demoProject")

  if (dir.exists(demoRoot)) {
    warning(
      'A folder called "demoProject" already exists and is left as is.',
      call. = FALSE
    )
  } else {
    dir.create(demoRoot)
    testRoot <- system.file("tinytest", "testProject", package = "RPhosFate")

    file.copy(file.path(testRoot, "Input"), demoRoot, recursive = TRUE)
    file.copy(file.path(testRoot, "cdt.txt"), demoRoot)
    file.copy(file.path(testRoot, "parameters.yaml"), demoRoot)

    writeLines(
      c(
        "These _[RPhosFate](https://gisler.github.io/RPhosFate/)_ demonstration project data are a derivative of the",
        dataLicense
      ),
      file.path(demoRoot, "LICENSE.md")
    )
  }

  normalizePath(demoRoot, winslash = .Platform$file.sep)
}

#' DEM related input
#'
#' @description
#' Clips, pre-processes and calculates or determines all input data related to
#' the digital elevation model (DEM) in the broader sense: \emph{acc, acc_wtd,
#' cha, dem, dir, rds, slp,} and _wsh._
#'
#' Requires _[TauDEM](http://hydrology.usu.edu/taudem/taudem5/downloads.html)_
#' 5.3.7 and the
#' _[WhiteboxTools](https://www.whiteboxgeo.com/download-whiteboxtools/)_ binary
#' ([`whitebox::install_whitebox`]) to be installed on your computer.
#'
#' @param cv_dir A character vector specifying the desired project root
#'   directory (first position).
#' @param cs_dem A character string specifying a path to a potentially large
#'   raster digital elevation model.
#' @param cs_cha A character string specifying a path to a potentially large
#'   raster providing channels.
#' @param sp_msk An [`sp::SpatialPolygonsDataFrame-class`] providing a somewhat
#'   oversized catchment mask used to clip the potentially large input rasters
#'   for further processing.
#' @param sp_olp An [`sp::SpatialPointsDataFrame-class`] providing the desired
#'   catchment outlet.
#' @param sp_sds An [`sp::SpatialPointsDataFrame-class`] providing channel
#'   sources.
#' @param cs_rds An optional character string specifying a path to a potentially
#'   large raster providing roads.
#' @param cs_wgs An optional character string specifying a path to a potentially
#'   large raster providing flow accumulation weights.
#' @param cs_dir An optional character string specifying a path to a potentially
#'   large raster providing D8 flow directions using _ArcGIS_ codes.
#' @param ns_brn A numeric scalar specifying the stream burning step size in m.
#' @param is_adj A numeric scalar specifying how many cells adjacent to channels
#'   shall be burnt.
#' @param is_ths An integer scalar specifying the number of threads to use
#'   during computation (no effect in case _OpenMP_ is not supported by the
#'   toolchain and/or platform).
#' @param ls_tmp A logical scalar specifying if the temporary files created
#'   during computation shall be kept.
#' @param cs_fex A character string specifying the file extension of the created
#'   raster files (either the default `"tif"` or `"img"` for backward
#'   compatibility).
#'
#' @details
#' This function applies the following (pre-processing) steps to ensure
#' hydrologic consistency of the generated input data:
#'
#' * Stream burning and orientation of cells adjacent to channel cells
#' approximately into the direction of channel cells (no effect with `ns_brn =
#' 0`).
#' * Depression breaching.
#' * Tracing of downslope flowpaths from the provided channel sources.
#'
#' When roads are provided, they are considered as flow obstacles breaking the
#' continuity of the calculated flow accumulation.
#'
#' In case no flow accumulation weights are provided, _acc_ and \emph{acc_wtd}
#' are identical.
#'
#' Providing existing flow directions prevents calculating them, which, for
#' example, may be useful in case the effect of tillage directions has been
#' enforced on topographic flow directions in advance. Please note that doing so
#' renders stream burning and depression breaching without effect.
#'
#' _slp_ is calculated from the breached DEM (stream burning is undone
#' beforehand) and represents D8 slopes.
#'
#' @return A numeric [`matrix`] specifying the catchment outlet coordinates.
#'
#' @references
#' \cite{Lindsay, J.B., 2016. Efficient hybrid breaching-filling sink removal
#' methods for flow path enforcement in digital elevation models. Hydrological
#' Processes 30, 846–857.}
#'
#' @seealso [`RPhosFate`], [`catchment`]
#'
#' @examples
#' \dontrun{
#' # obtain temporary project root directory
#' cv_dir <- normalizePath(
#'   tempfile("cmt"),
#'   winslash = .Platform$file.sep,
#'   mustWork = FALSE
#' )
#' # obtain directory holding "large" rasters and other required data sets
#' cs_dir_lrg <- system.file("tinytest", "largeData", package = "RPhosFate")
#'
#' nm_olc <- DEMrelatedInput(
#'   cv_dir = cv_dir,
#'   cs_dem = file.path(cs_dir_lrg, "dem_lrg.tif"),
#'   cs_cha = file.path(cs_dir_lrg, "cha_lrg.tif"),
#'   sp_msk = raster::shapefile(file.path(cs_dir_lrg, "msk.shp")),
#'   sp_olp = raster::shapefile(file.path(cs_dir_lrg, "olp.shp")),
#'   sp_sds = raster::shapefile(file.path(cs_dir_lrg, "sds.shp")),
#'   cs_rds = file.path(cs_dir_lrg, "rds_lrg.tif"),
#'   cs_wgs = file.path(cs_dir_lrg, "wgs_lrg.tif"),
#'   ls_tmp = TRUE
#' )
#' }
#'
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
  cs_dir = NULL,
  ns_brn = 50,
  is_adj = 1L,
  is_ths = 1L,
  ls_tmp = FALSE,
  cs_fex = c("tif", "img")
) {
  if (!requireNamespace("whitebox", quietly = TRUE) ||
      packageVersion("whitebox") < package_version("2.0.0")) {
    stop(paste(
      'Package "whitebox" (v2.0.0 or higher) must be installed for this',
      'functionality.',
      ), call. = FALSE)
  }
  if (!whitebox::check_whitebox_binary()) {
    stop(paste(
      'The "WhiteboxTools" binary must be installed for this functionality.',
      'Consider calling "whitebox::install_whitebox()" first.'
      ), call. = FALSE)
  }
  if (Sys.which("mpiexec") == "" || Sys.which("AreaD8") == "") {
    stop(paste(
      '"TauDEM" must be installed and added to the "PATH" environment variable',
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
  if (!is.null(cs_dir)) {
    qassert(cs_dir, "S1")
  }
  qassert(ns_brn, "N1[0,)")
  qassert(is_adj, "X1[0,)")
  qassert(is_ths, "X1[1,)")
  qassert(ls_tmp, "B1")
  cs_fex <- match.arg(cs_fex)

  dir.create(
    file.path(cv_dir[1L], "Input", "temp"),
    showWarnings = FALSE,
    recursive = TRUE
  )

  cs_dir_old <- setwd(file.path(cv_dir[1L], "Input", "temp"))
  on.exit(setwd(cs_dir_old))

  # Extract oversized DEM by mask
  rl_dem_ovr <- raster(cs_dem)
  rl_dem_ovr <- adjustExtent(rl_dem_ovr, sp_msk)
  rl_dem_ovr <- mask(
    rl_dem_ovr,
    sp_msk,
    filename = "dem_ovr.tif",
    datatype = "FLT8S",
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
    datatype = "FLT8S",
    overwrite = TRUE
  )
  rl_dem_bnt <- raster("dem_bnt.tif")
  rm(rl_cha_map)

  # Breach depressions (oversized DEM)
  whitebox::wbt_breach_depressions(
    dem = file.path(normalizePath("."), "dem_bnt.tif"),
    output = file.path(normalizePath("."), "dem_bnt_brd.tif")
  )

  # Calculate or extract D8 flow directions (oversized DEM)
  if (is.null(cs_dir)) {
    whitebox::wbt_d8_pointer(
      dem = file.path(normalizePath("."), "dem_bnt_brd.tif"),
      output = file.path(normalizePath("."), "dir_ovr.tif"),
      esri_pntr = TRUE
    )
  } else {
    rl_dir_ovr <- raster(cs_dir)
    rl_dir_ovr <- adjustExtent(rl_dir_ovr, sp_msk)
    rl_dir_ovr <- mask(
      rl_dir_ovr,
      sp_msk,
      filename = "dir_ovr.tif",
      datatype = "INT4S",
      overwrite = TRUE
    )
  }

  # Identify watershed
  shapefile(sp_olp, "olp.shp", overwrite = TRUE)
  whitebox::wbt_watershed(
    d8_pntr = file.path(normalizePath("."), "dir_ovr.tif"),
    pour_pts = file.path(normalizePath("."), "olp.shp"),
    output = file.path(normalizePath("."), "wsh_ovr.tif"),
    esri_pntr = TRUE
  )

  rl_wsh <- trim(
    raster("wsh_ovr.tif"),
    filename = "wsh.tif",
    datatype = "INT1U",
    overwrite = TRUE
  )

  # Extract DEM by watershed
  rl_dem <- mask(
    crop(raster("dem_ovr.tif"), rl_wsh),
    rl_wsh,
    filename = "dem.tif",
    datatype = "FLT8S",
    overwrite = TRUE
  )

  # Extract flow directions by watershed
  rl_dir <- mask(
    crop(raster("dir_ovr.tif"), rl_wsh),
    rl_wsh,
    filename = "dir.tif",
    datatype = "INT4S",
    overwrite = TRUE
  )

  # Determine channel cells
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
    overwrite = TRUE
  )
  rl_cha <- raster("cha.tif")

  # Determine road cells
  if (!is.null(cs_rds)) {
    rl_rds <- raster(cs_rds)
    rl_rds <- adjustExtent(rl_rds, rl_wsh)
    rl_rds[!rl_rds %in% c(0L, 1L)] <- NA_integer_
    rl_rds <- mask(
      rl_rds,
      rl_wsh,
      filename = "rds.tif",
      datatype = "INT1U",
      overwrite = TRUE
    )
  } else {
    rl_rds <- rl_wsh
    rl_rds[] <- NA_integer_
    writeRaster(
      rl_rds,
      filename = "rds.tif",
      datatype = "INT1U",
      overwrite = TRUE
    )
    rl_rds <- raster("rds.tif")
  }

  # Calculate flow accumulations
  rl_dir_tau <- subs(
    rl_dir,
    data.frame(
      by    = c(1L, 2L, 4L, 8L, 16L, 32L, 64L, 128L),
      which = c(1L, 8L, 7L, 6L,  5L,  4L,  3L,   2L)
    ),
    filename = "dir_tau.tif",
    datatype = "INT1U",
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
      datatype = "FLT8S",
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
      datatype = "FLT8S",
      overwrite = TRUE
    )
  } else {
    rl_acc <- raster("acc.tif")
    rl_acc_wtd <- raster("acc_wtd.tif")
  }

  # Undo stream burning (breached DEM)
  rl_cha_map <- raster(cs_cha)
  rl_cha_map <- adjustExtent(rl_cha_map, sp_msk)

  rl_dem_brd <- raster("dem_bnt_brd.tif")
  rl_dem_brd <- overlay(
    x = rl_dem_brd,
    y = rl_cha_map,
    fun = function(x, y) {
      ifelse(is.na(y), x, x + ns_brn)
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

    rl_dem_brd <- overlay(
      x = rl_dem_brd,
      y = rl_cha_map,
      fun = function(x, y) {
        ifelse(is.na(y), x, x + ns_brn)
      }
    )
  }

  if (ls_tmp) {
    writeRaster(
      rl_dem_brd,
      filename = "dem_brd.tif",
      datatype = "FLT8S",
      overwrite = TRUE
    )
  }

  # Calculate D8 slopes (oversized DEM)
  nm_slp_ovr <- D8slope(
    im_dir = as.matrix(raster("dir_ovr.tif")),
    nm_dem = as.matrix(rl_dem_brd),
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
    datatype = "FLT8S",
    overwrite = TRUE
  )
  rm(nm_slp_ovr)

  # Copy data to "Input" directory
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
    toInput, file.path("..", sprintf("%s.%s", names(toInput), cs_fex))
  )

  # Determine outlet coordinates
  nm_olc <- xyFromCell(
    rl_acc,
    Which(rl_acc == cellStats(rl_acc, max), cells = TRUE)
  )

  # Clean up temporary files
  if (!ls_tmp) {
    setwd("..")
    unlink("temp", recursive = TRUE)
  }

  nm_olc
}

dataLicense <- readLines("inst/tinytest/LICENSE.md")[-1]
#' Demonstration project
#'
#' @description
#' Copies a demonstration project to an existing or a temporary directory.
#'
#' The demonstration project data are a derivative of the
#' `r paste(dataLicense, collapse = "\n")`
#'
#' @param cs_dir An optional character string specifying an existing directory.
#'
#' @return A character string containing the demonstration project root
#'   directory.
#'
#' @seealso [`RPhosFate`], [`catchment`]
#'
#' @examples
#'
#' demoProject()
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
#' Requires _[TauDEM](https://hydrology.usu.edu/taudem/taudem5/downloads.html)_
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
#' @param sp_msk A [`terra::SpatVector-class`] providing a somewhat oversized
#'   catchment polygon mask used to clip the potentially large input rasters for
#'   further processing.
#' @param sp_olp A [`terra::SpatVector-class`] providing the desired catchment
#'   outlet point(s).
#' @param sp_sds A [`terra::SpatVector-class`] providing channel source points.
#' @param cs_rds An optional character string specifying a path to a potentially
#'   large raster providing roads.
#' @param cs_wgs An optional character string specifying a path to a potentially
#'   large raster providing flow accumulation weights.
#' @param cs_dir An optional character string specifying a path to a potentially
#'   large raster providing D8 flow directions using _ArcGIS_ codes.
#' @param ns_cha An optional numeric scalar specifying the minimum (weighted)
#'   flow accumulation determining a channel.
#' @param ns_brn A numeric scalar specifying the stream burning step size in m.
#' @param is_adj A numeric scalar specifying how many cells adjacent to channels
#'   shall be burnt.
#' @param is_ths An integer scalar specifying the number of threads to use for
#'   processing, where applicable.
#' @param ls_tmp A logical scalar specifying if the temporary files created
#'   during computation shall be kept.
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
#' continuity of the calculated flow accumulations.
#'
#' In case no flow accumulation weights are provided, _acc_ and \emph{acc_wtd}
#' are identical.
#'
#' Providing existing flow directions prevents calculating them, which, for
#' example, may be useful in case the effect of tillage directions has been
#' enforced on topographic flow directions in advance. Please note that doing so
#' renders stream burning and depression breaching without effect.
#'
#' `ns_cha` can be used to enhance the channel network obtained by the tracing
#' of downslope flowpaths from the provided channel sources.
#'
#' _dem_ represents the breached DEM with reversed stream burning if applicable.
#' This processed DEM also serves as the basis for the calculation of the D8
#' slopes provided by _slp._
#'
#' @return A two column numeric [`matrix`] specifying one or more catchment
#'   outlet coordinates and side effects in the form of raster files.
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
#'   sp_msk = terra::vect(file.path(cs_dir_lrg, "msk.shp")),
#'   sp_olp = terra::vect(file.path(cs_dir_lrg, "olp.shp")),
#'   sp_sds = terra::vect(file.path(cs_dir_lrg, "sds.shp")),
#'   cs_rds = file.path(cs_dir_lrg, "rds_lrg.tif"),
#'   cs_wgs = file.path(cs_dir_lrg, "wgs_lrg.tif"),
#'   ls_tmp = TRUE
#' )}
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
  ns_cha = NULL,
  ns_brn = 50,
  is_adj = 1L,
  is_ths = 1L,
  ls_tmp = FALSE
) {
  if (!requireNamespace("whitebox", quietly = TRUE) ||
        packageVersion("whitebox") < package_version("2.0.0")) {
    stop(
      'Package "whitebox" (v2.0.0 or higher) must be installed for this ',
      "functionality.",
      call. = FALSE
    )
  }
  if (!whitebox::check_whitebox_binary()) {
    stop(
      'The "WhiteboxTools" binary must be installed for this functionality. ',
      'Consider calling "whitebox::install_whitebox()" first.',
      call. = FALSE
    )
  }
  if (Sys.which("mpiexec") == "" || Sys.which("AreaD8") == "") {
    stop(
      '"TauDEM" must be installed and added to the "PATH" environment ',
      "variable for this functionality.",
      call. = FALSE
    )
  }
  qassert(cv_dir, "S+")
  qassert(cs_dem, "S1")
  qassert(cs_cha, "S1")
  assertClass(sp_msk, "SpatVector")
  assertTRUE(is.polygons(sp_msk))
  assertClass(sp_olp, "SpatVector")
  assertTRUE(is.points(sp_olp))
  assertClass(sp_sds, "SpatVector")
  assertTRUE(is.points(sp_sds))
  if (!is.null(cs_rds)) {
    qassert(cs_rds, "S1")
  }
  if (!is.null(cs_wgs)) {
    qassert(cs_wgs, "S1")
  }
  if (!is.null(cs_dir)) {
    qassert(cs_dir, "S1")
  }
  if (!is.null(ns_cha)) {
    qassert(ns_cha, "N1[1,)")
  }
  qassert(ns_brn, "N1[0,)")
  qassert(is_adj, "X1[0,)")
  is_ths <- assertCount(is_ths, positive = TRUE, coerce = TRUE)
  qassert(ls_tmp, "B1")

  dir.create(
    file.path(cv_dir[1L], "Input", "temp"),
    showWarnings = FALSE,
    recursive = TRUE
  )

  cs_dir_old <- setwd(file.path(cv_dir[1L], "Input", "temp"))
  on.exit(setwd(cs_dir_old))

  # Extract oversized DEM by mask
  rl_dem_ovr <- rast(cs_dem)
  rl_dem_ovr <- adjustExtent(rl_dem_ovr, sp_msk)
  rl_dem_ovr <- mask(
    rl_dem_ovr,
    sp_msk,
    filename = "dem_ovr.tif",
    datatype = "FLT8S",
    overwrite = TRUE
  )

  # Burn streams (oversized DEM)
  rl_cha_map <- rast(cs_cha)
  rl_cha_map <- adjustExtent(rl_cha_map, sp_msk)

  rl_dem_bnt <- lapp(
    c(x = rl_dem_ovr, y = rl_cha_map),
    fun = function(x, y) {
      ifelse(is.na(y), x, x - ns_brn)
    },
    cores = is_ths
  )

  for (i in seq_len(is_adj)) {
    rl_cha_map[as.integer(adjacent(
      rl_cha_map,
      cells(rl_cha_map),
      directions = "queen",
      include = TRUE
    ))] <- 1L

    rl_dem_bnt <- lapp(
      c(x = rl_dem_bnt, y = rl_cha_map),
      fun = function(x, y) {
        ifelse(is.na(y), x, x - ns_brn)
      },
      cores = is_ths
    )
  }

  writeRaster(
    rl_dem_bnt,
    filename = "dem_bnt.tif",
    datatype = "FLT8S",
    overwrite = TRUE
  )
  rl_dem_bnt <- rast("dem_bnt.tif")
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
    rl_dir_ovr <- rast(cs_dir)
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
  writeVector(sp_olp, "olp.shp", overwrite = TRUE)
  whitebox::wbt_watershed(
    d8_pntr = file.path(normalizePath("."), "dir_ovr.tif"),
    pour_pts = file.path(normalizePath("."), "olp.shp"),
    output = file.path(normalizePath("."), "wsh_ovr.tif"),
    esri_pntr = TRUE
  )

  rl_wsh <- trim(
    rast("wsh_ovr.tif"),
    filename = "wsh.tif",
    datatype = "INT1U",
    overwrite = TRUE
  )

  # Extract flow directions by watershed
  rl_dir <- mask(
    crop(rast("dir_ovr.tif"), rl_wsh),
    rl_wsh,
    filename = "dir.tif",
    datatype = "INT4S",
    overwrite = TRUE
  )

  # Determine channel cells
  writeVector(sp_sds, "sds.shp", overwrite = TRUE)
  whitebox::wbt_trace_downslope_flowpaths(
    seed_pts = file.path(normalizePath("."), "sds.shp"),
    d8_pntr = file.path(normalizePath("."), "dir.tif"),
    output = file.path(normalizePath("."), "cha.tif"),
    esri_pntr = TRUE
  )

  rl_cha <- rast("cha.tif")
  rl_cha[!is.na(rl_cha)] <- 1L
  writeRaster(
    rl_cha,
    filename = "cha.tif",
    datatype = "INT1U",
    overwrite = TRUE
  )
  rl_cha <- rast("cha.tif")

  # Determine road cells
  if (!is.null(cs_rds)) {
    rl_rds <- rast(cs_rds)
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
    rl_rds <- rast("rds.tif")
  }

  # Calculate flow accumulations
  rl_dir_tau <- subst(
    rl_dir,
    from = c(1L, 2L, 4L, 8L, 16L, 32L, 64L, 128L),
    to   = c(1L, 8L, 7L, 6L,  5L,  4L,  3L,   2L),
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
    rl_wgs <- rast(cs_wgs)
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

  # Enhance channels and "backup" traced ones
  if (!is.null(ns_cha)) {
    writeRaster(
      rl_cha,
      filename = "cha_trc.tif",
      datatype = "INT1U",
      overwrite = TRUE
    )

    rl_cha[rast("acc_wtd.tif") >= ns_cha] <- 1L

    writeRaster(
      rl_cha,
      filename = "cha.tif",
      datatype = "INT1U",
      overwrite = TRUE
    )
    rl_cha <- rast("cha.tif")
  }

  # Calculate flow accumulations considering roads
  if (!is.null(cs_rds)) {
    lapp(
      c(
        x = rl_dir_tau,
        y = rl_cha,
        z = rl_rds
      ),
      fun = function(x, y, z) {
        ifelse(is.na(y), ifelse(is.na(z), x, NA_integer_), x) # nolint
      },
      cores = is_ths,
      filename = "dir_tau_rds.tif",
      overwrite = TRUE,
      wopt = list(datatype = "INT1U")
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

    rl_acc <- lapp(
      c(
        x = rl_cha,
        y = rast("acc_rds.tif"),
        z = rast("acc.tif")
      ),
      fun = function(x, y, z) {
        ifelse(is.na(x), y, z)
      },
      cores = is_ths
    )
    writeRaster(
      rl_acc,
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

    rl_acc_wtd <- lapp(
      c(
        x = rl_cha,
        y = rast("acc_wtd_rds.tif"),
        z = rast("acc_wtd.tif")
      ),
      fun = function(x, y, z) {
        ifelse(is.na(x), y, z)
      },
      cores = is_ths
    )
    writeRaster(
      rl_acc_wtd,
      filename = "acc_wtd.tif",
      datatype = "FLT8S",
      overwrite = TRUE
    )
  }

  rl_acc <- rast("acc.tif")
  rl_acc_wtd <- rast("acc_wtd.tif")

  # Undo stream burning (oversized DEM)
  rl_cha_map <- rast(cs_cha)
  rl_cha_map <- adjustExtent(rl_cha_map, sp_msk)

  rl_cha_map_cha <- rl_cha_map
  if (!is.null(ns_cha)) {
    rl_cha <- rast("cha_trc.tif")
  } else {
    rl_cha <- rast("cha.tif")
  }
  rl_cha_map_cha[extend(rl_cha, rl_cha_map_cha) == 1L] <- 1L

  rl_dem_brd <- rast("dem_bnt_brd.tif")
  rl_dem_brd <- lapp(
    c(x = rl_dem_brd, y = rl_cha_map_cha),
    fun = function(x, y) {
      ifelse(is.na(y), x, x + ns_brn)
    },
    cores = is_ths
  )

  for (i in seq_len(is_adj)) {
    rl_cha_map[union(as.integer(adjacent(
      rl_cha_map,
      cells(rl_cha_map),
      directions = "queen",
      include = TRUE
    )), unlist(cells(rl_cha_map_cha, 1L)))] <- 1L

    rl_dem_brd <- lapp(
      c(x = rl_dem_brd, y = rl_cha_map),
      fun = function(x, y) {
        ifelse(is.na(y), x, x + ns_brn)
      },
      cores = is_ths
    )
  }

  # Extract DEM by watershed
  rl_dem <- mask(
    crop(rl_dem_brd, rl_wsh),
    rl_wsh,
    filename = "dem.tif",
    datatype = "FLT8S",
    overwrite = TRUE
  )

  # Calculate D8 slopes (oversized DEM)
  nm_slp_ovr <- D8slope(
    im_dir = as.matrix(rast("dir_ovr.tif"), wide = TRUE),
    nm_dem = as.matrix(rl_dem_brd, wide = TRUE),
    im_fDo = matrix(
      c(32L, 16L, 8L, 64L, 0L, 4L, 128L, 1L, 2L),
      3L
    ),
    ns_fpl = xres(rl_dem),
    is_ths = is_ths
  )

  rl_slp <- mask(
    crop(rast(nm_slp_ovr, crs = crs(rl_dem_ovr), extent = ext(rl_dem_ovr)), rl_wsh),
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
    dem     = rl_dem    ,
    dir     = rl_dir    ,
    rds     = rl_rds    ,
    slp     = rl_slp    ,
    wsh     = rl_wsh
  )
  for (item in names(toInput)) {
    set.names(toInput[[item]], item)

    writeRaster(
      toInput[[item]],
      filename = file.path("..", paste0(item, ".tif")),
      datatype = datatype(toInput[[item]]),
      overwrite = TRUE
    )
  }

  # Determine outlet coordinates
  nm_slp <- D8slope(
    im_dir = as.matrix(rl_dir, wide = TRUE),
    nm_dem = as.matrix(rl_dem, wide = TRUE),
    im_fDo = matrix(
      c(32L, 16L, 8L, 64L, 0L, 4L, 128L, 1L, 2L),
      3L
    ),
    ns_fpl = xres(rl_dem),
    is_ths = is_ths
  )
  rl_slp <- rast(nm_slp, crs = crs(rl_dem), extent = ext(rl_dem))
  rm(nm_slp)

  nm_olc <- xyFromCell(
    rl_acc,
    unlist(cells(is.na(rl_slp) & !is.na(rl_cha), 1L))
  )

  # Clean up temporary files
  if (!ls_tmp) {
    setwd("..")
    unlink("temp", recursive = TRUE)
  }

  nm_olc
}

#' Convert _ERDAS IMAGINE_ to _GeoTIFF_ raster files
#'
#' @description
#' Converts all _ERDAS IMAGINE_ raster files in a directory and its
#' subdirectories into _GeoTIFF_ raster files.
#'
#' @param cs_dir A character string specifying an existing directory.
#' @param cs_crs An optional character string used to set the coordinate
#'   reference system of all output raster files. See [`terra::crs`] for further
#'   information.
#'
#' @return A character vector containing the paths to the processed _ERDAS
#'   IMAGINE_ raster files.
#'
#' @export
img2tif <- function(cs_dir, cs_crs = NULL) {
  assertDirectoryExists(cs_dir, access = "w")

  files <- list.files(
    cs_dir,
    pattern = "\\.img$",
    full.names = TRUE,
    recursive = TRUE
  )

  for (file in files) {
    rl <- rast(file)
    datatype <- datatype(rl)

    if (!is.null(cs_crs)) {
      qassert(cs_crs, "S1")

      crs(rl) <- cs_crs
    }

    writeRaster(
      rl,
      paste0(tools::file_path_sans_ext(file), ".tif"),
      datatype = datatype
    )
  }

  files
}

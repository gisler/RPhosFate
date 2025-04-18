adjustExtent <- function(rl, ex) {
  extend(crop(rl, ex), ex)
}

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
    warning('A folder called "demoProject" already exists and is left as is.')
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
#' the digital elevation model (DEM) in the broader sense: \emph{acc_inf, cha,
#' dem, dir_inf, rds, slp_inf,} and _wsh._
#'
#' Requires the
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
#' @param ns_cha An optional numeric scalar specifying the minimum D8 flow
#'   accumulation in number of upslope grid cells determining a channel.
#' @param ns_brn A numeric scalar specifying the stream burning step size in m.
#' @param is_adj A numeric scalar specifying how many cells adjacent to channels
#'   shall be burnt.
#' @param is_ths An integer scalar specifying the number of threads to use for
#'   processing, where applicable.
#' @param ls_mD8 A logical scalar specifying if D8 flow directions shall be
#'   mimicked, i.e. the D-infinity flow directions are rounded to the nearest
#'   multiple of 45 degrees. Please note that this treatment is always applied
#'   to channel cells independently of this argument.
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
#' `ns_cha` can be used to enhance the channel network obtained by the tracing
#' of downslope flowpaths from the provided channel sources.
#'
#' _dem_ represents the breached DEM with reversed stream burning if applicable.
#' The basis for the calculation of the D-infinity slopes provided by
#' \emph{slp_inf,} however, is the original DEM.
#'
#' @return A two column numeric [`matrix`] specifying one or more catchment
#'   outlet coordinates and side effects in the form of raster files.
#'
#' @references
#' \cite{Lindsay, J.B., 2016. Efficient hybrid breaching-filling sink removal
#' methods for flow path enforcement in digital elevation models. Hydrological
#' Processes 30, 846–857. https://doi.org/10.1002/hyp.10648}
#'
#' \cite{Tarboton, D.G., 1997. A new method for the determination of flow
#' directions and upslope areas in grid digital elevation models. Water Resour.
#' Res. 33, 309–319. https://doi.org/10.1029/96WR03137}
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
  ns_cha = NULL,
  ns_brn = 50,
  is_adj = 1L,
  is_ths = 1L,
  ls_mD8 = FALSE,
  ls_tmp = FALSE
) {
  if (!requireNamespace("whitebox", quietly = TRUE) ||
        packageVersion("whitebox") < package_version("2.0.0")) {
    stop(
      'Package "whitebox" (v2.0.0 or higher) must be installed for this ',
      "functionality."
    )
  }
  if (!whitebox::check_whitebox_binary()) {
    stop(
      'The "WhiteboxTools" binary must be installed for this functionality. ',
      'Consider calling "whitebox::install_whitebox()" first.'
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
  rl_cha_ovr <- rast(cs_cha)
  rl_cha_ovr <- adjustExtent(rl_cha_ovr, sp_msk)

  rl_dem_ovr_bnt <- lapp(
    c(x = rl_dem_ovr, y = rl_cha_ovr),
    fun = function(x, y) {
      ifelse(is.na(y), x, x - ns_brn)
    },
    cores = is_ths
  )

  for (i in seq_len(is_adj)) {
    rl_cha_ovr[as.integer(adjacent(
      rl_cha_ovr,
      cells(rl_cha_ovr),
      directions = "queen",
      include = TRUE
    ))] <- 1L

    rl_dem_ovr_bnt <- lapp(
      c(x = rl_dem_ovr_bnt, y = rl_cha_ovr),
      fun = function(x, y) {
        ifelse(is.na(y), x, x - ns_brn)
      },
      cores = is_ths
    )
  }
  rm(rl_cha_ovr)

  writeRaster(
    rl_dem_ovr_bnt,
    filename = "dem_ovr_bnt.tif",
    datatype = "FLT8S",
    overwrite = TRUE
  )
  rl_dem_ovr_bnt <- rast("dem_ovr_bnt.tif")

  # Breach depressions (oversized DEM)
  whitebox::wbt_breach_depressions(
    dem = "dem_ovr_bnt.tif",
    output = "dem_ovr_bnt_brd.tif"
  )

  # Calculate D8 flow directions (oversized DEM)
  whitebox::wbt_d8_pointer(
    dem = "dem_ovr_bnt_brd.tif",
    output = "dir_ovr.tif",
    esri_pntr = TRUE
  )

  # Identify watershed
  writeVector(sp_olp, "olp.shp", overwrite = TRUE)
  whitebox::wbt_watershed(
    d8_pntr = "dir_ovr.tif",
    pour_pts = "olp.shp",
    output = "wsh_ovr.tif",
    esri_pntr = TRUE
  )

  rl_wsh <- trim(
    rast("wsh_ovr.tif"),
    filename = "wsh.tif",
    datatype = "INT1U",
    overwrite = TRUE
  )

  # Extract D8 flow directions by watershed
  mask(
    crop(rast("dir_ovr.tif"), rl_wsh),
    rl_wsh,
    filename = "dir.tif",
    datatype = "INT4S",
    overwrite = TRUE
  )

  # Trace channel cells
  writeVector(sp_sds, "sds.shp", overwrite = TRUE)
  whitebox::wbt_trace_downslope_flowpaths(
    seed_pts = "sds.shp",
    d8_pntr = "dir.tif",
    output = "cha_trc.tif",
    esri_pntr = TRUE
  )

  rl_cha <- rast("cha_trc.tif")
  rl_cha[!is.na(rl_cha)] <- 1L

  writeRaster(
    rl_cha,
    filename = "cha_trc.tif",
    datatype = "INT1U",
    overwrite = TRUE
  )
  rl_cha <- rast("cha_trc.tif")

  # Enhance channels
  if (!is.null(ns_cha)) {
    whitebox::wbt_d8_flow_accumulation(
      input = "dir.tif",
      output = "acc.tif",
      pntr = TRUE,
      esri_pntr = TRUE
    )

    rl_cha[rast("acc.tif") >= ns_cha] <- 1L

    writeRaster(
      rl_cha,
      filename = "cha.tif",
      datatype = "INT1U",
      overwrite = TRUE
    )
    rl_cha <- rast("cha.tif")
  }

  # Extract burnt and breached DEM by watershed
  mask(
    crop(rast("dem_ovr_bnt_brd.tif"), rl_wsh),
    rl_wsh,
    filename = "dem_bnt_brd.tif",
    datatype = "FLT8S",
    overwrite = TRUE
  )

  # Calculate DInf flow directions
  whitebox::wbt_d_inf_pointer(
    dem = "dem_bnt_brd.tif",
    output = "dir_inf.tif"
  )
  rl_dir_inf <- rast("dir_inf.tif")

  rl_dir_inf[rl_dir_inf == -1] <- NA_real_

  if (ls_mD8) {
    rl_dir_inf <- app(
      rl_dir_inf,
      function(x) {
        round(x / 45) * 45
      },
      cores = is_ths
    )
  } else {
    rl_dir_inf <- lapp(
      c(x = rl_dir_inf, y = rl_cha),
      function(x, y) {
        ifelse(is.na(y), x, round(x / 45) * 45)
      },
      cores = is_ths
    )
  }

  # Determine outlet coordinates
  nm_olc <- xyFromCell(
    rl_cha,
    unlist(cells(is.na(rl_dir_inf) & !is.na(rl_cha), 1L))
  )

  # Complement DInf at outlet coordinates
  DInfInsteadD8 <- c(
      "1" =  90,
      "2" = 135,
      "4" = 180,
      "8" = 225,
     "16" = 270,
     "32" = 315,
     "64" =   0,
    "128" =  45
  )

  df_dir_olc <- rast("dir.tif")[cellFromXY(rl_dir_inf, nm_olc)]
  df_dir_olc[[1L]] <- DInfInsteadD8[as.character(df_dir_olc[[1L]])]
  rl_dir_inf[cellFromXY(rl_dir_inf, nm_olc)] <- df_dir_olc

  writeRaster(
    rl_dir_inf,
    filename = "dir_inf.tif",
    datatype = "FLT8S",
    overwrite = TRUE
  )
  rl_dir_inf <- rast("dir_inf.tif")

  # Calculate DInf flow accumulations
  whitebox::wbt_d_inf_flow_accumulation(
    input = "dir_inf.tif",
    output = "acc_inf.tif",
    out_type = "cells",
    pntr = TRUE
  )

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

  # Calculate DInf flow accumulations considering roads
  if (!is.null(cs_rds)) {
    lapp(
      c(
        x = rl_dir_inf,
        y = rl_cha,
        z = rl_rds
      ),
      fun = function(x, y, z) {
        # nolint start
        ifelse(
          is.na(y),
          ifelse(
            is.na(z),
            x,
            NA_integer_
          ),
          x
        )
        # nolint end
      },
      cores = is_ths,
      filename = "dir_inf_rds.tif",
      overwrite = TRUE,
      wopt = list(datatype = "FLT8S")
    )

    whitebox::wbt_d_inf_flow_accumulation(
      input = "dir_inf_rds.tif",
      output = "acc_inf_rds.tif",
      out_type = "cells",
      pntr = TRUE
    )

    rl_acc_inf <- lapp(
      c(
        x = rl_cha,
        y = rast("acc_inf_rds.tif"),
        z = rast("acc_inf.tif")
      ),
      fun = function(x, y, z) {
        ifelse(is.na(x), y, z)
      },
      cores = is_ths
    )

    writeRaster(
      rl_acc_inf,
      filename = "acc_inf.tif",
      datatype = "FLT8S",
      overwrite = TRUE
    )
  }

  rl_acc_inf <- rast("acc_inf.tif")

  # Undo stream burning (oversized DEM)
  rl_cha_ovr <- rast(cs_cha)
  rl_cha_ovr <- adjustExtent(rl_cha_ovr, sp_msk)

  rl_cha_ovr_cha <- rl_cha_ovr
  rl_cha_ovr_cha[extend(rast("cha_trc.tif"), rl_cha_ovr_cha) == 1L] <- 1L

  rl_dem_ovr_brd <- lapp(
    c(x = rast("dem_ovr_bnt_brd.tif"), y = rl_cha_ovr_cha),
    fun = function(x, y) {
      ifelse(is.na(y), x, x + ns_brn)
    },
    cores = is_ths
  )

  for (i in seq_len(is_adj)) {
    rl_cha_ovr[union(as.integer(adjacent(
      rl_cha_ovr,
      cells(rl_cha_ovr),
      directions = "queen",
      include = TRUE
    )), unlist(cells(rl_cha_ovr_cha, 1L)))] <- 1L

    rl_dem_ovr_brd <- lapp(
      c(x = rl_dem_ovr_brd, y = rl_cha_ovr),
      fun = function(x, y) {
        ifelse(is.na(y), x, x + ns_brn)
      },
      cores = is_ths
    )
  }
  rm(rl_cha_ovr, rl_cha_ovr_cha)

  # Extract breached DEM by watershed
  rl_dem_brd <- mask(
    crop(rl_dem_ovr_brd, rl_wsh),
    rl_wsh,
    filename = "dem_brd.tif",
    datatype = "FLT8S",
    overwrite = TRUE
  )

  # Calculate (original oversized DEM) and extract DInf slopes by watershed
  nm_slp_inf_ovr <- dinfSlopeCpp(
    nm_dir_inf = as.matrix(extend(rl_dir_inf, rl_dem_ovr), wide = TRUE),
    nm_dem = as.matrix(rl_dem_ovr, wide = TRUE),
    ns_res = xres(rl_dem_ovr),
    is_ths = is_ths
  )

  rl_slp_inf <- crop(
    rast(
      nm_slp_inf_ovr,
      crs = crs(rl_dem_ovr),
      extent = ext(rl_dem_ovr)
    ),
    rl_wsh,
    filename = "slp_inf.tif",
    datatype = "FLT8S",
    overwrite = TRUE
  )
  rm(nm_slp_inf_ovr)

  # Copy data to "Input" directory
  toInput <- list(
    acc_inf = rl_acc_inf,
    cha     = rl_cha    ,
    dem     = rl_dem_brd,
    dir_inf = rl_dir_inf,
    rds     = rl_rds    ,
    slp_inf = rl_slp_inf,
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

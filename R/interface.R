#' @include aaa.R
NULL

#### Constructors ####
#' @rdname catchment
#' @export
RPhosFate <- function(...) {
  new("RPhosFate", list(...))
}
#' Initialise project
#'
#' @description
#' Initialises a project from scratch or loads the state of an existing one
#' utilising _GeoTIFF_ (*.tif) raster files from, by convention, the following
#' three project root subdirectories:
#'
#' * _Input_
#' * _Intermediate_
#' * _Result_
#'
#' See subdirectory sections for further information.
#'
#' `catchment` is an alias for `RPhosFate`.
#'
#' @param \dots Arguments used to initialise the project. See argument sections
#'   for further information.
#'
#' @section _Input_ subdirectory:
#' This directory holds all possible user input raster data (flow obstacles like
#' roads must be considered during generation of the flow accumulation layers
#' and also be cut out from them in order to be properly respected):
#'
#' * _acc:_ Flow accumulations required for [`transportCalcOrder`].
#' * \emph{acc_wtd:} Weighted flow accumulations required for everything (can be
#' equal to _acc_).
#' * _CFa:_ (R)USLE C-factors required for [`erosion`].
#' * _cha:_ Channel cells required for everything (`1`: channel cell, `NA`: no
#' channel cell).
#' * _clc:_ Clay contents of top soils in % required for substance
#' [`emission`]s.
#' * _dem:_ Digital elevation model in m a.s.l. (optional).
#' * _dir:_ D8 flow directions required for [`transportPrerequisites`] and
#' substance [`transport`].
#' * _fid:_ Field IDs (optional).
#' * _KFa:_ (R)USLE K-factors required for [`erosion`].
#' * _lue:_ Land use classes (optional).
#' * _man:_ Manning's roughness coefficients required for substance
#' [`transport`].
#' * _xxc:_ Substance contents of top soils in mg/kg required for substance
#' [`emission`]s, for example, _ppc_ for PP top soil contents.
#' * _rds:_ Road cells required for [`transportPrerequisites`] (`0`: road cell
#' without subsurface drainage, `1`: road cell with subsurface drainage, `NA`:
#' no road cell).
#' * _RFa:_ (R)USLE R-factors required for [`erosion`].
#' * _slp:_ Slopes in % required for everything.
#' * _wsh:_ Watershed (optional).
#'
#' @section _Intermediate_ subdirectory:
#' This directory holds intermediate calculations:
#'
#' * _inl:_ Cells representing inlets at roads (storm drains).
#' * _LFa:_ L-factors.
#' * _rhy:_ Hydraulic radii in m.
#' * _rip:_ Cells representing the riparian zones within channel cells.
#' * _SFa:_ RUSLE S-factors.
#' * \emph{slp_cap:} Capped slopes in %.
#'
#' @section _Result_ subdirectory:
#' This directory holds the model results:
#'
#' * _ero:_ Erosion in t/cell/yr.
#' * _xxe:_ Substance emissions in kg/cell/yr, for example, _ppe_ for PP
#' emissions.
#' * _xxr:_ Substance retentions in t/cell/yr (SS) or kg/cell/yr, for example,
#' _ppr_ for PP retentions.
#' * _xxt:_ Substance transports in t/cell/yr (SS) or kg/cell/yr, for example,
#' _ppt_ for PP transports.
#' * \emph{xxt_cld:} Substance cell loads in t/cell/yr (SS) or kg/cell/yr, for
#' example, \emph{ppt_cld} for PP cell loads.
#' * \emph{xxt_ctf:} Substance cell transfers in t/cell/yr (SS) or kg/cell/yr,
#' for example, \emph{ppt_ctf} for PP transfers.
#' * \emph{xxt_inp:} Substance inputs into surface waters in t/cell/yr (SS) or
#' kg/cell/yr, for example, \emph{ppt_inp} for PP inputs into surface waters.
#' * \emph{xxt_out:} Substance outlet loads of subsurface drainages in t/cell/yr
#' (SS) or kg/cell/yr, for example, \emph{ppt_out} for PP outlet loads.
#'
#' @section Data management arguments:
#' * `cv_dir`: A character vector specifying the project root (first position)
#' and optionally the Monte Carlo input data directory (second position).
#' * `ls_ini`: A logical scalar specifying if the state of an existing project
#' shall be loaded from disk (defaults to `FALSE`). Parameters or substance
#' parameter values specified via the `...` argument take precedence over loaded
#' ones.
#' * `is_MCi`: An integer scalar specifying the current Monte Carlo iteration if
#' applicable (defaults to `integer()`, which means Monte Carlo simulation mode
#' is disabled).
#' * `cv_MCl`: A character vector specifying the names of the layers, which
#' shall be written to disk with the associated Monte Carlo iteration in their
#' filenames upon calling the appropriate methods (defaults to `"xxt"`; no
#' effect in case Monte Carlo simulation mode is disabled).
#'
#' @section Model parameter arguments:
#' * `ns_slp_min`: A numeric scalar specifying the minimum bounding slope in %
#' (defaults to `0.001`).
#' * `ns_slp_max`: A numeric scalar specifying the maximum bounding slope in %
#' (defaults to `999.0`).
#' * `ns_rhy_a`: A numeric scalar specifying a network constant depending on the
#' discharge frequency needed for the calculation of the hydraulic radius, which
#' in turn is a prerequisite for substance transport (defaults to `0.09`
#' representing a discharge frequency of approximately six years).
#' * `ns_rhy_b`: A numeric scalar specifying a geometry scaling exponent
#' depending on the discharge frequency needed for the calculation of the
#' hydraulic radius, which in turn is a prerequisite for substance transport
#' (defaults to `0.50` representing a discharge frequency of approximately six
#' years).
#' * `ns_cha_rto`: A numeric scalar specifying the ratio of the channel to the
#' cell width determining the widths of the riparian zones required for
#' substance [`transport`] (defaults to `0.5`).
#' * `ns_man_rip`: A numeric scalar specifying Manning's roughness coefficient
#' of the riparian zones within channel cells required for substance
#' [`transport`] (defaults to `0.32`).
#' * `ns_man_cha`: A numeric scalar specifying Manning's roughness coefficient
#' of the channel within channel cells required for substance [`transport`]
#' (defaults to `0.04`).
#' * `ns_dep_ovl`: A numeric scalar specifying the overland deposition rate per
#' second required for substance [`transport`] (calibration parameter; no
#' default).
#' * `ns_dep_cha`: A numeric scalar specifying the channel deposition rate per
#' second required for substance [`transport`] (calibration parameter; no
#' default).
#' * `nv_tfc_inl`: A named numeric vector specifying the inlet transfer
#' coefficients required for substance [`transport`], for example, `c(SS = 0.6,
#' PP = 0.6)` (no default).
#' * `nv_enr_rto` A named numeric vector specifying the substance enrichment
#' ratios required for substance except SS [`transport`], for example, `c(PP =
#' 2.0)` (calibration parameter; no default).
#' * `iv_fDo`: An integer vector specifying the outflow direction vector
#' required for substance [`transport`] (defaults to _ArcGIS_ codes).
#' * `nm_olc`: A two column numeric [`matrix`] specifying one or more catchment
#' outlet coordinates required for the in-channel retention ratio of
#' [`calibrationQuality`] (no default).
#' * `df_cdt`: A [`data.frame`] with calibration data, which must have at least
#' the following three columns and one or more columns with substance river
#' loads in t/yr (no default):
#'
#'   * _ID:_ ID(s) of the gauge(s)
#'   * _x:_ x-coordinate(s) of the gauge(s)
#'   * _y:_ y-coordinate(s) of the gauge(s)
#'
#' @section Monte Carlo simulation mode:
#' This mode can make use of repeated random samples, i.e. raster data, of
#' distributions of about all input data. The filenames of the Monte Carlo input
#' raster data must contain the specified iteration, for example, _CFa12.tif_
#' for the twelfth iteration of the C-factors input data, and can reside in a
#' separate directory. In case no Monte Carlo raster file is found for a certain
#' layer in the designated directory, the respective project root subdirectory
#' is searched for one and finally the \dQuote{normal} project input raster data
#' is utilised.
#'
#' @return An S4 [`RPhosFate-class`] river catchment object.
#'
#' @seealso [`saveState`], [`demoProject`]
#'
#' @examples
#'
#' # temporary demonstration project copy
#' cv_dir <- demoProject()
#'
#' # initialise project from scratch
#' x <- RPhosFate(
#'   cv_dir = cv_dir,
#'   ns_dep_ovl = 25e-4,
#'   ns_dep_cha = 0.0,
#'   nv_tfc_inl = c(SS = 0.6, PP = 0.6),
#'   nv_enr_rto = c(PP = 2.0),
#'   nm_olc = matrix(c(4704255, 2795195), ncol = 2L),
#'   df_cdt = read.table(
#'     file.path(cv_dir, "cdt.txt"),
#'     header = TRUE,
#'     stringsAsFactors = FALSE
#'   )
#' )
#'
#' # load state of existing project in Monte Carlo simulation mode
#' x <- RPhosFate(
#'   cv_dir = c(
#'     cv_dir,
#'     system.file("tinytest", "testProject", package = "RPhosFate")
#'   ),
#'   ls_ini = TRUE,
#'   is_MCi = 1L,
#'   cv_MCl = c("xxt", "xxt_cld")
#' )
#'
#' @export
catchment <- function(...) {
  new("RPhosFate", list(...))
}

#### firstRun ####
#' @export
setGeneric(
  "firstRun",
  function(x, ...) standardGeneric("firstRun")
)
#' First run
#'
#' Calls [`erosionPrerequisites`], [`erosion`], [`emission`],
#' [`transportPrerequisites`], [`transportCalcOrder`] and [`transport`] in the
#' mentioned order. While [`transport`] is called for the specified substance
#' only, [`emission`] is called for all substances whose top soil concentrations
#' have been provided.
#'
#' @inheritParams emission,RPhosFate-method
#'
#' @inherit erosionPrerequisites,RPhosFate-method return
#'
#' @seealso [`subsequentRun`]
#'
#' @examples
#' \donttest{
#' # temporary demonstration project copy
#' cv_dir <- demoProject()
#' # load temporary demonstration project
#' x <- RPhosFate(
#'   cv_dir = cv_dir,
#'   ls_ini = TRUE
#' )
#'
#' x <- firstRun(x, "SS")}
#'
#' @aliases firstRun
#'
#' @export
setMethod(
  "firstRun",
  "RPhosFate",
  function(x, substance = "PP") {
    assertChoice(substance, slotNames(x@substances))

    x <- erosionPrerequisites(x)
    x <- erosion(x)
    for (emissiveSubstance in setdiff(slotNames(x@substances), "SS")) {
      if (compareRaster(
        x@topo@rl_acc_wtd,
        slot(x@substances, emissiveSubstance)@rl_xxc,
        stopiffalse = FALSE
      )) {
        x <- emission(x, emissiveSubstance)
      }
    }
    x <- transportPrerequisites(x)
    x <- transportCalcOrder(x)
    transport(x, substance)
  }
)

#### subsequentRun ####
#' @export
setGeneric(
  "subsequentRun",
  function(x, ...) standardGeneric("subsequentRun")
)
#' Subsequent run
#'
#' Calls [`transport`] for the specified substance and optionally
#' [`erosionPrerequisites`], [`erosion`], [`emission`],
#' [`transportPrerequisites`] and/or [`transportCalcOrder`] beforehand.
#'
#' @inheritParams emission,RPhosFate-method
#' @param erosionPrerequisites A logical scalar specifying if
#'   [`erosionPrerequisites`] is called.
#' @param erosion A logical scalar specifying if [`erosion`] is called.
#' @param emission A logical scalar specifying if [`emission`] is called. It is
#'   never called with `substance = "SS"` though.
#' @param transportPrerequisites A logical scalar specifying if
#'   [`transportPrerequisites`] is called.
#' @param transportCalcOrder A logical scalar specifying if
#'   [`transportCalcOrder`] is called.
#'
#' @inherit erosionPrerequisites,RPhosFate-method return
#'
#' @seealso [`firstRun`]
#'
#' @examples
#' \donttest{
#' # temporary demonstration project copy
#' cv_dir <- demoProject()
#' # load temporary demonstration project
#' x <- RPhosFate(
#'   cv_dir = cv_dir,
#'   ls_ini = TRUE
#' )
#' # presupposed method call
#' x <- firstRun(x, "SS")
#'
#' x <- subsequentRun(x, "PP")}
#'
#' @aliases subsequentRun
#'
#' @export
setMethod(
  "subsequentRun",
  "RPhosFate",
  function(
    x,
    substance = "PP",
    erosionPrerequisites = FALSE,
    erosion = FALSE,
    emission = FALSE,
    transportPrerequisites = FALSE,
    transportCalcOrder = FALSE
  ) {
    assertChoice(substance, slotNames(x@substances))

    if (erosionPrerequisites) {
      x <- erosionPrerequisites(x)
    }
    if (erosion) {
      x <- erosion(x)
    }
    if (emission && substance != "SS") {
      x <- emission(x, substance)
    }
    if (transportPrerequisites) {
      x <- transportPrerequisites(x)
    }
    if (transportCalcOrder || length(x@helpers@order@iv_ord_row) == 0L) {
      x <- transportCalcOrder(x)
    }

    transport(x, substance)
  }
)

#### snapGauges ####
#' @export
setGeneric(
  "snapGauges",
  function(x, ...) standardGeneric("snapGauges")
)
#' Snap gauge(s)
#'
#' Snaps the coordinates of the provided calibration gauges to the respective
#' midpoint of the nearest channel cell.
#'
#' @inheritParams erosionPrerequisites,RPhosFate-method
#'
#' @inherit catchment return
#'
#' @seealso [`calibrationQuality`], [`autoCalibrate`], [`autoCalibrate2`]
#'
#' @examples
#'
#' # temporary demonstration project copy
#' cv_dir <- demoProject()
#' # load temporary demonstration project
#' x <- RPhosFate(
#'   cv_dir = cv_dir,
#'   ls_ini = TRUE
#' )
#'
#' x <- snapGauges(x)
#'
#' @aliases snapGauges
#'
#' @export
setMethod(
  "snapGauges",
  "RPhosFate",
  function(x) {
    assertdf_cdt(x)

    df_ggs <- findNearestNeighbour(
      x@parameters@df_cdt[, c("x", "y", "ID")],
      rasterToPoints(x@topo@rl_cha),
      x@helpers@ex_cmt
    )
    x@parameters@df_cdt[, c("x", "y")] <- df_ggs[, c("Y.x", "Y.y")]

    x
  }
)

#### calibrationQuality ####
#' @export
setGeneric(
  "calibrationQuality",
  function(x, ...) standardGeneric("calibrationQuality")
)
#' Calibration quality
#'
#' @description
#' Assesses the model's calibration quality via the following metrics:
#'
#' * _NSE:_ Nash-Sutcliffe Efficiency
#' * _mNSE:_ Modified Nash-Sutcliffe Efficiency (`j = 1`)
#' * _KGE:_ Kling-Gupta Efficiency (`s = c(1, 1, 1), method = "2012"`)
#' * _RMSE:_ Root Mean Square Error
#' * _PBIAS:_ Percent Bias
#' * _RSR:_ Ratio of the RMSE to the standard deviation of the observations
#' * _RCV:_ Ratio of the coefficients of variation
#' * _GMRAE:_ Geometric Mean Relative Absolute Error
#' * _MdRAE:_ Median Relative Absolute Error
#'
#' In addition, a scatter plot with the observed river loads on the x- and the
#' modelled river loads on the y-axis is displayed and provides a visual
#' impression of the model performance. Other elements of this plot are an
#' identity line (solid) and plus/minus 30% deviation lines (dashed).
#'
#' @inheritParams emission,RPhosFate-method
#' @param col A character string specifying the calibration data column with the
#'   respective substance river loads.
#'
#' @return A named numeric vector containing the assessed metrics along with the
#'   in-channel retention ratio (one minus sum of _xxt_ at catchment outlet(s)
#'   divided by sum of \emph{xxt_inp}).
#'
#' @seealso [`snapGauges`], [`autoCalibrate`], [`autoCalibrate2`],
#'   [`hydroGOF::NSE`], [`hydroGOF::mNSE`], [`hydroGOF::KGE`],
#'   [`hydroGOF::rmse`], [`hydroGOF::pbias`], [`hydroGOF::rsr`]
#'
#' @examples
#' \donttest{
#' # temporary demonstration project copy
#' cv_dir <- demoProject()
#' # load temporary demonstration project
#' x <- RPhosFate(
#'   cv_dir = cv_dir,
#'   ls_ini = TRUE
#' )
#' # presupposed method calls
#' x <- firstRun(x, "SS")
#' x <- snapGauges(x)
#'
#' calibrationQuality(x, "SS", "SS_load")}
#'
#' @aliases calibrationQuality
#'
#' @export
setMethod(
  "calibrationQuality",
  "RPhosFate",
  function(x, substance, col) {
    assertChoice(substance, slotNames(x@substances))
    assertCol(x, col)
    assertMatrix(
      x@parameters@nm_olc,
      "numeric",
      any.missing = FALSE,
      min.rows = 1L,
      ncols = 2L,
      .var.name = "nm_olc"
    )

    nv_mld <- extract(
      slot(x@substances, substance)@rl_xxt,
      as.matrix(x@parameters@df_cdt[, c("x", "y")])
    )
    nv_old <- x@parameters@df_cdt[[col]]

    assertNumeric(
      nv_mld,
      lower = 0,
      finite = TRUE,
      all.missing = FALSE,
      len = length(nv_old),
      .var.name = "Modelled load(s)"
    )

    if (substance != "SS") {
      nv_mld <- nv_mld * 1e-3
    }
    nv_rae <- abs(nv_old - nv_mld) / abs(nv_old - mean(nv_old, na.rm = TRUE))

    metrics <- c(
      tryCatch(NSE(  nv_mld, nv_old                 ), error = function(e) NA_real_),
      tryCatch(mNSE( nv_mld, nv_old                 ), error = function(e) NA_real_),
      tryCatch(KGE(  nv_mld, nv_old, method = "2012"), error = function(e) NA_real_),
      tryCatch(rmse( nv_mld, nv_old                 ), error = function(e) NA_real_),
      tryCatch(pbias(nv_mld, nv_old                 ), error = function(e) NA_real_),
      tryCatch(rsr(  nv_mld, nv_old                 ), error = function(e) NA_real_),
      (sd(nv_mld, na.rm = TRUE) / mean(nv_mld, na.rm = TRUE)) /
        (sd(nv_old, na.rm = TRUE) / mean(nv_old, na.rm = TRUE))    ,
      exp(mean(log(nv_rae), na.rm = TRUE))                         ,
      median(nv_rae, na.rm = TRUE)                                 ,
      1 - (sum(extract(
        slot(x@substances, substance)@rl_xxt,
        x@parameters@nm_olc
      )) / cellStats(slot(x@substances, substance)@rl_xxt_inp, sum))
    )
    names(metrics) <- c(x@helpers@cv_met, "inChannelRetentionRatio")

    cat("NSE:   ", metrics["NSE"  ], "\n", sep = "")
    cat("mNSE:  ", metrics["mNSE" ], "\n", sep = "")
    cat("KGE:   ", metrics["KGE"  ], "\n", sep = "")
    cat("RMSE:  ", metrics["RMSE" ], "\n", sep = "")
    cat("PBIAS: ", metrics["PBIAS"], "\n", sep = "")
    cat("RSR:   ", metrics["RSR"  ], "\n", sep = "")
    cat("RCV:   ", metrics["RCV"  ], "\n", sep = "")
    cat("GMRAE: ", metrics["GMRAE"], "\n", sep = "")
    cat("MdRAE: ", metrics["MdRAE"], "\n", sep = "")
    cat(
      "\nIn-channel retention ratio: ", metrics["inChannelRetentionRatio"],
      "\n\n", sep = ""
    )

    clippingRectangle <- c(
      0,
      max(nv_old, na.rm = TRUE),
      0,
      max(nv_mld, na.rm = TRUE)
    )

    li_par_old <- par(no.readonly = TRUE)
    on.exit(par(li_par_old))

    par(ask = FALSE)
    plot(
      NULL,
      NULL,
      xlab = sprintf("Observed load(s) in t/yr: %s", col),
      ylab = "Modelled load(s) in t/yr",
      xlim = clippingRectangle[1:2],
      ylim = clippingRectangle[3:4],
      asp = 1
    )
    clip(
      clippingRectangle[1L],
      clippingRectangle[2L],
      clippingRectangle[3L],
      clippingRectangle[4L]
    )
    abline(0, 0.7, col = "grey50", lty = "longdash")
    abline(0, 1.3, col = "grey50", lty = "longdash")
    abline(0, 1.0)
    par(xpd = NA)
    points(
      nv_old,
      nv_mld,
      pch = 21L,
      col = "black",
      bg = "#e69800",
      cex = 1.2
    )

    metrics
  }
)

#### autoCalibrate ####
#' @export
setGeneric(
  "autoCalibrate",
  function(x, ...) standardGeneric("autoCalibrate")
)
#' One dimensional automatic model calibration
#'
#' Automatically calibrates the model with the help of a combination of golden
#' section search and successive parabolic interpolation.
#'
#' @inheritParams calibrationQuality,RPhosFate-method
#' @param interval A numeric vector specifying the end-points of the interval to
#'   be searched.
#' @param metric A character string specifying the metric to optimise. See
#'   [`calibrationQuality`] for available metrics.
#' @param tol A numeric scalar specifying the desired accuracy of the parameter
#'   used for optimisation (not the metric).
#' @param parameter By default, SS are calibrated utilising the overland
#'   deposition rate and all other substances are calibrated utilising their
#'   respective enrichment ratio. This argument can be used to specify a
#'   dedicated parameter utilised for calibration via a character string:
#'   `"ns_dep_ovl"` for overland or `"ns_dep_cha"` for channel deposition rate.
#'
#' @inherit erosionPrerequisites,RPhosFate-method return
#'
#' @seealso [`snapGauges`], [`optimize`]
#'
#' @examples
#' \donttest{
#' # temporary demonstration project copy
#' cv_dir <- demoProject()
#' # load temporary demonstration project
#' x <- RPhosFate(
#'   cv_dir = cv_dir,
#'   ls_ini = TRUE
#' )
#' # presupposed method calls
#' x <- firstRun(x, "SS")
#' x <- snapGauges(x)
#'
#' x <- autoCalibrate(
#'   x,
#'   "SS",
#'   col = "SS_load",
#'   interval = c(10e-4, 20e-4),
#'   metric = "NSE"
#' )}
#'
#' @aliases autoCalibrate
#'
#' @export
setMethod(
  "autoCalibrate",
  "RPhosFate",
  function(
    x,
    substance,
    col,
    interval,
    metric,
    tol = min(interval) * 0.1,
    parameter = NULL
  ) {
    assertChoice(substance, slotNames(x@substances))
    assertCol(x, col)
    qassert(interval, "N2[0,)")
    assertChoice(metric, x@helpers@cv_met)
    qassert(tol, "N1(0,)")
    assertChoice(parameter, c("ns_dep_ovl", "ns_dep_cha"), null.ok = TRUE)

    value <- optimize(
      calibrate,
      interval,
      cmt = x,
      substance = substance,
      col = col,
      metric = metric,
      parameter = parameter,
      maximum = if (metric %in% c("NSE", "mNSE", "KGE")) TRUE else FALSE,
      tol = tol
    )

    print(value)

    if (!is.null(parameter)) {
      slot(x@parameters, parameter) <- value[[1L]]
    } else if (substance == "SS") {
      x@parameters@ns_dep_ovl <- value[[1L]]
    } else {
      x@parameters@nv_enr_rto[substance] <- value[[1L]]
    }

    if (any(abs(interval - value[[1L]]) <= tol)) {
      warning(
        "Parameter approximately within tolerance of interval end-point. ",
        "Optimum may not have been found.",
        call. = FALSE
      )
    }

    x
  }
)

#### autoCalibrate2 ####
#' @export
setGeneric(
  "autoCalibrate2",
  function(x, ...) standardGeneric("autoCalibrate2")
)
#' Two dimensional automatic model calibration
#'
#' Automatically calibrates the model with the help of a general-purpose
#' optimisation function. In contrast to [`autoCalibrate`], this method always
#' utilises the overland and channel deposition rate at the same time and never
#' the respective enrichment ratio for calibration. Beware of local optima and
#' parameters approximately within the convergence tolerance of interval
#' end-points.
#'
#' @inheritParams autoCalibrate,RPhosFate-method
#' @param method A character string specifying the utilised optimisation method.
#'   See [`optim`] for further information (use [`autoCalibrate`] instead of
#'   method `"Brent"`).
#' @param lower A numeric scalar or vector specifying the lower end-point(s) of
#'   the interval(s) to be searched.
#' @param upper A numeric scalar or vector specifying the upper end-point(s) of
#'   the interval(s) to be searched.
#' @param control A [`list`] of control parameters passed on to [`optim`]. See
#'   [`optim`] for further information.
#'
#' @inherit erosionPrerequisites,RPhosFate-method return
#'
#' @seealso [`snapGauges`]
#'
#' @examples
#' \donttest{
#' # temporary demonstration project copy
#' cv_dir <- demoProject()
#' # load temporary demonstration project
#' x <- RPhosFate(
#'   cv_dir = cv_dir,
#'   ls_ini = TRUE
#' )
#' # presupposed method calls
#' x <- firstRun(x, "SS")
#' x <- snapGauges(x)
#'
#' x <- autoCalibrate2(
#'   x,
#'   "SS",
#'   col = "SS_load",
#'   metric = "NSE",
#'   method = "L-BFGS-B",
#'   lower = c(10e-4, 0),
#'   upper = c(20e-4, 20e-4),
#'   control = list(fnscale = -1, parscale = c(1e-3, 1e-3), factr = 1e12)
#' )}
#'
#' @aliases autoCalibrate2
#'
#' @export
setMethod(
  "autoCalibrate2",
  "RPhosFate",
  function(
    x,
    substance,
    col,
    metric,
    method = "Nelder-Mead",
    lower = 0,
    upper = 0.1,
    control = list(fnscale = if (metric %in% c("NSE", "mNSE", "KGE")) -1 else 1)
  ) {
    assertChoice(substance, slotNames(x@substances))
    assertCol(x, col)
    assertChoice(metric, x@helpers@cv_met)
    qassert(method, "S1")
    assertDisjunct(method, "Brent")
    qassert(lower, "N+[0,)")
    qassert(upper, "N+[0,)")
    qassert(control, "L*")

    values <- optim(
      c(x@parameters@ns_dep_ovl, x@parameters@ns_dep_cha),
      calibrate2,
      method = method,
      cmt = x,
      substance = substance,
      col = col,
      metric = metric,
      lower = lower,
      upper = upper,
      control = control
    )

    print(values)

    x@parameters@ns_dep_ovl <- values$par[1L]
    x@parameters@ns_dep_cha <- values$par[2L]

    x
  }
)

#### saveState ####
#' @export
setGeneric(
  "saveState",
  function(x, ...) standardGeneric("saveState")
)
#' Save state
#'
#' Saves parameters _(parameters.yaml)_ and transport calculation order
#' _(order.rds)_ to disk.
#'
#' @inheritParams erosionPrerequisites,RPhosFate-method
#'
#' @return `NULL` invisibly and side effects in the form of files.
#'
#' @seealso [`RPhosFate`], [`catchment`]
#'
#' @examples
#'
#' # temporary demonstration project copy
#' cv_dir <- demoProject()
#' # load temporary demonstration project
#' x <- RPhosFate(
#'   cv_dir = cv_dir,
#'   ls_ini = TRUE
#' )
#'
#' saveState(x)
#'
#' @aliases saveState
#'
#' @export
setMethod(
  "saveState",
  "RPhosFate",
  function(x) {
    cs_dir_old <- setwd(x@cv_dir[1L])
    on.exit(setwd(cs_dir_old))

    writeParameters(x@parameters)
    saveRDS(x@helpers@order, "order.rds")
  }
)

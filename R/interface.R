#' @include aaa.R
NULL

#### Constructors ####
#' @rdname catchment
#' @export
RPhosFate <- function(...) {
  new("RPhosFate", list(...))
}
#' Initialise Project
#'
#' @description
#' Creates a project from scratch or loads an existing one utilising _Erdas
#' Imagine_ raster files (*.img) from, by convention, the following three
#' project root subdirectories:
#' * _Input_
#' * _Intermediate_
#' * _Result_
#'
#' See subdirectory sections for further information.
#'
#' @param \dots Arguments used to initialise the project. See argument sections
#'   for further information.
#'
#' @section _Input_ subdirectory:
#' This directory holds all possible user input raster data:
#' * _acc:_ Flow accumulations required for transport (the top-most cells must
#' have a value of one).
#' * \emph{acc_wtd:} Weighted flow accumulations (can be equal to _acc_).
#' * _CFa:_ (R)USLE C-factors.
#' * _cha:_ Channel cells (1: channel cell, `NA`: no channel cell).
#' * _clc:_ Clay content of top soils in \% required for substance emissions.
#' * _dem:_ Digital elevation model in m a.s.l. (optional).
#' * _dir:_ D8 flow directions required for transport.
#' * _fid:_ Field IDs (optional).
#' * _KFa:_ (R)USLE K-factors.
#' * _lue:_ Land use classes (optional).
#' * _man:_ Manning's roughness coefficients required for transport.
#' * _xxc:_ Substance content of top soils in mg/kg required for substance
#' emissions, for example, _ppc_ (PP).
#' * _rds:_ Road cells required for transport (0: road cell without subsurface
#' drainage, 1: road cell with subsurface drainage, `NA`: no road cell).
#' * _RFa:_ (R)USLE R-factors.
#' * _slp:_ Slopes in \%.
#' * _wsh:_ Watershed (optional).
#'
#' @section _Intermediate_ subdirectory:
#' This directory holds intermediate results:
#' * _inl:_ Cells representing inlets at roads (storm drains).
#' * _LFa:_ L-factors.
#' * _rhy:_ Hydraulic radii in m.
#' * _rip:_ Cells representing the riparian zones within channel cells.
#' * _SFa:_ RUSLE S-factors.
#' * \emph{slp_cap:} Capped slopes in \%.
#'
#' @section _Result_ subdirectory:
#' This directory holds the model results:
#' * _ero:_ Erosion in t/cell/yr.
#' * _xxe:_ Substance emissions in kg/cell/yr, for example, _ppe_ (PP).
#' * _xxr:_ Substance retentions in t/cell/yr (SS) or kg/cell/yr, for example,
#' _ppr_ (PP).
#' * _xxt:_ Substance transports in t/cell/yr (SS) or kg/cell/yr, for example,
#' _ppt_ (PP).
#' * \emph{xxt_cld:} Substance cell loads in t/cell/yr (SS) or kg/cell/yr, for
#' example, \emph{ppt_cld} (PP).
#' * \emph{xxt_ctf:} Substance cell transfers in t/cell/yr (SS) or kg/cell/yr,
#' for example, \emph{ppt_ctf} (PP).
#' * \emph{xxt_inp:} Substance inputs into surface waters in t/cell/yr (SS) or
#' kg/cell/yr, for example, \emph{ppt_inp} (PP).
#' * \emph{xxt_out:} Substance outlet loads in t/cell/yr (SS) or kg/cell/yr, for
#' example, \emph{ppt_out} (PP).
#'
#' @section Data management arguments:
#' * `cv_dir`: A character vector specifying the project root (first position)
#' and optionally the Monte Carlo input data directory (second position).
#' * `ls_ini`: A logical scalar specifying if an existing project shall be
#' loaded from disk (defaults to `FALSE`). Specified parameters or substance
#' parameter values via the `\dots` argument take precedence over saved ones.
#' * `is_MCi`: An integer scalar specifying the current Monte Carlo iteration if
#' applicable (defaults to `integer()`, which means Monte Carlo simulation mode
#' is disabled).
#'
#' @section Model parameter arguments:
#' * `ns_slp_min`: A numeric scalar specifying the minimum bounding slope in \%
#' (defaults to 0.001).
#' * `ns_slp_max`: A numeric scalar specifying the maximum bounding slope in \%
#' (defaults to 999.0).
#' * `ns_rhy_a`: A numeric scalar specifying a network constant depending on the
#' discharge frequency needed for the calculation of the hydraulic radius, which
#' in turn is a prerequisite for substance transport (defaults to 0.09
#' representing a discharge frequency of approximately six years).
#' * `ns_rhy_b`: A numeric scalar specifying a geometry scaling exponent
#' depending on the discharge frequency needed for the calculation of the
#' hydraulic radius, which in turn is a prerequisite for substance transport
#' (defaults to 0.50 representing a discharge frequency of approximately six
#' years).
#' * `ns_cha_rto`: A numeric scalar specifying the ratio of the channel to the
#' cell width determining the widths of the riparian zones required for
#' substance transport (defaults to 0.5).
#' * `ns_man_rip`: A numeric scalar specifying Manning's roughness coefficient
#' of the riparian zones within channel cells required for substance transport
#' (defaults to 0.32).
#' * `ns_man_cha`: A numeric scalar specifying Manning's roughness coefficient
#' of the channel within channel cells required for substance transport
#' (defaults to 0.04).
#' * `ns_dep_ovl`: A numeric scalar specifying the overland deposition
#' coefficient in \eqn{s^{-1}}{s^(-1)} required for substance transport (no
#' default).
#' * `ns_dep_cha`: A numeric scalar specifying the channel deposition
#' coefficient in \eqn{s^{-1}}{s^(-1)} required for substance transport (no
#' default). * `nv_enr_rto` A named numeric vector specifying the substance
#' enrichment ratios required for substance except SS transport, for example,
#' `c(PP = 2.0)` (no default).
#' * `nv_tfc_inl`: A named numeric vector specifying the inlet transfer
#' coefficients required for substance transport, for example, `c(SS = 0.6, PP =
#' 0.6)` (no default).
#' * `iv_fDo`: An integer vector specifying the outflow direction vector
#' required for substance transport (defaults to _ArcGIS_ codes).
#' * `nm_olc`: A numeric [`matrix`] specifying the catchment outlet coordinates
#' required for calibration (no default).
#' * `df_cdt`: A [`data.frame`] with calibration data, which must have at least
#' the following three columns and one or more columns with substance loads in
#' t/yr (no default):
#'   * _ID:_ ID(s) of the gauge(s)
#'   * _x:_ x-coordinate(s) of the gauge(s)
#'   * _y:_ y-coordinate(s) of the gauge(s)
#'
#' @return An S4 [`RPhosFate-class`] river catchment object.
#'
#' @export
catchment <- function(...) {
  new("RPhosFate", list(...))
}

#### firstRun ####
setGeneric(
  "firstRun",
  function(x, ...) standardGeneric("firstRun")
)
#' @export
setMethod(
  "firstRun",
  "RPhosFate",
  function(x, substance = "PP") {
    assertSubstance(x, substance)

    x <- erosionPrerequisites(x)
    x <- erosion(x)
    for (emmisiveSubstance in setdiff(slotNames(x@substances), "SS")) {
      if (compareRaster(
        x@topo@rl_acc_wtd,
        slot(x@substances, emmisiveSubstance)@rl_xxc,
        stopiffalse = FALSE
      )) {
        x <- emission(x, emmisiveSubstance)
      }
    }
    x <- transportPrerequisites(x)
    x <- transportCalcOrder(x)
    x <- transport(x, substance)

    x
  }
)

#### subsequentRun ####
setGeneric(
  "subsequentRun",
  function(x, ...) standardGeneric("subsequentRun")
)
#' @export
setMethod(
  "subsequentRun",
  "RPhosFate",
  function(x, substance = "PP") {
    assertSubstance(x, substance)

    if (length(x@is_MCi) == 1L) {
      x <- erosion(x)

      if (substance != "SS") {
        x <- emission(x, substance)
      }
    }
    if (length(x@helpers@order@iv_ord_row) == 0L) {
      x <- transportCalcOrder(x)
    }
    x <- transport(x, substance)

    x
  }
)

#### snapGauges ####
setGeneric(
  "snapGauges",
  function(x, ...) standardGeneric("snapGauges")
)
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
setGeneric(
  "calibrationQuality",
  function(x, ...) standardGeneric("calibrationQuality")
)
#' @export
setMethod(
  "calibrationQuality",
  "RPhosFate",
  function(x, substance = "PP", col) {
    assertSubstance(x, substance)
    assertCol(x, col)
    assertMatrix(
      x@parameters@nm_olc,
      "numeric",
      any.missing = FALSE,
      nrows = 1L,
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
      tryCatch(NSE(  nv_mld, nv_old), error = function(e) NA_real_),
      tryCatch(mNSE( nv_mld, nv_old), error = function(e) NA_real_),
      tryCatch(rmse( nv_mld, nv_old), error = function(e) NA_real_),
      tryCatch(nrmse(nv_mld, nv_old), error = function(e) NA_real_),
      tryCatch(pbias(nv_mld, nv_old), error = function(e) NA_real_),
      tryCatch(rsr(  nv_mld, nv_old), error = function(e) NA_real_),
      exp(mean(log(nv_rae), na.rm = TRUE))                         ,
      median(nv_rae, na.rm = TRUE)                                 ,
      1 - (
        extract(slot(x@substances, substance)@rl_xxt, x@parameters@nm_olc) /
          cellStats(slot(x@substances, substance)@rl_xxt_inp, sum)
      )
    )
    names(metrics) <- c(x@helpers@cv_met, "inChannelRetention")

    cat("NSE:   ", metrics["NSE"  ], "\n", sep = "")
    cat("mNSE:  ", metrics["mNSE" ], "\n", sep = "")
    cat("RMSE:  ", metrics["RMSE" ], "\n", sep = "")
    cat("NRMSE: ", metrics["NRMSE"], "\n", sep = "")
    cat("PBIAS: ", metrics["PBIAS"], "\n", sep = "")
    cat("RSR:   ", metrics["RSR"  ], "\n", sep = "")
    cat("GMRAE: ", metrics["GMRAE"], "\n", sep = "")
    cat("MdRAE: ", metrics["MdRAE"], "\n", sep = "")
    cat(
      "\nIn-channel retention: ", metrics["inChannelRetention"],
      "\n\n", sep = ""
    )

    plot(
      nv_old,
      nv_mld,
      pch = 16L,
      xlab = col,
      ylab = "Modelled load(s)",
      xlim = c(0, max(nv_old, na.rm = TRUE)),
      ylim = c(0, max(nv_mld, na.rm = TRUE))
    )
    clip(0, max(nv_old, na.rm = TRUE), 0, max(nv_mld, na.rm = TRUE))
    abline(0, 1.3, col = "grey50", lty = 2L)
    abline(0, 1.0)
    abline(0, 0.7, col = "grey50", lty = 2L)

    metrics
  }
)

#### autoCalibrate ####
setGeneric(
  "autoCalibrate",
  function(x, ...) standardGeneric("autoCalibrate")
)
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
    assertSubstance(x, substance)
    assertCol(x, col)
    qassert(interval, "N2(0,)")
    qassert(metric, "S1")
    assertSubset(metric, x@helpers@cv_met)
    qassert(tol, "N1(0,)")
    if (!is.null(parameter)) {
      qassert(parameter, "S1")
      assertSubset(parameter, c("ns_dep_ovl", "ns_dep_cha"))
    }

    value <- optimize(
      calibrate,
      interval,
      cmt = x,
      substance = substance,
      col = col,
      metric = metric,
      parameter = parameter,
      maximum = if (metric %in% c("NSE", "mNSE")) TRUE else FALSE,
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
      warning(paste(
        "Parameter approximately within tolerance of interval end-point.",
        "Optimum may not have been found."
      ), call. = FALSE)
    }

    x
  }
)

#### saveState ####
setGeneric(
  "saveState",
  function(x, ...) standardGeneric("saveState")
)
#' Save State
#'
#' Saves parameters (_parameters.yaml_) and transport calculation order
#' (_order.rds_) to disk.
#'
#' @inheritParams erosionPrerequisites,RPhosFate-method
#'
#' @return `NULL` invisibly.
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

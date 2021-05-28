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
#' * _cha:_ Channel cells (0: channel cell, `NA`: no channel cell).
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
  function(cmt, ...) standardGeneric("firstRun")
)
#' @export
setMethod(
  "firstRun",
  "RPhosFate",
  function(cmt, substance = "PP") {
    assertSubstance(cmt, substance)

    cmt <- erosionPrerequisites(cmt)
    cmt <- erosion(cmt)
    for (emmisiveSubstance in setdiff(slotNames(cmt@substances), "SS")) {
      if (compareRaster(
        cmt@topo@rl_acc_wtd,
        slot(cmt@substances, emmisiveSubstance)@rl_xxc,
        stopiffalse = FALSE
      )) {
        cmt <- emission(cmt, emmisiveSubstance)
      }
    }
    cmt <- transportPrerequisites(cmt)
    cmt <- transportCalcOrder(cmt)
    cmt <- transport(cmt, substance)

    cmt
  }
)

#### subsequentRun ####
setGeneric(
  "subsequentRun",
  function(cmt, ...) standardGeneric("subsequentRun")
)
#' @export
setMethod(
  "subsequentRun",
  "RPhosFate",
  function(cmt, substance = "PP") {
    assertSubstance(cmt, substance)

    if (length(cmt@is_MCi) == 1L) {
      cmt <- erosion(cmt)

      if (substance != "SS") {
        cmt <- emission(cmt, substance)
      }
    }
    if (length(cmt@helpers@order@iv_ord_row) == 0L) {
      cmt <- transportCalcOrder(cmt)
    }
    cmt <- transport(cmt, substance)

    cmt
  }
)

#### snapGauges ####
setGeneric(
  "snapGauges",
  function(cmt, ...) standardGeneric("snapGauges")
)
#' @export
setMethod(
  "snapGauges",
  "RPhosFate",
  function(cmt) {
    assertdf_cdt(cmt)

    df_ggs <- findNearestNeighbour(
      cmt@parameters@df_cdt[, c("x", "y", "ID")],
      rasterToPoints(cmt@topo@rl_cha),
      cmt@helpers@ex_cmt
    )
    cmt@parameters@df_cdt[, c("x", "y")] <- df_ggs[, c("Y.x", "Y.y")]

    cmt
  }
)

#### calibrationQuality ####
setGeneric(
  "calibrationQuality",
  function(cmt, ...) standardGeneric("calibrationQuality")
)
#' @export
setMethod(
  "calibrationQuality",
  "RPhosFate",
  function(cmt, substance = "PP", col) {
    assertSubstance(cmt, substance)
    assertCol(cmt, col)
    assertMatrix(
      cmt@parameters@nm_olc,
      "numeric",
      any.missing = FALSE,
      nrows = 1L,
      ncols = 2L,
      .var.name = "nm_olc"
    )

    nv_mld <- extract(
      slot(cmt@substances, substance)@rl_xxt,
      as.matrix(cmt@parameters@df_cdt[, c("x", "y")])
    )
    nv_old <- cmt@parameters@df_cdt[[col]]

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
        extract(slot(cmt@substances, substance)@rl_xxt, cmt@parameters@nm_olc) /
          cellStats(slot(cmt@substances, substance)@rl_xxt_inp, sum)
      )
    )
    names(metrics) <- c(cmt@helpers@cv_met, "inChannelRetention")

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
  function(cmt, ...) standardGeneric("autoCalibrate")
)
#' @export
setMethod(
  "autoCalibrate",
  "RPhosFate",
  function(
    cmt,
    substance,
    col,
    interval,
    metric,
    tol = min(interval) * 0.1,
    parameter = NULL
  ) {
    assertSubstance(cmt, substance)
    assertCol(cmt, col)
    qassert(interval, "N2(0,)")
    qassert(metric, "S1")
    assertSubset(metric, cmt@helpers@cv_met)
    qassert(tol, "N1(0,)")
    if (!is.null(parameter)) {
      qassert(parameter, "S1")
      assertSubset(parameter, c("ns_dep_ovl", "ns_dep_cha"))
    }

    value <- optimize(
      calibrate,
      interval,
      cmt = cmt,
      substance = substance,
      col = col,
      metric = metric,
      parameter = parameter,
      maximum = if (metric %in% c("NSE", "mNSE")) TRUE else FALSE,
      tol = tol
    )

    print(value)

    if (!is.null(parameter)) {
      slot(cmt@parameters, parameter) <- value[[1L]]
    } else if (substance == "SS") {
      cmt@parameters@ns_dep_ovl <- value[[1L]]
    } else {
      cmt@parameters@nv_enr_rto[substance] <- value[[1L]]
    }

    if (any(abs(interval - value[[1L]]) <= tol)) {
      warning(paste(
        "Parameter approximately within tolerance of interval end-point.",
        "Optimum may not have been found."
      ), call. = FALSE)
    }

    cmt
  }
)

#### saveState ####
setGeneric(
  "saveState",
  function(cmt, ...) standardGeneric("saveState")
)
#' @export
setMethod(
  "saveState",
  "RPhosFate",
  function(cmt) {
    cs_dir_old <- setwd(cmt@cv_dir[1L])
    on.exit(setwd(cs_dir_old))

    writeParameters(cmt@parameters)
    saveRDS(cmt@helpers@order, "order.rds")
  }
)

#' @include aaa.R
NULL

#### Constructors ####
#' @export
RPhosFate <- function(...) {
  new("RPhosFate", list(...))
}
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

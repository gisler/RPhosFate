#' @include aaa.R
NULL

#### Constructors ####
#' @export
RPhosFate <- function(...) {
  arguments <- list(...)
  new("RPhosFate", arguments)
}
#' @export
catchment <- function(...) {
  arguments <- list(...)
  new("RPhosFate", arguments)
}

#### firstRun ####
setGeneric(
  "firstRun",
  function(cmt, substance = "PP", ...) standardGeneric("firstRun")
)
#' @export
setMethod(
  "firstRun",
  "RPhosFate",
  function(cmt, substance) {
    cmt <- erosionPrerequisites(cmt)
    cmt <- erosion(cmt)
    if (substance != "SS") {
      cmt <- emission(cmt, substance)
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
  function(cmt, substance = "PP", ...) standardGeneric("subsequentRun")
)
#' @export
setMethod(
  "subsequentRun",
  "RPhosFate",
  function(cmt, substance) {
    if (length(cmt@is_MCi) == 1L) {
      cmt <- erosion(cmt)
    }

    if (substance != "SS") {
      cmt <- emission(cmt, substance)
    }

    if (length(cmt@helper@order@iv_ord_row) == 0L) {
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
    df_ggs <- findNearestNeighbour(
      data.frame(
        x  = cmt@parameters@df_cdt$x,
        y  = cmt@parameters@df_cdt$y,
        ID = cmt@parameters@df_cdt$ID
      ),
      rasterToPoints(cmt@topo@rl_cha),
      cmt@helper@ex_cmt
    )

    cmt@parameters@df_cdt$x <- df_ggs$Y.x
    cmt@parameters@df_cdt$y <- df_ggs$Y.y

    cmt
  }
)

#### calibrationQuality ####
setGeneric(
  "calibrationQuality",
  function(cmt, substance = "PP", ...) standardGeneric("calibrationQuality")
)
#' @export
setMethod(
  "calibrationQuality",
  "RPhosFate",
  function(cmt, substance, col) {
    nv_mld <- extract(
      slot(cmt@substance, substance)@rl_xxt,
      as.matrix(cmt@parameters@df_cdt[, c("x", "y")])
    )
    if (substance != "SS") {
      nv_mld <- nv_mld / 1000
    }
    nv_old <- cmt@parameters@df_cdt[[col]]
    nv_rae <- abs(nv_old - nv_mld) / abs(nv_old - mean(nv_old, na.rm = TRUE))

    if (length(nv_old) > 1L) {
      metrics <- c(
        NSE( nv_mld, nv_old),
        mNSE(nv_mld, nv_old),
        rsr( nv_mld, nv_old)
      )
    } else {
      metrics <- c(rep(NA_real_, 3L))
    }
    names(metrics) <- c("NSE", "mNSE", "RSR")
    metrics <- c(
      metrics,
      PBIAS = pbias(nv_mld, nv_old),
      GMRAE = exp(mean(log(nv_rae), na.rm = TRUE)),
      MdRAE = median(nv_rae, na.rm = TRUE),
      inChannelRetention = 1 - (
        extract(slot(cmt@substance, substance)@rl_xxt, cmt@parameters@nm_olc) /
          cellStats(slot(cmt@substance, substance)@rl_xxt_inp, sum)
      )
    )

    cat("NSE:   ", metrics["NSE"  ], "\n", sep = "")
    cat("mNSE:  ", metrics["mNSE" ], "\n", sep = "")
    cat("RSR:   ", metrics["RSR"  ], "\n", sep = "")
    cat("PBIAS: ", metrics["PBIAS"], "\n", sep = "")
    cat("GMRAE: ", metrics["GMRAE"], "\n", sep = "")
    cat("MdRAE: ", metrics["MdRAE"], "\n", sep = "")
    cat("\nIn-channel retention: ", metrics["inChannelRetention"], "\n", sep = "")

    plot(
      nv_old,
      nv_mld,
      pch = 16L,
      xlim = c(0, max(nv_old, na.rm = TRUE)),
      ylim = c(0, max(nv_mld, na.rm = TRUE))
    )
    clip(0, max(nv_old, na.rm = TRUE), 0, max(nv_mld, na.rm = TRUE))
    abline(0, 1.3, lty = 2L)
    abline(0, 1.0)
    abline(0, 0.7, lty = 2L)

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
    value <- optimize(
      calibrate,
      interval,
      cmt = cmt,
      substance = substance,
      col = col,
      metric = metric,
      parameter = parameter,
      maximum = if (metric %in% c("NSE", "mNSE")) {TRUE} else {FALSE},
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
    saveRDS(cmt@helper@order, "order.rds")
  }
)

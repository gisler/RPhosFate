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
    if (!requireNamespace("hydroGOF", quietly = TRUE)) {
      stop(
        "Package \"hydroGOF\" must be installed for this functionality.",
        call. = FALSE
      )
    }

    nv_mld <- extract(
      slot(cmt@substance, substance)@rl_xxt,
      cbind(cmt@parameters@df_cdt$x, cmt@parameters@df_cdt$y)
    )
    if (substance != "SS") {
      nv_mld <- nv_mld / 1000
    }
    nv_old <- cmt@parameters@df_cdt[[col]]
    nv_rae <- abs(nv_old - nv_mld) / abs(nv_old - mean(nv_old, na.rm = TRUE))

    if (length(nv_old) > 1L) {
      cat("NSE:   ", hydroGOF::NSE.default(  nv_mld, nv_old), "\n", sep = "")
      cat("mNSE:  ", hydroGOF::mNSE.default( nv_mld, nv_old), "\n", sep = "")
      cat("RSR:   ", hydroGOF::rsr.default(  nv_mld, nv_old), "\n", sep = "")
    }
    {
      cat("PBIAS: ", hydroGOF::pbias.default(nv_mld, nv_old), "\n", sep = "")
      cat("GMRAE: ", exp(mean(log(nv_rae), na.rm = TRUE)),    "\n", sep = "")
      cat("MdRAE: ", median(nv_rae, na.rm = TRUE),            "\n", sep = "")
      cat(
        "\nIn-channel retention: ",
        1 - (extract(slot(cmt@substance, substance)@rl_xxt, cmt@parameters@nm_olc) /
          cellStats(slot(cmt@substance, substance)@rl_xxt_inp, sum)),
        "\n",
        sep = ""
      )
    }

    plot(
      nv_old, nv_mld,
      pch = 16L,
      xlim = c(0, max(nv_old, na.rm = TRUE)),
      ylim = c(0, max(nv_mld, na.rm = TRUE))
    )
    graphics::clip(0, max(nv_old, na.rm = TRUE), 0, max(nv_mld, na.rm = TRUE))
    graphics::abline(0, 1.3, lty = 2L)
    graphics::abline(0, 1.0)
    graphics::abline(0, 0.7, lty = 2L)
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

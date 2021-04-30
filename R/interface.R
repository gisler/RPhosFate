#' @include aaa.R
NULL

#### Constructors ####
#' @export
RPhosFate <- function(...) {
  arguments <- list(...)
  return(new("RPhosFate", arguments))
}
#' @export
catchment <- function(...) {
  arguments <- list(...)
  return(new("RPhosFate", arguments))
}

#### firstRun ####
setGeneric(
  name = "firstRun",
  def = function(cmt) {standardGeneric("firstRun")}
)
#' @export
setMethod(
  f = "firstRun",
  signature = "RPhosFate",
  definition = function(cmt) {
    cmt <- erosionPrerequisites(cmt)
    cmt <- erosion(cmt)
    cmt <- emission(cmt, cmt@PP)
    cmt <- transportPrerequisites(cmt)
    cmt <- transportCalcOrder(cmt)
    cmt <- transport(cmt, cmt@PP)

    return(cmt)
  }
)

#### subsequentRun ####
setGeneric(
  name = "subsequentRun",
  def = function(cmt) {standardGeneric("subsequentRun")}
)
#' @export
setMethod(
  f = "subsequentRun",
  signature = "RPhosFate",
  definition = function(cmt) {
    if (length(cmt@is_MCi) == 1L) {
      cmt <- erosion(cmt)
      cmt <- emission(cmt, cmt@PP)
    }

    if (length(cmt@helper@order@iv_ord_row) == 0L) {
      cmt <- transportCalcOrder(cmt)
    }

    cmt <- transport(cmt, cmt@PP)

    return(cmt)
  }
)

#### snapGauges ####
setGeneric(
  name = "snapGauges",
  def = function(cmt) {standardGeneric("snapGauges")}
)
#' @export
setMethod(
  f = "snapGauges",
  signature = "RPhosFate",
  definition = function(cmt) {
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

    return(cmt)
  }
)

#### calibrationQuality ####
setGeneric(
  name = "calibrationQuality",
  def = function(cmt, col) {standardGeneric("calibrationQuality")}
)
#' @export
setMethod(
  f = "calibrationQuality",
  signature = c("RPhosFate"),
  definition = function(cmt, col) {
    if (!requireNamespace("hydroGOF", quietly = TRUE)) {
      stop("Package \"hydroGOF\" must be installed for this functionality.", call. = FALSE)
    }

    nv_mld <- extract(cmt@PP@rl_ppt, cbind(cmt@parameters@df_cdt$x, cmt@parameters@df_cdt$y)) / 1000
    nv_old <- cmt@parameters@df_cdt[[col]]
    nv_rae <- abs(nv_old - nv_mld) / abs(nv_old - mean(nv_old, na.rm = TRUE))

    if (length(nv_old) > 1) {
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
        1 - (extract(cmt@PP@rl_ppt, cmt@parameters@nm_olc) / cellStats(cmt@PP@rl_ppt_inp, sum)),
        "\n",
        sep = ""
      )
    }

    plot(
      nv_old, nv_mld,
      pch = 16,
      xlim = c(0, max(nv_old, na.rm = TRUE)),
      ylim = c(0, max(nv_mld, na.rm = TRUE))
    )
    graphics::clip(0, max(nv_old, na.rm = TRUE), 0, max(nv_mld, na.rm = TRUE))
    graphics::abline(0, 1.3, lty = 2)
    graphics::abline(0, 1)
    graphics::abline(0, 0.7, lty = 2)
  }
)

#### saveState ####
setGeneric(
  name = "saveState",
  def = function(cmt) {standardGeneric("saveState")}
)
#' @export
setMethod(
  f = "saveState",
  signature = "RPhosFate",
  definition = function(cmt) {
    cs_dir_old <- setwd(cmt@cv_dir[1])
    on.exit(setwd(cs_dir_old))

    writeParameters(cmt@parameters)
    saveRDS(cmt@helper@order, "order.rds")
  }
)

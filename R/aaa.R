#' @import methods
#' @import raster
#' @importFrom Rcpp sourceCpp
#' @importFrom stats median
#' @useDynLib RPhosFate
NULL

#### Class RPhosFateParameters ####
setClassRPhosFateParameters <- function(env) {
  setClass(
    Class = "RPhosFateParameters",
    slots = c(
      ns_slp_min = "numeric",
      ns_slp_max = "numeric",
      ns_rhy_a   = "numeric",
      ns_rhy_b   = "numeric",
      ns_cha_rto = "numeric",
      ns_man_rip = "numeric",
      ns_man_cha = "numeric",
      ns_dep_ovl = "numeric",
      ns_dep_cha = "numeric",
      ns_tfc_inl = "numeric",
      nv_enr_rto = "numeric",
      iv_fDo     = "integer",
      nm_olc     = "matrix",
      df_cdt     = "data.frame"
    ),
    where = env
  )
}
setClassRPhosFateParameters(environment())
setMethod(
  f = "initialize",
  signature = "RPhosFateParameters",
  definition = function(.Object, arguments) {
    # Min slope cap in %
    if (!is.null(arguments$ns_slp_min       )) {.Object@ns_slp_min <- arguments$ns_slp_min} else {.Object@ns_slp_min <- 0.001}

    # Max slope cap in %
    if (!is.null(arguments$ns_slp_max       )) {.Object@ns_slp_max <- arguments$ns_slp_max} else {.Object@ns_slp_max <- 999}

    # Parameter related to the discharge recurrence interval (WetSpa, T = 6)
    if (!is.null(arguments$ns_rhy_a         )) {.Object@ns_rhy_a   <- arguments$ns_rhy_a  } else {.Object@ns_rhy_a   <- 0.09}

    # Parameter related to the discharge recurrence interval (WetSpa, T = 6)
    if (!is.null(arguments$ns_rhy_b         )) {.Object@ns_rhy_b   <- arguments$ns_rhy_b  } else {.Object@ns_rhy_b   <- 0.5}

    # Ratio of channel width to cell length determining the riparian zone
    if (!is.null(arguments$ns_cha_rto       )) {.Object@ns_cha_rto <- arguments$ns_cha_rto} else {.Object@ns_cha_rto <- 0.5}

    # Riparian zone manning n
    if (!is.null(arguments$ns_man_rip       )) {.Object@ns_man_rip <- arguments$ns_man_rip} else {.Object@ns_man_rip <- 0.32}

    # Channel manning n
    if (!is.null(arguments$ns_man_cha       )) {.Object@ns_man_cha <- arguments$ns_man_cha} else {.Object@ns_man_cha <- 0.04}

    # Overland deposition coefficient
    if (!is.null(arguments$ns_dep_ovl       )) {.Object@ns_dep_ovl <- arguments$ns_dep_ovl} else {stop("\"ns_dep_ovl\" must be supplied.")}

    # Channel deposition coefficient
    if (!is.null(arguments$ns_dep_cha       )) {.Object@ns_dep_cha <- arguments$ns_dep_cha} else {stop("\"ns_dep_cha\" must be supplied.")}

    # Inlet transfer coefficient
    if (!is.null(arguments$ns_tfc_inl       )) {.Object@ns_tfc_inl <- arguments$ns_tfc_inl} else {stop("\"ns_tfc_inl\" must be supplied.")}

    # Enrichment ratios
    if (!is.null(names(arguments$nv_enr_rto))) {.Object@nv_enr_rto <- arguments$nv_enr_rto} else {stop("\"nv_enr_rto\" must be supplied as named vector.")}

    # Outflow direction vector (ArcGIS coded)
    if (!is.null(arguments$iv_fDo           )) {.Object@iv_fDo     <- arguments$iv_fDo    } else {.Object@iv_fDo     <- as.integer(c(32, 16, 8, 64, 0, 4, 128, 1, 2))}

    # Catchment outlet coordinates
    if (!is.null(arguments$nm_olc           )) {.Object@nm_olc     <- arguments$nm_olc    } else {.Object@nm_olc     <- matrix()}

    # Calibration data
    if (!is.null(arguments$df_cdt           )) {.Object@df_cdt     <- arguments$df_cdt    } else {.Object@df_cdt     <- data.frame()}

    return(.Object)
  }
)

#### Class RPhosFateTopo ####
setClass(
  Class = "RPhosFateTopo",
  slots = c(
    rl_acc     = "RasterLayer", # Flow accumulation for transport calculation order
    rl_acc_wtd = "RasterLayer", # Weighted flow accumulation for erosion and transport
    rl_cha     = "RasterLayer", # Channel cells
    rl_clc     = "RasterLayer", # Clay content of topsoil in %
    rl_dem     = "RasterLayer", # DEM
    rl_dir     = "RasterLayer", # Flow direction
    rl_fid     = "RasterLayer", # Field plot IDs
    rl_inl     = "RasterLayer", # Inlet cells
    rl_lue     = "RasterLayer", # Land use
    rl_rds     = "RasterLayer", # Road cells
    rl_rip     = "RasterLayer", # Riparian zone cells
    rl_slp     = "RasterLayer", # Slope in %
    rl_slp_cap = "RasterLayer", # Capped slope in %
    rl_wsh     = "RasterLayer"  # Watershed
  )
)
setMethod(
  f = "initialize",
  signature = "RPhosFateTopo",
  definition = function(.Object, cmt) {
    cs_dir_old <- setwd(cmt@cv_dir[1])
    on.exit(setwd(cs_dir_old))

    setwd("Input")
    .Object@rl_acc     <- readLayer(cmt, "acc"          )
    .Object@rl_acc_wtd <- readLayer(cmt, "acc_wtd", TRUE)
    .Object@rl_cha     <- readLayer(cmt, "cha"    , TRUE)
    .Object@rl_clc     <- readLayer(cmt, "clc"          )
    .Object@rl_dem     <- readLayer(cmt, "dem"          )
    .Object@rl_dir     <- readLayer(cmt, "dir"          )
    .Object@rl_fid     <- readLayer(cmt, "fid"          )
    .Object@rl_lue     <- readLayer(cmt, "lue"          )
    .Object@rl_rds     <- readLayer(cmt, "rds"          )
    .Object@rl_slp     <- readLayer(cmt, "slp"    , TRUE)
    .Object@rl_wsh     <- readLayer(cmt, "wsh"          )

    setwd("../Intermediate")
    .Object@rl_inl     <- readLayer(cmt, "inl"    )
    .Object@rl_rip     <- readLayer(cmt, "rip"    )
    .Object@rl_slp_cap <- readLayer(cmt, "slp_cap")

    return(.Object)
  }
)

#### Class RPhosFateErosion ####
setClass(
  Class = "RPhosFateErosion",
  slots = c(
    rl_RFa = "RasterLayer", # R-factor
    rl_KFa = "RasterLayer", # K-factor
    rl_LFa = "RasterLayer", # L-factor
    rl_SFa = "RasterLayer", # S-factor
    rl_CFa = "RasterLayer", # C-factor
    rl_ero = "RasterLayer"  # Erosion in t/cell/yr
  )
)
setMethod(
  f = "initialize",
  signature = "RPhosFateErosion",
  definition = function(.Object, cmt) {
    cs_dir_old <- setwd(cmt@cv_dir[1])
    on.exit(setwd(cs_dir_old))

    setwd("Input")
    .Object@rl_RFa <- readLayer(cmt, "RFa", TRUE      )
    .Object@rl_KFa <- readLayer(cmt, "KFa", TRUE      )
    .Object@rl_CFa <- readLayer(cmt, "CFa", TRUE, TRUE)

    setwd("../Intermediate")
    .Object@rl_LFa <- readLayer(cmt, "LFa")
    .Object@rl_SFa <- readLayer(cmt, "SFa")

    setwd("../Result")
    .Object@rl_ero <- readLayer(cmt, "ero")

    return(.Object)
  }
)

#### Class RPhosFateTransport ####
setClass(
  Class = "RPhosFateTransport",
  slots = c(
    rl_man = "RasterLayer", # Manning n
    rl_rhy = "RasterLayer"  # Hydraulic radius in m
  )
)
setMethod(
  f = "initialize",
  signature = "RPhosFateTransport",
  definition = function(.Object, cmt) {
    cs_dir_old <- setwd(cmt@cv_dir[1])
    on.exit(setwd(cs_dir_old))

    setwd("Input")
    .Object@rl_man <- readLayer(cmt, "man", FALSE, TRUE)

    setwd("../Intermediate")
    .Object@rl_rhy <- readLayer(cmt, "rhy")

    return(.Object)
  }
)

#### Class RPhosFateSS ####
setClass(
  Class = "RPhosFateSS",
  slots = c(
    rl_ssr     = "RasterLayer", # SS retention     in kg/cell/yr
    rl_sst     = "RasterLayer", # SS transport     in kg/cell/yr
    rl_sst_inp = "RasterLayer", # SS input load    in kg/cell/yr
    rl_sst_out = "RasterLayer", # SS outlet load   in kg/cell/yr
    rl_sst_cld = "RasterLayer", # SS cell load     in kg/cell/yr
    rl_sst_ctf = "RasterLayer"  # SS cell transfer in kg/cell/yr
  )
)
setMethod(
  f = "initialize",
  signature = "RPhosFateSS",
  definition = function(.Object, cmt) {
    cs_dir_old <- setwd(cmt@cv_dir[1])
    on.exit(setwd(cs_dir_old))

    setwd("Result")
    .Object@rl_ssr     <- readLayer(cmt, "ssr"    )
    .Object@rl_sst     <- readLayer(cmt, "sst"    )
    .Object@rl_sst_inp <- readLayer(cmt, "sst_inp")
    .Object@rl_sst_out <- readLayer(cmt, "sst_out")
    .Object@rl_sst_cld <- readLayer(cmt, "sst_cld")
    .Object@rl_sst_ctf <- readLayer(cmt, "sst_ctf")

    return(.Object)
  }
)

#### Class RPhosFatePP ####
setClass(
  Class = "RPhosFatePP",
  slots = c(
    rl_ppc     = "RasterLayer", # PP content of topsoil in mg P/kg
    rl_ppe     = "RasterLayer", # PP emission      in kg P/cell/yr
    rl_ppr     = "RasterLayer", # PP retention     in kg P/cell/yr
    rl_ppt     = "RasterLayer", # PP transport     in kg P/cell/yr
    rl_ppt_inp = "RasterLayer", # PP input load    in kg P/cell/yr
    rl_ppt_out = "RasterLayer", # PP outlet load   in kg P/cell/yr
    rl_ppt_cld = "RasterLayer", # PP cell load     in kg P/cell/yr
    rl_ppt_ctf = "RasterLayer"  # PP cell transfer in kg P/cell/yr
  )
)
setMethod(
  f = "initialize",
  signature = "RPhosFatePP",
  definition = function(.Object, cmt) {
    cs_dir_old <- setwd(cmt@cv_dir[1])
    on.exit(setwd(cs_dir_old))

    setwd("Input")
    .Object@rl_ppc     <- readLayer(cmt, "ppc")

    setwd("../Result")
    .Object@rl_ppe     <- readLayer(cmt, "ppe"    )
    .Object@rl_ppr     <- readLayer(cmt, "ppr"    )
    .Object@rl_ppt     <- readLayer(cmt, "ppt"    )
    .Object@rl_ppt_inp <- readLayer(cmt, "ppt_inp")
    .Object@rl_ppt_out <- readLayer(cmt, "ppt_out")
    .Object@rl_ppt_cld <- readLayer(cmt, "ppt_cld")
    .Object@rl_ppt_ctf <- readLayer(cmt, "ppt_ctf")

    return(.Object)
  }
)

#### Class RPhosFateOrder ####
setClass(
  Class = "RPhosFateOrder",
  slots = c(
    iv_ord_row = "integer",
    iv_ord_col = "integer",
    iv_ord_ovl_row_rev = "integer",
    iv_ord_ovl_col_rev = "integer"
  )
)

#### Class RPhosFateHelper ####
setClass(
  Class = "RPhosFateHelper",
  slots = c(
    ex_cmt     = "Extent",
    is_res     = "integer",
    is_siz     = "integer",
    is_rws     = "integer",
    is_cls     = "integer",
    iv_fDo_dgl = "integer",
    im_fDo     = "matrix",
    im_fDi     = "matrix",
    order      = "RPhosFateOrder"
  )
)
setMethod(
  f = "initialize",
  signature = "RPhosFateHelper",
  definition = function(.Object, cmt) {
    cs_dir_old <- setwd(cmt@cv_dir[1])
    on.exit(setwd(cs_dir_old))

    # Extent of catchment area
    .Object@ex_cmt <- extent(cmt@topo@rl_acc_wtd)

    # Cell length in m
    .Object@is_res <- as.integer(xres(cmt@topo@rl_acc_wtd))

    # Cell area in m^2
    .Object@is_siz <- as.integer(.Object@is_res^2)

    # Number of rows
    .Object@is_rws <- nrow(cmt@topo@rl_acc_wtd)

    # Number of columns
    .Object@is_cls <- ncol(cmt@topo@rl_acc_wtd)

    # Diagonal outflow direction vector
    .Object@iv_fDo_dgl <- cmt@parameters@iv_fDo[c(1, 3, 7, 9)]

    # Outflow direction matrix
    .Object@im_fDo <- matrix(cmt@parameters@iv_fDo, nrow = 3, ncol = 3)

    # Inflow direction matrix
    .Object@im_fDi <- matrix(rev(cmt@parameters@iv_fDo), nrow = 3, ncol = 3)

    # Transport calculation order
    if (cmt@ls_ini && file.exists("order.rds")) {
      .Object@order <- readRDS("order.rds")
    }

    return(.Object)
  }
)

#### Class RPhosFate ####
setClass(
  Class = "RPhosFate",
  slots = c(
    cv_dir     = "character",
    ls_ini     = "logical",
    is_MCi     = "integer",
    parameters = "RPhosFateParameters",
    topo       = "RPhosFateTopo",
    erosion    = "RPhosFateErosion",
    transport  = "RPhosFateTransport",
    SS         = "RPhosFateSS",
    PP         = "RPhosFatePP",
    helper     = "RPhosFateHelper"
  )
)
setMethod(
  f = "initialize",
  signature = "RPhosFate",
  definition = function(.Object, arguments) {
    # Project directory
    if (!is.null(arguments$cv_dir)) .Object@cv_dir <- arguments$cv_dir else stop("\"cv_dir\" must be supplied.")

    # Load parameters from disc?
    if (!is.null(arguments$ls_ini)) .Object@ls_ini <- arguments$ls_ini else .Object@ls_ini <- FALSE

    # Monte Carlo iteration
    if (!is.null(arguments$is_MCi)) .Object@is_MCi <- arguments$is_MCi

    cs_dir_old <- setwd(.Object@cv_dir[1])
    on.exit(setwd(cs_dir_old))

    if (!dir.exists("Intermediate") || !dir.exists("Result")) {
      dir.create("Intermediate", showWarnings = FALSE)
      dir.create("Result", showWarnings = FALSE)
    }

    if (.Object@ls_ini && file.exists("parameters.rds")) {
      .Object@parameters <- tryCatch(
        readParameters(),
        error = updateRPhosFateParameters
      )
    } else {
      .Object@parameters <- new("RPhosFateParameters", arguments)
    }
    .Object@topo      <- new("RPhosFateTopo", .Object)
    .Object@erosion   <- new("RPhosFateErosion", .Object)
    .Object@transport <- new("RPhosFateTransport", .Object)
    .Object@SS        <- new("RPhosFateSS", .Object)
    .Object@PP        <- new("RPhosFatePP", .Object)
    .Object@helper    <- new("RPhosFateHelper", .Object)

    return(.Object)
  }
)

#### readLayer ####
setGeneric(
  name = "readLayer",
  def = function(cmt, layer, isRequiredInputLayer = FALSE, isMCinputLayer = FALSE) {standardGeneric("readLayer")}
)
setMethod(
  f = "readLayer",
  signature = c("RPhosFate", "character"),
  definition = function(cmt, layer, isRequiredInputLayer, isMCinputLayer) {
    if (length(cmt@is_MCi) == 1L && isMCinputLayer) {
      cs_dir_old <- setwd(cmt@cv_dir[2])
      on.exit(setwd(cs_dir_old))
    }

    filename <- paste0(layer, cmt@is_MCi, ".img")
    if (isRequiredInputLayer) {
      rasterLayer <- raster(filename)
    } else if (file.exists(filename)) {
      rasterLayer <- raster(filename)
    } else {
      rasterLayer <- new("RasterLayer")
    }

    return(rasterLayer)
  }
)

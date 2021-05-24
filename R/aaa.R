#' @import checkmate
#' @import methods
#' @import raster
#' @importFrom graphics abline clip
#' @importFrom hydroGOF mNSE NSE pbias rmse nrmse rsr
#' @importFrom Rcpp sourceCpp
#' @importFrom spatstat as.owin as.ppp nncross
#' @importFrom stats median optimize setNames
#' @importFrom yaml read_yaml write_yaml
#' @importFrom utils modifyList packageVersion
#' @useDynLib RPhosFate
NULL

#### Class RPhosFateParameters2 ####
setClass(
  "RPhosFateParameters2",
  slots = c(
    ns_slp_min = "numeric",   # Min slope cap in %
    ns_slp_max = "numeric",   # Max slope cap in %
    ns_rhy_a   = "numeric",   # Parameter related to the discharge recurrence interval (WetSpa, T = 6)
    ns_rhy_b   = "numeric",   # Parameter related to the discharge recurrence interval (WetSpa, T = 6)
    ns_cha_rto = "numeric",   # Ratio of channel width to cell length determining the riparian zone
    ns_man_rip = "numeric",   # Riparian zone manning n
    ns_man_cha = "numeric",   # Channel manning n
    ns_dep_ovl = "numeric",   # Overland deposition coefficient
    ns_dep_cha = "numeric",   # Channel deposition coefficient
    nv_enr_rto = "numeric",   # Enrichment ratios
    nv_tfc_inl = "numeric",   # Inlet transfer coefficients
    iv_fDo     = "integer",   # Outflow direction vector (ArcGIS coded)
    nm_olc     = "matrix",    # Catchment outlet coordinates
    df_cdt     = "data.frame" # Calibration data
  ),
  prototype = list(
    ns_slp_min = 0.001,
    ns_slp_max = 999.0,
    ns_rhy_a   = 0.09,
    ns_rhy_b   = 0.50,
    ns_cha_rto = 0.5,
    ns_man_rip = 0.32,
    ns_man_cha = 0.04,
    nv_enr_rto = numeric(),
    nv_tfc_inl = numeric(),
    iv_fDo     = c(32L, 16L, 8L, 64L, 0L, 4L, 128L, 1L, 2L),
    nm_olc     = matrix(NA_real_),
    df_cdt     = data.frame()
  )
)
setMethod(
  "initialize",
  "RPhosFateParameters2",
  function(.Object, arguments) {
    populateParameterSlots(.Object, arguments)
  }
)
setValidity(
  "RPhosFateParameters2",
  function(object) {
    qassert(object@ns_slp_min, "N1[0,)", .var.name = "ns_slp_min")
    assertNumber(
      object@ns_slp_max,
      lower = object@ns_slp_min,
      finite = TRUE,
      .var.name = "ns_slp_max"
    )
    qassert(object@ns_rhy_a  , "N1(0,)", .var.name = "ns_rhy_a"  )
    qassert(object@ns_rhy_b  , "N1[0,)", .var.name = "ns_rhy_b"  )
    qassert(object@ns_cha_rto, "N1(0,)", .var.name = "ns_cha_rto")
    qassert(object@ns_man_rip, "N1(0,)", .var.name = "ns_man_rip")
    qassert(object@ns_man_cha, "N1(0,)", .var.name = "ns_man_cha")
    qassert(object@ns_dep_ovl, "N1(0,)", .var.name = "ns_dep_ovl")
    qassert(object@ns_dep_cha, "N1[0,)", .var.name = "ns_dep_cha")
    qassert(object@iv_fDo    , "I9[0,)", .var.name = "iv_fDo"    )

    TRUE
  }
)

#### Class RPhosFateTopo ####
setClass(
  "RPhosFateTopo",
  slots = c(
    rl_acc     = "RasterLayer", # Flow accumulation for transport calculation order
    rl_acc_wtd = "RasterLayer", # Weighted flow accumulation for erosion and transport
    rl_cha     = "RasterLayer", # Channel cells
    rl_clc     = "RasterLayer", # Clay content of topsoil in %
    rl_dem     = "RasterLayer", # Digital elevation model
    rl_dir     = "RasterLayer", # Flow direction
    rl_fid     = "RasterLayer", # Field IDs
    rl_inl     = "RasterLayer", # Inlet cells at roads
    rl_lue     = "RasterLayer", # Land use
    rl_rds     = "RasterLayer", # Road cells
    rl_rip     = "RasterLayer", # Riparian zone cells at channels
    rl_slp     = "RasterLayer", # Slope in %
    rl_slp_cap = "RasterLayer", # Capped slope in %
    rl_wsh     = "RasterLayer"  # Extent of watershed
  )
)
setMethod(
  "initialize",
  "RPhosFateTopo",
  function(.Object, cmt) {
    cs_dir_old <- setwd(cmt@cv_dir[1L])
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

    .Object
  }
)

#### Class RPhosFateErosion ####
setClass(
  "RPhosFateErosion",
  slots = c(
    rl_RFa = "RasterLayer", # R factor
    rl_KFa = "RasterLayer", # K factor
    rl_LFa = "RasterLayer", # L factor
    rl_SFa = "RasterLayer", # S factor
    rl_CFa = "RasterLayer", # C factor
    rl_ero = "RasterLayer"  # Erosion in t/cell/yr
  )
)
setMethod(
  "initialize",
  "RPhosFateErosion",
  function(.Object, cmt) {
    cs_dir_old <- setwd(cmt@cv_dir[1L])
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

    .Object
  }
)

#### Class RPhosFateTransport ####
setClass(
  "RPhosFateTransport",
  slots = c(
    rl_man = "RasterLayer", # Manning n
    rl_rhy = "RasterLayer"  # Hydraulic radius in m
  )
)
setMethod(
  "initialize",
  "RPhosFateTransport",
  function(.Object, cmt) {
    cs_dir_old <- setwd(cmt@cv_dir[1L])
    on.exit(setwd(cs_dir_old))

    setwd("Input")
    .Object@rl_man <- readLayer(cmt, "man", FALSE, TRUE)

    setwd("../Intermediate")
    .Object@rl_rhy <- readLayer(cmt, "rhy")

    .Object
  }
)

#### Class RPhosFateBare ####
setClass(
  "RPhosFateBare",
  slots = c(
    rl_xxr     = "RasterLayer", # Substance retention                  in t (SS) or kg/cell/yr
    rl_xxt     = "RasterLayer", # Substance transport                  in t (SS) or kg/cell/yr
    rl_xxt_inp = "RasterLayer", # Substance inputs into surface waters in t (SS) or kg/cell/yr
    rl_xxt_out = "RasterLayer", # Substance outlet loads               in t (SS) or kg/cell/yr
    rl_xxt_cld = "RasterLayer", # Substance cell loads                 in t (SS) or kg/cell/yr
    rl_xxt_ctf = "RasterLayer"  # Substance cell transfers             in t (SS) or kg/cell/yr
  ),
  contains = "VIRTUAL"
)
#### Class RPhosFateConc ####
setClass(
  "RPhosFateConc",
  slots = c(
    rl_xxc = "RasterLayer", # Substance content of topsoil in mg/kg
    rl_xxe = "RasterLayer"  # Substance emission           in kg/cell/yr
  ),
  contains = c("VIRTUAL", "RPhosFateBare")
)

#### Class RPhosFateSS ####
setClass(
  "RPhosFateSS",
  contains = "RPhosFateBare"
)
setMethod(
  "initialize",
  "RPhosFateSS",
  function(.Object, cmt) {
    slots <- slotNames(.Object)
    layers <- file.path(
      rep("Result", length(slots)),
      sub("^rl_xx", "ss", slots)
    )

    populateLayerSlots(cmt, .Object, slots, layers)
  }
)

#### Class RPhosFatePP ####
setClass(
  "RPhosFatePP",
  contains = "RPhosFateConc"
)
setMethod(
  "initialize",
  "RPhosFatePP",
  function(.Object, cmt) {
    slots <- slotNames(.Object)
    layers <- file.path(
      c("Input", rep("Result", length(slots) - 1L)),
      sub("^rl_xx", "pp", slots)
    )

    populateLayerSlots(cmt, .Object, slots, layers)
  }
)

#### Class RPhosFateSubstance ####
setClass(
  "RPhosFateSubstance",
  slots = c(
    SS = "RPhosFateSS",
    PP = "RPhosFatePP"
  )
)
setMethod(
  "initialize",
  "RPhosFateSubstance",
  function(.Object, cmt) {
    .Object@SS <- new("RPhosFateSS", cmt)
    .Object@PP <- new("RPhosFatePP", cmt)

    .Object
  }
)

#### Class RPhosFateOrder ####
setClass(
  "RPhosFateOrder",
  slots = c(
    iv_ord_row = "integer",
    iv_ord_col = "integer",
    iv_ord_ovl_row_rev = "integer",
    iv_ord_ovl_col_rev = "integer"
  )
)

#### Class RPhosFateHelper ####
setClass(
  "RPhosFateHelper",
  slots = c(
    ex_cmt     = "Extent",        # Extent of catchment area
    is_res     = "integer",       # Cell length in m
    is_siz     = "integer",       # Cell area in m^2
    is_rws     = "integer",       # Number of rows
    is_cls     = "integer",       # Number of columns
    iv_fDo_dgl = "integer",       # Diagonal outflow direction vector
    im_fDo     = "matrix",        # Outflow direction matrix
    im_fDi     = "matrix",        # Inflow direction matrix
    cv_rls     = "character",     # Objects holding raster layers
    cv_met     = "character",     # Implemented calibration quality metrics
    order      = "RPhosFateOrder" # Transport calculation order
  )
)
setMethod(
  "initialize",
  "RPhosFateHelper",
  function(.Object, cmt) {
    cs_dir_old <- setwd(cmt@cv_dir[1L])
    on.exit(setwd(cs_dir_old))

    .Object@ex_cmt     <- extent(cmt@topo@rl_acc_wtd)
    .Object@is_res     <- as.integer(xres(cmt@topo@rl_acc_wtd))
    .Object@is_siz     <- as.integer(.Object@is_res^2)
    .Object@is_rws     <- nrow(cmt@topo@rl_acc_wtd)
    .Object@is_cls     <- ncol(cmt@topo@rl_acc_wtd)
    .Object@iv_fDo_dgl <- cmt@parameters@iv_fDo[c(1L, 3L, 7L, 9L)]
    .Object@im_fDo     <- matrix(cmt@parameters@iv_fDo, 3L)
    .Object@im_fDi     <- matrix(rev(cmt@parameters@iv_fDo), 3L)
    .Object@cv_rls     <- c("topo", "erosion", "transport")

    .Object@cv_met <- c(
      "NSE", "mNSE", "RMSE", "NRMSE",
      "PBIAS", "RSR", "GMRAE", "MdRAE"
    )
    if (cmt@ls_ini && file.exists("order.rds")) {
      .Object@order <- readRDS("order.rds")
    }

    .Object
  }
)

#### Class RPhosFate ####
setClass(
  "RPhosFate",
  slots = c(
    cv_dir     = "character", # Project directory
    ls_ini     = "logical",   # Load parameters from disc?
    is_MCi     = "integer",   # Monte Carlo iteration
    parameters = "RPhosFateParameters2",
    topo       = "RPhosFateTopo",
    erosion    = "RPhosFateErosion",
    transport  = "RPhosFateTransport",
    substance  = "RPhosFateSubstance",
    helper     = "RPhosFateHelper"
  ),
  prototype = list(
    ls_ini = FALSE,
    is_MCi = integer()
  )
)
setMethod(
  "initialize",
  "RPhosFate",
  function(.Object, arguments) {
    argumentNames <- names(arguments)

    .Object@cv_dir <- arguments$cv_dir
    if ("ls_ini" %in% argumentNames) {
      .Object@ls_ini <- arguments$ls_ini
    }
    if ("is_MCi" %in% argumentNames) {
      .Object@is_MCi <- arguments$is_MCi
    }
    validObject(.Object)

    cs_dir_old <- setwd(.Object@cv_dir[1L])
    on.exit(setwd(cs_dir_old))

    if (!dir.exists("Intermediate")) {
      dir.create("Intermediate")
    }
    if (!dir.exists("Result")) {
      dir.create("Result")
    }

    .Object@substance  <- new("RPhosFateSubstance", .Object)

    if (.Object@ls_ini && file.exists("parameters.yaml")) {
      arguments <- readParameters(arguments)
    } else if (.Object@ls_ini && file.exists("parameters.rds")) {
      arguments <- parametersRDS2YAML(slotNames(.Object@substance))
    }
    arguments <- arguments[setdiff(
      names(arguments),
      c("RPhosFate", "cv_dir", "ls_ini", "is_MCi")
    )]

    .Object@parameters <- new("RPhosFateParameters2", arguments)
    .Object@topo       <- new("RPhosFateTopo", .Object)
    .Object@erosion    <- new("RPhosFateErosion", .Object)
    .Object@transport  <- new("RPhosFateTransport", .Object)
    .Object@helper     <- new("RPhosFateHelper", .Object)

    .Object
  }
)
setValidity(
  "RPhosFate",
  function(object) {
    qassert(object@cv_dir, "S+"    , .var.name = "cv_dir")
    qassert(object@ls_ini, "B1"    , .var.name = "ls_ini")
    qassert(object@is_MCi, "I?[0,)", .var.name = "is_MCi")

    TRUE
  }
)

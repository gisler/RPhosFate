#' @import methods
#' @import raster
#' @importFrom graphics abline clip
#' @importFrom hydroGOF mNSE NSE pbias rsr
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
    ns_slp_min = "numeric",
    ns_slp_max = "numeric",
    ns_rhy_a   = "numeric",
    ns_rhy_b   = "numeric",
    ns_cha_rto = "numeric",
    ns_man_rip = "numeric",
    ns_man_cha = "numeric",
    ns_dep_ovl = "numeric",
    ns_dep_cha = "numeric",
    nv_enr_rto = "numeric",
    nv_tfc_inl = "numeric",
    iv_fDo     = "integer",
    nm_olc     = "matrix",
    df_cdt     = "data.frame"
  )
)
setMethod(
  "initialize",
  "RPhosFateParameters2",
  function(.Object, arguments) {
    # Min slope cap in %
    if (!is.null(arguments$ns_slp_min       )) {.Object@ns_slp_min <- arguments$ns_slp_min} else {.Object@ns_slp_min <- 0.001}

    # Max slope cap in %
    if (!is.null(arguments$ns_slp_max       )) {.Object@ns_slp_max <- arguments$ns_slp_max} else {.Object@ns_slp_max <- 999.0}

    # Parameter related to the discharge recurrence interval (WetSpa, T = 6)
    if (!is.null(arguments$ns_rhy_a         )) {.Object@ns_rhy_a   <- arguments$ns_rhy_a  } else {.Object@ns_rhy_a   <- 0.09}

    # Parameter related to the discharge recurrence interval (WetSpa, T = 6)
    if (!is.null(arguments$ns_rhy_b         )) {.Object@ns_rhy_b   <- arguments$ns_rhy_b  } else {.Object@ns_rhy_b   <- 0.50}

    # Ratio of channel width to cell length determining the riparian zone
    if (!is.null(arguments$ns_cha_rto       )) {.Object@ns_cha_rto <- arguments$ns_cha_rto} else {.Object@ns_cha_rto <- 0.5}

    # Riparian zone manning n
    if (!is.null(arguments$ns_man_rip       )) {.Object@ns_man_rip <- arguments$ns_man_rip} else {.Object@ns_man_rip <- 0.32}

    # Channel manning n
    if (!is.null(arguments$ns_man_cha       )) {.Object@ns_man_cha <- arguments$ns_man_cha} else {.Object@ns_man_cha <- 0.04}

    # Overland deposition coefficient
    if (!is.null(arguments$ns_dep_ovl       )) {.Object@ns_dep_ovl <- arguments$ns_dep_ovl} else {stop('"ns_dep_ovl" must be supplied.')}

    # Channel deposition coefficient
    if (!is.null(arguments$ns_dep_cha       )) {.Object@ns_dep_cha <- arguments$ns_dep_cha} else {stop('"ns_dep_cha" must be supplied.')}

    # Enrichment ratios
    if (!is.null(names(arguments$nv_enr_rto))) {.Object@nv_enr_rto <- arguments$nv_enr_rto} else {stop('"nv_enr_rto" must be supplied as a named vector.')}

    # Inlet transfer coefficients
    if (!is.null(names(arguments$nv_tfc_inl))) {.Object@nv_tfc_inl <- arguments$nv_tfc_inl} else {stop('"nv_tfc_inl" must be supplied as a named vector.')}

    # Outflow direction vector (ArcGIS coded)
    if (!is.null(arguments$iv_fDo           )) {.Object@iv_fDo     <- arguments$iv_fDo    } else {.Object@iv_fDo     <- c(32L, 16L, 8L, 64L, 0L, 4L, 128L, 1L, 2L)}

    # Catchment outlet coordinates
    if (!is.null(arguments$nm_olc           )) {.Object@nm_olc     <- arguments$nm_olc    } else {.Object@nm_olc     <- matrix()}

    # Calibration data
    if (!is.null(arguments$df_cdt           )) {.Object@df_cdt     <- arguments$df_cdt    } else {.Object@df_cdt     <- data.frame()}

    .Object
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
  "initialize",
  "RPhosFateHelper",
  function(.Object, cmt) {
    cs_dir_old <- setwd(cmt@cv_dir[1L])
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
    .Object@iv_fDo_dgl <- cmt@parameters@iv_fDo[c(1L, 3L, 7L, 9L)]

    # Outflow direction matrix
    .Object@im_fDo <- matrix(cmt@parameters@iv_fDo, 3L)

    # Inflow direction matrix
    .Object@im_fDi <- matrix(rev(cmt@parameters@iv_fDo), 3L)

    # Transport calculation order
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
    cv_dir     = "character",
    ls_ini     = "logical",
    is_MCi     = "integer",
    parameters = "RPhosFateParameters2",
    topo       = "RPhosFateTopo",
    erosion    = "RPhosFateErosion",
    transport  = "RPhosFateTransport",
    substance  = "RPhosFateSubstance",
    helper     = "RPhosFateHelper"
  )
)
setMethod(
  "initialize",
  "RPhosFate",
  function(.Object, arguments) {
    # Project directory
    if (!is.null(arguments$cv_dir)) .Object@cv_dir <- arguments$cv_dir else stop('"cv_dir" must be supplied.')

    # Load parameters from disc?
    if (!is.null(arguments$ls_ini)) .Object@ls_ini <- arguments$ls_ini else .Object@ls_ini <- FALSE

    # Monte Carlo iteration
    if (!is.null(arguments$is_MCi)) .Object@is_MCi <- arguments$is_MCi

    cs_dir_old <- setwd(.Object@cv_dir[1L])
    on.exit(setwd(cs_dir_old))

    if (!dir.exists("Intermediate") || !dir.exists("Result")) {
      dir.create("Intermediate", showWarnings = FALSE)
      dir.create("Result", showWarnings = FALSE)
    }

    .Object@substance  <- new("RPhosFateSubstance", .Object)

    if (.Object@ls_ini && file.exists("parameters.yaml")) {
      arguments <- readParameters(arguments)
    } else if (.Object@ls_ini && file.exists("parameters.rds")) {
      arguments <- parametersRDS2YAML(slotNames(.Object@substance))
    }

    .Object@parameters <- new("RPhosFateParameters2", arguments)
    .Object@topo       <- new("RPhosFateTopo", .Object)
    .Object@erosion    <- new("RPhosFateErosion", .Object)
    .Object@transport  <- new("RPhosFateTransport", .Object)
    .Object@helper     <- new("RPhosFateHelper", .Object)

    .Object
  }
)

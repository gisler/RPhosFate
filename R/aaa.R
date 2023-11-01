#' @import checkmate
#' @import methods
#' @import raster
#' @import terra
#' @importFrom graphics abline clip par
#' @importFrom Rcpp sourceCpp
#' @importFrom spatstat.geom as.owin as.ppp nncross
#' @importFrom stats cor optim optimize sd setNames
#' @importFrom yaml read_yaml write_yaml
#' @importFrom utils modifyList packageVersion
#' @useDynLib RPhosFate, .registration = TRUE
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
    nv_tfc_inl = "numeric",
    nv_enr_rto = "numeric",
    iv_fDo     = "integer",
    nm_olc     = "matrix",
    df_cdt     = "data.frame"
  ),
  prototype = list(
    ns_slp_min = 0.001,
    ns_slp_max = 999.0,
    ns_rhy_a   = 0.09,
    ns_rhy_b   = 0.50,
    ns_cha_rto = 0.5,
    ns_man_rip = 0.32,
    ns_man_cha = 0.04,
    ns_dep_ovl = numeric(),
    ns_dep_cha = numeric(),
    nv_tfc_inl = numeric(),
    nv_enr_rto = numeric(),
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
    qassert(object@ns_rhy_a  , "N1(0,)" , .var.name = "ns_rhy_a"  )
    qassert(object@ns_rhy_b  , "N1[0,)" , .var.name = "ns_rhy_b"  )
    qassert(object@ns_cha_rto, "N1(0,1]", .var.name = "ns_cha_rto")
    qassert(object@ns_man_rip, "N1(0,)" , .var.name = "ns_man_rip")
    qassert(object@ns_man_cha, "N1(0,)" , .var.name = "ns_man_cha")
    qassert(object@iv_fDo    , "I9[0,)" , .var.name = "iv_fDo"    )

    TRUE
  }
)

#### Class RPhosFateTopo ####
setClass(
  "RPhosFateTopo",
  slots = c(
    rl_acc     = "SpatRaster",
    rl_acc_wtd = "SpatRaster",
    rl_cha     = "SpatRaster",
    rl_clc     = "SpatRaster",
    rl_dem     = "SpatRaster",
    rl_dir     = "SpatRaster",
    rl_fid     = "SpatRaster",
    rl_inl     = "SpatRaster",
    rl_lue     = "SpatRaster",
    rl_rds     = "SpatRaster",
    rl_rip     = "SpatRaster",
    rl_slp     = "SpatRaster",
    rl_slp_cap = "SpatRaster",
    rl_wsh     = "SpatRaster"
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
    rl_RFa = "SpatRaster",
    rl_KFa = "SpatRaster",
    rl_LFa = "SpatRaster",
    rl_SFa = "SpatRaster",
    rl_CFa = "SpatRaster",
    rl_ero = "SpatRaster"
  )
)
setMethod(
  "initialize",
  "RPhosFateErosion",
  function(.Object, cmt) {
    cs_dir_old <- setwd(cmt@cv_dir[1L])
    on.exit(setwd(cs_dir_old))

    setwd("Input")
    .Object@rl_RFa <- readLayer(cmt, "RFa", TRUE)
    .Object@rl_KFa <- readLayer(cmt, "KFa", TRUE)
    .Object@rl_CFa <- readLayer(cmt, "CFa", TRUE)

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
    rl_man = "SpatRaster",
    rl_rhy = "SpatRaster"
  )
)
setMethod(
  "initialize",
  "RPhosFateTransport",
  function(.Object, cmt) {
    cs_dir_old <- setwd(cmt@cv_dir[1L])
    on.exit(setwd(cs_dir_old))

    setwd("Input")
    .Object@rl_man <- readLayer(cmt, "man")

    setwd("../Intermediate")
    .Object@rl_rhy <- readLayer(cmt, "rhy")

    .Object
  }
)

#### Class RPhosFateBare ####
setClass(
  "RPhosFateBare",
  slots = c(
    rl_xxr     = "SpatRaster",
    rl_xxt     = "SpatRaster",
    rl_xxt_inp = "SpatRaster",
    rl_xxt_out = "SpatRaster",
    rl_xxt_cld = "SpatRaster",
    rl_xxt_ctf = "SpatRaster"
  ),
  contains = "VIRTUAL"
)
#### Class RPhosFateConc ####
setClass(
  "RPhosFateConc",
  slots = c(
    rl_xxc = "SpatRaster",
    rl_xxe = "SpatRaster"
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
      "Result",
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

#### Class RPhosFateSubstances ####
setClass(
  "RPhosFateSubstances",
  slots = c(
    SS = "RPhosFateSS",
    PP = "RPhosFatePP"
  )
)
setMethod(
  "initialize",
  "RPhosFateSubstances",
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

#### Class RPhosFateHelpers ####
setClass(
  "RPhosFateHelpers",
  slots = c(
    cs_cmt     = "character",     # Coordinate reference system of river catchment
    ex_cmt     = "SpatExtent",    # Extent of river catchment
    is_res     = "integer",       # Cell resolution in m
    is_siz     = "integer",       # Cell area in m^2
    is_rws     = "integer",       # Number of rows
    is_cls     = "integer",       # Number of columns
    iv_fDo_dgl = "integer",       # Diagonal outflow direction vector
    im_fDo     = "matrix",        # Outflow direction matrix
    im_fDi     = "matrix",        # Inflow direction matrix
    cv_met     = "character",     # Implemented calibration quality metrics
    order      = "RPhosFateOrder" # Transport calculation order
  )
)
setMethod(
  "initialize",
  "RPhosFateHelpers",
  function(.Object, cmt) {
    cs_dir_old <- setwd(cmt@cv_dir[1L])
    on.exit(setwd(cs_dir_old))

    .Object@cs_cmt     <- crs(cmt@topo@rl_acc_wtd)
    .Object@ex_cmt     <- ext(cmt@topo@rl_acc_wtd)
    .Object@is_res     <- as.integer(xres(cmt@topo@rl_acc_wtd))
    .Object@is_siz     <- as.integer(.Object@is_res^2)
    .Object@is_rws     <- as.integer(nrow(cmt@topo@rl_acc_wtd))
    .Object@is_cls     <- as.integer(ncol(cmt@topo@rl_acc_wtd))
    .Object@iv_fDo_dgl <- cmt@parameters@iv_fDo[c(1L, 3L, 7L, 9L)]
    .Object@im_fDo     <- matrix(cmt@parameters@iv_fDo, 3L)
    .Object@im_fDi     <- matrix(rev(cmt@parameters@iv_fDo), 3L)

    .Object@cv_met <- c(
      "NSE", "mNSE", "KGE", "RMSE",
      "PBIAS", "RSR", "RCV",
      "GMRAE", "MdRAE"
    )
    if (cmt@ls_ini && file.exists("order.rds")) {
      .Object@order <- readRDS("order.rds")
    }

    .Object
  }
)

#### Class RPhosFate ####
#' RPhosFate class
#'
#' An S4 object representing a river catchment.
#'
#' @slot cv_dir A character vector holding the project root (first position) and
#'   optionally the Monte Carlo input data directory (second position).
#' @slot ls_ini A logical scalar specifying if the state of an existing project
#'   was loaded from disk.
#' @slot is_ths An integer scalar holding the number of threads to use for
#'   processing.
#' @slot is_MCi An integer scalar holding the current Monte Carlo iteration if
#'   applicable.
#' @slot cv_MCl A character vector holding the names of the layers, which shall
#'   be written to disk with the associated Monte Carlo iteration in their
#'   filenames upon calling the appropriate methods.
#' @slot parameters An S4 object holding the model parameters.
#' @slot topo An S4 object holding the raster layers related to topography in
#'   the broader sense.
#' @slot erosion An S4 object holding the raster layers related to erosion.
#' @slot transport An S4 object holding raster layers required for modelling
#'   transport.
#' @slot substances An S4 object holding the substance raster layer containers.
#' @slot helpers An S4 object holding helper data.
#'
#' @seealso [`RPhosFate`], [`catchment`]
#'
#' @export
setClass(
  "RPhosFate",
  slots = c(
    cv_dir     = "character",
    ls_ini     = "logical",
    is_ths     = "integer",
    is_MCi     = "integer",
    cv_MCl     = "character",
    parameters = "RPhosFateParameters2",
    topo       = "RPhosFateTopo",
    erosion    = "RPhosFateErosion",
    transport  = "RPhosFateTransport",
    substances = "RPhosFateSubstances",
    helpers    = "RPhosFateHelpers"
  ),
  prototype = list(
    cv_dir = character(),
    ls_ini = FALSE,
    is_ths = 1L,
    is_MCi = integer(),
    cv_MCl = "xxt"
  )
)
setMethod(
  "initialize",
  "RPhosFate",
  function(.Object, arguments) {
    argumentNames <- names(arguments)

    .Object@cv_dir <- normalizePath(
      arguments$cv_dir,
      winslash = .Platform$file.sep
    )
    if ("ls_ini" %in% argumentNames) {
      .Object@ls_ini <- arguments$ls_ini
    }
    if ("is_ths" %in% argumentNames) {
      .Object@is_ths <- arguments$is_ths
    }
    if ("is_MCi" %in% argumentNames) {
      .Object@is_MCi <- arguments$is_MCi
    }
    if ("cv_MCl" %in% argumentNames) {
      .Object@cv_MCl <- arguments$cv_MCl
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

    .Object@substances <- new("RPhosFateSubstances", .Object)

    if (.Object@ls_ini && file.exists("parameters.yaml")) {
      arguments <- readParameters(arguments)
      argumentNames <- names(arguments)
    } else if (.Object@ls_ini) {
      stop('"parameters.yaml" file does not exist.', call. = FALSE)
    }
    arguments <- arguments[setdiff(
      argumentNames,
      c("RPhosFate", "cv_dir", "ls_ini", "is_MCi", "cv_MCl")
    )]

    .Object@parameters <- new("RPhosFateParameters2", arguments)
    .Object@topo       <- new("RPhosFateTopo", .Object)
    .Object@erosion    <- new("RPhosFateErosion", .Object)
    .Object@transport  <- new("RPhosFateTransport", .Object)
    .Object@helpers    <- new("RPhosFateHelpers", .Object)

    .Object
  }
)
setValidity(
  "RPhosFate",
  function(object) {
    if (length(object@cv_dir) <= 2L) {
      assertDirectoryExists(object@cv_dir[1L], .var.name = "cv_dir[1L]")
    }
    if (length(object@cv_dir) > 1L) {
      qassert(object@cv_dir, "S2", .var.name = "cv_dir")
      assertDirectoryExists(object@cv_dir[2L], .var.name = "cv_dir[2L]")
    }
    qassert(object@ls_ini, "B1"    , .var.name = "ls_ini")
    qassert(object@is_ths, "I1(0,)", .var.name = "is_ths")
    qassert(object@is_MCi, "I?[0,)", .var.name = "is_MCi")
    qassert(object@cv_MCl, "S+"    , .var.name = "cv_MCl")

    TRUE
  }
)

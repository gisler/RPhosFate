#' @include aaa.R
NULL

#### erosionPrerequisites ####
#' @export
setGeneric(
  "erosionPrerequisites",
  function(x, ...) standardGeneric("erosionPrerequisites")
)
#' Erosion prerequisites
#'
#' Calculates capped slopes, L- and RUSLE S-factors. Weighted flow accumulations
#' less than one are set to one for the calculation of the L-factors.
#'
#' @param x An S4 [`RPhosFate-class`] river catchment object.
#'
#' @inherit catchment return
#'
#' @references
#' \cite{Desmet, P.J.J., Govers, G., 1996. A GIS procedure for automatically
#' calculating the USLE LS factor on topographically complex landscape units.
#' Journal of Soil and Water Conservation 51, 427–433.}
#'
#' \cite{Renard, K.G., Foster, G.R., Weesies, G.A., McCool, D.K., Yoder, D.C.,
#' 1997. Predicting soil erosion by water: a guide to conservation planning with
#' the Revised Universal Soil Loss Equation (RUSLE), Agriculture Handbook. U.S.
#' Government Printing Office, Washington, DC.}
#'
#' @examples
#' \dontrun{
#' # create temporary demonstration project
#' cv_dir <- demoProject()
#' # load temporary demonstration project
#' x <- RPhosFate(
#'   cv_dir = cv_dir,
#'   ls_ini = TRUE
#' )
#'
#' x <- erosionPrerequisites(x)
#' }
#'
#' @aliases erosionPrerequisites
#'
#' @export
setMethod(
  "erosionPrerequisites",
  "RPhosFate",
  function(x) {
    compareRaster(x@topo@rl_acc_wtd, x@topo@rl_slp, x@topo@rl_cha)

    cs_dir_old <- setwd(file.path(x@cv_dir[1L], "Intermediate"))
    on.exit(setwd(cs_dir_old))

    # Capped slope in % (also relevant for transport)
    x@topo@rl_slp_cap <- x@topo@rl_slp
    x@topo@rl_slp_cap[x@topo@rl_slp_cap < x@parameters@ns_slp_min] <- x@parameters@ns_slp_min
    x@topo@rl_slp_cap[x@topo@rl_slp_cap > x@parameters@ns_slp_max] <- x@parameters@ns_slp_max

    # Capped slope in radian
    rl_slp_cap_rad <- calc(
      x@topo@rl_slp_cap,
      function(x) {
        atan(x * 1e-2)
      }
    )

    # Weighted overland flow accumulation
    rl_acc_wtd_ovl <- x@topo@rl_acc_wtd
    rl_acc_wtd_ovl[!is.na(x@topo@rl_cha)] <- NA_real_
    rl_acc_wtd_ovl[rl_acc_wtd_ovl < 1] <- 1

    # Ratio of rill to interrill erosion
    rl_LFa_b <- calc(
      rl_slp_cap_rad,
      function(x) {
        sin(x) / (0.0896 * (3 * sin(x)^0.8 + 0.56))
      }
    )
    # Rill erodibility parameter
    rl_LFa_m <- calc(
      rl_LFa_b,
      function(x) {
        x / (1 + x)
      }
    )

    # L factor
    is_res = x@helpers@is_res
    x@erosion@rl_LFa <- overlay(
      x = rl_acc_wtd_ovl,
      y = rl_LFa_m,
      fun = function(x, y) {
        ((x * is_res)^(1 + y) - ((x - 1) * is_res)^(1 + y)) / # nolint
          (is_res * 22.13^y)
      },
      filename = sprintf("LFa%s", x@cs_fex),
      datatype = "FLT8S",
      options = "COMPRESSED=YES",
      overwrite = TRUE
    )

    # S factor
    x@erosion@rl_SFa <- overlay(
      x = rl_slp_cap_rad,
      y = x@topo@rl_slp_cap,
      fun = function(x, y) {
        ifelse(y < 9, 10.8 * sin(x) + 0.03, 16.8 * sin(x) - 0.5)
      },
      filename = sprintf("SFa%s", x@cs_fex),
      datatype = "FLT8S",
      options = "COMPRESSED=YES",
      overwrite = TRUE
    )

    filename <- sprintf("slp_cap%s", x@cs_fex)
    writeRaster(
      x@topo@rl_slp_cap,
      filename = filename,
      datatype = "FLT8S",
      options = "COMPRESSED=YES",
      overwrite = TRUE
    )
    x@topo@rl_slp_cap <- raster(filename)

    x
  }
)

#### erosion ####
#' @export
setGeneric(
  "erosion",
  function(x, ...) standardGeneric("erosion")
)
#' Erosion
#'
#' Calculates (R)USLE erosion.
#'
#' @inheritParams erosionPrerequisites,RPhosFate-method
#'
#' @inherit catchment return
#'
#' @references
#' \cite{Renard, K.G., Foster, G.R., Weesies, G.A., McCool, D.K., Yoder, D.C.,
#' 1997. Predicting soil erosion by water: a guide to conservation planning with
#' the Revised Universal Soil Loss Equation (RUSLE), Agriculture Handbook. U.S.
#' Government Printing Office, Washington, DC.}
#'
#' \cite{Wischmeier, W.H., Smith, D.D., 1978. Predicting rainfall erosion
#' losses. A guide to conservation planning, Agriculture Handbook. U.S.
#' Government Printing Office, Washington, DC.}
#'
#' @examples
#' \dontrun{
#' # create temporary demonstration project
#' cv_dir <- demoProject()
#' # load temporary demonstration project
#' x <- RPhosFate(
#'   cv_dir = cv_dir,
#'   ls_ini = TRUE
#' )
#' # presupposed function call
#' x <- erosionPrerequisites(x)
#'
#' x <- erosion(x)
#' }
#'
#' @aliases erosion
#'
#' @export
setMethod(
  "erosion",
  "RPhosFate",
  function(x) {
    compareRaster(
      x@topo@rl_acc_wtd,
      x@erosion@rl_RFa,
      x@erosion@rl_KFa,
      x@erosion@rl_LFa,
      x@erosion@rl_SFa,
      x@erosion@rl_CFa
    )

    cs_dir_old <- setwd(file.path(x@cv_dir[1L], "Result"))
    on.exit(setwd(cs_dir_old))

    # Erosion in t/cell/yr
    x@erosion@rl_ero <-
      x@erosion@rl_RFa *
      x@erosion@rl_KFa *
      x@erosion@rl_LFa *
      x@erosion@rl_SFa *
      x@erosion@rl_CFa *
      x@helpers@is_siz * 1e-4

    filename <- paste0("ero", x@is_MCi, x@cs_fex)
    writeRaster(
      x@erosion@rl_ero,
      filename = filename,
      datatype = "FLT8S",
      options = "COMPRESSED=YES",
      overwrite = TRUE
    )
    x@erosion@rl_ero <- raster(filename)

    x
  }
)

#### emission ####
#' @export
setGeneric(
  "emission",
  function(x, ...) standardGeneric("emission")
)
#' Emission
#'
#' Calculates substance emissions.
#'
#' @inheritParams erosionPrerequisites,RPhosFate-method
#' @param substance A character string specifying the substance to calculate.
#'
#' @inherit catchment return
#'
#' @examples
#' \dontrun{
#' # create temporary demonstration project
#' cv_dir <- demoProject()
#' # load temporary demonstration project
#' x <- RPhosFate(
#'   cv_dir = cv_dir,
#'   ls_ini = TRUE
#' )
#' # presupposed function calls
#' x <- erosionPrerequisites(x)
#' x <- erosion(x)
#'
#' x <- emission(x, "PP")
#' }
#'
#' @aliases emission
#'
#' @export
setMethod(
  "emission",
  "RPhosFate",
  function(x, substance = "PP") {
    assertSubstance(x, substance)
    compareRaster(
      x@topo@rl_acc_wtd,
      x@erosion@rl_ero,
      slot(x@substances, substance)@rl_xxc,
      x@topo@rl_clc
    )

    cs_dir_old <- setwd(file.path(x@cv_dir[1L], "Result"))
    on.exit(setwd(cs_dir_old))

    filename <- paste0(tolower(substance), "e", x@is_MCi, x@cs_fex)
    # Emission in kg/cell/yr
    slot(x@substances, substance)@rl_xxe <- overlay(
      x = x@erosion@rl_ero,
      y = slot(x@substances, substance)@rl_xxc,
      z = x@topo@rl_clc,
      fun = function(x, y, z) {
        x * y * (1 + z * 1e-2) * 1e-3
      },
      filename = filename,
      datatype = "FLT8S",
      options = "COMPRESSED=YES",
      overwrite = TRUE
    )

    x
  }
)

#### transportPrerequisites ####
#' @export
setGeneric(
  "transportPrerequisites",
  function(x, ...) standardGeneric("transportPrerequisites")
)
#' Transport prerequisites
#'
#' Calculates hydraulic radii and determines cells representing inlets as well
#' as riparian zones.
#'
#' @inheritParams erosionPrerequisites,RPhosFate-method
#'
#' @inherit catchment return
#'
#' @references
#' \cite{Molnár, P., Ramírez, J.A., 1998. Energy dissipation theories and
#' optimal channel characteristics of river networks. Water Resources Research
#' 34, 1809–1818.}
#'
#' @examples
#' \dontrun{
#' # create temporary demonstration project
#' cv_dir <- demoProject()
#' # load temporary demonstration project
#' x <- RPhosFate(
#'   cv_dir = cv_dir,
#'   ls_ini = TRUE
#' )
#'
#' x <- transportPrerequisites(x)
#' }
#'
#' @aliases transportPrerequisites
#'
#' @export
setMethod(
  "transportPrerequisites",
  "RPhosFate",
  function(x) {
    compareRaster(
      x@topo@rl_acc_wtd,
      x@topo@rl_dir,
      x@topo@rl_cha,
      x@topo@rl_rds
    )

    cs_dir_old <- setwd(file.path(x@cv_dir[1L], "Intermediate"))
    on.exit(setwd(cs_dir_old))

    # Hydraulic radius in m
    ns_rhy_a <- x@parameters@ns_rhy_a
    ns_rhy_b <- x@parameters@ns_rhy_b
    is_siz <- x@helpers@is_siz
    x@transport@rl_rhy <- calc(
      x@topo@rl_acc_wtd,
      function(x) {
        ns_rhy_a * (x * is_siz * 1e-6)^ns_rhy_b
      },
      filename = sprintf("rhy%s", x@cs_fex),
      datatype = "FLT8S",
      options = "COMPRESSED=YES",
      overwrite = TRUE
    )

    # Riparian zone cells
    x@topo@rl_rip <- raster(
      dir_sth(
        im_dir = as.matrix(x@topo@rl_dir),
        im_sth = as.matrix(x@topo@rl_cha),
        im_fDo = x@helpers@im_fDo
      ),
      template = x@topo@rl_acc_wtd
    )

    # Inlet cells
    x@topo@rl_inl <- raster(
      dir_sth(
        im_dir = as.matrix(x@topo@rl_dir),
        im_sth = as.matrix(x@topo@rl_rds),
        im_fDo = x@helpers@im_fDo
      ),
      template = x@topo@rl_acc_wtd
    )

    # No inlet cells at channel cells
    x@topo@rl_inl[!is.na(x@topo@rl_cha)] <- NA_integer_
    # No riparian zone cells at road cells
    x@topo@rl_rip[!is.na(x@topo@rl_rds)] <- NA_integer_
    # No inlet cells at riparian zone cells
    x@topo@rl_inl[!is.na(x@topo@rl_rip)] <- NA_integer_

    # Nearest channel cells for inlet cells
    df_out <- findNearestNeighbour(
      rasterToPoints(x@topo@rl_inl),
      rasterToPoints(x@topo@rl_cha),
      x@helpers@ex_cmt
    )

    # X-coordinates of nearest channel cells to column numbers
    df_out$Y.x <- (df_out$Y.x + x@helpers@is_res / 2 -
      x@helpers@ex_cmt[1L]) / x@helpers@is_res
    # Y-coordinates of nearest channel cells to row numbers
    df_out$Y.y <- x@helpers@is_rws - ((df_out$Y.y - x@helpers@is_res / 2 -
      x@helpers@ex_cmt[3L]) / x@helpers@is_res)

    # Substituting inlet values with integer codes identifying nearest channel
    # cells
    df_out$code <- as.integer(df_out$Y.y * x@helpers@is_cls + df_out$Y.x)
    # Bug in subs() {raster}: use default by and which
    x@topo@rl_inl <- subs(x@topo@rl_inl, y = df_out[c(3L, 8L)])

    filename <- sprintf("rip%s", x@cs_fex)
    writeRaster(
      x@topo@rl_rip,
      filename = filename,
      datatype = "INT4S",
      options = "COMPRESSED=YES",
      overwrite = TRUE
    )
    x@topo@rl_rip <- raster(filename)

    filename <- sprintf("inl%s", x@cs_fex)
    writeRaster(
      x@topo@rl_inl,
      filename = filename,
      datatype = "INT4S",
      options = "COMPRESSED=YES",
      overwrite = TRUE
    )
    x@topo@rl_inl <- raster(filename)

    x
  }
)

#### transportCalcOrder ####
#' @export
setGeneric(
  "transportCalcOrder",
  function(x, ...) standardGeneric("transportCalcOrder")
)
#' Transport calculation order
#'
#' Determines the cell transport calculation order.
#'
#' @inheritParams erosionPrerequisites,RPhosFate-method
#'
#' @inherit catchment return
#'
#' @examples
#' \dontrun{
#' # create temporary demonstration project
#' cv_dir <- demoProject()
#' # load temporary demonstration project
#' x <- RPhosFate(
#'   cv_dir = cv_dir,
#'   ls_ini = TRUE
#' )
#'
#' x <- transportCalcOrder(x)
#' }
#'
#' @aliases transportCalcOrder
#'
#' @export
setMethod(
  "transportCalcOrder",
  "RPhosFate",
  function(x) {
    compareRaster(
      x@topo@rl_acc_wtd,
      x@topo@rl_acc,
      x@topo@rl_cha
    )

    # Overland flow accumulation
    rl_acc_ovl <- x@topo@rl_acc
    rl_acc_ovl[!is.na(x@topo@rl_cha)] <- NA_integer_

    # Channel flow accumulation
    rl_acc_cha <- x@topo@rl_acc
    rl_acc_cha[ is.na(x@topo@rl_cha)] <- NA_integer_

    # Transport calculation order as column-major index
    im_acc_ovl <- as.matrix(rl_acc_ovl)
    im_acc_cha <- as.matrix(rl_acc_cha)
    ar_ord_ovl <- tapply(seq_along(im_acc_ovl), im_acc_ovl, identity)
    ar_ord_cha <- tapply(seq_along(im_acc_cha), im_acc_cha, identity)

    # Row order from index
    is_rws <- x@helpers@is_rws
    fun <- function(x) {
      (x + is_rws) - ceiling(x / is_rws) * is_rws
    }
    iv_ord_ovl_row <- as.integer(unlist(lapply(ar_ord_ovl, fun)))
    iv_ord_cha_row <- as.integer(unlist(lapply(ar_ord_cha, fun)))

    # Column order from index
    fun <- function(x) {
      ceiling(x / is_rws)
    }
    iv_ord_ovl_col <- as.integer(unlist(lapply(ar_ord_ovl, fun)))
    iv_ord_cha_col <- as.integer(unlist(lapply(ar_ord_cha, fun)))

    # Overland as well as channel row and column numbers for top-down
    # computation (C++ has zero-based numbering)
    x@helpers@order@iv_ord_row <- c(iv_ord_ovl_row, iv_ord_cha_row) - 1L
    x@helpers@order@iv_ord_col <- c(iv_ord_ovl_col, iv_ord_cha_col) - 1L

    # Reverse overland row and column numbers for bottom-up computation
    # (C++ has zero-based numbering)
    x@helpers@order@iv_ord_ovl_row_rev <- rev(iv_ord_ovl_row) - 1L
    x@helpers@order@iv_ord_ovl_col_rev <- rev(iv_ord_ovl_col) - 1L

    x
  }
)

#### transport ####
#' @export
setGeneric(
  "transport",
  function(x, ...) standardGeneric("transport")
)
#' Transport
#'
#' Calculates substance retentions, transports and cell loads as well as
#' transfers among others.
#'
#' @inheritParams emission,RPhosFate-method
#'
#' @inherit catchment return
#'
#' @references
#' \cite{Engman, E.T., 1986. Roughness coefficients for routing surface runoff.
#' Journal of Irrigation and Drainage Engineering 112, 39–53.}
#'
#' @examples
#' \dontrun{
#' # create temporary demonstration project
#' cv_dir <- demoProject()
#' # load temporary demonstration project
#' x <- RPhosFate(
#'   cv_dir = cv_dir,
#'   ls_ini = TRUE
#' )
#' # presupposed function calls
#' x <- erosionPrerequisites(x)
#' x <- erosion(x)
#' x <- emission(x, "PP")
#' x <- transportPrerequisites(x)
#' x <- transportCalcOrder(x)
#'
#' x <- transport(x, "PP")
#' }
#'
#' @aliases transport
#'
#' @export
setMethod(
  "transport",
  "RPhosFate",
  function(x, substance = "PP") {
    assertSubstance(x, substance)
    compareRaster(
      x@topo@rl_acc_wtd,
      x@topo@rl_cha,
      x@topo@rl_dir,
      x@topo@rl_inl,
      x@topo@rl_rip,
      x@transport@rl_man,
      if (substance == "SS") {
        x@erosion@rl_ero
      } else {
        slot(x@substances, substance)@rl_xxe
      },
      x@transport@rl_rhy,
      x@topo@rl_slp_cap
    )
    qassert(x@parameters@ns_dep_ovl, "N1(0,)", .var.name = "ns_dep_ovl")
    qassert(x@parameters@ns_dep_cha, "N1[0,)", .var.name = "ns_dep_cha")
    qassert(
      x@parameters@nv_tfc_inl[substance],
      "N1[0,1]",
      .var.name = "nv_tfc_inl[substance]"
    )
    if (substance != "SS") {
      qassert(
        x@parameters@nv_enr_rto[substance],
        "N1[1,)",
        .var.name = "nv_enr_rto[substance]"
      )
    }

    cs_dir_old <- setwd(file.path(x@cv_dir[1L], "Result"))
    on.exit(setwd(cs_dir_old))

    li_tpt <- transportCpp(
      parameters = x@parameters,
      ns_dep_ovl = if (substance == "SS") {
        x@parameters@ns_dep_ovl
      } else {
        x@parameters@ns_dep_ovl / x@parameters@nv_enr_rto[substance]
      },
      ns_tfc_inl = x@parameters@nv_tfc_inl[substance],
      helpers    = x@helpers,
      order      = x@helpers@order,
      im_cha     = as.matrix(x@topo@rl_cha),
      im_dir     = as.matrix(x@topo@rl_dir),
      im_inl     = as.matrix(x@topo@rl_inl),
      im_rip     = as.matrix(x@topo@rl_rip),
      nm_man     = as.matrix(x@transport@rl_man),
      nm_xxe     = if (substance == "SS") {
        as.matrix(x@erosion@rl_ero)
      } else {
        as.matrix(slot(x@substances, substance)@rl_xxe)
      },
      nm_rhy     = as.matrix(x@transport@rl_rhy),
      nm_slp     = as.matrix(x@topo@rl_slp_cap)
    )

    layers <- paste0("xx", c("r", "t_inp", "t_out", "t_cld", "t_ctf", "t"))
    filenames <- setNames(
      paste0(
        sub("^xx", tolower(substance), layers),
        x@is_MCi,
        x@cs_fex
      ),
      layers
    )

    if (length(x@is_MCi) == 0L) {
      writeRaster(raster(li_tpt$nm_xxr    , template = x@topo@rl_acc_wtd), filename = filenames["xxr"    ], datatype = "FLT8S", options = "COMPRESSED=YES", overwrite = TRUE)
      writeRaster(raster(li_tpt$nm_xxt_inp, template = x@topo@rl_acc_wtd), filename = filenames["xxt_inp"], datatype = "FLT8S", options = "COMPRESSED=YES", overwrite = TRUE)
      writeRaster(raster(li_tpt$nm_xxt_out, template = x@topo@rl_acc_wtd), filename = filenames["xxt_out"], datatype = "FLT8S", options = "COMPRESSED=YES", overwrite = TRUE)
      writeRaster(raster(li_tpt$nm_xxt_cld, template = x@topo@rl_acc_wtd), filename = filenames["xxt_cld"], datatype = "FLT8S", options = "COMPRESSED=YES", overwrite = TRUE)
      writeRaster(raster(li_tpt$nm_xxt_ctf, template = x@topo@rl_acc_wtd), filename = filenames["xxt_ctf"], datatype = "FLT8S", options = "COMPRESSED=YES", overwrite = TRUE)
      slot(x@substances, substance)@rl_xxr     <- raster(filenames["xxr"    ])
      slot(x@substances, substance)@rl_xxt_inp <- raster(filenames["xxt_inp"])
      slot(x@substances, substance)@rl_xxt_out <- raster(filenames["xxt_out"])
      slot(x@substances, substance)@rl_xxt_cld <- raster(filenames["xxt_cld"])
      slot(x@substances, substance)@rl_xxt_ctf <- raster(filenames["xxt_ctf"])
    }
    writeRaster(
      raster(li_tpt$nm_xxt, template = x@topo@rl_acc_wtd),
      filename = filenames["xxt"],
      datatype = "FLT8S",
      options = "COMPRESSED=YES",
      overwrite = TRUE
    )
    slot(x@substances, substance)@rl_xxt <- raster(filenames["xxt"])

    x
  }
)

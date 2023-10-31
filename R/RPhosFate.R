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
#' Calculates and writes capped slopes, L- and RUSLE S-factors (equations for
#' summer conditions and slopes \eqn{\geq}{≥} 15 ft) to disk. Weighted flow
#' accumulations less than one are set to one for the calculation of the
#' L-factors.
#'
#' @param x An S4 [`RPhosFate-class`] river catchment object.
#'
#' @return An S4 [`RPhosFate-class`] river catchment object and side effects in
#'   the form of raster files.
#'
#' @references
#' \cite{Renard, K.G., Foster, G.R., Weesies, G.A., McCool, D.K., Yoder, D.C.,
#' 1997. Predicting soil erosion by water: a guide to conservation planning with
#' the Revised Universal Soil Loss Equation (RUSLE), Agriculture Handbook. U.S.
#' Government Printing Office, Washington, DC.}
#'
#' @seealso [`firstRun`], [`subsequentRun`]
#'
#' @examples
#' \donttest{
#' # temporary demonstration project copy
#' cv_dir <- demoProject()
#' # load temporary demonstration project
#' x <- RPhosFate(
#'   cv_dir = cv_dir,
#'   ls_ini = TRUE
#' )
#'
#' x <- erosionPrerequisites(x)}
#'
#' @aliases erosionPrerequisites
#'
#' @export
setMethod(
  "erosionPrerequisites",
  "RPhosFate",
  function(x) {
    compareGeom(x@topo@rl_acc_wtd, x@topo@rl_slp, x@topo@rl_cha)

    cs_dir_old <- setwd(file.path(x@cv_dir[1L], "Intermediate"))
    on.exit(setwd(cs_dir_old))

    # Capped slope in % (also relevant for transport)
    x@topo@rl_slp_cap <- x@topo@rl_slp
    x@topo@rl_slp_cap[x@topo@rl_slp_cap < x@parameters@ns_slp_min] <- x@parameters@ns_slp_min
    x@topo@rl_slp_cap[x@topo@rl_slp_cap > x@parameters@ns_slp_max] <- x@parameters@ns_slp_max

    # Capped slope in radian
    rl_slp_cap_rad <- atan(x@topo@rl_slp_cap * 1e-2)

    # Weighted overland flow accumulation
    rl_acc_wtd_ovl <- x@topo@rl_acc_wtd
    rl_acc_wtd_ovl[!is.na(x@topo@rl_cha)] <- NA_real_
    rl_acc_wtd_ovl[rl_acc_wtd_ovl < 1] <- 1

    # Ratio of rill to interrill erosion
    rl_LFa_b <- app(
      rl_slp_cap_rad,
      function(x) {
        sin(x) / (0.0896 * (3 * sin(x)^0.8 + 0.56))
      }
    )
    # Rill erodibility parameter
    rl_LFa_m <- rl_LFa_b / (1 + rl_LFa_b)

    # L factor
    is_res <- x@helpers@is_res
    x@erosion@rl_LFa <- lapp(
      c(x = rl_acc_wtd_ovl, y = rl_LFa_m),
      fun = function(x, y) {
        ((x * is_res)^(1 + y) - ((x - 1) * is_res)^(1 + y)) / # nolint
          (is_res * 22.13^y)
      }
    )

    # S factor
    x@erosion@rl_SFa <- lapp(
      c(x = rl_slp_cap_rad, y = x@topo@rl_slp_cap),
      fun = function(x, y) {
        ifelse(y < 9, 10.8 * sin(x) + 0.03, 16.8 * sin(x) - 0.5)
      }
    )

    x@topo@rl_slp_cap <- writeLayer(x, "slp_cap", x@topo@rl_slp_cap, "FLT8S")
    x@erosion@rl_LFa  <- writeLayer(x, "LFa"    , x@erosion@rl_LFa , "FLT8S")
    x@erosion@rl_SFa  <- writeLayer(x, "SFa"    , x@erosion@rl_SFa , "FLT8S")

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
#' Calculates and writes (R)USLE erosion to disk.
#'
#' @inheritParams erosionPrerequisites,RPhosFate-method
#'
#' @inherit erosionPrerequisites,RPhosFate-method return
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
#' @seealso [`firstRun`], [`subsequentRun`]
#'
#' @examples
#' \donttest{
#' # temporary demonstration project copy
#' cv_dir <- demoProject()
#' # load temporary demonstration project
#' x <- RPhosFate(
#'   cv_dir = cv_dir,
#'   ls_ini = TRUE
#' )
#' # presupposed method call
#' x <- erosionPrerequisites(x)
#'
#' x <- erosion(x)}
#'
#' @aliases erosion
#'
#' @export
setMethod(
  "erosion",
  "RPhosFate",
  function(x) {
    compareGeom(
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

    x@erosion@rl_ero <- writeLayer(x, "ero", x@erosion@rl_ero, "FLT8S")

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
#' Calculates and writes substance emissions to disk.
#'
#' @inheritParams erosionPrerequisites,RPhosFate-method
#' @param substance A character string specifying the substance to calculate.
#'
#' @inherit erosionPrerequisites,RPhosFate-method return
#'
#' @seealso [`firstRun`], [`subsequentRun`]
#'
#' @examples
#' \donttest{
#' # temporary demonstration project copy
#' cv_dir <- demoProject()
#' # load temporary demonstration project
#' x <- RPhosFate(
#'   cv_dir = cv_dir,
#'   ls_ini = TRUE
#' )
#' # presupposed method calls
#' x <- erosionPrerequisites(x)
#' x <- erosion(x)
#'
#' x <- emission(x, "PP")}
#'
#' @aliases emission
#'
#' @export
setMethod(
  "emission",
  "RPhosFate",
  function(x, substance = "PP") {
    assertChoice(substance, slotNames(x@substances))
    compareGeom(
      x@topo@rl_acc_wtd,
      x@erosion@rl_ero,
      slot(x@substances, substance)@rl_xxc,
      x@topo@rl_clc
    )

    cs_dir_old <- setwd(file.path(x@cv_dir[1L], "Result"))
    on.exit(setwd(cs_dir_old))

    # Emission in kg/cell/yr
    slot(x@substances, substance)@rl_xxe <- lapp(
      c(
        x = x@erosion@rl_ero,
        y = slot(x@substances, substance)@rl_xxc,
        z = x@topo@rl_clc
      ),
      fun = function(x, y, z) {
        x * y * (1 + z * 1e-2) * 1e-3
      }
    )

    slot(x@substances, substance)@rl_xxe <- writeLayer(
      x,
      "xxe",
      slot(x@substances, substance)@rl_xxe,
      "FLT8S",
      substance
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
#' as riparian zones before writing them to disk.
#'
#' @inheritParams erosionPrerequisites,RPhosFate-method
#'
#' @inherit erosionPrerequisites,RPhosFate-method return
#'
#' @references
#' \cite{Molnár, P., Ramírez, J.A., 1998. Energy dissipation theories and
#' optimal channel characteristics of river networks. Water Resources Research
#' 34, 1809–1818.}
#'
#' @seealso [`firstRun`], [`subsequentRun`]
#'
#' @examples
#' \donttest{
#' # temporary demonstration project copy
#' cv_dir <- demoProject()
#' # load temporary demonstration project
#' x <- RPhosFate(
#'   cv_dir = cv_dir,
#'   ls_ini = TRUE
#' )
#'
#' x <- transportPrerequisites(x)}
#'
#' @aliases transportPrerequisites
#'
#' @export
setMethod(
  "transportPrerequisites",
  "RPhosFate",
  function(x) {
    compareGeom(
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
    x@transport@rl_rhy <- app(
      x@topo@rl_acc_wtd,
      function(x) {
        ns_rhy_a * (x * is_siz * 1e-6)^ns_rhy_b
      }
    )

    # Riparian zone cells
    x@topo@rl_rip <- rast(
      dir_sth(
        im_dir = as.matrix(x@topo@rl_dir, wide = TRUE),
        im_sth = as.matrix(x@topo@rl_cha, wide = TRUE),
        im_fDo = x@helpers@im_fDo
      ),
      crs = x@helpers@cs_cmt,
      extent = x@helpers@ex_cmt
    )

    # Inlet cells
    x@topo@rl_inl <- rast(
      dir_sth(
        im_dir = as.matrix(x@topo@rl_dir, wide = TRUE),
        im_sth = as.matrix(x@topo@rl_rds, wide = TRUE),
        im_fDo = x@helpers@im_fDo
      ),
      crs = x@helpers@cs_cmt,
      extent = x@helpers@ex_cmt
    )
    set.names(x@topo@rl_inl, "inl")

    # No inlet cells at channel cells
    x@topo@rl_inl[!is.na(x@topo@rl_cha)] <- NA_integer_
    # No riparian zone cells at road cells
    x@topo@rl_rip[!is.na(x@topo@rl_rds)] <- NA_integer_
    # No inlet cells at riparian zone cells
    x@topo@rl_inl[!is.na(x@topo@rl_rip)] <- NA_integer_

    # Nearest channel cells for inlet cells
    df_out <- findNearestNeighbour(
      cbind(crds(x@topo@rl_inl), values(x@topo@rl_inl, na.rm = TRUE)),
      cbind(crds(x@topo@rl_cha), values(x@topo@rl_cha, na.rm = TRUE)),
      x@helpers@ex_cmt
    )

    # X-coordinates of nearest channel cells to column numbers
    df_out$Y.x <- colFromX(x@topo@rl_inl, df_out$Y.x)
    # Y-coordinates of nearest channel cells to row numbers
    df_out$Y.y <- rowFromY(x@topo@rl_inl, df_out$Y.y)

    # Substituting inlet values with integer codes identifying nearest channel
    # cells
    df_out$code <- as.integer(df_out$Y.y * x@helpers@is_cls + df_out$Y.x)
    x@topo@rl_inl <- subst(x@topo@rl_inl, df_out[["inl"]], df_out[["code"]])

    x@transport@rl_rhy <- writeLayer(x, "rhy", x@transport@rl_rhy, "FLT8S")
    x@topo@rl_rip      <- writeLayer(x, "rip", x@topo@rl_rip     , "INT4S")
    x@topo@rl_inl      <- writeLayer(x, "inl", x@topo@rl_inl     , "INT4S")

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
#' @seealso [`firstRun`], [`subsequentRun`]
#'
#' @examples
#' \donttest{
#' # temporary demonstration project copy
#' cv_dir <- demoProject()
#' # load temporary demonstration project
#' x <- RPhosFate(
#'   cv_dir = cv_dir,
#'   ls_ini = TRUE
#' )
#'
#' x <- transportCalcOrder(x)}
#'
#' @aliases transportCalcOrder
#'
#' @export
setMethod(
  "transportCalcOrder",
  "RPhosFate",
  function(x) {
    compareGeom(
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

    # Transport calculation order as row-major index
    im_acc_ovl <- as.matrix(rl_acc_ovl)
    im_acc_cha <- as.matrix(rl_acc_cha)
    ar_ord_ovl <- tapply(seq_along(im_acc_ovl), im_acc_ovl, identity, simplify = FALSE)
    ar_ord_cha <- tapply(seq_along(im_acc_cha), im_acc_cha, identity, simplify = FALSE)

    # Row and column numbers from index (C++ has zero-based indexing)
    iv_ord_ovl <- rowColFromCell(rl_acc_ovl, unlist(ar_ord_ovl)) - 1
    iv_ord_cha <- rowColFromCell(rl_acc_cha, unlist(ar_ord_cha)) - 1
    storage.mode(iv_ord_ovl) <- "integer"
    storage.mode(iv_ord_cha) <- "integer"

    # Overland as well as channel row and column numbers for top-down
    # computation
    x@helpers@order@iv_ord_row <- c(iv_ord_ovl[, 1L], iv_ord_cha[, 1L])
    x@helpers@order@iv_ord_col <- c(iv_ord_ovl[, 2L], iv_ord_cha[, 2L])

    # Reverse overland row and column numbers for bottom-up computation
    x@helpers@order@iv_ord_ovl_row_rev <- rev(iv_ord_ovl[, 1L])
    x@helpers@order@iv_ord_ovl_col_rev <- rev(iv_ord_ovl[, 2L])

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
#' Calculates and writes substance retentions, transports and cell loads as well
#' as transfers to disk.
#'
#' @inheritParams emission,RPhosFate-method
#'
#' @inherit erosionPrerequisites,RPhosFate-method return
#'
#' @references
#' \cite{Engman, E.T., 1986. Roughness coefficients for routing surface runoff.
#' Journal of Irrigation and Drainage Engineering 112, 39–53.}
#'
#' @seealso [`firstRun`], [`subsequentRun`]
#'
#' @examples
#' \donttest{
#' # temporary demonstration project copy
#' cv_dir <- demoProject()
#' # load temporary demonstration project
#' x <- RPhosFate(
#'   cv_dir = cv_dir,
#'   ls_ini = TRUE
#' )
#' # presupposed method calls
#' x <- erosionPrerequisites(x)
#' x <- erosion(x)
#' x <- emission(x, "PP")
#' x <- transportPrerequisites(x)
#' x <- transportCalcOrder(x)
#'
#' x <- transport(x, "PP")}
#'
#' @aliases transport
#'
#' @export
setMethod(
  "transport",
  "RPhosFate",
  function(x, substance = "PP") {
    assertChoice(substance, slotNames(x@substances))
    compareGeom(
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
    qassert(x@parameters@ns_dep_ovl, "N1[0,)", .var.name = "ns_dep_ovl")
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
      im_cha     = as.matrix(x@topo@rl_cha, wide = TRUE),
      im_dir     = as.matrix(x@topo@rl_dir, wide = TRUE),
      im_inl     = as.matrix(x@topo@rl_inl, wide = TRUE),
      im_rip     = as.matrix(x@topo@rl_rip, wide = TRUE),
      nm_man     = as.matrix(x@transport@rl_man, wide = TRUE),
      nm_xxe     = if (substance == "SS") {
        as.matrix(x@erosion@rl_ero, wide = TRUE)
      } else {
        as.matrix(slot(x@substances, substance)@rl_xxe, wide = TRUE)
      },
      nm_rhy     = as.matrix(x@transport@rl_rhy, wide = TRUE),
      nm_slp     = as.matrix(x@topo@rl_slp_cap, wide = TRUE)
    )

    slot(x@substances, substance)@rl_xxr     <- writeLayer(x, "xxr"    , rast(li_tpt$nm_xxr    , crs = x@helpers@cs_cmt, extent = x@helpers@ex_cmt), "FLT8S", substance)
    slot(x@substances, substance)@rl_xxt_inp <- writeLayer(x, "xxt_inp", rast(li_tpt$nm_xxt_inp, crs = x@helpers@cs_cmt, extent = x@helpers@ex_cmt), "FLT8S", substance)
    slot(x@substances, substance)@rl_xxt_out <- writeLayer(x, "xxt_out", rast(li_tpt$nm_xxt_out, crs = x@helpers@cs_cmt, extent = x@helpers@ex_cmt), "FLT8S", substance)
    slot(x@substances, substance)@rl_xxt     <- writeLayer(x, "xxt"    , rast(li_tpt$nm_xxt    , crs = x@helpers@cs_cmt, extent = x@helpers@ex_cmt), "FLT8S", substance)
    slot(x@substances, substance)@rl_xxt_ctf <- writeLayer(x, "xxt_ctf", rast(li_tpt$nm_xxt_ctf, crs = x@helpers@cs_cmt, extent = x@helpers@ex_cmt), "FLT8S", substance)
    slot(x@substances, substance)@rl_xxt_cld <- writeLayer(x, "xxt_cld", rast(li_tpt$nm_xxt_cld, crs = x@helpers@cs_cmt, extent = x@helpers@ex_cmt), "FLT8S", substance)

    x
  }
)

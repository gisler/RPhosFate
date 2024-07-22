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
#' summer conditions and slopes \eqn{\geq}{≥} 15 ft) to disk.
#'
#' @param x An S4 [`RPhosFate-class`] river catchment object.
#'
#' @return An S4 [`RPhosFate-class`] river catchment object and side effects in
#'   the form of raster files.
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
    compareGeom(
      x@topo@rl_acc_inf,
      x@topo@rl_dir_inf,
      x@topo@rl_slp_inf,
      x@topo@rl_cha
    )

    cs_dir_old <- setwd(file.path(x@cv_dir[1L], "Intermediate"))
    on.exit(setwd(cs_dir_old))

    # Capped slope in % (also relevant for transport)
    x@topo@rl_slp_cap <- x@topo@rl_slp_inf
    x@topo@rl_slp_cap[x@topo@rl_slp_cap < x@parameters@ns_slp_min] <- x@parameters@ns_slp_min
    x@topo@rl_slp_cap[x@topo@rl_slp_cap > x@parameters@ns_slp_max] <- x@parameters@ns_slp_max

    # Capped slope in radian
    rl_slp_cap_rad <- app(
      x@topo@rl_slp_cap,
      function(x) {
        atan(x * 1e-2)
      },
      cores = x@is_ths
    )

    # Overland flow accumulation
    rl_acc_inf_ovl <- x@topo@rl_acc_inf
    rl_acc_inf_ovl[!is.na(x@topo@rl_cha)] <- NA_real_

    # L factor
    ns_siz <- x@helpers@ns_siz
    ns_res <- x@helpers@ns_res

    x@erosion@rl_LFa <- lapp(
      c(x = rl_acc_inf_ovl, y = rl_slp_cap_rad, z = x@topo@rl_dir_inf),
      fun = function(x, y, z) {
        # Ratio of rill to interrill erosion
        b <- sin(y) / (0.0896 * (3 * sin(y)^0.8 + 0.56))
        # Rill erodibility parameter
        m <- b / (1 + b)

        z_rad = z * pi / 180

        ((x * ns_siz)^(m + 1) - ((x - 1) * ns_siz)^(m + 1)) /
          (ns_res^(m + 2) * (abs(sin(z_rad)) + abs(cos(z_rad)))^m * 22.13^m)
      },
      cores = x@is_ths
    )

    # S factor
    x@erosion@rl_SFa <- lapp(
      c(x = rl_slp_cap_rad, y = x@topo@rl_slp_cap),
      fun = function(x, y) {
        ifelse(
          y < 9,
          10.8 * sin(x) + 0.03,
          16.8 * sin(x) - 0.5
        )
      },
      cores = x@is_ths
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
      x@topo@rl_acc_inf,
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
      x@helpers@ns_siz * 1e-4

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
      x@topo@rl_acc_inf,
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
      },
      cores = x@is_ths
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
#' Determines cells representing inlets as well as riparian zones before writing
#' them to disk.
#'
#' @inheritParams erosionPrerequisites,RPhosFate-method
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
      x@topo@rl_acc_inf,
      x@topo@rl_dir_inf,
      x@topo@rl_cha,
      x@topo@rl_rds
    )

    cs_dir_old <- setwd(file.path(x@cv_dir[1L], "Intermediate"))
    on.exit(setwd(cs_dir_old))

    # Riparian zone and inlet cells
    li_rip_inl <- ripInlCpp(
      nm_dir_inf = as.matrix(x@topo@rl_dir_inf, wide = TRUE),
      im_cha = as.matrix(x@topo@rl_cha, wide = TRUE),
      im_rds = as.matrix(x@topo@rl_rds, wide = TRUE),
      is_ths = x@is_ths
    )

    x@topo@rl_rip <- rast(li_rip_inl$im_rip, crs = x@helpers@cs_cmt, extent = x@helpers@ex_cmt)
    x@topo@rl_inl <- rast(li_rip_inl$im_inl, crs = x@helpers@cs_cmt, extent = x@helpers@ex_cmt)
    set.names(x@topo@rl_inl, "inl")

    # Nearest channel cells for inlet cells
    df_out <- findNearestNeighbour(
      sort(as.points(x@topo@rl_inl), "inl"),
      as.points(x@topo@rl_cha)
    )

    # X-coordinates of nearest channel cells to column numbers
    df_out$x <- colFromX(x@topo@rl_inl, df_out$x)
    # Y-coordinates of nearest channel cells to row numbers
    df_out$y <- rowFromY(x@topo@rl_inl, df_out$y)

    # Substituting inlet values with integer codes identifying nearest channel
    # cells
    df_out$code <- as.integer(df_out$y * x@helpers@is_cls + df_out$x)
    x@topo@rl_inl <- subst(x@topo@rl_inl, df_out[["from_id"]], df_out[["code"]])

    x@topo@rl_rip      <- writeLayer(x, "rip", x@topo@rl_rip     , "INT4S")
    x@topo@rl_inl      <- writeLayer(x, "inl", x@topo@rl_inl     , "INT4S")

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
#' \cite{Molnár, P., Ramírez, J.A., 1998. Energy dissipation theories and
#' optimal channel characteristics of river networks. Water Resources Research
#' 34, 1809–1818. https://doi.org/10.1029/98WR00983}
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
      x@topo@rl_acc_inf,
      x@topo@rl_dir_inf,
      x@topo@rl_slp_cap,
      x@transport@rl_man,
      if (substance == "SS") {
        x@erosion@rl_ero
      } else {
        slot(x@substances, substance)@rl_xxe
      },
      x@topo@rl_cha,
      x@topo@rl_rds,
      x@topo@rl_rip,
      x@topo@rl_inl
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
      nm_acc_inf = as.matrix(x@topo@rl_acc_inf, wide = TRUE),
      nm_dir_inf = as.matrix(x@topo@rl_dir_inf, wide = TRUE),
      nm_slp_cap = as.matrix(x@topo@rl_slp_cap, wide = TRUE),
      nm_man = as.matrix(x@transport@rl_man, wide = TRUE),
      nm_xxe = if (substance == "SS") {
        as.matrix(x@erosion@rl_ero, wide = TRUE)
      } else {
        as.matrix(slot(x@substances, substance)@rl_xxe, wide = TRUE)
      },
      im_cha = as.matrix(x@topo@rl_cha, wide = TRUE),
      im_rds = as.matrix(x@topo@rl_rds, wide = TRUE),
      im_rip = as.matrix(x@topo@rl_rip, wide = TRUE),
      im_inl = as.matrix(x@topo@rl_inl, wide = TRUE),
      substance = substance,
      parameters = x@parameters,
      helpers = x@helpers,
      is_ths = x@is_ths
    )

    slot(x@substances, substance)@rl_xxr     <- writeLayer(x, "xxr"    , rast(li_tpt$nm_xxr    , crs = x@helpers@cs_cmt, extent = x@helpers@ex_cmt), "FLT8S", substance)
    slot(x@substances, substance)@rl_xxt     <- writeLayer(x, "xxt"    , rast(li_tpt$nm_xxt    , crs = x@helpers@cs_cmt, extent = x@helpers@ex_cmt), "FLT8S", substance)
    slot(x@substances, substance)@rl_xxt_inp <- writeLayer(x, "xxt_inp", rast(li_tpt$nm_xxt_inp, crs = x@helpers@cs_cmt, extent = x@helpers@ex_cmt), "FLT8S", substance)
    slot(x@substances, substance)@rl_xxt_out <- writeLayer(x, "xxt_out", rast(li_tpt$nm_xxt_out, crs = x@helpers@cs_cmt, extent = x@helpers@ex_cmt), "FLT8S", substance)
    slot(x@substances, substance)@rl_xxt_ctf <- writeLayer(x, "xxt_ctf", rast(li_tpt$nm_xxt_ctf, crs = x@helpers@cs_cmt, extent = x@helpers@ex_cmt), "FLT8S", substance)
    slot(x@substances, substance)@rl_xxt_cld <- writeLayer(x, "xxt_cld", rast(li_tpt$nm_xxt_cld, crs = x@helpers@cs_cmt, extent = x@helpers@ex_cmt), "FLT8S", substance)

    writeRaster(
      rast(li_tpt$im_ifl, crs = x@helpers@cs_cmt, extent = x@helpers@ex_cmt),
      filename = "ifl.tif",
      datatype = "INT4S",
      overwrite = TRUE
    )
    # writeRaster(
    #   rast(li_tpt$im_ord, crs = x@helpers@cs_cmt, extent = x@helpers@ex_cmt),
    #   filename = "ord.tif",
    #   datatype = "INT4S",
    #   overwrite = TRUE
    # )

    x
  }
)

#' @include aaa.R
NULL

#### erosionPrerequisites ####
setGeneric(
  "erosionPrerequisites",
  function(cmt, ...) standardGeneric("erosionPrerequisites")
)
#' @export
setMethod(
  "erosionPrerequisites",
  "RPhosFate",
  function(cmt) {
    cs_dir_old <- setwd(file.path(cmt@cv_dir[1L], "Intermediate"))
    on.exit(setwd(cs_dir_old))

    # Capped slope in % (also relevant for transport)
    cmt@topo@rl_slp_cap <- cmt@topo@rl_slp
    cmt@topo@rl_slp_cap[cmt@topo@rl_slp_cap < cmt@parameters@ns_slp_min] <- cmt@parameters@ns_slp_min
    cmt@topo@rl_slp_cap[cmt@topo@rl_slp_cap > cmt@parameters@ns_slp_max] <- cmt@parameters@ns_slp_max

    # Capped slope in radian
    rl_slp_cap_rad <- calc(
      cmt@topo@rl_slp_cap,
      function(x) {
        atan(x / 1e2)
      }
    )

    # Weighted overland flow accumulation
    rl_acc_wtd_ovl <- cmt@topo@rl_acc_wtd
    rl_acc_wtd_ovl[!is.na(cmt@topo@rl_cha)] <- NA_real_
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
    cmt@erosion@rl_LFa <- overlay(
      x = rl_acc_wtd_ovl,
      y = rl_LFa_m,
      fun = function(x, y) {
        ((x * cmt@helper@is_res)^(1 + y) -
          ((x - 1) * cmt@helper@is_res)^(1 + y)) / # nolint
          (cmt@helper@is_res * 22.13^y)
      },
      filename = "LFa.img",
      datatype = "FLT4S",
      options = "COMPRESSED=YES",
      overwrite = TRUE
    )

    # S factor
    cmt@erosion@rl_SFa <- overlay(
      x = rl_slp_cap_rad,
      y = cmt@topo@rl_slp_cap,
      fun = function(x, y) {
        ifelse(y < 9, 10.8 * sin(x) + 0.03, 16.8 * sin(x) - 0.5)
      },
      filename = "SFa.img",
      datatype = "FLT4S",
      options = "COMPRESSED=YES",
      overwrite = TRUE
    )

    writeRaster(
      cmt@topo@rl_slp_cap,
      "slp_cap.img",
      datatype = "FLT4S",
      options = "COMPRESSED=YES",
      overwrite = TRUE
    )
    cmt@topo@rl_slp_cap <- raster("slp_cap.img")

    cmt
  }
)

#### erosion ####
setGeneric(
  "erosion",
  function(cmt, ...) standardGeneric("erosion")
)
#' @export
setMethod(
  "erosion",
  "RPhosFate",
  function(cmt) {
    cs_dir_old <- setwd(file.path(cmt@cv_dir[1L], "Result"))
    on.exit(setwd(cs_dir_old))

    # Erosion in t/cell/yr
    cmt@erosion@rl_ero <-
      cmt@erosion@rl_RFa *
      cmt@erosion@rl_KFa *
      cmt@erosion@rl_LFa *
      cmt@erosion@rl_SFa *
      cmt@erosion@rl_CFa *
      cmt@helper@is_siz / 1e4

    filename <- paste0("ero", cmt@is_MCi, ".img")
    writeRaster(
      cmt@erosion@rl_ero,
      filename,
      datatype = "FLT4S",
      options = "COMPRESSED=YES",
      overwrite = TRUE
    )
    cmt@erosion@rl_ero <- raster(filename)

    cmt
  }
)

#### emission ####
setGeneric(
  "emission",
  function(cmt, ...) standardGeneric("emission")
)
#' @export
setMethod(
  "emission",
  "RPhosFate",
  function(cmt, substance = "PP") {
    assertSubstance(cmt, substance)

    cs_dir_old <- setwd(file.path(cmt@cv_dir[1L], "Result"))
    on.exit(setwd(cs_dir_old))

    filename <- paste0(tolower(substance), "e", cmt@is_MCi, ".img")
    # Emission in kg/cell/yr
    slot(cmt@substance, substance)@rl_xxe <- overlay(
      x = cmt@erosion@rl_ero,
      y = slot(cmt@substance, substance)@rl_xxc,
      z = cmt@topo@rl_clc,
      fun = function(x, y, z) {
        x * y * (1 + z / 1e2) * 1e-3
      },
      filename = filename,
      datatype = "FLT4S",
      options = "COMPRESSED=YES",
      overwrite = TRUE
    )

    cmt
  }
)

#### transportPrerequisites ####
setGeneric(
  "transportPrerequisites",
  function(cmt, ...) standardGeneric("transportPrerequisites")
)
#' @export
setMethod(
  "transportPrerequisites",
  "RPhosFate",
  function(cmt) {
    cs_dir_old <- setwd(file.path(cmt@cv_dir[1L], "Intermediate"))
    on.exit(setwd(cs_dir_old))

    # Hydraulic radius in m
    cmt@transport@rl_rhy <- calc(
      cmt@topo@rl_acc_wtd,
      function(x) {
        cmt@parameters@ns_rhy_a *
          (x * cmt@helper@is_siz / 1e6)^cmt@parameters@ns_rhy_b
      },
      filename = "rhy.img",
      datatype = "FLT4S",
      options = "COMPRESSED=YES",
      overwrite = TRUE
    )

    # Riparian zone cells
    cmt@topo@rl_rip <- raster(
      dir_sth(
        im_dir = as.matrix(cmt@topo@rl_dir),
        im_sth = as.matrix(cmt@topo@rl_cha),
        im_fDo = cmt@helper@im_fDo
      ),
      template = cmt@topo@rl_acc_wtd
    )

    # Inlet cells
    cmt@topo@rl_inl <- raster(
      dir_sth(
        im_dir = as.matrix(cmt@topo@rl_dir),
        im_sth = as.matrix(cmt@topo@rl_rds),
        im_fDo = cmt@helper@im_fDo
      ),
      template = cmt@topo@rl_acc_wtd
    )

    # No inlet cells at channel cells
    cmt@topo@rl_inl[!is.na(cmt@topo@rl_cha)] <- NA_integer_
    # No riparian zone cells at road cells
    cmt@topo@rl_rip[!is.na(cmt@topo@rl_rds)] <- NA_integer_
    # No inlet cells at riparian zone cells
    cmt@topo@rl_inl[!is.na(cmt@topo@rl_rip)] <- NA_integer_

    # Nearest channel cells for inlet cells
    df_out <- findNearestNeighbour(
      rasterToPoints(cmt@topo@rl_inl),
      rasterToPoints(cmt@topo@rl_cha),
      cmt@helper@ex_cmt
    )

    # X-coordinates of nearest channel cells to column numbers
    df_out$Y.x <- (df_out$Y.x + cmt@helper@is_res / 2 -
      cmt@helper@ex_cmt[1L]) / cmt@helper@is_res
    # Y-coordinates of nearest channel cells to row numbers
    df_out$Y.y <- cmt@helper@is_rws - ((df_out$Y.y - cmt@helper@is_res / 2 -
      cmt@helper@ex_cmt[3L]) / cmt@helper@is_res)

    # Substituting inlet values with integer codes identifying nearest channel
    # cells
    df_out$code <- as.integer(df_out$Y.y * cmt@helper@is_cls + df_out$Y.x)
    # Bug in subs() {raster}: use default by and which
    cmt@topo@rl_inl <- subs(cmt@topo@rl_inl, y = df_out[c(3L, 8L)])

    writeRaster(cmt@topo@rl_rip, "rip.img", datatype = "INT4S", options = "COMPRESSED=YES", overwrite = TRUE)
    writeRaster(cmt@topo@rl_inl, "inl.img", datatype = "INT4S", options = "COMPRESSED=YES", overwrite = TRUE)
    cmt@topo@rl_rip <- raster("rip.img")
    cmt@topo@rl_inl <- raster("inl.img")

    cmt
  }
)

#### transportCalcOrder ####
setGeneric(
  "transportCalcOrder",
  function(cmt, ...) standardGeneric("transportCalcOrder")
)
#' @export
setMethod(
  "transportCalcOrder",
  "RPhosFate",
  function(cmt) {
    # Overland flow accumulation
    rl_acc_ovl <- cmt@topo@rl_acc
    rl_acc_ovl[!is.na(cmt@topo@rl_cha)] <- NA_integer_

    # Channel flow accumulation
    rl_acc_cha <- cmt@topo@rl_acc
    rl_acc_cha[ is.na(cmt@topo@rl_cha)] <- NA_integer_

    # Transport calculation order as column-major index
    im_acc_ovl <- as.matrix(rl_acc_ovl)
    im_acc_cha <- as.matrix(rl_acc_cha)
    ar_ord_ovl <- tapply(seq_along(im_acc_ovl), im_acc_ovl, identity)
    ar_ord_cha <- tapply(seq_along(im_acc_cha), im_acc_cha, identity)

    # Row order from index
    iv_ord_ovl_row <- as.integer(unlist(lapply(
      ar_ord_ovl,
      function(x) {
        (x + cmt@helper@is_rws) - ceiling(x / cmt@helper@is_rws) *
          cmt@helper@is_rws
      }
    )))
    iv_ord_cha_row <- as.integer(unlist(lapply(
      ar_ord_cha,
      function(x) {
        (x + cmt@helper@is_rws) - ceiling(x / cmt@helper@is_rws) *
          cmt@helper@is_rws
      }
    )))

    # Column order from index
    iv_ord_ovl_col <- as.integer(unlist(lapply(
      ar_ord_ovl,
      function(x) {
        ceiling(x / cmt@helper@is_rws)
      }
    )))
    iv_ord_cha_col <- as.integer(unlist(lapply(
      ar_ord_cha,
      function(x) {
        ceiling(x / cmt@helper@is_rws)
      }
    )))

    # Overland as well as channel row and column numbers for top-down
    # computation (C++ has zero-based numbering)
    cmt@helper@order@iv_ord_row <- c(iv_ord_ovl_row, iv_ord_cha_row) - 1L
    cmt@helper@order@iv_ord_col <- c(iv_ord_ovl_col, iv_ord_cha_col) - 1L

    # Reverse overland row and column numbers for bottom-up computation
    # (C++ has zero-based numbering)
    cmt@helper@order@iv_ord_ovl_row_rev <- rev(iv_ord_ovl_row) - 1L
    cmt@helper@order@iv_ord_ovl_col_rev <- rev(iv_ord_ovl_col) - 1L

    cmt
  }
)

#### transport ####
setGeneric(
  "transport",
  function(cmt, ...) standardGeneric("transport")
)
#' @export
setMethod(
  "transport",
  "RPhosFate",
  function(cmt, substance = "PP") {
    assertSubstance(cmt, substance)

    cs_dir_old <- setwd(file.path(cmt@cv_dir[1L], "Result"))
    on.exit(setwd(cs_dir_old))

    li_tpt <- transportCpp(
      parameters = cmt@parameters,
      ns_dep_ovl = if (substance == "SS") {
        cmt@parameters@ns_dep_ovl
      } else {
        cmt@parameters@ns_dep_ovl / cmt@parameters@nv_enr_rto[substance]
      },
      ns_tfc_inl = cmt@parameters@nv_tfc_inl[substance],
      helper     = cmt@helper,
      order      = cmt@helper@order,
      im_cha     = as.matrix(cmt@topo@rl_cha),
      im_dir     = as.matrix(cmt@topo@rl_dir),
      im_inl     = as.matrix(cmt@topo@rl_inl),
      im_rip     = as.matrix(cmt@topo@rl_rip),
      nm_man     = as.matrix(cmt@transport@rl_man),
      nm_xxe     = if (substance == "SS") {
        as.matrix(cmt@erosion@rl_ero)
      } else {
        as.matrix(slot(cmt@substance, substance)@rl_xxe)
      },
      nm_rhy     = as.matrix(cmt@transport@rl_rhy),
      nm_slp     = as.matrix(cmt@topo@rl_slp_cap)
    )

    layers <- paste0("xx", c("r", "t_inp", "t_out", "t_cld", "t_ctf", "t"))
    filenames <- setNames(
      paste0(sub("^xx", tolower(substance), layers), cmt@is_MCi, ".img"),
      layers
    )

    if (length(cmt@is_MCi) == 0L) {
      writeRaster(raster(li_tpt$nm_xxr    , template = cmt@topo@rl_acc_wtd), filenames["xxr"    ], datatype = "FLT4S", options = "COMPRESSED=YES", overwrite = TRUE)
      writeRaster(raster(li_tpt$nm_xxt_inp, template = cmt@topo@rl_acc_wtd), filenames["xxt_inp"], datatype = "FLT4S", options = "COMPRESSED=YES", overwrite = TRUE)
      writeRaster(raster(li_tpt$nm_xxt_out, template = cmt@topo@rl_acc_wtd), filenames["xxt_out"], datatype = "FLT4S", options = "COMPRESSED=YES", overwrite = TRUE)
      writeRaster(raster(li_tpt$nm_xxt_cld, template = cmt@topo@rl_acc_wtd), filenames["xxt_cld"], datatype = "FLT4S", options = "COMPRESSED=YES", overwrite = TRUE)
      writeRaster(raster(li_tpt$nm_xxt_ctf, template = cmt@topo@rl_acc_wtd), filenames["xxt_ctf"], datatype = "FLT4S", options = "COMPRESSED=YES", overwrite = TRUE)
      slot(cmt@substance, substance)@rl_xxr     <- raster(filenames["xxr"    ])
      slot(cmt@substance, substance)@rl_xxt_inp <- raster(filenames["xxt_inp"])
      slot(cmt@substance, substance)@rl_xxt_out <- raster(filenames["xxt_out"])
      slot(cmt@substance, substance)@rl_xxt_cld <- raster(filenames["xxt_cld"])
      slot(cmt@substance, substance)@rl_xxt_ctf <- raster(filenames["xxt_ctf"])
    }
    writeRaster(
      raster(li_tpt$nm_xxt, template = cmt@topo@rl_acc_wtd),
      filenames["xxt"],
      datatype = "FLT4S",
      options = "COMPRESSED=YES",
      overwrite = TRUE
    )
    slot(cmt@substance, substance)@rl_xxt <- raster(filenames["xxt"])

    cmt
  }
)

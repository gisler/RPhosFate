#' @include aaa.R
NULL

#### erosionPrerequisites ####
setGeneric(
  name = "erosionPrerequisites",
  def = function(cmt) {standardGeneric("erosionPrerequisites")}
)
#' @export
setMethod(
  f = "erosionPrerequisites",
  signature = "RPhosFate",
  definition = function(cmt) {
    cs_dir_old <- setwd(file.path(cmt@cv_dir[1], "Intermediate"))
    on.exit(setwd(cs_dir_old))

    # Capped slope in % (also relevant for transport() {RPhosFate})
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

    # Overland weighted flow accumulation
    rl_acc_wtd_ovl <- cmt@topo@rl_acc_wtd
    rl_acc_wtd_ovl[!is.na(cmt@topo@rl_cha)] <- NA
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

    # L-factor
    cmt@erosion@rl_LFa <- overlay(
      x = rl_acc_wtd_ovl,
      y = rl_LFa_m,
      fun = function(x, y) {
        ((x * cmt@helper@is_res)^(1 + y) - ((x - 1) * cmt@helper@is_res)^(1 + y)) / # nolint
          (cmt@helper@is_res * 22.13^y)
      }
    )

    # S-factor
    cmt@erosion@rl_SFa <- overlay(
      x = rl_slp_cap_rad,
      y = cmt@topo@rl_slp_cap,
      fun = function(x, y) {
        ifelse(y < 9, 10.8 * sin(x) + 0.03, 16.8 * sin(x) - 0.5)
      }
    )

    writeRaster(cmt@topo@rl_slp_cap, "slp_cap.img", datatype = "FLT4S", options = "COMPRESSED=YES", overwrite = TRUE)
    writeRaster(cmt@erosion@rl_LFa , "LFa.img"    , datatype = "FLT4S", options = "COMPRESSED=YES", overwrite = TRUE)
    writeRaster(cmt@erosion@rl_SFa , "SFa.img"    , datatype = "FLT4S", options = "COMPRESSED=YES", overwrite = TRUE)
    cmt@topo@rl_slp_cap <- raster("slp_cap.img")
    cmt@erosion@rl_LFa  <- raster("LFa.img")
    cmt@erosion@rl_SFa  <- raster("SFa.img")

    return(cmt)
  }
)

#### erosion ####
setGeneric(
  name = "erosion",
  def = function(cmt) {standardGeneric("erosion")}
)
#' @export
setMethod(
  f = "erosion",
  signature = "RPhosFate",
  definition = function(cmt) {
    cs_dir_old <- setwd(file.path(cmt@cv_dir[1], "Result"))
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

    return(cmt)
  }
)

#### emission ####
setGeneric(
  name = "emission",
  def = function(cmt, ety) {standardGeneric("emission")}
)
#' @export
setMethod(
  f = "emission",
  signature = c("RPhosFate", "RPhosFatePP"),
  definition = function(cmt, ety) {
    cs_dir_old <- setwd(file.path(cmt@cv_dir[1], "Result"))
    on.exit(setwd(cs_dir_old))

    # PP emission in kg P/cell/yr
    cmt@PP@rl_ppe <- overlay(
      x = cmt@erosion@rl_ero,
      y = cmt@PP@rl_ppc,
      z = cmt@topo@rl_clc,
      fun = function(x, y, z) {
        x * y * (1 + z / 1e2) * 1e-3
      }
    )

    filename <- paste0("ppe", cmt@is_MCi, ".img")
    writeRaster(
      cmt@PP@rl_ppe,
      filename,
      datatype = "FLT4S",
      options = "COMPRESSED=YES",
      overwrite = TRUE
    )
    cmt@PP@rl_ppe <- raster(filename)

    return(cmt)
  }
)

#### transportPrerequisites ####
setGeneric(
  name = "transportPrerequisites",
  def = function(cmt) {standardGeneric("transportPrerequisites")}
)
#' @export
setMethod(
  f = "transportPrerequisites",
  signature = "RPhosFate",
  definition = function(cmt) {
    cs_dir_old <- setwd(file.path(cmt@cv_dir[1], "Intermediate"))
    on.exit(setwd(cs_dir_old))

    # Hydraulic radius in m
    cmt@transport@rl_rhy <- calc(
      cmt@topo@rl_acc_wtd,
      function(x) {
        cmt@parameters@ns_rhy_a * (x * cmt@helper@is_siz / 1e6)^cmt@parameters@ns_rhy_b
      }
    )

    # Riparian zone cells
    cmt@topo@rl_rip <- raster(
      dir_sth(
        im_dir = as.matrix(cmt@topo@rl_dir),
        im_sth = as.matrix(cmt@topo@rl_cha),
        im_fDo = cmt@helper@im_fDo
      ),
      template = cmt@topo@rl_dir
    )

    # Inlet cells
    cmt@topo@rl_inl <- raster(
      dir_sth(
        im_dir = as.matrix(cmt@topo@rl_dir),
        im_sth = as.matrix(cmt@topo@rl_rds),
        im_fDo = cmt@helper@im_fDo
      ),
      template = cmt@topo@rl_dir
    )

    # No inlet cells at channel cells
    cmt@topo@rl_inl[!is.na(cmt@topo@rl_cha)] <- NA

    # No riparian zone cells at road cells
    cmt@topo@rl_rip[!is.na(cmt@topo@rl_rds)] <- NA

    # No inlet cells at riparian zone cells
    cmt@topo@rl_inl[!is.na(cmt@topo@rl_rip)] <- NA

    # Nearest channel cells for inlet cells
    df_out <- findNearestNeighbour(
      rasterToPoints(cmt@topo@rl_inl),
      rasterToPoints(cmt@topo@rl_cha),
      cmt@helper@ex_cmt
    )

    # X-coordinates of nearest channel cells to column numbers
    df_out$Y.x <- (df_out$Y.x + cmt@helper@is_res / 2 - cmt@helper@ex_cmt[1]) / cmt@helper@is_res

    # Y-coordinates of nearest channel cells to row numbers
    df_out$Y.y <- cmt@helper@is_rws - ((df_out$Y.y - cmt@helper@is_res / 2 - cmt@helper@ex_cmt[3]) /
      cmt@helper@is_res)

    # Substituting inlet values with integer codes identifying nearest channel cells
    df_out$code <- as.integer(df_out$Y.y * cmt@helper@is_cls + df_out$Y.x)
    # Bug in subs() {raster}: use default by and which
    cmt@topo@rl_inl <- subs(cmt@topo@rl_inl, y = df_out[c(3, 8)])

    writeRaster(cmt@transport@rl_rhy, "rhy.img", datatype = "FLT4S", options = "COMPRESSED=YES", overwrite = TRUE)
    writeRaster(cmt@topo@rl_rip     , "rip.img", datatype = "INT4S", options = "COMPRESSED=YES", overwrite = TRUE)
    writeRaster(cmt@topo@rl_inl     , "inl.img", datatype = "INT4S", options = "COMPRESSED=YES", overwrite = TRUE)
    cmt@transport@rl_rhy <- raster("rhy.img")
    cmt@topo@rl_rip      <- raster("rip.img")
    cmt@topo@rl_inl      <- raster("inl.img")

    return(cmt)
  }
)

#### transportCalcOrder ####
setGeneric(
  name = "transportCalcOrder",
  def = function(cmt) {standardGeneric("transportCalcOrder")}
)
#' @export
setMethod(
  f = "transportCalcOrder",
  signature = "RPhosFate",
  definition = function(cmt) {
    # Overland flow accumulation
    rl_acc_ovl <- cmt@topo@rl_acc
    rl_acc_ovl[!is.na(cmt@topo@rl_cha)] <- NA

    # Channel flow accumulation
    rl_acc_cha <- cmt@topo@rl_acc
    rl_acc_cha[is.na(cmt@topo@rl_cha)] <- NA

    # Transport calculation order as column-major index
    im_acc_ovl <- as.matrix(rl_acc_ovl)
    im_acc_cha <- as.matrix(rl_acc_cha)
    ar_ord_ovl <- tapply(seq_along(im_acc_ovl), im_acc_ovl, identity)
    ar_ord_cha <- tapply(seq_along(im_acc_cha), im_acc_cha, identity)

    # Row numbers from index
    iv_ord_ovl_row <- as.integer(unlist(lapply(
      ar_ord_ovl,
      function(x) {
        (x + cmt@helper@is_rws) - ceiling(x / cmt@helper@is_rws) * cmt@helper@is_rws
      }
    )))
    iv_ord_cha_row <- as.integer(unlist(lapply(
      ar_ord_cha,
      function(x) {
        (x + cmt@helper@is_rws) - ceiling(x / cmt@helper@is_rws) * cmt@helper@is_rws
      }
    )))

    # Column numbers from index
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

    # Reverse overland row and column numbers for bottom-up calculations (C++ is zero-based)
    cmt@helper@order@iv_ord_ovl_row_rev <- rev(iv_ord_ovl_row) - 1L
    cmt@helper@order@iv_ord_ovl_col_rev <- rev(iv_ord_ovl_col) - 1L

    # Overland and channel row and column numbers for top-down calculations (C++ is zero-based)
    cmt@helper@order@iv_ord_row <- c(iv_ord_ovl_row, iv_ord_cha_row) - 1L
    cmt@helper@order@iv_ord_col <- c(iv_ord_ovl_col, iv_ord_cha_col) - 1L

    return(cmt)
  }
)

#### transport ####
setGeneric(
  name = "transport",
  def = function(cmt, tty) {standardGeneric("transport")}
)
#' @export
setMethod(
  f = "transport",
  signature = c("RPhosFate", "RPhosFateSS"),
  definition = function(cmt, tty) {
    cs_dir_old <- setwd(file.path(cmt@cv_dir[1], "Result"))
    on.exit(setwd(cs_dir_old))

    parameters <- cmt@parameters
    parameters@nv_tfc_inl <- parameters@nv_tfc_inl["SS"]
    li_tpt <- transportCpp(
      parameters = parameters,
      helper     = cmt@helper,
      order      = cmt@helper@order,
      im_cha     = as.matrix(cmt@topo@rl_cha),
      im_dir     = as.matrix(cmt@topo@rl_dir),
      im_inl     = as.matrix(cmt@topo@rl_inl),
      im_rip     = as.matrix(cmt@topo@rl_rip),
      nm_man     = as.matrix(cmt@transport@rl_man),
      nm_xxe     = as.matrix(cmt@erosion@rl_ero),
      nm_rhy     = as.matrix(cmt@transport@rl_rhy),
      nm_slp     = as.matrix(cmt@topo@rl_slp_cap)
    )

    if (length(cmt@is_MCi) < 1) {
      writeRaster(raster(li_tpt$nm_xxr    , template = cmt@topo@rl_dir), "ssr.img"    , datatype = "FLT4S", options = "COMPRESSED=YES", overwrite = TRUE)
      writeRaster(raster(li_tpt$nm_xxt_inp, template = cmt@topo@rl_dir), "sst_inp.img", datatype = "FLT4S", options = "COMPRESSED=YES", overwrite = TRUE)
      writeRaster(raster(li_tpt$nm_xxt_out, template = cmt@topo@rl_dir), "sst_out.img", datatype = "FLT4S", options = "COMPRESSED=YES", overwrite = TRUE)
      writeRaster(raster(li_tpt$nm_xxt_cld, template = cmt@topo@rl_dir), "sst_cld.img", datatype = "FLT4S", options = "COMPRESSED=YES", overwrite = TRUE)
      writeRaster(raster(li_tpt$nm_xxt_ctf, template = cmt@topo@rl_dir), "sst_ctf.img", datatype = "FLT4S", options = "COMPRESSED=YES", overwrite = TRUE)
      cmt@SS@rl_ssr     <- raster("ssr.img")
      cmt@SS@rl_sst_inp <- raster("sst_inp.img")
      cmt@SS@rl_sst_out <- raster("sst_out.img")
      cmt@SS@rl_sst_cld <- raster("sst_cld.img")
      cmt@SS@rl_sst_ctf <- raster("sst_ctf.img")
    }
    filename <- paste0("sst", cmt@is_MCi, ".img")
    writeRaster(raster(li_tpt$nm_xxt, template = cmt@topo@rl_dir), filename, datatype = "FLT4S", options = "COMPRESSED=YES", overwrite = TRUE)
    cmt@SS@rl_sst <- raster(filename)

    return(cmt)
  }
)
setMethod(
  f = "transport",
  signature = c("RPhosFate", "RPhosFatePP"),
  definition = function(cmt, tty) {
    cs_dir_old <- setwd(file.path(cmt@cv_dir[1], "Result"))
    on.exit(setwd(cs_dir_old))

    parameters <- cmt@parameters
    parameters@ns_dep_ovl <- parameters@ns_dep_ovl / parameters@nv_enr_rto["PP"]
    parameters@nv_tfc_inl <- parameters@nv_tfc_inl["PP"]
    li_tpt <- transportCpp(
      parameters = parameters,
      helper     = cmt@helper,
      order      = cmt@helper@order,
      im_cha     = as.matrix(cmt@topo@rl_cha),
      im_dir     = as.matrix(cmt@topo@rl_dir),
      im_inl     = as.matrix(cmt@topo@rl_inl),
      im_rip     = as.matrix(cmt@topo@rl_rip),
      nm_man     = as.matrix(cmt@transport@rl_man),
      nm_xxe     = as.matrix(cmt@PP@rl_ppe),
      nm_rhy     = as.matrix(cmt@transport@rl_rhy),
      nm_slp     = as.matrix(cmt@topo@rl_slp_cap)
    )

    if (length(cmt@is_MCi) < 1) {
      writeRaster(raster(li_tpt$nm_xxr    , template = cmt@topo@rl_dir), "ppr.img"    , datatype = "FLT4S", options = "COMPRESSED=YES", overwrite = TRUE)
      writeRaster(raster(li_tpt$nm_xxt_inp, template = cmt@topo@rl_dir), "ppt_inp.img", datatype = "FLT4S", options = "COMPRESSED=YES", overwrite = TRUE)
      writeRaster(raster(li_tpt$nm_xxt_out, template = cmt@topo@rl_dir), "ppt_out.img", datatype = "FLT4S", options = "COMPRESSED=YES", overwrite = TRUE)
      writeRaster(raster(li_tpt$nm_xxt_cld, template = cmt@topo@rl_dir), "ppt_cld.img", datatype = "FLT4S", options = "COMPRESSED=YES", overwrite = TRUE)
      writeRaster(raster(li_tpt$nm_xxt_ctf, template = cmt@topo@rl_dir), "ppt_ctf.img", datatype = "FLT4S", options = "COMPRESSED=YES", overwrite = TRUE)
      cmt@PP@rl_ppr     <- raster("ppr.img")
      cmt@PP@rl_ppt_inp <- raster("ppt_inp.img")
      cmt@PP@rl_ppt_out <- raster("ppt_out.img")
      cmt@PP@rl_ppt_cld <- raster("ppt_cld.img")
      cmt@PP@rl_ppt_ctf <- raster("ppt_ctf.img")
    }
    filename <- paste0("ppt", cmt@is_MCi, ".img")
    writeRaster(raster(li_tpt$nm_xxt, template = cmt@topo@rl_dir), filename, datatype = "FLT4S", options = "COMPRESSED=YES", overwrite = TRUE)
    cmt@PP@rl_ppt <- raster(filename)

    return(cmt)
  }
)

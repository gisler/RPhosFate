populateLayerSlots <- function(
  cmt,
  .Object,
  slots,
  layers,
  areRequiredInputLayers = rep(FALSE, length(slots)),
  areMCinputLayers = rep(FALSE, length(slots))
) {
  cs_dir_old <- setwd(cmt@cv_dir[1L])
  on.exit(setwd(cs_dir_old))

  for (i in seq_along(slots)) {
    slot(.Object, slots[i]) <- readLayer(
      cmt,
      layers[i],
      areRequiredInputLayers[i],
      areMCinputLayers[i]
    )
  }

  .Object
}

readLayer <- function(
  cmt,
  layer,
  isRequiredInputLayer = FALSE,
  isMCinputLayer = FALSE
) {
  if (length(cmt@is_MCi) == 1L && isMCinputLayer) {
    cs_dir_old <- setwd(cmt@cv_dir[2L])
    on.exit(setwd(cs_dir_old))
  }

  filename <- paste0(layer, cmt@is_MCi, ".img")

  if (isRequiredInputLayer) {
    raster(filename)
  } else if (file.exists(filename)) {
    raster(filename)
  } else {
    new("RasterLayer")
  }
}

readParameters <- function(arguments) {
  parameters <- read_yaml("parameters.yaml")

  parameters[["nv_enr_rto"]] <- unlist(parameters[["nv_enr_rto"]])
  parameters[["nv_tfc_inl"]] <- unlist(parameters[["nv_tfc_inl"]])
  parameters[["nm_olc"]] <- matrix(parameters[["nm_olc"]], 1L)
  parameters[["df_cdt"]] <- as.data.frame(parameters[["df_cdt"]])

  modifyList(parameters, arguments)
}

writeParameters <- function(parameters) {
  parameters <- slots2list(parameters)

  parameters[["nv_enr_rto"]] <- as.list(parameters[["nv_enr_rto"]])
  parameters[["nv_tfc_inl"]] <- as.list(parameters[["nv_tfc_inl"]])

  write_yaml(
    c(list(RPhosFate = as.character(packageVersion("RPhosFate"))), parameters),
    "parameters.yaml",
    indent.mapping.sequence = TRUE
  )
}

adjustExtent <- function(rl, ex) {
  extend(crop(rl, ex), ex)
}

calibrate <- function(value, cmt, substance, col, metric, parameter) {
  if (!is.null(parameter)) {
    slot(cmt@parameters, parameter) <- value
  } else if (substance == "SS") {
    cmt@parameters@ns_dep_ovl <- value
  } else {
    cmt@parameters@nv_enr_rto[substance] <- value
  }

  subsequentRun(cmt, substance)
  metrics <- calibrationQuality(cmt, substance, col)

  if (metric == "PBIAS") {
    abs(metrics[metric])
  } else {
    metrics[metric]
  }
}

populateLayerSlots <- function(
  cmt,
  object,
  slots,
  layers,
  areRequiredInputLayers = rep(FALSE, length(slots)),
  areMCinputLayers = rep(FALSE, length(slots))
) {
  cs_dir_old <- setwd(cmt@cv_dir[1L])
  on.exit(setwd(cs_dir_old))

  for (i in seq_along(slots)) {
    slot(object, slots[i]) <- readLayer(
      cmt,
      layers[i],
      areRequiredInputLayers[i],
      areMCinputLayers[i]
    )
  }

  object
}

populateParameterSlots <- function(parameters, arguments) {
  argumentNames <- names(arguments)

  assertSubset(argumentNames, slotNames(parameters))

  for (i in seq_along(arguments)) {
    if (argumentNames[i] %in% c("nv_tfc_inl", "nv_enr_rto")) {
      assertCharacter(
        names(arguments[[i]]),
        min.chars = 1L,
        any.missing = FALSE,
        unique = TRUE,
        .var.name = sprintf("names(%s)", argumentNames[i])
      )

      slot(parameters, argumentNames[i]) <- unlist(modifyList(
        as.list(slot(parameters, argumentNames[i])),
        as.list(arguments[[i]])
      ))
    } else {
      slot(parameters, argumentNames[i]) <- arguments[[i]]
    }
  }
  validObject(parameters)

  parameters
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

  filename <- paste0(layer, cmt@is_MCi, cmt@cs_fex)

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

  parameters[["nv_tfc_inl"]] <- unlist(parameters[["nv_tfc_inl"]])
  parameters[["nv_enr_rto"]] <- unlist(parameters[["nv_enr_rto"]])
  parameters[["nm_olc"]] <- matrix(parameters[["nm_olc"]], 1L)
  parameters[["df_cdt"]] <- as.data.frame(
    parameters[["df_cdt"]],
    stringsAsFactors = FALSE
  )

  modifyList(parameters, arguments)
}

slots2list <- function(parameters) {
  parameterNames <- slotNames(parameters)

  setNames(lapply(
    parameterNames,
    function(name, parameters) {slot(parameters, name)},
    parameters = parameters
  ), parameterNames)
}

writeParameters <- function(parameters) {
  parameters <- slots2list(parameters)

  parameters[["nv_enr_rto"]] <- as.list(parameters[["nv_enr_rto"]])
  parameters[["nv_tfc_inl"]] <- as.list(parameters[["nv_tfc_inl"]])

  write_yaml(
    c(list(RPhosFate = as.character(packageVersion("RPhosFate"))), parameters),
    "parameters.yaml",
    precision = 15L,
    indent.mapping.sequence = TRUE
  )
}

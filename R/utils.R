adjustExtent <- function(rl, ex) {
  extend(crop(rl, ex), ex)
}

adjustMetric <- function(metric, metrics) {
  if (metric == "PBIAS") {
    abs(metrics[metric])
  } else if (metric == "RCV") {
    abs(metrics[metric] - 1)
  } else {
    metrics[metric]
  }
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

  adjustMetric(metric, metrics)
}

calibrate2 <- function(values, cmt, substance, col, metric) {
  cmt@parameters@ns_dep_ovl <- values[1L]
  cmt@parameters@ns_dep_cha <- values[2L]

  subsequentRun(cmt, substance)
  metrics <- calibrationQuality(cmt, substance, col)

  adjustMetric(metric, metrics)
}

determineMCfilename <- function(cmt, layer) {
  if (length(cmt@cv_dir) == 2L) {
    MCfilename <- file.path(
      cmt@cv_dir[2L],
      paste0(basename(layer), cmt@is_MCi, ".tif")
    )
    if (file.exists(MCfilename)) {
      return(MCfilename)
    }
  }

  MCfilename <- paste0(layer, cmt@is_MCi, ".tif")
  if (file.exists(MCfilename)) {
    MCfilename
  } else {
    paste0(layer, ".tif")
  }
}

findNearestNeighbour <- function(X, Y) {
  nn <- nearest(X, Y)
  Y <- Y[nn$to_id, ]

  cbind(
    values(X),
    values(nn),
    crds(Y, df = TRUE)
  )
}

populateLayerSlots <- function(
  cmt,
  object,
  slots,
  layers,
  areRequiredInputLayers = rep(FALSE, length(slots))
) {
  cs_dir_old <- setwd(cmt@cv_dir[1L])
  on.exit(setwd(cs_dir_old))

  for (i in seq_along(slots)) {
    slot(object, slots[i]) <- readLayer(
      cmt,
      layers[i],
      areRequiredInputLayers[i]
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

readLayer <- function(cmt, layer, isRequiredInputLayer = FALSE) {
  if (length(cmt@is_MCi) == 1L) {
    filename <- determineMCfilename(cmt, layer)
  } else {
    filename <- paste0(layer, ".tif")
  }

  if (isRequiredInputLayer) {
    rast(filename)
  } else if (file.exists(filename)) {
    rast(filename)
  } else {
    rast()
  }
}

readParameters <- function(arguments) {
  parameters <- read_yaml("parameters.yaml")

  parameters[["nv_tfc_inl"]] <- unlist(parameters[["nv_tfc_inl"]])
  parameters[["nv_enr_rto"]] <- unlist(parameters[["nv_enr_rto"]])
  parameters[["nm_olc"]] <- matrix(parameters[["nm_olc"]], ncol = 2L)
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
    function(name, parameters) slot(parameters, name),
    parameters = parameters
  ), parameterNames)
}

writeLayer <- function(cmt, layer, rl, datatype, substance = NULL) {
  if (length(cmt@is_MCi) == 0L || (length(cmt@is_MCi) == 1L &&
      layer %in% cmt@cv_MCl)) {
    if (!is.null(substance)) {
      layer <- sub("^xx", tolower(substance), layer)
    }
    name <- paste0(layer, cmt@is_MCi)
    set.names(rl, name)
    filename <- paste0(name, ".tif")

    writeRaster(
      rl,
      filename = filename,
      datatype = datatype,
      overwrite = TRUE
    )

    rast(filename)
  } else {
    rl
  }
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

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
      paste0(basename(layer), cmt@is_MCi, cmt@cs_fex)
    )
    if (file.exists(MCfilename)) {
      return(MCfilename)
    }
  }

  MCfilename <- paste0(layer, cmt@is_MCi, cmt@cs_fex)
  if (file.exists(MCfilename)) {
    MCfilename
  } else {
    paste0(layer, cmt@cs_fex)
  }
}

findNearestNeighbour <- function(X, Y, Extent) {
  win <- as.owin(Extent[1:4])
  pppX <- as.ppp(X, W = win)
  pppY <- as.ppp(Y, W = win)

  nn <- data.frame(
    X.x = X[, 1L],
    X.y = X[, 2L],
    X[, 3L, drop = FALSE],
    nncross(pppX, pppY)
  )
  dfY <- data.frame(
    Y.x = Y[, 1L],
    Y.y = Y[, 2L],
    Y[, 3L, drop = FALSE],
    index = seq_len(nrow(Y))
  )
  nn <- merge(nn, dfY, by.x = "which", by.y = "index", sort = FALSE)

  nn[2:8]
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
    filename <- paste0(layer, cmt@cs_fex)
  }

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
    filename <- paste0(layer, cmt@is_MCi, cmt@cs_fex)

    writeRaster(
      rl,
      filename = filename,
      datatype = datatype,
      options = if (cmt@cs_fex == ".img") "COMPRESSED=YES" else "COMPRESS=LZW",
      overwrite = TRUE
    )

    raster(filename)
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

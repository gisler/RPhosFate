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

#' Demonstration project
#'
#' Copies a demonstration project to an existing or a temporary directory.
#'
#' @param cs_dir An optional character string specifying an existing directory.
#'
#' @return A character string containing the demonstration project root
#'   directory.
#'
#' @export
demoProject <- function(cs_dir = tempdir(TRUE)) {
  assertDirectoryExists(cs_dir, access = "w")

  file.copy(
    system.file("demoData", "demoProject", package = "RPhosFate"),
    cs_dir,
    overwrite = FALSE,
    recursive = TRUE
  )

  file.path(cs_dir, "demoProject")
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
  parameterNames <- slotNames(parameters)
  argumentNames <- names(arguments)

  assertSubset(argumentNames, parameterNames)

  for (i in seq_along(arguments)) {
    if (argumentNames[i] %in% c("nv_enr_rto", "nv_tfc_inl")) {
      slot(parameters, argumentNames[i])[intersect(
        names(slot(parameters, argumentNames[i])),
        names(arguments[[i]])
      )] <- arguments[[i]][intersect(
        names(slot(parameters, argumentNames[i])),
        names(arguments[[i]])
      )]
      slot(parameters, argumentNames[i]) <- c(
        slot(parameters, argumentNames[i]),
        arguments[[i]][setdiff(
          names(arguments[[i]]),
          names(slot(parameters, argumentNames[i]))
        )]
      )
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

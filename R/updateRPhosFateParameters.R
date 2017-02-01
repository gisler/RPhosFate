readParameters <- function() {
  parameters <- readRDS("parameters.rds")
  lapply(slotNames(parameters), function(x) {slot(parameters, x)})
  parameters
}

slots2list <- function(parameters) {
  parameterNames <- slotNames(parameters)
  parameters <- lapply(parameterNames, function(name) {slot(parameters, name)})
  names(parameters) <- parameterNames
  return(parameters)
}

firstClassUpdate <- function() {
  setClass(
    Class = "RPhosFateParameters",
    slots = c(
      ns_slp_min = "numeric",
      ns_slp_max = "numeric",
      ns_rhy_a   = "numeric",
      ns_rhy_b   = "numeric",
      ns_cha_rto = "numeric",
      ns_man_rip = "numeric",
      ns_man_cha = "numeric",
      ns_dep_ovl = "numeric",
      ns_dep_cha = "numeric",
      ns_tfc_inl = "numeric",
      iv_fDo     = "integer",
      nm_olc     = "matrix",
      df_cdt     = "data.frame"
    ),
    where = as.environment("package:RPhosFate")
  )

  result <- tryCatch(
    readParameters(),
    error = function(e) {NULL}
  )

  return(result)
}

firstParametersUpdate <- function(parameters) {
  parameters$nv_enr_rto <- c(PP = 1)
  return(parameters)
}

updateRPhosFateParameters <- function(e) {
  unlockEnvironment(as.environment("package:RPhosFate"))

  if (!is.null(firstClassUpdate())) {
    parameters <- firstClassUpdate()
    parameters <- slots2list(parameters)
    parameters <- firstParametersUpdate(parameters)
  # nolint start
    # parameters <- secondParametersUpdate(parameters)
  } # else if (!is.null(secondClassUpdate())) {
    # parameters <- secondClassUpdate()
    # parameters <- slots2list(parameters)
    # parameters <- secondParametersUpdate(parameters)
  # }
  # nolint end

  setClassRPhosFateParameters(as.environment("package:RPhosFate"))

  lockEnvironment(as.environment("package:RPhosFate"))

  parameters <- new("RPhosFateParameters", parameters)

  message("New parameters available: Use \"saveState(cmt)\" to save the updated project.")

  return(parameters)
}

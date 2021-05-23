slots2list <- function(parameters) {
  parameterNames <- slotNames(parameters)

  setNames(lapply(
    parameterNames,
    function(name, parameters) {slot(parameters, name)},
    parameters = parameters
  ), parameterNames)
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

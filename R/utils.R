slots2list <- function(parameters) {
  parameterNames <- slotNames(parameters)

  setNames(lapply(
    parameterNames,
    function(name, parameters) {slot(parameters, name)},
    parameters = parameters
  ), parameterNames)
}

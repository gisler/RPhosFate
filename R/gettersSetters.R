#### getLayer ####
setGeneric(
  "getLayer",
  function(cmt, ...) standardGeneric("getLayer")
)
#' @export
setMethod(
  "getLayer",
  "RPhosFate",
  function(cmt, layer, substance = NULL) {
    qassert(layer, "S1")

    if (!is.null(substance)) {
      assertSubstance(cmt, substance)
      assertSubset(
        layer,
        sub("^rl_", "", slotNames(slot(cmt@substances, substance)))
      )

      return(slot(slot(cmt@substances, substance), sprintf("rl_%s", layer)))
    } else {
      for (object in cmt@helpers@cv_rlo) {
        if (layer %in% sub("^rl_", "", slotNames(slot(cmt, object)))) {
          return(slot(slot(cmt, object), sprintf("rl_%s", layer)))
        }
      }
    }

    stop(sprintf("Layer %s was not found.", deparse(layer)), call. = FALSE)
  }
)

#### getParameter ####
setGeneric(
  "getParameter",
  function(cmt, ...) standardGeneric("getParameter")
)
#' @export
setMethod(
  "getParameter",
  "RPhosFate",
  function(cmt, parameter) {
    qassert(parameter, "S1")
    assertSubset(parameter, slotNames(cmt@parameters))

    slot(cmt@parameters, parameter)
  }
)

#### setParameter ####
setGeneric(
  "setParameter",
  function(cmt, ...) standardGeneric("setParameter")
)
#' @export
setMethod(
  "setParameter",
  "RPhosFate",
  function(cmt, ...) {
    cmt@parameters <- populateParameterSlots(cmt@parameters, list(...))

    cmt
  }
)

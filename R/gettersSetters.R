#### getLayer ####
setGeneric(
  "getLayer",
  function(x, ...) standardGeneric("getLayer")
)
#' @export
setMethod(
  "getLayer",
  "RPhosFate",
  function(x, layer, substance = NULL) {
    qassert(layer, "S1")

    if (!is.null(substance)) {
      assertSubstance(x, substance)
      assertSubset(
        layer,
        sub("^rl_", "", slotNames(slot(x@substances, substance)))
      )

      return(slot(slot(x@substances, substance), sprintf("rl_%s", layer)))
    } else {
      for (object in x@helpers@cv_rlo) {
        if (layer %in% sub("^rl_", "", slotNames(slot(x, object)))) {
          return(slot(slot(x, object), sprintf("rl_%s", layer)))
        }
      }
    }

    stop(sprintf("Layer %s was not found.", deparse(layer)), call. = FALSE)
  }
)
#' @export
setMethod(
  "[",
  "RPhosFate",
  function(x, i, j) {
    if (missing(j)) {
      j <- NULL
    }

    getLayer(x, i, j)
  }
)

#### getParameter ####
setGeneric(
  "getParameter",
  function(x, ...) standardGeneric("getParameter")
)
#' @export
setMethod(
  "getParameter",
  "RPhosFate",
  function(x, parameter) {
    qassert(parameter, "S1")
    assertSubset(parameter, slotNames(x@parameters))

    slot(x@parameters, parameter)
  }
)

#### setParameter ####
setGeneric(
  "setParameter",
  function(x, ...) standardGeneric("setParameter")
)
#' @export
setMethod(
  "setParameter",
  "RPhosFate",
  function(x, ...) {
    x@parameters <- populateParameterSlots(x@parameters, list(...))

    x
  }
)

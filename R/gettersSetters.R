#### getLayer ####
setGeneric(
  "getLayer",
  function(x, ...) standardGeneric("getLayer")
)
#' Get Layer
#'
#' Obtains a project raster layer for further analysis.
#'
#' @inheritParams erosionPrerequisites,RPhosFate-method
#' @param i A character string specifying a layer name. Substance related layers
#'   whose names start with _xx_ are treated differently. They have to be
#'   queried by their name (not filename), for example, `"xxc"` together with
#'   the respective substance specified in argument `j`, for example, `"PP"`,
#'   which queries the particulate phosphorus concentration in top soil. See
#'   subdirectory sections for further information.
#' @param j A character string specifying a substance if applicable.
#'
#' @inheritSection catchment _Input_ subdirectory
#'
#' @inheritSection catchment _Intermediate_ subdirectory
#'
#' @inheritSection catchment _Result_ subdirectory
#'
#' @return A [`raster::RasterLayer-class`] object.
#'
#' @aliases getLayer
#'
#' @export
setMethod(
  "getLayer",
  "RPhosFate",
  function(x, i, j = NULL) {
    qassert(i, "S1")

    if (!is.null(j)) {
      assertSubstance(x, j)
      assertSubset(i, sub("^rl_", "", slotNames(slot(x@substances, j))))

      return(slot(slot(x@substances, j), sprintf("rl_%s", i)))
    } else {
      for (object in x@helpers@cv_rlo) {
        if (i %in% sub("^rl_", "", slotNames(slot(x, object)))) {
          return(slot(slot(x, object), sprintf("rl_%s", i)))
        }
      }
    }

    stop(sprintf("Layer %s was not found.", deparse(i)), call. = FALSE)
  }
)
#' @rdname getLayer-RPhosFate-method
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

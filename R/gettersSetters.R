#' @include aaa.R
NULL

#### getLayer ####
#' @export
setGeneric(
  "getLayer",
  function(x, ...) standardGeneric("getLayer")
)
#' Get layer
#'
#' Obtains a project raster layer for further analysis.
#'
#' @inheritParams erosionPrerequisites,RPhosFate-method
#' @param i A character string specifying a layer name. Substance related layers
#'   whose names start with _xx_ are treated differently. They have to be
#'   queried by their name (not filename), for example, `"xxc"` in combination
#'   with `"PP"` in argument `j` queries the particulate phosphorus
#'   concentrations in top soils. See subdirectory sections for further
#'   information.
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
#' @examples
#' \dontrun{
#' # create temporary demonstration project
#' cv_dir <- demoProject()
#' # load temporary demonstration project
#' x <- RPhosFate(
#'   cv_dir = cv_dir,
#'   ls_ini = TRUE
#' )
#' # presupposed function call
#' x <- firstRun(x, "SS")
#'
#' getLayer(x, "dir")
#' getLayer(x, "xxt", "SS")
#' getLayer(x, "xxe", "PP")
#' }
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
      for (object in c("topo", "erosion", "transport")) {
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
#' @export
setGeneric(
  "getParameter",
  function(x, ...) standardGeneric("getParameter")
)
#' Get parameter(s)
#'
#' Obtains a single model parameter or all model parameters at once.
#'
#' @inheritParams erosionPrerequisites,RPhosFate-method
#' @param parameter A character string specifying a parameter name or `NULL` for
#'   a [`list`] of all parameters. See model parameter arguments section for
#'   further information.
#'
#' @inheritSection catchment Model parameter arguments
#'
#' @return Depends on the queried parameter or a [`list`] in case of all
#'   parameters. See model parameter arguments section for further information.
#'
#' @seealso [`setParameter`]
#'
#' @examples
#' \dontrun{
#' # create temporary demonstration project
#' cv_dir <- demoProject()
#' # load temporary demonstration project
#' x <- RPhosFate(
#'   cv_dir = cv_dir,
#'   ls_ini = TRUE
#' )
#'
#' getParameter(x)
#' getParameter(x, "ns_dep_ovl")
#' }
#'
#' @aliases getParameter
#'
#' @export
setMethod(
  "getParameter",
  "RPhosFate",
  function(x, parameter = NULL) {
    if (!is.null(parameter)) {
      qassert(parameter, "S1")
      assertSubset(parameter, slotNames(x@parameters))

      slot(x@parameters, parameter)
    } else {
      slots2list(x@parameters)
    }
  }
)

#### setParameter ####
#' @export
setGeneric(
  "setParameter",
  function(x, ...) standardGeneric("setParameter")
)
#' Set parameter(s)
#'
#' Sets one or more model parameters or substance parameter values.
#'
#' @inheritParams erosionPrerequisites,RPhosFate-method
#' @param \dots Names and values of the parameters to set. See model parameter
#'   arguments section for further information.
#'
#' @inheritSection catchment Model parameter arguments
#'
#' @inherit catchment return
#'
#' @seealso [`getParameter`]
#'
#' @examples
#' \dontrun{
#' # create temporary demonstration project
#' cv_dir <- demoProject()
#' # load temporary demonstration project
#' x <- RPhosFate(
#'   cv_dir = cv_dir,
#'   ls_ini = TRUE
#' )
#'
#' x <- setParameter(x, ns_dep_ovl = 15e-4)
#' x <- setParameter(
#'   x,
#'   nv_tfc_inl = c(SS = 0.6, PP = 0.6),
#'   nv_enr_rto = c(PP = 1.4)
#' )
#' }
#'
#' @aliases setParameter
#'
#' @export
setMethod(
  "setParameter",
  "RPhosFate",
  function(x, ...) {
    x@parameters <- populateParameterSlots(x@parameters, list(...))

    x
  }
)

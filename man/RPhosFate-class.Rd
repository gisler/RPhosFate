% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aaa.R
\docType{class}
\name{RPhosFate-class}
\alias{RPhosFate-class}
\title{RPhosFate class}
\description{
An S4 object representing a river catchment.
}
\section{Slots}{

\describe{
\item{\code{cv_dir}}{A character vector holding the project root (first position) and
optionally the Monte Carlo input data directory (second position).}

\item{\code{ls_ini}}{A logical scalar specifying if the state of an existing project
was loaded from disk.}

\item{\code{is_ths}}{An integer scalar holding the number of threads to use for
processing, where applicable.}

\item{\code{is_MCi}}{An integer scalar holding the current Monte Carlo iteration if
applicable.}

\item{\code{cv_MCl}}{A character vector holding the names of the layers, which shall
be written to disk with the associated Monte Carlo iteration in their
filenames upon calling the appropriate methods.}

\item{\code{parameters}}{An S4 object holding the model parameters.}

\item{\code{topo}}{An S4 object holding the raster layers related to topography in
the broader sense.}

\item{\code{erosion}}{An S4 object holding the raster layers related to erosion.}

\item{\code{transport}}{An S4 object holding raster layers required for modelling
transport.}

\item{\code{substances}}{An S4 object holding the substance raster layer containers.}

\item{\code{helpers}}{An S4 object holding helper data.}
}}

\seealso{
\code{\link{RPhosFate}}, \code{\link{catchment}}
}

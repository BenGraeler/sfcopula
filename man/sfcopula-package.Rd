\name{sfcopula-package}
\alias{sfcopula-package}
\alias{sfcopula}
\docType{package}
\title{Copula Driven Analysis of Spatial and Spatio-Temporal Data Based on 'sf' and 'stars'}
\description{
A framework to analyse multivariate, spatial and spatio-temporal data via copulas, vine copulas and hierarchical Kendall copulas. Spatial data is provided as simple features (\code{\link[sf]{sf}} objects with POINT geometries) and spatio-temporal data as \code{stars} vector data cubes with one geometry and one time dimension (see \code{\link[stars]{st_as_stars}} and \code{\link{as_spacetime_cube}}). Additionally, support for calculating different multivariate return periods based on copulas and vine copulas is implemented.
}
\details{
\pkg{sfcopula} is the successor of the package \pkg{spcopula}, which was based on the classes of the \pkg{sp} package and its spatio-temporal extension. The copula functionality is unchanged; the spatial and spatio-temporal data handling has been ported to \pkg{sf} and \pkg{stars}:
\itemize{
  \item spatial data (formerly \code{Spatial*DataFrame}) are \code{sf} objects, or bare \code{sfc} geometries where only locations are needed,
  \item spatio-temporal data (formerly full space-time grid data frames) are \code{stars} vector data cubes; spatio-temporal prediction targets (formerly full space-time grids) are either such cubes or \code{sf} objects with a time column,
  \item distances are Euclidean in the units of the coordinate reference system for projected (or unknown) CRS and great circle distances in metres (\pkg{spcopula}: kilometres) for geographic coordinates.
}

The package imports \pkg{sf}, \pkg{stars}, \pkg{VineCopula} and \pkg{methods} and depends on \pkg{copula} and \pkg{VC2copula}.

Central functions are \code{\link{calc_bins}}, \code{\link{fit_cor_fun}}, \code{\link{loglik_by_lags}}, \code{\link{spatial_copula}}, \code{\link{spacetime_copula}}, \code{\link{distance_vine_copula}}, \code{\link{neighbours}}, \code{\link{cond_vine}} and \code{\link[=predict,distance_vine_copula-method]{predict}}. See \code{demo(package = "sfcopula")} for worked examples.
}
\author{
Benedikt Graeler, with contributions by Marius Appel

Maintainer: Benedikt Graeler <b.graeler@52north.org>
}
\references{
Graeler, B. and E. Pebesma (2011): The pair-copula construction for spatial data: a new approach to model spatial dependency. Procedia Environmental Sciences (Vol. 7, pp. 206 - 211), Elsevier.
}
\keyword{ package }

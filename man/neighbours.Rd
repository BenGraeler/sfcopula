\name{neighbours}
\alias{neighbours}

\title{
Creating Local (Spatio-Temporal) Neighbourhoods
}
\description{
This function calculates local spatial or spatio-temporal neighbourhoods to be used for fitting of spatial/spatio-temporal vine copulas and for prediction using spatial/spatio-temporal vine copulas. Spatial neighbourhoods are derived from \code{\link[sf]{sf}} objects, spatio-temporal ones from \code{stars} vector data cubes.
}
\usage{
neighbours(data, target = NULL, size = NULL, var = NULL, covar = character(),
           prediction = !is.null(target), min_dist = 0.01,
           tlags = -(0:2), time_steps = NA, time_col = "time", ...)
}
\arguments{
  \item{data}{
Spatial case: an \code{\link[sf]{sf}} object with POINT geometries holding the data used for estimation/prediction; a bare \code{sfc} geometry is accepted as well, in which case the neighbourhood carries no data.
Spatio-temporal case: a \code{\link[stars]{st_as_stars}} vector data cube with exactly two dimensions, one geometry dimension (\code{sfc} with POINT geometries) and one time dimension (see also \code{\link{as_spacetime_cube}}).
}
  \item{target}{
the prediction targets, \code{NULL} (default) if the neighbourhood is used for fitting.
Spatial case: an \code{\link[sf]{sf}} or \code{sfc} object with POINT geometries defining the prediction locations.
Spatio-temporal case: either a \code{stars} vector data cube as \code{data} (all combinations of its locations and time instances are used; the neighbourhoods are ordered by location first and by time within each location) or an \code{\link[sf]{sf}} object with POINT geometries holding one time stamp per feature in the column \code{time_col}. All target time instances need to be part of the time dimension of \code{data}.
}
  \item{size}{The (spatial) size of the neighbourhood including the location of interest (for fitting as well for prediction). Defaults to 5 in the spatial and to 4 in the spatio-temporal case. In the spatio-temporal case, each of the \code{size - 1} spatial neighbours is used at every temporal lag in \code{tlags}.}
  \item{var}{the variable name of interest, by default (\code{NULL}) the first attribute (non-geometry) column of \code{data} is used}
  \item{covar}{the name of the covariate, if any. Its values at the central location are added as last data column(s): taken from \code{data} for fitting and from \code{target} for prediction (the target then needs to be an \code{sf} object or \code{stars} cube holding the covariate).}
  \item{prediction}{whether the neighbourhood should be used for prediction (TRUE) or spatial/spatio-temporal vine copula fitting; defaults to \code{TRUE} if \code{target} is provided.}
  \item{min_dist}{the minimal distance for a location to be included. Must be larger than 0 for fitting purposes and might be 0 for prediction}
  \item{tlags}{spatio-temporal case only: the temporal lags to be used in the spatio-temporal neighbourhood.}
  \item{time_steps}{spatio-temporal case only: the number of time instances that should randomly be selected from \code{data}. The default, \code{NA}, selects all time instances.}
  \item{time_col}{spatio-temporal case only: the name of the time column if \code{target} is an \code{sf} object (ignored for \code{sftime} objects, whose time column is used).}
  \item{\dots}{deprecated argument names of \pkg{spcopula} (\code{coVar}, \code{min.dist}, \code{timeSteps}, \code{timeCol}, \code{spSize}) are accepted with a warning.}
}
\details{
Distances are Euclidean distances in the units of the coordinate reference system for projected (or unknown) CRS. For geographic (longitude/latitude) coordinates, great circle distances in metres are calculated via \code{\link[sf]{st_distance}}. Temporal lags are counted in steps of the time dimension of \code{data}.
}
\value{An object of \code{\linkS4class{neighbourhood}}: spatial neighbourhoods hold distance and index matrices, spatio-temporal ones arrays with the spatial and temporal component in the third dimension.}
\author{Benedikt Graeler}

\seealso{
See \code{\link{neighbourhood}} for the native constructor of a \code{\linkS4class{neighbourhood}} class; \code{\link{as_spacetime_cube}} builds a vector data cube from a long \code{sf} table and \code{\link{reduce_neighbours}} selects the strongest correlated spatio-temporal neighbours.
}
\examples{
library("sf")
spdf <- st_as_sf(data.frame(x=c(112,154,212,289,345), y=c(124,198,85,168,346),
                            measure=rlnorm(5)),
                 coords=c("x","y"))

# spatial neighbourhoods for fitting and prediction
neighbours(spdf, size=4)
neighbours(spdf[1:4,], spdf[5,], size=3)

library("stars")

pts <- st_sfc(st_point(c(181000,333500)), st_point(c(181100,333600)),
              crs = 28992)
time <- as.Date("2014-03-18") + 0:2

stData <- st_as_stars(list(measure = matrix(runif(6), 2, 3)),
                      dimensions = st_dimensions(geometry = pts, time = time))

# prediction targets as vector data cube (all location/time combinations)
stQuerry <- st_as_stars(list(dummy = matrix(NA_real_, 2, 2)),
                        dimensions = st_dimensions(
                          geometry = st_sfc(st_point(c(181000,333600)),
                                            st_point(c(181200,333600)),
                                            crs = 28992),
                          time = time[2:3]))

neighbours(stData, stQuerry, var="measure", size=3, tlags=-(0:1))

# the same targets as sf object with a time column
stTargets <- st_sf(time = rep(time[2:3], 2),
                   geometry = st_sfc(st_point(c(181000,333600)),
                                     st_point(c(181000,333600)),
                                     st_point(c(181200,333600)),
                                     st_point(c(181200,333600)),
                                     crs = 28992))

neighbours(stData, stTargets, var="measure", size=3, tlags=-(0:1))

# spatio-temporal neighbourhoods for fitting
neighbours(stData, var="measure", size=2, tlags=-(0:1))
}
\keyword{ spatial }
\keyword{ spatio-temporal }

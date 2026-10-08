\name{as_spacetime_cube}
\alias{as_spacetime_cube}

\title{
Building a Spatio-Temporal Vector Data Cube from a Long Table
}
\description{
Converts a long \code{\link[sf]{sf}} table holding one row per location and time instance into a \code{stars} vector data cube with a geometry and a time dimension, the spatio-temporal data structure used throughout \pkg{sfcopula}.
}
\usage{
as_spacetime_cube(x, time = "time", vars = NULL)
}
\arguments{
  \item{x}{an \code{\link[sf]{sf}} object with POINT geometries and one row per combination of location and time instance.}
  \item{time}{the name of the column holding the time stamps (e.g. \code{Date} or \code{POSIXct}). For \code{sftime} objects, the time column of the object is used.}
  \item{vars}{the names of the columns to be stored as attributes of the cube; by default all attribute columns besides \code{time}.}
}
\details{
Locations are identified by their coordinates and time instances by their unique (sorted) time stamps. Each combination of location and time instance may occur at most once; combinations missing in \code{x} are filled with \code{NA}.
}
\value{
A \code{stars} object (see \code{\link[stars]{st_as_stars}}) with the dimensions \code{geometry} (an \code{sfc} of the unique locations in order of their first appearance) and \code{time} (the sorted unique time stamps), and one attribute per entry in \code{vars}.
}
\author{
Benedikt Graeler
}

\seealso{
\code{\link{neighbours}}, \code{\link{calc_bins}}
}
\examples{
library("sf")

long <- data.frame(x = rep(c(181000, 181100), each = 3),
                   y = rep(c(333500, 333600), each = 3),
                   time = rep(as.Date("2014-03-18") + 0:2, 2),
                   measure = runif(6))
long <- st_as_sf(long[-4, ], coords = c("x", "y"), crs = 28992)

stData <- as_spacetime_cube(long)
stData
stData[["measure"]]
}
\keyword{ spatio-temporal }

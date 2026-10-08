\name{neighbourhood}
\alias{neighbourhood}
\title{
Constructor of the \code{\linkS4class{neighbourhood}} class
}
\description{An object that stores the nearest k (spatial or spatio-temporal) neighbours for a set of locations.}

\usage{neighbourhood(data, distances, index, var, coVar = character(), prediction = FALSE)}

\arguments{
  \item{data}{A data.frame (or matrix) that contains the variables of interest for each local neighbourhood.}
  \item{distances}{Spatial neighbourhoods: a matrix \code{[n, k]} with all the separating distances per local neighbourhood. Spatio-temporal neighbourhoods: a 3-dimensional array \code{[n, k, 2]} with all the separating spatial (\code{[,,1]}) and temporal (\code{[,,2]}) distances per local neighbourhood.}
  \item{index}{Spatial neighbourhoods: a matrix \code{[n, k+1]} pointing from each neighbourhood to the corresponding row indices of the underlying \code{\link[sf]{sf}} object. Spatio-temporal neighbourhoods: a 3-dimensional array \code{[n, k+1, 2]} pointing from each neighbourhood to the corresponding indices of the geometry (\code{[,,1]}) and time (\code{[,,2]}) dimension of the underlying vector data cube.}
  \item{var}{The variable's name.}
  \item{coVar}{The covariate's name, if any.}
  \item{prediction}{whether this neighbourhood is intended to be used for prediction or fitting purpose.}
}

\details{
The first column of \code{index} refers to the location of interest (the prediction location in case of prediction). In the spatio-temporal prediction case, the data of the central location are not part of \code{data}. Spatio-temporal data columns are renamed as \code{N<i>.<var>}.
}

\value{A \code{\linkS4class{neighbourhood}} object.}

\author{Benedikt Graeler}

\seealso{
\code{\linkS4class{neighbourhood}}, \code{\link{neighbours}}
}
\examples{
library("sf")
spdf <- st_as_sf(data.frame(x=c(112,154,212,289), y=c(124,198,85,168),
                            measure=rlnorm(4)),
                 coords=c("x","y"))

neigh <- neighbours(spdf, size=4)
neigh

# rebuilding neigh
neighbourhood(neigh@data, neigh@distances, neigh@index, neigh@var,
              neigh@coVar, neigh@prediction)

# a spatio-temporal neighbourhood
library("stars")

pts <- st_sfc(st_point(c(181000,333500)), st_point(c(181100,333600)),
              crs = 28992)
time <- as.Date("2014-03-18") + 0:2

stData <- st_as_stars(list(var1 = matrix(runif(6), 2, 3)),
                      dimensions = st_dimensions(geometry = pts, time = time))

stNeigh <- neighbours(stData, size=2, tlags=-(0:1))
stNeigh

# rebuilding stNeigh
neighbourhood(stNeigh@data, stNeigh@distances, stNeigh@index, stNeigh@var,
              stNeigh@coVar, stNeigh@prediction)
}

\keyword{spatial}
\keyword{spatio-temporal}

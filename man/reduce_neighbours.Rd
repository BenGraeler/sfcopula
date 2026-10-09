\name{reduce_neighbours}
\alias{reduce_neighbours}

\title{
Selecting the strongest correlated neighbours
}
\description{
A function selecting the strongest correlated neighbours from a larger set of neighbours
}
\usage{
reduce_neighbours(neigh, dep_fun, n,
                 prediction=neigh@prediction, drop_empty=!prediction, ...)
}

\arguments{
  \item{neigh}{
the proxy spatio-temporal \code{\linkS4class{neighbourhood}} to be investigated
}
  \item{dep_fun}{
a spatio-temporal dependence function that return correlation estimates based on a spatial and temporal distance
}
  \item{n}{the number of neighbours to be selected.}
  \item{prediction}{whether the neighbourhood is used for prediction (the data slot does not provide the central location's data)}
  \item{drop_empty}{whether empty neighbourhoods (i.e. neighbourhoods with \code{NA}'s at the central location) shall be dropped.}
  \item{\dots}{deprecated argument names of \pkg{spcopula} (\code{stNeigh}, \code{stDepFun}, \code{dropEmpty}) are accepted with a warning.}
}
\value{A spatio-temporal neighbourhood \code{\linkS4class{neighbourhood}} with fewer neighbours.}

\author{
Benedikt Graeler
}


\seealso{
\code{\link{neighbours}}
}

\examples{
library("sf")
library("stars")

pts <- st_sfc(st_point(c(181000,333500)), st_point(c(181100,333600)),
              crs = 28992)
time <- as.Date("2014-03-18") + 0:4

stData <- st_as_stars(list(var1 = matrix(runif(10), 2, 5)),
                      dimensions = st_dimensions(geometry = pts, time = time))

neigh <- neighbours(stData, size=2, tlags=-(0:2))

reduce_neighbours(neigh, function(h,delta) return(1/h/delta), 2)
}
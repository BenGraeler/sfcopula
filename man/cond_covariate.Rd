\name{cond_covariate}
\alias{cond_covariate}

\title{
Conditioning of a Covariate
}
\description{
Conditions the covariate (last data column) on the central location (first data column) of a spatial or spatio-temporal neighbourhood, e.g. for the second tree of a \code{\linkS4class{covariate_vine_copula}}.
}
\usage{
cond_covariate(neigh, covar_cop, ...)
}

\arguments{
  \item{neigh}{a spatial or spatio-temporal \code{\linkS4class{neighbourhood}} with one covariate (see \code{\link{neighbours}}).}
  \item{covar_cop}{a function returning a \code{\linkS4class{copula}} object describing the dependence between the variable of interest and the covariate. It receives the index of the central location: the pair (location, time) for spatio-temporal and the location index for spatial neighbourhoods.}
  \item{\dots}{deprecated argument names of \pkg{spcopula} (\code{stNeigh}, \code{coVarCop}) are accepted with a warning.}
}

\value{A vector of conditioned data, i.e. covariate|variable of interest}

\author{
Benedikt Graeler
}

\seealso{\code{\linkS4class{neighbourhood}}, \code{\link{cond_vine}}}

\examples{
library("sf")
library("stars")

pts <- st_sfc(st_point(c(181000,333500)), st_point(c(181100,333600)),
              crs = 28992)
time <- as.Date("2014-03-18") + 0:2
var <- runif(6)
covar <- invdduCopula(var, gumbelCopula(7), runif(6))

stData <- st_as_stars(list(var = matrix(var, 2, 3),
                           covar = matrix(covar, 2, 3)),
                      dimensions = st_dimensions(geometry = pts, time = time))

# for prediction, the covariate is taken from the targets
stQuerry <- st_as_stars(list(covar = matrix(runif(4), 2, 2)),
                        dimensions = st_dimensions(
                          geometry = st_sfc(st_point(c(181000,333600)),
                                            st_point(c(181200,333600)),
                                            crs = 28992),
                          time = time[2:3]))

neigh <- neighbours(stData, stQuerry, size=3, tlags=-(0:1),
                      var="var", covar="covar")

cond_covariate(neigh, function(x) gumbelCopula(7))
}

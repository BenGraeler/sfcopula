\name{neighbourhood-class}
\Rdversion{1.1}
\docType{class}
\alias{neighbourhood-class}
\alias{names,neighbourhood-method}
\alias{show,neighbourhood-method}
\alias{[,neighbourhood,numeric,ANY,ANY-method}
\alias{calc_bins,neighbourhood-method}

\title{Class \code{"neighbourhood"}}
\description{A class representing a set of local spatial or spatio-temporal neighbourhoods including data, indices and distances.}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("neighbourhood", ...)}, from the constructor \code{\link{neighbourhood}} or through the easier function \code{\link{neighbours}}.
}
\section{Slots}{
  \describe{
    \item{\code{data}:}{Object of class \code{"data.frame"} holding the data for each local neighbourhood.}
    \item{\code{distances}:}{Object of class \code{"array"} providing the separating distances per neighbourhood: a matrix \code{[n, k]} for spatial neighbourhoods, or an array \code{[n, k, 2]} holding the spatial and temporal distances for spatio-temporal neighbourhoods.}
    \item{\code{index}:}{Object of class \code{"array"} pointing to the row indices of the underlying \code{\link[sf]{sf}} object (a matrix \code{[n, k+1]}), or to the geometry and time indices of the underlying vector data cube (an array \code{[n, k+1, 2]}). The first column refers to the location of interest, i.e. the prediction location in case of prediction.}
    \item{\code{var}:}{Object of class \code{"character"} giving the variable's name.}
    \item{\code{coVar}:}{Object of class \code{"character"} giving the covariate's name.}
    \item{\code{prediction}:}{Object of class \code{"logical"}, whether this neighbourhood is suited for fitting or prediction purposes.}
  }
}
\section{Methods}{
  \describe{
    \item{names}{\code{signature(x = "neighbourhood")}: provides the variable names of the neighbourhood. }
    \item{show}{\code{signature(object = "neighbourhood")}: a brief description of the characteristics of the neighbourhood.}
    \item{[}{\code{signature(x = "neighbourhood", i = "numeric")}: subsets the selection of neighbourhoods returning a subset of these (row wise).}
    \item{calc_bins}{\code{signature(data = "neighbourhood")}: calculates bins from an existing (conditional) spatial neighbourhood for repeated application of spatial trees in a spatial vine copula; with the arguments \code{var = data@var}, \code{nbins = 9}, \code{boundaries}, \code{cutoff}, \code{cor.method = "kendall"} and \code{plot} (see \code{\link{calc_bins}}). Only available for spatial neighbourhoods.}
	 }
}
\author{
Benedikt Graeler
}

\seealso{
\code{\link{neighbours}}, \code{\link{drop_tree}}
}
\examples{
library("sf")
spdf <- st_as_sf(data.frame(x=c(112,154,212,289), y=c(124,198,85,168),
                            measure=rlnorm(4)),
                 coords=c("x","y"))

neigh <- neighbours(spdf, size=3)
neigh
names(neigh)
neigh[1:2]@distances
}
\keyword{classes}

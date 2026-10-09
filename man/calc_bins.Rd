\name{calc_bins}
\alias{calc_bins}
\alias{calc_bins-methods}
\alias{calc_bins,sf-method}
\alias{calc_bins,stars-method}

\title{
A function calculating the spatial/spatio-temporal bins
}
\description{
The (spatio-temporal) space is subdivided into pairs of observations that belong to certain spatial/spatio-temporal distance classes. For each distance class, the mean separating distance of all pairs involved is calculated alongside a correlation measure. The spatial/spatio-temporal correlogram is plotted by default.
}
\usage{
calc_bins(data, var, nbins = 15, boundaries = NA, cutoff = NA,
         ..., cor_method="fasttau", plot=TRUE)
}

\arguments{
  \item{data}{
Either an \code{\link[sf]{sf}} object with POINT geometries (spatial case) or a \code{stars} vector data cube with one geometry and one time dimension (spatio-temporal case, see \code{\link{as_spacetime_cube}} and \code{\link[stars]{st_as_stars}}). A method for spatial \code{\linkS4class{neighbourhood}} objects exists as well.
}
  \item{var}{
The name of the variable (attribute column) of interest.
}
  \item{nbins}{
the number of spatial bins.
}
  \item{boundaries}{
The right hand-side limits of the lag classes to be investigated, if set to \code{NA}, the \code{cutoff} will be split in equal parts according to \code{nbins}.
}
  \item{cutoff}{
The largest spatial distance to be investigated, if set to \code{NA}, one third of the bounding box's diagonal will be used.
}
\item{\dots}{Additional arguments for the spatio-temporal case:
\describe{
    \item{instances}{To reduce the data size or circumvent unwanted autocorrelation effects, one might provide a single number of randomly selected time instances from the vector data cube. If this parameter is set to \code{NA} (default), the complete time series will be used. Otherwise, \code{instances} is interpreted as a vector of time indices or time stamps (\code{Date}/\code{POSIXct} values that are part of the time dimension) to be used; \code{tlags} is set to 0 in this case.}
    \item{tlags}{a vector indicating the time lags to be investigated}
}
}
  \item{cor_method}{
  defining the correlation method that should be used. Possible choices are \code{kendall}, \code{pearson}, \code{spearman} and \code{fasttau} (default) that re-uses a very fast implementation of Kendall's tau from \code{\link{VineCopula-package}}.
  }
  \item{plot}{Whether the correlogram should be plotted (default \code{TRUE}; the \code{stars} method defaults to \code{FALSE}).}
}
\details{
Distances are Euclidean distances in the units of the coordinate reference system for projected (or unknown) CRS, and great circle distances in metres for geographic (longitude/latitude) coordinates.
}
\value{A list containing (a subset of) the following entries
\item{np}{The number of pairs per lag (spatial and \code{neighbourhood} case).}
\item{meanDists}{A vector holding the mean separating distance for each lag.}
\item{lagCor}{The calculated correlation measure for each lag; a matrix with one row per temporal lag in the spatio-temporal case.}
\item{lagNp}{A matrix with the number of complete pairs per temporal (rows) and spatial lag (columns), only in the spatio-temporal case.}
\item{lagData}{The data pairs for each lag (\code{neighbourhood} case).}
\item{lags}{In the spatial case, a list holding for each lag a matrix of the indices of the locations of each pair and their distance. In the spatio-temporal case, a list with the entries \code{sp} (the spatial lags as in the spatial case) and \code{time} (a matrix of temporal indices).}
The attributes \code{cor_method} and \code{variable} store the correlation method and the variable's name.
}
\author{
Benedikt Graeler
}
\examples{
library("sf")
data("meuse", package = "sfcopula")
meuse$rtZinc <- rank(meuse$zinc)/(nrow(meuse)+1)

## lag classes ##
bins <- calc_bins(meuse, var="rtZinc", nbins=10, cutoff=800)

## spatio-temporal lag classes ##
library("stars")
data("EU_RB", package = "sfcopula")
stBins <- calc_bins(EU_RB, var="PM10", nbins=10, cutoff=5e5,
                   instances=10, tlags=-(0:2))
stBins$lagCor
}

\keyword{spatial preparation}
\keyword{spatio-temporal preparation}

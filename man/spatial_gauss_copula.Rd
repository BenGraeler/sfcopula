\name{spatial_gauss_copula}
\alias{spatial_gauss_copula}
\title{
Spatial Gaussian Copula
}
\description{
Constructor of a spatial Gaussian copula \code{\linkS4class{spatial_gauss_copula}}, i.e. a Gaussian copula whose correlation matrix is derived from the distances between the locations of a neighbourhood through a correlation function.
}
\usage{
spatial_gauss_copula(cor_fun, ...)
}
\arguments{
  \item{cor_fun}{
A valid correlogram (i.e. producing a valid correlation matrix when applied to a distance matrix; e.g. based on a variogram).
}
  \item{\dots}{deprecated argument names of \pkg{spcopula} (\code{corFun}) are accepted with a warning.}
}
\value{
An instance of the \code{\linkS4class{spatial_gauss_copula}} class.
}
\author{
Benedikt Graeler
}
\seealso{
\code{\link[=predict,spatial_gauss_copula-method]{predict}} and \code{\link{spatial_gauss_loglik}}
}
\examples{
spGaussCop <- spatial_gauss_copula(function(h) exp(-h/400))
spGaussCop
}
\keyword{ distribution }

\name{compose_spatial_copula}
\alias{compose_spatial_copula}
\title{
Composing a bivariate Spatial Copula
}
\description{
A bivariate spatial copula is composed out of a set of bivariate copulas. These are combined using a convex linear combination with weights based on distances where for copulas with a 1-1 correspondence of Kendall's tau or Spearman's rho a dependence function providing measures of association based on distances might be used.
}
\usage{
compose_spatial_copula(best_fit, families, bins, calc_cor, range=max(bins$meanDists), ...)
}
\arguments{
  \item{best_fit}{
A vector indicating the best fitting copula family index (as in \code{families}) per lag-class. If this vector is shorter than the actual number of bins, only these first bins will be used. It is an alternative way of setting a range.
}
  \item{families}{
The list of investigated families.
}
  \item{bins}{
The spatial bins used.
}
  \item{calc_cor}{
Optional, some spatially dependent function providing values of Kendall's tau for distances.
}
  \item{range}{
  The maximal distance to be used for the spatial copula. By default, the maximum distance provided by bins will be used. 
  }
  \item{\dots}{deprecated argument names of \pkg{spcopula} (\code{bestFit}, \code{calcCor}) are accepted with a warning.}
}
\value{
A (\code{\linkS4class{spatial_copula}}) object.
}
\author{
Benedikt Graeler
}

\examples{
compose_spatial_copula(c(1,1,2,3),families=list(frankCopula(.4), gumbelCopula(1.6),gumbelCopula(1.4)),
                bins=data.frame(meanDists=c(0,500,1000,1500,2000)),range=1750)
}

\keyword{ spatial }
\keyword{ multivariate }
\keyword{ distribution }

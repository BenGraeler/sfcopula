\name{spatial_gauss_loglik}
\alias{spatial_gauss_loglik}

\title{
Density evaluation for a spatial Gaussian Copula
}
\description{
Evaluates the density for a spatial Gaussian Copula.
}
\usage{
spatial_gauss_loglik(corFun, neigh, dataLocs, log = TRUE)
}

\arguments{
  \item{corFun}{A valid correlogram (i.e. producing a valid correlation matrix; e.g. based on a variogram) or a \code{\linkS4class{spatial_gauss_copula}}.}
  \item{neigh}{A \code{\linkS4class{neighbourhood}} object to be evaluated.}
  \item{dataLocs}{The same \code{\link[sf]{sf}} (or \code{sfc}) object used to generate \code{neigh}.}
  \item{log}{Should the log-likelihood be returned?}
}
\details{
Evaluates the density for all neighbourhoods in \code{neigh} and returns the (log)-likelihood.
}
\value{
The (log)-likelihood value.
}
\author{
Benedikt Graeler
}

\seealso{
\code{\link{spatial_gauss_copula}}, \code{\link[=predict,spatial_gauss_copula-method]{predict}}
}

\examples{
library("sfcopula")

# load data from the Meuse demo
data("spCopDemo", package = "sfcopula")

# calculate the correlation function based on Kendall's tau
calcKTauPol <- fit_cor_fun(bins, degree=1)

# translate Kendall's tau correlation function into Gaussian Copula parameters 
# using a linear variogram
meuseGaussCorFun <- function(h) {
  res <- pmax(iTau(normalCopula(0), calcKTauPol(0))/658*(658-h),0)
  res[h ==0] <- 1
  return(res)
}

# get the neighbours
library("sf")
data("meuse", package = "sfcopula")

meuse$rtZinc <- rank(meuse$zinc)/(nrow(meuse)+1)

meuseNeigh <- neighbours(meuse, size=5L, var="rtZinc")

# calculate the log-likelihood
spatial_gauss_loglik(meuseGaussCorFun, meuseNeigh, meuse)

# the same based on a spatial Gaussian copula
spatial_gauss_loglik(spatial_gauss_copula(meuseGaussCorFun), meuseNeigh, meuse)
}

\keyword{ distribution }

\name{spacetime_copula}
\alias{spacetime_copula}
\title{
Spatio-Temporal Copula
}
\description{
Constructor of a bivariate spatio-temporal copula \code{\linkS4class{spacetime_copula}}.
}
\usage{
spacetime_copula(components, tlags, distances=NA, stDepFun, unit="m", tres="day")
}

\arguments{
  \item{components}{A list of bivariate spatial copulas (\code{\linkS4class{spatial_copula}}) to be used at each temporal lag. Or a list of with lists of the spatial components per temporal lag together with the argument \code{distances}.}
  \item{tlags}{The temporal lags used in the spatio-temporal copula.}
  \item{distances}{This and the following 2 arguments are only necessary when the provided \code{components} argument is not yet a list of \code{\linkS4class{spatial_copula}}s: A vector of the mean distances of the spatial lag classes.}
  \item{stDepFun}{A list of spatial dependence functions; one per temporal lag. This argument is only needed when components is not yet a list of \code{\linkS4class{spatial_copula}}s.}
  \item{unit}{The spatial unit, default: m (meters). This argument is only needed when components is not yet a list of \code{\linkS4class{spatial_copula}}s.}
  \item{tres}{The temporal resolution, default: day}
}
\value{
An instance of the spatio-temporal Copula class \code{\linkS4class{spacetime_copula}}.
}
\author{
Benedikt Graeler
}

\seealso{
\code{\link{spatial_copula}}
}
\examples{
spCopT0 <- spatial_copula(components=list(claytonCopula(8), claytonCopula(4), 
                                    claytonCopula(2), claytonCopula(1),
                                    claytonCopula(0.5), indepCopula()),
                    distances=c(0,100,200,300,400,500),
                    unit="km")
spCopT1 <- spatial_copula(components=list(claytonCopula(4), claytonCopula(2), 
                                    claytonCopula(1), claytonCopula(0.5),
                                    indepCopula()),
                    distances=c(0,100,200,300,400),
                    unit="km")
spCopT2 <- spatial_copula(components=list(claytonCopula(2), claytonCopula(1), 
                                    claytonCopula(0.5), indepCopula()),
                    distances=c(0,100,200,300),
                    unit="km")

stCop <- spacetime_copula(components=list(spCopT0, spCopT1, spCopT2),
                  tlags=-(0:2))
}
\keyword{sfcopula}
\keyword{copula}

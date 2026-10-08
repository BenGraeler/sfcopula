\name{spatial_copula}
\alias{spatial_copula}
\title{
Spatial Copula
}
\description{
Constructor of a bivariate spatial copula \code{\linkS4class{spatial_copula}}. The bivariate copulas provided for a set of distances are combined into a distance dependent copula, either through a convex combination or through geometric means.
}
\usage{
spatial_copula(components, distances, spDepFun, unit = "m",
               combination = c("convex", "geometric"))
}

\arguments{
  \item{components}{
A list of bivariate copulas to be used at each lag. Typically, the last one is the independence copula used beyond the range.
}
  \item{distances}{
A strictly increasing vector of the mean distances of the lag classes PLUS the range which will define the use of the last copula in \code{components} (typically the independence copula). For the geometric combination, \code{distances} must contain 0.
}
  \item{spDepFun}{
The spatial dependence function. Typically, the result of \code{\link{fit_cor_fun}}. If no spatial dependence function is provided, the spatial copula will be a pure linear convex combination of the provided components (i.e. no further fit of their parameters will take place). Only supported for the convex combination.
}
  \item{unit}{
The spatial unit, default: \code{"m"}
}
  \item{combination}{
How neighbouring components are combined between two distances: \code{"convex"} (default) for a distance dependent convex combination, or \code{"geometric"} for a distance dependent geometric mean of their distribution functions.
}
}
\details{
For the convex combination, any \code{\linkS4class{indepCopula}} in \code{components} is internally replaced by a Gaussian copula with a correlation of 0.01. For a pair at distance \eqn{h} between \eqn{h_i} and \eqn{h_{i+1}}, the geometric combination evaluates \eqn{C_i(u,v)^{p} C_{i+1}(u,v)^{1-p}} with \eqn{p = (h_{i+1}-h)/(h_{i+1}-h_i)}. Beyond the last distance, the last component is used.
}
\value{
An instance of the \code{\linkS4class{spatial_copula}} class.
}
\references{
Graeler, B. and E. Pebesma (2011): The pair-copula construction for spatial data: a new approach to model spatial dependency. Procedia Environmental Sciences (Vol. 7, pp. 206 - 211), Elsevier. 
}
\author{
Benedikt Graeler
}

\seealso{
See also \code{\link{calc_bins}} or \code{\link{loglik_by_lags}} and \code{demo("spatial_copula", package = "sfcopula")}. Spatial copulas form the trees of a \code{\link{distance_vine_copula}}.
}
\examples{
data("spCopDemo")

calcKTauPol <- fit_cor_fun(bins, degree=3)

spCop <- spatial_copula(components=list(normalCopula(), tCopula(),
                                  frankCopula(), normalCopula(), claytonCopula(),
                                  claytonCopula(), claytonCopula(), claytonCopula(),
                                  claytonCopula(), indepCopula()),
                  distances=bins$meanDists,
                  spDepFun=calcKTauPol, unit="m")

# a spatial copula based on geometric means
spGeomCop <- spatial_copula(list(normalCopula(0.2), gumbelCopula(1.2), claytonCopula(2)),
                            c(0, 120, 300), combination = "geometric")
spGeomCop
dCopula(c(0.3, 0.7), spGeomCop, h = 200)
}
\keyword{sfcopula}
\keyword{copula}

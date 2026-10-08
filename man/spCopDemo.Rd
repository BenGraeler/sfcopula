\name{spCopDemo}
\alias{spCopDemo}
\alias{bestFitTau}
\alias{bins}
\alias{calcKTauPol}
\alias{spCop}
\alias{meuseSpVine}
\docType{data}
\title{
workspace produced in \code{demo("spatial_copula", package = "sfcopula")}
}
\description{
This workspace is produced following \code{demo("spatial_copula", package = "sfcopula")} (see the script \code{data-raw/spCopDemo.R} in the package sources) and used in a couple of examples throughout the package. It is based on the zinc concentrations of the \code{\link{meuse}} data set, transformed to the unit interval with a fitted log-normal distribution (\code{marZinc}).
}
\usage{data("spCopDemo", package = "sfcopula")}
\format{
The workspace contains the following objects:
\describe{
  \item{\code{bins}}{the spatial bins of \code{marZinc} (10 lag classes up to 800 m) as produced by \code{\link{calc_bins}}.}
  \item{\code{calcKTauPol}}{a polynomial of degree 3 describing Kendall's tau as function of distance, as produced by \code{\link{fit_cor_fun}}.}
  \item{\code{bestFitTau}}{the indices of the best fitting copula families per lag class, as selected via \code{\link{loglik_by_lags}}.}
  \item{\code{spCop}}{the resulting bivariate spatial copula (\code{\linkS4class{spatial_copula}}).}
  \item{\code{meuseSpVine}}{a spatial vine copula (\code{\linkS4class{distance_vine_copula}}) fitted to neighbourhoods of size 5.}
}
}
\references{
Graeler, B. & E. Pebesma (2011): The pair-copula construction for spatial data: a new approach to model spatial dependency. Poster at: Spatial Statistics 2011 - Mapping global change. Enschede, The Netherlands, 23-25 March 2011. DOI: 10.1016/j.proenv.2011.07.036.}
\examples{
data("spCopDemo", package = "sfcopula")
ls()
bins$meanDists
spCop
}
\keyword{datasets}

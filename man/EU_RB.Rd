\name{EU_RB}
\alias{EU_RB}
\docType{data}
\title{
Daily mean PM10 concentrations over Europe in June and July 2005
}
\description{
Daily mean PM10 concentrations over Europe in June and July 2005
}
\usage{data("EU_RB")}
\format{
A \code{stars} vector data cube (see \code{\link[stars]{st_as_stars}}) with two dimensions:
\describe{
  \item{\code{geometry}}{194 rural background stations as POINT geometries (\code{sfc}) in a Lambert azimuthal equal-area projection (\code{+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m}; units: metres),}
  \item{\code{time}}{61 days (class \code{Date}) from 2005-06-01 to 2005-07-31,}
}
and three attributes, each a 194 x 61 array (stations x days):
\describe{
  \item{\code{PM10}}{daily mean PM10 concentration (numeric); missing observations (474 out of 11834) are \code{NA},}
  \item{\code{EMEP}}{daily mean PM10 concentration as modelled by the EMEP chemical transport model at the station locations (numeric),}
  \item{\code{station_altitude}}{the altitude of the station (integer), constant over time.}
}
}
\source{
Obtained from the european Air Qualtiy airbase http://acm.eionet.europa.eu/databases/airbase/.
}
\references{
http://acm.eionet.europa.eu/databases/airbase/

Graeler, B., L. E. Gerharz, & E. Pebesma (2012): Spatio-temporal analysis and interpolation of PM10 measurements in Europe. ETC/ACM Technical Paper 2011/10, January 2012.
http://acm.eionet.europa.eu/reports/ETCACM_TP_2011_10_spatio-temp_AQinterpolation
}
\examples{
library("stars")
data("EU_RB", package = "sfcopula")
EU_RB
dim(EU_RB[["PM10"]])
}
\keyword{datasets}

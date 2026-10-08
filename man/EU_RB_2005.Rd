\name{EU_RB_2005}
\alias{EU_RB_2005}
\docType{data}
\title{
Daily mean PM10 concentrations over Europe in 2005 as used in the JSS manuscript
}
\description{
Daily mean PM10 concentrations over Europe in 2005 as used in the JSS manuscript
}
\usage{data("EU_RB_2005")}
\format{
A \code{stars} vector data cube (see \code{\link[stars]{st_as_stars}}) with two dimensions:
\describe{
  \item{\code{geometry}}{194 rural background stations as POINT geometries (\code{sfc}) in a Lambert azimuthal equal-area projection (\code{+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m}; units: metres),}
  \item{\code{time}}{365 days (class \code{Date}) from 2005-01-01 to 2005-12-31,}
}
and the following attributes, each a 194 x 365 array (stations x days):
\describe{
  \item{\code{PM10}}{daily mean PM10 concentration (numeric); missing observations are \code{NA},}
  \item{\code{EMEP}}{daily mean PM10 concentration as modelled by the EMEP chemical transport model at the station locations (numeric),}
  \item{\code{logResidKrige}}{PM10 predictions on the original scale obtained by residual kriging on the log scale (numeric; some values are \code{NA}), as used in the JSS manuscript for comparison,}
  \item{\code{station_altitude}}{the altitude of the station (integer),}
  \item{\code{station_european_code}, \code{country_iso_code}, \code{station_start_date}, \code{station_end_date}, \code{type_of_station}, \code{station_type_of_area}, \code{street_type}}{station metadata from AirBase (character).}
}
All station related attributes are constant over time.
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
data("EU_RB_2005", package = "sfcopula")
EU_RB_2005["PM10"]

}
\keyword{datasets}

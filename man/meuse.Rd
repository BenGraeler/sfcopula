\name{meuse}
\alias{meuse}
\docType{data}
\title{
Meuse river data set as simple features
}
\description{
This data set gives locations and topsoil heavy metal concentrations, along with a number of soil and landscape variables at the observation locations, collected in a flood plain of the river Meuse, near the village of Stein (NL). Heavy metal concentrations are from composite samples of an area of approximately 15 m x 15 m. It is the classic \code{meuse} data set of the \pkg{sp} package converted to an \code{\link[sf]{sf}} object.
}
\usage{data("meuse", package = "sfcopula")}
\format{
An \code{\link[sf]{sf}} object with 155 POINT features in the Dutch Rijksdriehoek coordinate reference system (Amersfoort / RD New, EPSG:28992; units: metres) and the following attribute columns:
\describe{
  \item{\code{cadmium}}{topsoil cadmium concentration, mg kg-1 soil ("ppm"); zero cadmium values in the original data set have been shifted to 0.2 (half the lowest non-zero value)}
  \item{\code{copper}}{topsoil copper concentration, mg kg-1 soil ("ppm")}
  \item{\code{lead}}{topsoil lead concentration, mg kg-1 soil ("ppm")}
  \item{\code{zinc}}{topsoil zinc concentration, mg kg-1 soil ("ppm")}
  \item{\code{elev}}{relative elevation above local river bed, m}
  \item{\code{dist}}{distance to the Meuse; obtained from the nearest cell of a grid derived by a spread (spatial distance) GIS operation, horizontal precision 20 metres; then normalized to [0,1]}
  \item{\code{om}}{organic matter, kg (100 kg)-1 soil (percent)}
  \item{\code{ffreq}}{flooding frequency class (factor): 1 = once in two years; 2 = once in ten years; 3 = one in 50 years}
  \item{\code{soil}}{soil type (factor) according to the 1:50 000 soil map of the Netherlands. 1 = Rd10A (Calcareous weakly-developed meadow soils, light sandy clay); 2 = Rd90C/VII (Non-calcareous weakly-developed meadow soils, heavy sandy clay to light clay); 3 = Bkd26/VII (Red Brick soil, fine-sandy, silty light clay)}
  \item{\code{lime}}{lime class (factor): 0 = absent, 1 = present by field test with 5\% HCl}
  \item{\code{landuse}}{landuse class (factor): Aa Agriculture/unspecified, Ab = Agr/sugar beetsm, Ag = Agr/small grains, Ah = Agr/??, Am = Agr/maize, B = woods, Bw = trees in pasture, DEN = ??, Fh = tall fruit trees, Fl = low fruit trees; Fw = fruit trees in pasture, Ga = home gardens, SPO = sport field, STA = stable yard, Tv = ?? , W = pasture}
  \item{\code{dist.m}}{distance to river Meuse in metres, as obtained during the field survey}
  \item{\code{geometry}}{the point locations (\code{sfc_POINT})}
}
}
\details{
Row names refer to the original sample number. As with \code{sp::coordinates()}, the coordinates are only part of the geometry (use \code{sf::st_coordinates(meuse)} to retrieve them); hence \code{cadmium} is the first attribute column and the default variable of \code{\link{neighbours}}.
}
\source{
The \code{meuse} data set of the \pkg{sp} package (Pebesma and Bivand). Field data were collected by Ruud van Rijn and Mathieu Rikken; compiled for R by Edzer Pebesma; description extended by David Rossiter.
}
\references{
M G J Rikken and R P G Van Rijn, 1993. Soil pollution with heavy metals - an inquiry into spatial variation, cost of mapping and the risk evaluation of copper, cadmium, lead and zinc in the floodplains of the Meuse west of Stein, the Netherlands. Doctoraalveldwerkverslag, Dept. of Physical Geography, Utrecht University

P.A. Burrough, R.A. McDonnell, 1998. Principles of Geographical Information Systems. Oxford University Press.

Stichting voor Bodemkartering (STIBOKA), 1970. Bodemkaart van Nederland : Blad 59 Peer, Blad 60 West en 60 Oost Sittard: schaal 1 : 50 000. Wageningen, STIBOKA.
}
\examples{
library("sf")
data("meuse", package = "sfcopula")
summary(meuse)
st_crs(meuse)$epsg
}
\keyword{datasets}

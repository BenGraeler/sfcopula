\name{predict-methods}
\docType{methods}
\alias{predict}
\alias{predict-methods}
\alias{predict,distance_vine_copula-method}
\alias{predict,spatial_gauss_copula-method}
\alias{predict,covariate_vine_copula-method}

\title{
Spatial and spatio-temporal prediction based on copulas
}
\description{
A spatial or spatio-temporal (vine) copula is used to predict values at unobserved (spatio-temporal) locations conditioned on observations of a local neighbourhood.
}
\usage{
\S4method{predict}{distance_vine_copula}(object, neigh, data, target, margin,
        method = c("quantile", "expectation"), p = 0.5, ...)

\S4method{predict}{spatial_gauss_copula}(object, neigh, data, target, margin,
        method = c("quantile", "expectation"), p = 0.5, ..., n = 1000)

\S4method{predict}{covariate_vine_copula}(object, neigh, data, target, margin,
        method = c("quantile", "expectation"), p = 0.5, ...)
}
\arguments{
  \item{object}{the copula describing the (spatio-temporal) dependence: a spatial or spatio-temporal \code{\linkS4class{distance_vine_copula}}, a \code{\linkS4class{spatial_gauss_copula}} or a \code{\linkS4class{covariate_vine_copula}}.}
  \item{neigh}{the \code{\linkS4class{neighbourhood}} used for prediction, its slot \code{prediction} must be \code{TRUE} (see \code{\link{neighbours}}). For a \code{\linkS4class{covariate_vine_copula}} (spatial or spatio-temporal tree), it needs to hold a single covariate which is stored in the last column of its data; for prediction, the covariate is taken from the target (see \code{\link{neighbours}}).}
  \item{data}{the data used in the generation of \code{neigh}: an \code{\link[sf]{sf}} (or \code{sfc}) object with POINT geometries in the spatial case, the \code{stars} vector data cube in the spatio-temporal case.}
  \item{target}{the same prediction targets as used in the generation of \code{neigh}: an \code{\link[sf]{sf}} or \code{sfc} object with POINT geometries in the spatial case; either a \code{stars} vector data cube or an \code{\link[sf]{sf}} object with a time column in the spatio-temporal case.}
  \item{margin}{the marginal distribution as a list with entries named "d" for the density function (PDF), "q" for the quantile function and "p" for cumulative distribution function (CDF). Currently only the entry \code{q} is required.}
  \item{method}{one of \code{"quantile"} or \code{"expectation"} denoting the type of predictor.}
  \item{p}{only used for the quantile predictor indicating the desired fraction the quantile should correspond to. With \code{p = "random"}, a random fraction is drawn for each location (e.g. for conditional simulation).}
  \item{...}{Further arguments passed to \code{\link{integrate}} for the expectation.}
  \item{n}{the approximate number of points used in the linear approximation of the conditional distribution function.}
}
\details{
For each neighbourhood, the conditional distribution of the central location is derived from the copula (see \code{\link{cond_vine}}) and transformed with the marginal quantile function through numerical integration/optimisation. 

For the \code{\linkS4class{spatial_gauss_copula}}, a Gaussian copula (\code{\linkS4class{normalCopula}}) is generated for each neighbourhood based on the correlation function applied to the distance matrix of the neighbourhood and conditioned under the data of the neighbouring locations.
}
\value{
The prediction targets \code{target} extended by the predictions named \code{"expect"} for the expectation or \code{"quantile.<p>"} (e.g. \code{"quantile.0.5"}) for the quantile predictor: an additional column for an \code{sf} object (a bare \code{sfc} \code{target} is turned into an \code{sf} object with this single column) or an additional attribute (location x time) for a \code{stars} vector data cube.
}
\references{
Graeler, B. and E. Pebesma (2011): The pair-copula construction for spatial data: a new approach to model spatial dependency. Procedia Environmental Sciences (Vol. 7, pp. 206 - 211), Elsevier.
}
\author{
Benedikt Graeler
}

\seealso{
\code{\link{cond_vine}}, \code{\link{neighbours}}
}
\examples{
library("sf")
data("meuse", package = "sfcopula")

# a few prediction locations close to the first observations
predLocs <- st_as_sf(expand.grid(x = c(180900, 181000, 181100, 181200),
                                 y = c(333400, 333500, 333600)),
                     coords = c("x", "y"), crs = 28992)

data("spCopDemo", package = "sfcopula")

calcKTauPol <- fit_cor_fun(bins, degree=3)

spCop <- spatial_copula(components=list(normalCopula(), tCopula(),
                                  frankCopula(), normalCopula(), claytonCopula(),
                                  claytonCopula(), claytonCopula(), claytonCopula(),
                                  claytonCopula(), indepCopula()),
                  distances=bins$meanDists,
                  spDepFun=calcKTauPol, unit="m")

spVineCop <- distance_vine_copula(spCop, vineCopula(4L))

meuse$rtZinc <- rank(meuse$zinc)/(nrow(meuse)+1)

dataLocs <- meuse[1:4,]
predMeuseNeigh <- neighbours(dataLocs, predLocs, 5, "rtZinc", min.dist=-1)

qMar <- function(x) {
  qlnorm(x,mean(log(meuse$zinc)),sd(log(meuse$zinc)))
}

predMedian <- predict(spVineCop, predMeuseNeigh, dataLocs, predLocs,
                      list(q=qMar), "quantile", p=0.5)
predMedian

## spatial Gaussian copula

# translate Kendall's tau correlation function into Gaussian Copula parameters 
# using a linear variogram
calcKTauLin <- fit_cor_fun(bins, degree=1)
meuseGaussCorFun <- function(h) {
  res <- pmax(iTau(normalCopula(0), calcKTauLin(0))/658*(658-h),0)
  res[h ==0] <- 1
  return(res)
}

predGauss <- predict(spatial_gauss_copula(meuseGaussCorFun), predMeuseNeigh,
                     dataLocs, predLocs, list(q=qMar))
predGauss$quantile.0.5

\dontrun{
  plot(predMedian["quantile.0.5"], pch = 15, cex = 3, reset = FALSE)
  plot(st_geometry(meuse), add = TRUE, pch = 19, col = "red")
}

## spatio-temporal vine copula
library("stars")

spCopT0 <- spatial_copula(components=list(claytonCopula(8), claytonCopula(4), 
                                    claytonCopula(2), claytonCopula(1),
                                    claytonCopula(0.5), indepCopula()),
                    distances=c(0,100,200,300,400,500))
spCopT1 <- spatial_copula(components=list(claytonCopula(4), claytonCopula(2), 
                                    claytonCopula(1), claytonCopula(0.5),
                                    indepCopula()),
                    distances=c(0,100,200,300,400))

stCop <- spacetime_copula(components=list(spCopT0, spCopT1),
                  tlags=-(0:1))

stVineCop <- distance_vine_copula(stCop, vineCopula(4L))

pts <- st_sfc(st_point(c(181000,333500)), st_point(c(181100,333600)),
              crs = 28992)
time <- as.Date("2014-03-18") + 0:2

stData <- st_as_stars(list(var1 = matrix(runif(6), 2, 3)),
                      dimensions = st_dimensions(geometry = pts, time = time))

# prediction targets: all combinations of two locations and two days
stQuerry <- st_as_stars(list(dummy = matrix(NA_real_, 2, 2)),
                        dimensions = st_dimensions(
                          geometry = st_sfc(st_point(c(181000,333600)),
                                            st_point(c(181200,333600)),
                                            crs = 28992),
                          time = time[2:3]))

stNeigh <- neighbours(stData, stQuerry, var="var1", size=3, tlags=-(0:1))

predST <- predict(stVineCop, stNeigh, stData, stQuerry, list(q=qunif),
                  "quantile", 0.5)
predST[["quantile.0.5"]]

## spatio-temporal covariate vine copula
# the covariate of the targets is taken from the target cube
stData$cov <- matrix(runif(6), 2, 3)
stQuerry$cov <- matrix(runif(4), 2, 2)
stCoNeigh <- neighbours(stData, stQuerry["cov"], var="var1", coVar="cov",
                        size=3, tlags=-(0:1))

stCVVC <- covariate_vine_copula(function(stInd) normalCopula(0.5), stCop,
                                      vineCopula(5L))
predict(stCVVC, stCoNeigh, stData, stQuerry["cov"], list(q=qunif))[["quantile.0.5"]]
}

\keyword{ distribution }
\keyword{ prediction }
\keyword{ spatial }

\name{cond_vine}
\alias{cond_vine}
\alias{cond_vine-methods}
\alias{cond_vine,ANY,ANY,distance_vine_copula-method}
\alias{cond_vine,ANY,ANY,covariate_vine_copula-method}

\title{
Conditions a spatial or spatio-temporal vine copula for conditional prediction
}
\description{
A spatial or spatio-temporal vine copula is conditioned under the observations of all but one neighbour generating a conditional univariate distribution used for prediction.
}
\usage{
cond_vine(cond_var, dists, vine, n = 1000, ...)
}
\arguments{
  \item{cond_var}{
the conditioning variables, i.e. the (uniform) data of the neighbours. For a \code{\linkS4class{covariate_vine_copula}}, the neighbours are followed by the covariate at the central location, i.e. the order of the columns of a \code{\linkS4class{neighbourhood}} (\pkg{spcopula}'s \code{condStCoVarVine} expected the covariate first).
}
  \item{dists}{
the distances between the neighbours to adjust the copula parameters. For a \code{\linkS4class{distance_vine_copula}}, a list with the distances per tree (see \code{\link{tree_dists}}): a vector (or a single row matrix) for spatial trees and an array \code{[1, pairs, 2]} holding the spatial and temporal distances for spatio-temporal trees. For a \code{\linkS4class{covariate_vine_copula}}, the distances of its tree: a vector for a spatial tree, an array \code{[1, k, 2]} for a spatio-temporal tree.
}
  \item{vine}{
the vine copula, either a \code{\linkS4class{distance_vine_copula}} or a \code{\linkS4class{covariate_vine_copula}}
}
  \item{n}{
a proxy to the number of intervals to be used in the numerical integration. The intervals will not be split uniform with shorter intervals at the corners of the copula.
}
  \item{\dots}{further arguments of the methods: the \code{\linkS4class{covariate_vine_copula}} method needs \code{centre}, the index of the central location passed to the covariate copula function (\code{covar_cop}) of the vine: the location index for a spatial tree, the pair (location, time) for a spatio-temporal tree.}
}
\section{Methods}{
\describe{
\item{\code{signature(vine = "distance_vine_copula")}}{conditions a spatial or spatio-temporal \code{\linkS4class{distance_vine_copula}}.}
\item{\code{signature(vine = "covariate_vine_copula")}}{conditions a \code{\linkS4class{covariate_vine_copula}}; requires the argument \code{centre}.}
}}
\value{
A function describing the conditional density over [0,1]. Its attribute \code{"xVals"} holds the points of the linear approximation.
}
\note{
The density is linearly approximated at a limited number (approximately \code{n}) of points.
}
\author{
Benedikt Graeler
}

\seealso{
\code{\linkS4class{distance_vine_copula}}, \code{\link[=predict,distance_vine_copula-method]{predict}}
}
\examples{
## spatial vine copula
data("spCopDemo")

calcKTauPol <- fit_cor_fun(bins, degree=3)

spCop <- spatial_copula(components=list(normalCopula(0.2), tCopula(0.2, dispstr = "un"),
                                  frankCopula(1.2), normalCopula(0.2), claytonCopula(0.2),
                                  claytonCopula(0.2), claytonCopula(0.2), claytonCopula(0.2),
                                  claytonCopula(0.2), indepCopula()),
                  distances=c(0, bins$meanDists[1:9]),
                  dep_fun=calcKTauPol, unit="m")

spVineCop <- distance_vine_copula(spCop, vineCopula(4L))

dists <- list(c(473, 124, 116, 649))
cond_var <- c(0.29, 0.55, 0.05, 0.41)
condDensity <- cond_vine(cond_var, dists, spVineCop)

curve(condDensity)
mtext(paste("Dists:",paste(round(dists[[1]],0),collapse=", ")),line=0)
mtext(paste("Cond.:",paste(round(cond_var,2),collapse=", ")),line=1)

## spatio-temporal vine copula
spCopT0 <- spatial_copula(components=list(claytonCopula(8), claytonCopula(4), 
                                    claytonCopula(2), claytonCopula(1),
                                    claytonCopula(0.5), indepCopula()),
                    distances=c(0,100,200,300,400,500))
spCopT1 <- spatial_copula(components=list(claytonCopula(4), claytonCopula(2), 
                                    claytonCopula(1), claytonCopula(0.5),
                                    indepCopula()),
                    distances=c(0,100,200,300,400))
spCopT2 <- spatial_copula(components=list(claytonCopula(2), claytonCopula(1), 
                                    claytonCopula(0.5), indepCopula()),
                    distances=c(0,100,200,300))

stCop <- spacetime_copula(components=list(spCopT0, spCopT1),
                  tlags=-(0:1))

stVineCop <- distance_vine_copula(stCop, vineCopula(4L))

dists <- array(c(150, 250, 150, 250,0,0,-1,-1),dim=c(1,4,2))

condDensity <- cond_vine(cond_var, list(dists), stVineCop)
curve(condDensity)

## spatio-temporal covariate vine copula
stCop3 <- spacetime_copula(components=list(spCopT0, spCopT1, spCopT2),
                   tlags=-(0:2))

# only a constant copula is used for the covariate
stCVVC <- covariate_vine_copula(function(x) gumbelCopula(7), stCop3, vineCopula(5L))

# the neighbours followed by the covariate
cond_var <- c(0.29, 0.55, 0.05, 0.41, 0.95)

condDensity <- cond_vine(cond_var, dists, stCVVC, centre = c(1,1))
curve(condDensity)
}

\keyword{ distribution }

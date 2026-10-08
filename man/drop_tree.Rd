\name{drop_tree}
\alias{drop_tree}

\title{
Drops a spatial or spatio-temporal tree
}
\description{
A spatial or spatio-temporal tree is dropped in order to fit/evaluate the following copula. The data of the neighbourhood are conditioned on the central location using the partial derivatives of the tree's bivariate copula.
}
\usage{drop_tree(neigh, data, copula)}

\arguments{
  \item{neigh}{the current \code{\linkS4class{neighbourhood}} (holding the data of the central location, i.e. for fitting)}
  \item{data}{the locations of the data already used to generate \code{neigh}: an \code{\link[sf]{sf}} or \code{sfc} object in the spatial case; the \code{stars} vector data cube itself or its \code{sfc} geometries in the spatio-temporal case.}
  \item{copula}{the current bivariate copula performing the conditioning: a \code{\linkS4class{spatial_copula}} for spatial and a \code{\linkS4class{spacetime_copula}} for spatio-temporal neighbourhoods.}
}

\details{
Covariates are not carried over to the conditioned neighbourhood; condition them on the central location with \code{\link{cond_covariate}}.

The first neighbour of \code{neigh} becomes the central location of the conditioned neighbourhood. Its distances are the (spatial and temporal) distances between this new central location and the remaining neighbours, i.e. the distances needed for the next tree (see \code{\link{tree_dists}}).
}

\value{
A conditioned \code{\linkS4class{neighbourhood}} of dimension 1 smaller than the current one.
}

\author{
Benedikt Graeler
}

\seealso{
\code{\linkS4class{neighbourhood}}, \code{\link{tree_dists}}, \code{\linkS4class{distance_vine_copula}}
}

\examples{
library("sf")
data("meuse", package = "sfcopula")
meuse$marZinc <- plnorm(meuse$zinc, mean(log(meuse$zinc)), sd(log(meuse$zinc)))

meuseNeigh <- neighbours(meuse, var = "marZinc", size = 4)
spCop <- spatial_copula(list(claytonCopula(4), gumbelCopula(2), indepCopula()),
                        c(0, 200, 500))

condNeigh <- drop_tree(meuseNeigh, meuse, spCop)
condNeigh
head(condNeigh@data)

# a spatio-temporal neighbourhood
library("stars")
pts <- st_sfc(st_point(c(0, 0)), st_point(c(100, 0)), st_point(c(0, 100)),
              st_point(c(100, 100)), crs = 28992)
cube <- st_as_stars(list(u = matrix(runif(16), 4, 4)),
                    dimensions = st_dimensions(geometry = pts,
                                               time = as.Date("2020-01-01") + 0:3))
stNeigh <- neighbours(cube, var = "u", size = 3, tlags = -(0:1))
stCop <- spacetime_copula(list(spatial_copula(list(claytonCopula(2), indepCopula()), c(0, 300)),
                               spatial_copula(list(claytonCopula(1), indepCopula()), c(0, 300))),
                          tlags = -(0:1))
drop_tree(stNeigh, cube, stCop)
}

\keyword{ spatial }
\keyword{ spatio-temporal }
\keyword{ multivariate }
\keyword{distribution}

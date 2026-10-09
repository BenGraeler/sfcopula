\name{covariate_vine_copula}
\alias{covariate_vine_copula}

\title{
Covariate vine copula with a spatial or spatio-temporal tree
}
\description{
Defines a \code{\linkS4class{covariate_vine_copula}}: a C-vine copula for a variable of interest at a central location, its spatial or spatio-temporal neighbours and a covariate observed at the central location.
}
\usage{
covariate_vine_copula(covar_cop, tree, top_cop, ...)
}

\arguments{
  \item{covar_cop}{
A function returning the bivariate copula between the variable of interest and the covariate at the central location. It receives the index of the central location: the location index for a spatial tree and the pair (location, time) for a spatio-temporal tree. Indices refer to the data when fitting and to the target when predicting (see \code{\link{neighbours}}).
}
  \item{tree}{
A \code{\linkS4class{spatial_copula}} or \code{\linkS4class{spacetime_copula}} modelling the dependence between the central location and its neighbours on the first tree.
}
  \item{top_cop}{
Object of class \code{"copula"} joining the conditioned covariate and neighbours to a full vine copula; its dimension is the number of neighbours plus one.
}
  \item{\dots}{deprecated argument names of \pkg{spcopula} (\code{coVarCop}, \code{topCop}, \code{stCop}) are accepted with a warning.}
}
\details{
For a random field Z with covariate Y a C-vine is assumed: the first tree couples the central value z_0 with its neighbours z_1, ..., z_n through \code{tree} and with the covariate y_0 through \code{covar_cop}; \code{top_cop} joins the conditioned values (y_0|z_0, z_1|z_0, ..., z_n|z_0).

The density is evaluated with \code{dCopula(u, copula, h, centre, log = FALSE)}, where the columns of \code{u} are ordered as in a \code{\linkS4class{neighbourhood}}: (z_0, z_1, ..., z_n, y_0). \code{h} holds the distances of the neighbours (spatial: matrix [n or 1, number of neighbours]; spatio-temporal: array [n or 1, number of neighbours, 2]) and \code{centre} the index of the central location passed to \code{covar_cop} (spatial: vector of length n or 1; spatio-temporal: matrix [n or 1, 2]).

\code{fitCopula(copula, neigh, method)} takes a \code{\linkS4class{neighbourhood}} with one covariate (see \code{\link{neighbours}}), drops the tree, conditions the covariate (see \code{\link{cond_covariate}}) and fits \code{top_cop} to the conditioned values; the returned log-likelihood covers all three parts.
}

\value{An object of class \code{\linkS4class{covariate_vine_copula}}.}

\author{
Benedikt Graeler
}

\seealso{
\code{\link{distance_vine_copula}}, \code{\link{cond_vine}}, \code{\link[=predict,covariate_vine_copula-method]{predict}}
}

\examples{
library("sf")
data("meuse", package = "sfcopula")
meuse$marZinc <- rank(meuse$zinc) / (nrow(meuse) + 1)
meuse$marCopper <- rank(meuse$copper) / (nrow(meuse) + 1)

# spatial tree, a constant copula for the covariate
tree <- spatial_copula(list(claytonCopula(2), indepCopula()), distances = c(0, 800))
cvc <- covariate_vine_copula(function(loc) normalCopula(0.8), tree, vineCopula(3L))

neigh <- neighbours(meuse, var = "marZinc", covar = "marCopper", size = 3)
fit <- fitCopula(cvc, neigh, method = list(familyset = 1:4))
fit@loglik

# prediction: the covariate is taken from the target
target <- meuse[150:155, ]
predNeigh <- neighbours(meuse[1:149, ], target, var = "marZinc", covar = "marCopper", size = 3)
pred <- predict(fit@copula, predNeigh, meuse[1:149, ], target, 
                margin = list(q = function(p) quantile(meuse$zinc, p)))
pred[c("zinc", "quantile.0.5")]

# a spatio-temporal tree
spCopT0 <- spatial_copula(list(claytonCopula(4), indepCopula()), distances = c(0, 500), unit = "km")
spCopT1 <- spatial_copula(list(claytonCopula(2), indepCopula()), distances = c(0, 500), unit = "km")
stCop <- spacetime_copula(components = list(spCopT0, spCopT1), tlags = -(0:1))
covariate_vine_copula(function(stInd) gumbelCopula(2), stCop, vineCopula(5L))
}
\keyword{ distribution }

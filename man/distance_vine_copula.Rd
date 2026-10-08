\name{distance_vine_copula}
\alias{distance_vine_copula}

\title{
Constructor of the class \code{\linkS4class{distance_vine_copula}}
}
\description{
Construct an instance of the \code{\linkS4class{distance_vine_copula}}: a C-vine copula whose lower trees are distance dependent bivariate copulas, either spatial (\code{\linkS4class{spatial_copula}}) or spatio-temporal (\code{\linkS4class{spacetime_copula}}), optionally joined by a non distance dependent copula at the top.
}
\usage{
distance_vine_copula(trees, topCop = NULL)
}
\arguments{
  \item{trees}{A list of objects of class \code{\linkS4class{spatial_copula}} or a list of objects of class \code{\linkS4class{spacetime_copula}} (mixing both is not supported) describing the distance dependent trees of the vine. A single bivariate copula is wrapped into a list.}
  \item{topCop}{Object of class \code{\linkS4class{copula}}. A bivariate, multivariate or R-Vine copula joining the conditioned components (neighbours) remaining after the distance dependent trees. In case all trees are distance dependent copulas, this argument is \code{NULL} (default).}
}
\details{
The dimension of the vine is the number of trees plus one if \code{topCop} is \code{NULL}, and the dimension of \code{topCop} plus the number of trees otherwise. The data of a neighbourhood are ordered as the central location followed by its neighbours (see \code{\link{neighbours}}).
}
\value{
An instance of the \code{\linkS4class{distance_vine_copula}} class.
}
\references{
Graeler, B. and E. Pebesma (2011): The pair-copula construction for spatial data: a new approach to model spatial dependency. Procedia Environmental Sciences (Vol. 7, pp. 206 - 211), Elsevier.
}
\author{
Benedikt Graeler
}
\seealso{
\code{\link{neighbours}}, \code{\link{cond_vine}} and \code{\link[=predict,distance_vine_copula-method]{predict}}
}
\examples{
# a spatial C-vine copula (with independent dummy copulas in the upper vine)
data("spCopDemo")

calcKTauPol <- fit_cor_fun(bins, degree=3)

spCop <- spatial_copula(components=list(normalCopula(0.2), normalCopula(0.2), 
                                  tCopula(0.2, dispstr = "un"),
                                  frankCopula(1.2), normalCopula(0.2), 
                                  claytonCopula(0.2), claytonCopula(0.2),
                                  claytonCopula(0.2), claytonCopula(0.2),
                                  claytonCopula(0.2), indepCopula()),
                  distances=c(0,bins$meanDists),
                  spDepFun=calcKTauPol, unit="m")

# VineCopula is not attached to avoid masking vineCopula() from VC2copula
RVM <- VineCopula::RVineMatrix(matrix(c(1,0,0,2,2,0,3,3,3),3,byrow=TRUE))
distance_vine_copula(spCop, vineCopula(RVM))

# a vine of two spatial trees only
tree1 <- spatial_copula(list(claytonCopula(3), indepCopula()), c(0, 800))
tree2 <- spatial_copula(list(gumbelCopula(1.5), indepCopula()), c(0, 800))
distance_vine_copula(list(tree1, tree2))

# a spatio-temporal C-vine copula (with independent copulas in the upper vine)
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

stCop <- spacetime_copula(components=list(spCopT0, spCopT1, spCopT2),
                  tlags=-(0:2))

distance_vine_copula(stCop, vineCopula(9L))
}
\keyword{ mulitvariate }
\keyword{ distribution }

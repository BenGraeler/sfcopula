\name{distance_vine_copula-class}
\Rdversion{1.1}
\docType{class}
\alias{distance_vine_copula-class}
\alias{optionalCopula-class}
\alias{show,distance_vine_copula-method}
\alias{dCopula,matrix,distance_vine_copula-method}
\alias{dCopula,numeric,distance_vine_copula-method}
\alias{dCopula,data.frame,distance_vine_copula-method}
\alias{fitCopula,distance_vine_copula,ANY-method}
\alias{rCopula,numeric,distance_vine_copula-method}

\title{Class \code{"distance_vine_copula"}}
\description{
A S4-class representing a spatial or spatio-temporal vine copula structure to model spatial or spatio-temporal processes. The lower trees of this C-vine are distance dependent bivariate copulas, either all spatial (\code{\linkS4class{spatial_copula}}) or all spatio-temporal (\code{\linkS4class{spacetime_copula}}). The vine might be composed purely out of distance dependent trees or being mixed as lower distance dependent trees with an upper non distance dependent joint copula.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("distance_vine_copula", ...)} or to the constructor \code{\link{distance_vine_copula}}.
}
\section{Slots}{
  \describe{
    \item{\code{trees}:}{Object of class \code{"list"} holding the bivariate distance dependent copulas of the lower trees, either \code{\linkS4class{spatial_copula}} or \code{\linkS4class{spacetime_copula}} objects. }
    \item{\code{topCop}:}{Object of class \code{"optionalCopula"}, i.e. a \code{\linkS4class{copula}} or \code{NULL}. The top copula being constant over space (and time). It is either a R-Vine (\code{\linkS4class{vineCopula}}) or a bivariate copula (\code{\linkS4class{copula}}) joining the conditioned components (neighbours) or all trees of the vine are distance dependent bivariate copulas (\code{NULL}).}
    \item{\code{dimension}:}{Object of class \code{"integer"} giving the dimension of the vine copula (size of the neighbourhood). }
    \item{\code{parameters}:}{Object of class \code{"numeric"} empty }
    \item{\code{param.names}:}{Object of class \code{"character"} empty }
    \item{\code{param.lowbnd}:}{Object of class \code{"numeric"} empty }
    \item{\code{param.upbnd}:}{Object of class \code{"numeric"} empty }
    \item{\code{fullname}:}{Object of class \code{"character"} A descriptive name of the this copula family.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{copula}"}, directly.
Class \code{"\linkS4class{Copula}"}, by class "copula", distance 2.
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "distance_vine_copula")}: a short description of the vine.}
    \item{dCopula}{\code{signature(u = "matrix", copula = "distance_vine_copula")}: the density of the vine, also for \code{u} being \code{"numeric"} or \code{"data.frame"}. The additional argument \code{h} is a list with the distances of each tree: a matrix with one row per row of \code{u} and one column per pair of the tree (or a vector for a single row) for spatial trees, and an array \code{[n, pairs, 2]} holding the spatial and temporal distances for spatio-temporal trees (see \code{\link{tree_dists}}). }
    \item{fitCopula}{\code{signature(copula = "distance_vine_copula")}: fits the top copula of the vine for given distance dependent trees. The argument \code{data} is a \code{\linkS4class{neighbourhood}} or, if the vine has more than one tree, a list holding the \code{\linkS4class{neighbourhood}} and the data locations used to create it (\code{list(neigh, dataLocs)}). After dropping all distance dependent trees, a bivariate copula is selected by \code{\link[VineCopula]{BiCopSelect}} if two dimensions remain; for more than two dimensions \code{topCop} is fitted (a \code{\linkS4class{vineCopula}} if \code{topCop} is \code{NULL}) passing on \code{method}, which defaults to \code{list(StructureSelect = FALSE, indeptest = FALSE, familyset = NA)}. Returns a \code{\linkS4class{fitCopula}} object whose copula is the completed \code{distance_vine_copula}.}
    \item{rCopula}{\code{signature(n = "numeric", copula = "distance_vine_copula")}: simulates \code{n} samples from a vine without top copula; the additional argument \code{h} is a list with the distances of each tree: vectors for spatial trees and matrices with columns space and time for spatio-temporal trees.}
    \item{cond_vine}{\code{signature(vine = "distance_vine_copula")}: see \code{\link{cond_vine}}.}
    \item{predict}{\code{signature(object = "distance_vine_copula")}: see \code{\link[=predict,distance_vine_copula-method]{predict}}.}
	 }
}
\references{
Graeler, B. and E. Pebesma (2011): The pair-copula construction for spatial data: a new approach to model spatial dependency. Procedia Environmental Sciences (Vol. 7, pp. 206 - 211), Elsevier.
}
\author{
Benedikt Graeler
}

\seealso{
\code{\linkS4class{spatial_copula}}, \code{\linkS4class{spacetime_copula}}, \code{\linkS4class{vineCopula}}
}
\examples{
showClass("distance_vine_copula")

library("sf")
data("meuse", package = "sfcopula")
meuse$marZinc <- plnorm(meuse$zinc, mean(log(meuse$zinc)), sd(log(meuse$zinc)))

tree1 <- spatial_copula(list(claytonCopula(3), indepCopula()), c(0, 800))
tree2 <- spatial_copula(list(gumbelCopula(1.5), indepCopula()), c(0, 800))

# a pure spatial vine of two trees
pureVine <- distance_vine_copula(list(tree1, tree2))
dCopula(c(0.3, 0.4, 0.5), pureVine, h = list(c(50, 100), 80))
rCopula(3, pureVine, h = list(c(50, 50), 50))

# fit the top copula of a spatial vine with a single tree
meuseNeigh <- neighbours(meuse, var = "marZinc", size = 3)
fit <- fitCopula(distance_vine_copula(tree1, normalCopula(0)), meuseNeigh)
fit@copula@topCop

# vines with more than one tree need the data locations as well
fitCopula(pureVine, list(meuseNeigh, meuse))@loglik
}
\keyword{classes}

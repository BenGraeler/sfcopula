\name{spatial_copula-class}
\Rdversion{1.1}
\docType{class}
\alias{spatial_copula-class}
\alias{show,spatial_copula-method}
\alias{pCopula,numeric,spatial_copula-method}
\alias{pCopula,matrix,spatial_copula-method}
\alias{dCopula,numeric,spatial_copula-method}
\alias{dCopula,matrix,spatial_copula-method}
\alias{dduCopula,matrix,spatial_copula-method}
\alias{ddvCopula,matrix,spatial_copula-method}
\alias{dduCopula,numeric,spatial_copula-method}
\alias{ddvCopula,numeric,spatial_copula-method}
\alias{rCopula,numeric,spatial_copula-method}

\title{Class \code{"spatial_copula"}}
\description{
A class representing a bivariate spatial copula, i.e. a distance dependent convex combination or geometric mean of bivariate copulas.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("spatial_copula", ...)} or by the function \code{\link{spatial_copula}}.
}

\section{Slots}{
  \describe{
    \item{\code{components}:}{Object of class \code{"list"} A list of copula components per lag.}
    \item{\code{distances}:}{Object of class \code{"numeric"} The mean distances of the lag classes. }
    \item{\code{calibMoa}:}{Object of class \code{"function"} A function providing the copulas' parameter based on a measure of association.}
    \item{\code{unit}:}{Object of class \code{"character"} The spatial unit. }
    \item{\code{combination}:}{Object of class \code{"character"} How the components are combined, either \code{"convex"} or \code{"geometric"}. }
    \item{\code{dimension}:}{Object of class \code{"integer"} The dimension of the copula (2). }
    \item{\code{parameters}:}{Object of class \code{"numeric"} The parameters of the components. }
    \item{\code{param.names}:}{Object of class \code{"character"} Names of the above parameters. }
    \item{\code{param.lowbnd}:}{Object of class \code{"numeric"} The lower bounds of the above parameters. }
    \item{\code{param.upbnd}:}{Object of class \code{"numeric"} The upper bounds of the above parameters. }
    \item{\code{fullname}:}{Object of class \code{"character"}A phrase describing the spatial copula.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{copula}"}, directly.
Class \code{"\linkS4class{Copula}"}, by class "copula", distance 2.
}
\section{Methods}{
All methods require the additional argument \code{h} providing the distance(s) of the pair(s).
  \describe{
    \item{show}{\code{signature(object = "spatial_copula")}: lists the components and their distances. }
    \item{pCopula}{\code{signature(u = "numeric", copula = "spatial_copula")}: the distribution function, also for \code{u = "matrix"}. }
    \item{dCopula}{\code{signature(u = "numeric", copula = "spatial_copula")}: the density, also for \code{u = "matrix"}. }
    \item{dduCopula}{\code{signature(u = "numeric", copula = "spatial_copula")}: the partial derivative with respect to \code{u}, also for \code{u = "matrix"}. }
    \item{ddvCopula}{\code{signature(u = "numeric", copula = "spatial_copula")}: the partial derivative with respect to \code{v}, also for \code{u = "matrix"}. }
    \item{invdduCopula}{\code{signature(u = "numeric", copula = "spatial_copula", y = "ANY")}: see \code{\link{invdduCopula}}. }
    \item{invddvCopula}{\code{signature(v = "numeric", copula = "spatial_copula", y = "ANY")}: see \code{\link{invddvCopula}}. }
    \item{rCopula}{\code{signature(n = "numeric", copula = "spatial_copula")}: random samples for the distances \code{h}. }
	 }
}
\references{
Graeler, B. and E. Pebesma (2011): The pair-copula construction for spatial data: a new approach to model spatial dependency. Procedia Environmental Sciences (Vol. 7, pp. 206 - 211), Elsevier.
}
\author{
Benedikt Graeler
}

\seealso{
\code{\link{spatial_copula}}, \code{\linkS4class{distance_vine_copula}}
}
\examples{
# data from demo("spatial_copula", package = "sfcopula")
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

dCopula(u=matrix(c(.3,.3,.7,.7),ncol=2),spCop,h=c(200,400))
pCopula(u=matrix(c(.3,.3,.7,.7),ncol=2),spCop,h=c(200,400))

# geometric means
spGeomCop <- spatial_copula(list(claytonCopula(4), gumbelCopula(2), indepCopula()),
                            c(0, 200, 500), combination = "geometric")
spGeomCop@combination
dduCopula(c(0.3, 0.7), spGeomCop, h = 300)
}
\keyword{classes}
\keyword{sfcopula}
\keyword{copula}

\name{covariate_vine_copula-class}
\Rdversion{1.1}
\docType{class}
\alias{covariate_vine_copula-class}
\alias{distance_copula-class}
\alias{show,covariate_vine_copula-method}
\alias{dCopula,matrix,covariate_vine_copula-method}
\alias{dCopula,numeric,covariate_vine_copula-method}
\alias{dCopula,data.frame,covariate_vine_copula-method}
\alias{fitCopula,covariate_vine_copula,ANY-method}

\title{Class \code{"covariate_vine_copula"}}
\description{
A S4-class representing a covariate vine copula with a spatial or spatio-temporal tree. The class union \code{"distance_copula"} combines \code{\linkS4class{spatial_copula}} and \code{\linkS4class{spacetime_copula}}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("covariate_vine_copula", ...)} or by calls to the constructor \code{\link{covariate_vine_copula}}.
}
\section{Slots}{
  \describe{
    \item{\code{coVarCop}:}{Object of class \code{"function"} returning a bivariate copula object for the index of the central location.}
    \item{\code{tree}:}{Object of class \code{"distance_copula"} (a \code{\linkS4class{spatial_copula}} or \code{\linkS4class{spacetime_copula}}) modelling the dependence on the first tree.}
    \item{\code{topCop}:}{Object of class \code{"copula"} joining the conditioned covariate and neighbours to a full vine copula.}
    \item{\code{dimension}:}{Object of class \code{"integer"} denoting the dimension of the overall copula. }
    \item{\code{parameters}:}{Object of class \code{"numeric"}, not used.}
    \item{\code{param.names}:}{Object of class \code{"character"}, not used.}
    \item{\code{param.lowbnd}:}{Object of class \code{"numeric"}, not used.}
    \item{\code{param.upbnd}:}{Object of class \code{"numeric"}, not used.}
    \item{\code{fullname}:}{Object of class \code{"character"} providing a textual short description of this class.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{copula}"}, directly.
Class \code{"\linkS4class{Copula}"}, by class "copula", distance 2.
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "covariate_vine_copula")}: a short description of the vine.}
    \item{dCopula}{\code{signature(u = "matrix", copula = "covariate_vine_copula")}: the density of the vine, also for \code{u} being \code{"numeric"} or \code{"data.frame"}; needs the additional arguments \code{h} (distances) and \code{centre} (index of the central location), see \code{\link{covariate_vine_copula}}.}
    \item{fitCopula}{\code{signature(copula = "covariate_vine_copula")}: fits the top copula for a \code{\linkS4class{neighbourhood}} with covariate, see \code{\link{covariate_vine_copula}}.}
    \item{cond_vine}{\code{signature(vine = "covariate_vine_copula")}: see \code{\link{cond_vine}}.}
    \item{predict}{\code{signature(object = "covariate_vine_copula")}: see \code{\link[=predict,covariate_vine_copula-method]{predict}}.}
  }
}
\author{
Benedikt Graeler
}

\seealso{\code{\linkS4class{distance_vine_copula}}}
\examples{
showClass("covariate_vine_copula")
}
\keyword{classes}

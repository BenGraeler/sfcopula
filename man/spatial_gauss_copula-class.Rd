\name{spatial_gauss_copula-class}
\Rdversion{1.1}
\docType{class}
\alias{spatial_gauss_copula-class}
\alias{show,spatial_gauss_copula-method}

\title{Class \code{"spatial_gauss_copula"}}
\description{
A class representing a spatial Gaussian copula through its correlation function.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("spatial_gauss_copula", ...)} or by the constructor \code{\link{spatial_gauss_copula}}.
}
\section{Slots}{
  \describe{
    \item{\code{corFun}:}{Object of class \code{"function"}: the correlation function turning a distance matrix into a correlation matrix.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "spatial_gauss_copula")}: prints the correlation function.}
    \item{predict}{\code{signature(object = "spatial_gauss_copula")}: see \code{\link[=predict,spatial_gauss_copula-method]{predict}}.}
	 }
}
\author{
Benedikt Graeler
}
\seealso{
\code{\link{spatial_gauss_loglik}}
}
\examples{
showClass("spatial_gauss_copula")
}
\keyword{classes}

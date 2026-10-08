\name{invdduCopula-methods}
\docType{methods}
\alias{invdduCopula}
\alias{invdduCopula-methods}
\alias{invdduCopula,numeric,claytonCopula,numeric-method}
\alias{invdduCopula,numeric,frankCopula,numeric-method}
\alias{invdduCopula,numeric,indepCopula,numeric-method}
\alias{invdduCopula,numeric,normalCopula,numeric-method}
\alias{invdduCopula,numeric,spatial_copula,ANY-method}
\alias{invdduCopula,numeric,spacetime_copula,ANY-method}

\title{Methods for Function \code{invdduCopula} in Package \pkg{sfcopula}}
\description{
Methods for function \code{invdduCopula} in package \pkg{sfcopula}
}
\section{Methods}{
\describe{
\item{\code{signature(u = "ANY", copula = "ANY", y = "ANY")}}{
The inverse of the partial derivative of a copula is evaluated. For a given \code{u} and \code{y} a \code{v} is returned such that C(u,v)=y. In case no closed form is known, the evaluation is done numerically using \code{\link{optimise}}.}
\item{\code{signature(u = "numeric", copula = "asCopula", y = "numeric")}}{
Closed form evaluation of \code{invdduCopula} for the \code{\linkS4class{asCopula}}.}
\item{\code{signature(u = "numeric", copula = "claytonCopula", y = "numeric")}}{
Closed form evaluation of \code{invdduCopula} for the \code{\linkS4class{claytonCopula}}.}
\item{\code{signature(u = "numeric", copula = "cqsCopula", y = "numeric")}}{
Closed form evaluation of \code{invdduCopula} for the \code{\linkS4class{cqsCopula}}.}
\item{\code{signature(u = "numeric", copula = "frankCopula", y = "numeric")}}{
Closed form evaluation of \code{invdduCopula} for the \code{\linkS4class{frankCopula}}.}
\item{\code{signature(u = "numeric", copula = "indepCopula", y = "numeric")}}{
Closed form evaluation of \code{invdduCopula} for the \code{\linkS4class{indepCopula}}.}
\item{\code{signature(u = "numeric", copula = "normalCopula", y = "numeric")}}{
Closed form evaluation of \code{invdduCopula} for the \code{\linkS4class{normalCopula}}.}
\item{\code{signature(u = "numeric", copula = "spatial_copula", y = "ANY")}}{
Numerical evaluation of \code{invdduCopula} for the \code{\linkS4class{spatial_copula}}.}
\item{\code{signature(u = "numeric", copula = "spacetime_copula", y = "ANY")}}{
Numerical evaluation of \code{invdduCopula} for the \code{\linkS4class{spacetime_copula}}.}
}}
\keyword{methods}
\keyword{inverse partial derivatives}
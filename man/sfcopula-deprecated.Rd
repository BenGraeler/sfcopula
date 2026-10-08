\name{sfcopula-deprecated}
\alias{sfcopula-deprecated}
\alias{stCube}
\alias{calcBins}
\alias{composeSpCopula}
\alias{condCovariate}
\alias{condStCoVarVine}
\alias{condStVine}
\alias{condSpVine}
\alias{dropStTree}
\alias{dropSpTree}
\alias{fitCorFun}
\alias{fitSpCopula}
\alias{loglikByCopulasLags}
\alias{loglikByCopulasStLags}
\alias{reduceNeighbours}
\alias{stCopula}
\alias{stCoVarVineCopula}
\alias{stNeighbourhood}
\alias{getStNeighbours}
\alias{stCopPredict}
\alias{stVineCopula}
\alias{spCopula}
\alias{spGaussLogLik}
\alias{spGaussCopPredict}
\alias{spGeomCopula}
\alias{getNeighbours}
\alias{spCopPredict}
\alias{calcSpTreeDists}
\alias{spVineCopula}
\title{Deprecated spcopula-style function names}
\description{
These function names stem from the \pkg{spcopula} package and are deprecated in favour of snake_case names and of functions that unify the spatial and spatio-temporal case. They issue a deprecation warning and call the new function: most of them forward all arguments unchanged, the others translate their \pkg{spcopula} argument list into a call of the new function.
}
\usage{
stCube(...)
calcBins(...)
composeSpCopula(...)
condCovariate(...)
condStCoVarVine(condVar, dists, stCVVC, stInd, n = 1000)
condStVine(...)
condSpVine(...)
dropStTree(...)
dropSpTree(...)
fitCorFun(...)
fitSpCopula(...)
loglikByCopulasLags(...)
loglikByCopulasStLags(...)
reduceNeighbours(...)
stCopula(...)
stCoVarVineCopula(coVarCop, stCop, topCop)
stNeighbourhood(...)
getStNeighbours(stData, ST, spSize = 4, tlags = -(0:2), var = names(stData)[1],
                coVar = character(), timeSteps = NA, prediction = FALSE,
                min.dist = 0.01, timeCol = "time")
stCopPredict(predNeigh, dataST, predST, stVine, margin,
             method = "quantile", p = 0.5, ...)
stVineCopula(stCop, topCop)
spCopula(...)
spGaussLogLik(...)
spGaussCopPredict(corFun, predNeigh, dataLocs, predLocs, margin, p = 0.5, ...,
                  n = 1000)
spGeomCopula(components, distances, unit = "m")
getNeighbours(dataLocs, predLocs, size = 5, var = NULL, coVar = character(),
              prediction = FALSE, min.dist = 0.01)
spCopPredict(predNeigh, dataLocs, predLocs, spVine, margin,
             method = "quantile", p = 0.5, ...)
calcSpTreeDists(...)
spVineCopula(spCop, topCop = NULL)
}
\arguments{
  \item{...}{arguments passed on to the new function.}
  \item{components, distances, unit}{see \code{\link{spatial_copula}}; \code{spGeomCopula} uses \code{combination = "geometric"}.}
  \item{spCop, stCop}{the bivariate spatial or spatio-temporal copula of the first tree, passed as \code{trees} to \code{\link{distance_vine_copula}}.}
  \item{topCop}{see \code{\link{distance_vine_copula}}.}
  \item{dataLocs, stData}{the data, passed as \code{data} to \code{\link{neighbours}} or \code{\link[=predict,distance_vine_copula-method]{predict}}.}
  \item{predLocs, ST}{the prediction targets, passed as \code{target} to \code{\link{neighbours}} or \code{\link[=predict,distance_vine_copula-method]{predict}}; might be missing for fitting.}
  \item{size, spSize}{the (spatial) size of the neighbourhood, passed as \code{size} to \code{\link{neighbours}}.}
  \item{var, coVar, prediction, min.dist, tlags, timeSteps, timeCol}{see \code{\link{neighbours}}.}
  \item{condVar, dists, n}{see \code{\link{cond_vine}}.}
  \item{coVarCop}{the covariate copula function, see \code{\link{covariate_vine_copula}}.}
  \item{stCVVC}{the \code{\linkS4class{covariate_vine_copula}}, see \code{\link{cond_vine}}.}
  \item{stInd}{the spatio-temporal index pair passed on to the covariate copula function, see \code{\link{cond_vine}}.}
  \item{predNeigh}{the prediction neighbourhood, passed as \code{neigh} to \code{\link[=predict,distance_vine_copula-method]{predict}}.}
  \item{dataST, predST}{the data and the targets in the spatio-temporal case, passed as \code{data} and \code{target} to \code{\link[=predict,distance_vine_copula-method]{predict}}.}
  \item{spVine, stVine}{the vine copula, passed as \code{object} to \code{\link[=predict,distance_vine_copula-method]{predict}}.}
  \item{margin, method, p}{see \code{\link[=predict,distance_vine_copula-method]{predict}}.}
  \item{corFun}{the correlation function, turned into a \code{\link{spatial_gauss_copula}} for \code{\link[=predict,spatial_gauss_copula-method]{predict}}; \code{spGaussCopPredict} always uses the quantile predictor.}
}
\details{
\tabular{ll}{
\bold{deprecated} \tab \bold{use instead} \cr
\code{stCube} \tab \code{\link{as_spacetime_cube}} \cr
\code{calcBins} \tab \code{\link{calc_bins}} \cr
\code{composeSpCopula} \tab \code{\link{compose_spatial_copula}} \cr
\code{condCovariate} \tab \code{\link{cond_covariate}} \cr
\code{condStCoVarVine} \tab \code{\link{cond_vine}} \cr
\code{condStVine} \tab \code{\link{cond_vine}} \cr
\code{condSpVine} \tab \code{\link{cond_vine}} \cr
\code{dropStTree} \tab \code{\link{drop_tree}} \cr
\code{dropSpTree} \tab \code{\link{drop_tree}} \cr
\code{fitCorFun} \tab \code{\link{fit_cor_fun}} \cr
\code{fitSpCopula} \tab \code{\link{fit_spatial_copula}} \cr
\code{loglikByCopulasLags} \tab \code{\link{loglik_by_lags}} \cr
\code{loglikByCopulasStLags} \tab \code{\link{loglik_by_lags}} \cr
\code{reduceNeighbours} \tab \code{\link{reduce_neighbours}} \cr
\code{stCopula} \tab \code{\link{spacetime_copula}} \cr
\code{stCoVarVineCopula} \tab \code{\link{covariate_vine_copula}} \cr
\code{stNeighbourhood} \tab \code{\link{neighbourhood}} \cr
\code{getStNeighbours} \tab \code{\link{neighbours}} \cr
\code{stCopPredict} \tab \code{\link[=predict,distance_vine_copula-method]{predict}} \cr
\code{stVineCopula} \tab \code{\link{distance_vine_copula}} \cr
\code{spCopula} \tab \code{\link{spatial_copula}} \cr
\code{spGaussLogLik} \tab \code{\link{spatial_gauss_loglik}} \cr
\code{spGaussCopPredict} \tab \code{\link[=predict,spatial_gauss_copula-method]{predict}} with \code{\link{spatial_gauss_copula}} \cr
\code{spGeomCopula} \tab \code{\link{spatial_copula}(..., combination = "geometric")} \cr
\code{getNeighbours} \tab \code{\link{neighbours}} \cr
\code{spCopPredict} \tab \code{\link[=predict,distance_vine_copula-method]{predict}} \cr
\code{calcSpTreeDists} \tab \code{\link{tree_dists}} \cr
\code{spVineCopula} \tab \code{\link{distance_vine_copula}} \cr
}
S4 classes have been renamed and merged accordingly, e.g. \code{spCopula} and \code{spGeomCopula} to \code{\linkS4class{spatial_copula}}, \code{mixedSpVineCopula}, \code{pureSpVineCopula} and \code{stVineCopula} to \code{\linkS4class{distance_vine_copula}}, and \code{neighbourhood} and \code{stNeighbourhood} to \code{\linkS4class{neighbourhood}}; objects of the old classes can be translated with \code{\link{upgrade_spcopula}}.
}
\keyword{internal}

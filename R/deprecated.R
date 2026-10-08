############################################################
## deprecated spcopula-style names, kept to ease migration ##
############################################################

# forwards all arguments to the new function
.deprecatedAlias <- function(old, new) {
  force(old); force(new)
  function(...) {
    .Deprecated(new, package = "sfcopula", old = old)
    get(new, envir = asNamespace("sfcopula"))(...)
  }
}

.deprecate <- function(old, new)
  .Deprecated(new, package = "sfcopula", old = old)

## copulas

spCopula <- .deprecatedAlias("spCopula", "spatial_copula")
stCopula <- .deprecatedAlias("stCopula", "spacetime_copula")
stCoVarVineCopula <- function(coVarCop, stCop, topCop) {
  .deprecate("stCoVarVineCopula", "covariate_vine_copula")
  covariate_vine_copula(coVarCop, stCop, topCop)
}

spGeomCopula <- function(components, distances, unit = "m") {
  .deprecate("spGeomCopula", "spatial_copula(..., combination = \"geometric\")")
  spatial_copula(components, distances, unit = unit, combination = "geometric")
}

spVineCopula <- function(spCop, topCop = NULL) {
  .deprecate("spVineCopula", "distance_vine_copula")
  distance_vine_copula(spCop, topCop)
}

stVineCopula <- function(stCop, topCop) {
  .deprecate("stVineCopula", "distance_vine_copula")
  distance_vine_copula(stCop, topCop)
}

## neighbourhoods

stNeighbourhood <- .deprecatedAlias("stNeighbourhood", "neighbourhood")

getNeighbours <- function(dataLocs, predLocs, size = 5, var = NULL, coVar = character(),
                          prediction = FALSE, min.dist = 0.01) {
  .deprecate("getNeighbours", "neighbours")
  neighbours(dataLocs, if (missing(predLocs)) NULL else predLocs, size = size, var = var,
             coVar = coVar, prediction = prediction, min.dist = min.dist)
}

getStNeighbours <- function(stData, ST, spSize = 4, tlags = -(0:2), var = names(stData)[1],
                            coVar = character(), timeSteps = NA, prediction = FALSE,
                            min.dist = 0.01, timeCol = "time") {
  .deprecate("getStNeighbours", "neighbours")
  neighbours(stData, if (missing(ST)) NULL else ST, size = spSize, var = var, coVar = coVar,
             prediction = prediction, min.dist = min.dist, tlags = tlags,
             timeSteps = timeSteps, timeCol = timeCol)
}

reduceNeighbours <- .deprecatedAlias("reduceNeighbours", "reduce_neighbours")
stCube <- .deprecatedAlias("stCube", "as_spacetime_cube")

## binning and fitting

calcBins <- .deprecatedAlias("calcBins", "calc_bins")
fitCorFun <- .deprecatedAlias("fitCorFun", "fit_cor_fun")
loglikByCopulasLags <- .deprecatedAlias("loglikByCopulasLags", "loglik_by_lags")
loglikByCopulasStLags <- .deprecatedAlias("loglikByCopulasStLags", "loglik_by_lags")
fitSpCopula <- .deprecatedAlias("fitSpCopula", "fit_spatial_copula")
composeSpCopula <- .deprecatedAlias("composeSpCopula", "compose_spatial_copula")

## trees and conditioning

dropSpTree <- .deprecatedAlias("dropSpTree", "drop_tree")
dropStTree <- .deprecatedAlias("dropStTree", "drop_tree")
calcSpTreeDists <- .deprecatedAlias("calcSpTreeDists", "tree_dists")
condSpVine <- .deprecatedAlias("condSpVine", "cond_vine")
condStVine <- .deprecatedAlias("condStVine", "cond_vine")
condCovariate <- .deprecatedAlias("condCovariate", "cond_covariate")

condStCoVarVine <- function(condVar, dists, stCVVC, stInd, n = 1000) {
  .deprecate("condStCoVarVine", "cond_vine")
  # spcopula expected the covariate first, cond_vine expects it last
  cond_vine(c(condVar[-1], condVar[1]), dists, stCVVC, n = n, centre = stInd)
}

## prediction

spCopPredict <- function(predNeigh, dataLocs, predLocs, spVine, margin,
                         method = "quantile", p = 0.5, ...) {
  .deprecate("spCopPredict", "predict")
  predict(spVine, predNeigh, dataLocs, predLocs, margin, method = method, p = p, ...)
}

stCopPredict <- function(predNeigh, dataST, predST, stVine, margin,
                         method = "quantile", p = 0.5, ...) {
  .deprecate("stCopPredict", "predict")
  predict(stVine, predNeigh, dataST, predST, margin, method = method, p = p, ...)
}

spGaussCopPredict <- function(corFun, predNeigh, dataLocs, predLocs, margin, p = 0.5, ..., n = 1000) {
  .deprecate("spGaussCopPredict", "predict")
  predict(spatial_gauss_copula(corFun), predNeigh, dataLocs, predLocs, margin,
          method = "quantile", p = p, ..., n = n)
}

spGaussLogLik <- .deprecatedAlias("spGaussLogLik", "spatial_gauss_loglik")

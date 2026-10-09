############################################################
## deprecated spcopula-style names, kept to ease migration ##
############################################################

# forwards all arguments to the new function
.deprecatedAlias <- function(old, new) {
  force(old); force(new)
  function(...) {
    .Deprecated(new, package = "sfcopula", old = old)
    .quietArgs(get(new, envir = asNamespace("sfcopula"))(...))
  }
}

.deprecate <- function(old, new)
  .Deprecated(new, package = "sfcopula", old = old)

## deprecated argument names (spcopula style), accepted with a warning

.deprecateArg <- function(old, new)
  .Deprecated(msg = sprintf("Argument '%s' is deprecated; use '%s' instead.", old, new))

# value of a renamed argument passed through '...' (for S4 methods)
.dotsArg <- function(dots, old, new, value) {
  if (!old %in% names(dots))
    return(value)
  .deprecateArg(old, new)
  dots[[old]]
}

# wraps f such that the old argument names in map (c(old = "new")) are accepted
.withOldArgs <- function(f, map) {
  force(f); force(map)
  wrapper <- function() {
    call <- match.call(expand.dots = TRUE)
    argNames <- names(call)
    isOld <- !is.na(argNames) & argNames %in% names(map)
    for (old in argNames[isOld])
      .deprecateArg(old, map[[old]])
    names(call)[isOld] <- map[argNames[isOld]]
    call[[1]] <- f
    eval(call, parent.frame())
  }
  args <- formals(f)
  if (!"..." %in% names(args))
    args <- c(args, alist(... = ))
  formals(wrapper) <- args
  wrapper
}

# keeps the warnings about deprecated argument names out of the deprecated functions
.quietArgs <- function(expr)
  withCallingHandlers(expr, deprecatedWarning = function(w) {
    if (startsWith(conditionMessage(w), "Argument '"))
      invokeRestart("muffleWarning")
  })

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
             covar = coVar, prediction = prediction, min_dist = min.dist)
}

getStNeighbours <- function(stData, ST, spSize = 4, tlags = -(0:2), var = names(stData)[1],
                            coVar = character(), timeSteps = NA, prediction = FALSE,
                            min.dist = 0.01, timeCol = "time") {
  .deprecate("getStNeighbours", "neighbours")
  neighbours(stData, if (missing(ST)) NULL else ST, size = spSize, var = var, covar = coVar,
             prediction = prediction, min_dist = min.dist, tlags = tlags,
             time_steps = timeSteps, time_col = timeCol)
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

## return periods, Kendall distribution and tail dependence

kendallRP <- .deprecatedAlias("kendallRP", "kendall_rp")
criticalLevel <- .deprecatedAlias("criticalLevel", "critical_level")
criticalPair <- .deprecatedAlias("criticalPair", "critical_pair")
criticalTriple <- .deprecatedAlias("criticalTriple", "critical_triple")
kendallDistribution <- .deprecatedAlias("kendallDistribution", "kendall_distribution")
getKendallDistr <- .deprecatedAlias("getKendallDistr", "get_kendall_distr")
genEmpKenFun <- .deprecatedAlias("genEmpKenFun", "gen_emp_ken_fun")
genInvKenFun <- .deprecatedAlias("genInvKenFun", "gen_inv_ken_fun")
bivJointDepFun <- .deprecatedAlias("bivJointDepFun", "biv_joint_dep_fun")
lowerBivJointDepFun <- .deprecatedAlias("lowerBivJointDepFun", "lower_biv_joint_dep_fun")
upperBivJointDepFun <- .deprecatedAlias("upperBivJointDepFun", "upper_biv_joint_dep_fun")
empBivJointDepFun <- .deprecatedAlias("empBivJointDepFun", "emp_biv_joint_dep_fun")
lowerEmpBivJointDepFun <- .deprecatedAlias("lowerEmpBivJointDepFun", "lower_emp_biv_joint_dep_fun")
upperEmpBivJointDepFun <- .deprecatedAlias("upperEmpBivJointDepFun", "upper_emp_biv_joint_dep_fun")

## utilities

rankTransform <- .deprecatedAlias("rankTransform", "rank_transform")
dependencePlot <- .deprecatedAlias("dependencePlot", "dependence_plot")
unitScatter <- .deprecatedAlias("unitScatter", "unit_scatter")
univScatter <- function(formula = NULL, smpl) {
  .deprecate("univScatter", "unit_scatter")
  unit_scatter(formula, smpl)
}

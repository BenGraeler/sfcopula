###################################################################
## prediction from conditional densities of the central location ##
###################################################################

setGeneric("predict")

# conditioning values of the i-th neighbourhood (the central location is
# either missing (spatial prediction) or not part of the data (spatio-temporal));
# covariates (last columns) are only kept on request
.condValues <- function(neigh, i, covar = FALSE) {
  x <- as.numeric(neigh@data[i,])
  if (!covar && length(neigh@coVar) > 0)
    x <- x[seq_len(length(x) - length(neigh@coVar))]
  if (.isStNeighbourhood(neigh) && neigh@prediction)
    return(x)
  x[-1]
}

# quantile of a conditional density function (see .condDensityFun)
.condQuantile <- function(condFun, p, margin) {
  xVals <- attr(condFun,"xVals")
  density <- condFun(xVals)
  nx <- length(xVals)
  int <- cumsum(c(0,diff(xVals)*(0.5*diff(density)+density[-nx])))
  pVal <- if (is.numeric(p)) p else runif(1)
  lower <- max(which(int <= pVal))
  m <- (density[lower+1]-density[lower])/(xVals[lower+1]-xVals[lower])
  b <- density[lower]
  xRes <- -b/m+sign(m)*sqrt(b^2/m^2+2*(pVal-int[lower])/m)
  margin$q(xVals[lower]+xRes)
}

# expected value of a conditional density function
.condExpectation <- function(condFun, margin, loc, ..., stop.on.error=FALSE) {
  condExp <-  function(x) margin$q(x)*condFun(x)
  ePred <- integrate(condExp, 0+.Machine$double.eps, 1-.Machine$double.eps,
                     subdivisions=10000L, stop.on.error=stop.on.error, ...)
  if(ePred$abs.error > 0.05)
    warning("Numerical integration in predExpectation performed at a level of absolute error of only ",
            ePred$abs.error, " for location ", loc, ".")
  ePred$value
}

# evaluates condFunAt(i) for all neighbourhoods and adds the predictions to the target
.predictAll <- function(neigh, target, condFunAt, margin, method, p, ...) {
  stopifnot(is(neigh, "neighbourhood"))
  stopifnot(is.list(margin), is.function(margin$q))
  if (method == "quantile") {
    if (is.numeric(p))
      stopifnot(0 < p & p < 1)
    else
      stopifnot(p == "random")
  }

  nLocs <- nrow(neigh@data)
  pred <- numeric(nLocs)
  pb <- txtProgressBar(0, nLocs, 0, width=getOption("width")-10, style=3)
  for(i in 1:nLocs) {
    setTxtProgressBar(pb, i)
    condFun <- condFunAt(i)
    pred[i] <- switch(method,
                      quantile = .condQuantile(condFun, p, margin),
                      expectation = .condExpectation(condFun, margin, i, ...))
  }
  close(pb)

  name <- switch(method, quantile = paste0("quantile.", p), expectation = "expect")
  if (.isStNeighbourhood(neigh))
    return(.addStPrediction(target, pred, name))
  .addPrediction(target, pred, name)
}

## distance vine copulas
#########################

predictDistanceVine <- function(object, neigh, data, target, margin,
                                method = c("quantile", "expectation"), p = 0.5, ...) {
  method <- match.arg(method)
  dists <- tree_dists(neigh, data, length(object@trees))
  rowDists <- function(i) lapply(dists, function(x) {
    if (length(dim(x)) == 3) x[i,,,drop=FALSE] else x[i,]
  })
  .predictAll(neigh, target,
              function(i) cond_vine(.condValues(neigh, i), rowDists(i), object),
              margin, method, p, ...)
}

setMethod("predict", signature("distance_vine_copula"), predictDistanceVine)

## spatial Gaussian copula
###########################

spatial_gauss_copula <- function(cor_fun) {
  stopifnot(is.function(cor_fun))
  new("spatial_gauss_copula", corFun = cor_fun)
}

setMethod("show", signature("spatial_gauss_copula"), function(object) {
  cat("Spatial Gaussian copula with correlation function:\n")
  print(object@corFun)
})

predictSpatialGauss <- function(object, neigh, data, target, margin,
                                method = c("quantile", "expectation"), p = 0.5, ..., n = 1000) {
  method <- match.arg(method)
  stopifnot(!.isStNeighbourhood(neigh))
  neighDim <- ncol(neigh@data) - length(neigh@coVar)
  allDataDists <- .spDistMat(.pointGeom(data))
  xVals <- .condGrid(n)

  condFunAt <- function(i) {
    tmpDataDists <- allDataDists[neigh@index[i,-1], neigh@index[i,-1]]
    tmpDists <- rbind(c(0,neigh@distances[i,]),
                      cbind(neigh@distances[i,], tmpDataDists))
    tmpCor <- object@corFun(tmpDists)
    tmpGaussCop <- normalCopula(tmpCor[lower.tri(tmpCor)], neighDim, dispstr="un")
    density <- dCopula(cbind(xVals, matrix(rep(.condValues(neigh, i), length(xVals)),
                                           ncol=neighDim-1, byrow=T)),
                       tmpGaussCop)
    .condDensityFun(xVals, density)
  }

  .predictAll(neigh, target, condFunAt, margin, method, p, ...)
}

setMethod("predict", signature("spatial_gauss_copula"), predictSpatialGauss)

## spatio-temporal covariate vine copula
#########################################

# condVar: the neighbours followed by the covariate (as in the neighbourhood data)
# centre: index of the central location passed to the covariate copula function
condCovariateVine <- function (cond_var, dists, vine, n = 1000, ..., centre) {
  if (missing(cond_var))
    cond_var <- .dotsArg(list(...), "condVar", "cond_var", NULL)
  if (missing(centre))
    centre <- .dotsArg(list(...), "stInd", "centre", NULL)
  xVals <- .condGrid(n)
  repCondVar <- matrix(cond_var, ncol = length(cond_var), nrow = length(xVals), byrow = T)
  density <- dCovariateVine(cbind(xVals, repCondVar), vine, h = dists, centre = centre)
  .condDensityFun(xVals, density)
}

setMethod("cond_vine", signature(vine = "covariate_vine_copula"), condCovariateVine)

# the covariate is stored in the last column of the neighbourhood's data
predictCovariateVine <- function(object, neigh, data, target, margin,
                                 method = c("quantile", "expectation"), p = 0.5, ...) {
  method <- match.arg(method)
  stopifnot(length(neigh@coVar) == 1)
  isSt <- .isStNeighbourhood(neigh)
  centre <- .centreIndex(neigh)
  condFunAt <- function(i) {
    dists <- if (isSt) neigh@distances[i,,,drop=FALSE] else neigh@distances[i,]
    cond_vine(.condValues(neigh, i, covar = TRUE), dists, object,
              centre = if (isSt) centre[i,] else centre[i])
  }
  .predictAll(neigh, target, condFunAt, margin, method, p, ...)
}

setMethod("predict", signature("covariate_vine_copula"), predictCovariateVine)

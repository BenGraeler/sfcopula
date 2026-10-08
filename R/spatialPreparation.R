########################################################
##                                                    ##
## functions based on sf preparing the use of copulas ##
##                                                    ##
########################################################

## neighbourhood constructor (spatial and spatio-temporal)
#############################################################

# spatial: distances [n, k] and index [n, k+1] matrices
# spatio-temporal: distances [n, k, 2] and index [n, k+1, 2] arrays holding
#                  space (1) and time (2) in the third dimension
neighbourhood <- function(data, distances, index, var, coVar=character(), prediction=FALSE) {
  data <- as.data.frame(data)
  
  if (length(dim(distances)) == 3) {
    sizeN <- nrow(data)
    dimDists <- dim(distances)
    dimInd <- dim(index)
    
    stopifnot(length(dimInd) == 3)
    stopifnot(dimDists[1] == sizeN)
    stopifnot(dimInd[1] == dimDists[1])
    stopifnot(((dimDists[2] + !prediction) + length(coVar)) == ncol(data))
    stopifnot(dimInd[2] == dimDists[2]+1)
    stopifnot(dimDists[3] == 2)
    stopifnot(dimInd[3] == dimDists[3])
    
    colnames(data) <- paste(paste("N", (0+prediction):dimDists[2], sep=""), var, sep=".")
    if(length(coVar)>0)
      colnames(data)[ncol(data) + 1 - (length(coVar):1)] <- paste("N0", coVar)
  }
  
  if (anyDuplicated(rownames(data))>0)
    rownames(data) <- 1:nrow(data)
  
  new("neighbourhood", data=data, distances=distances, index=index,
      var=var, coVar=coVar, prediction=prediction)
}

.isStNeighbourhood <- function(neigh) length(dim(neigh@distances)) == 3

## show
showNeighbourhood <- function(object){
  cat("A set of", ifelse(.isStNeighbourhood(object), "spatio-temporal", "spatial"),
      "neighbourhoods consisting of", dim(object@distances)[2]+1, "locations each \n")
  if (length(object@var)>0) {
    cat("with", nrow(object@data), "rows of observations for:\n")
    cat(object@var, "\n")
  } else {
    cat("without data \n")
  }
  if(length(object@coVar)>0)
    cat("with covariate", object@coVar, "\n")
}

setMethod(show, signature("neighbourhood"), showNeighbourhood)

## names
setMethod(names, signature("neighbourhood"), function(x) c(x@var,x@coVar))

selectFromNeighbourhood <- function(x, i) {
  if (.isStNeighbourhood(x))
    return(new("neighbourhood", data=x@data[i,,drop=F], 
               distances=x@distances[i,,,drop=F], index=x@index[i,,,drop=F], 
               var=x@var, coVar=x@coVar, prediction=x@prediction))
  new("neighbourhood", data=x@data[i,,drop=F], 
      distances=x@distances[i,,drop=F], index=x@index[i,,drop=F], 
      var=x@var, coVar=x@coVar, prediction=x@prediction)
}

setMethod("[", signature("neighbourhood","numeric"), selectFromNeighbourhood) 

## neighbourhoods from sf points or stars vector data cubes
############################################################

neighbours <- function(data, target = NULL, size = NULL, var = NULL, coVar = character(),
                       prediction = !is.null(target), min.dist = 0.01,
                       tlags = -(0:2), timeSteps = NA, timeCol = "time") {
  if (is.null(var))
    var <- .attrNames(data)[1]
  if (inherits(data, "stars"))
    return(.stNeighbours(data, target, spSize = if (is.null(size)) 4 else size, 
                         tlags = tlags, var = var, coVar = coVar, timeSteps = timeSteps,
                         prediction = prediction, min.dist = min.dist, timeCol = timeCol))
  .spNeighbours(data, target, size = if (is.null(size)) 5 else size, var = var, 
                coVar = coVar, prediction = prediction, min.dist = min.dist)
}

## calculate neighbourhood from sf points
.spNeighbours <- function (dataLocs, predLocs = NULL, size = 5, 
                           var = .attrNames(dataLocs)[1], coVar=character(),
                           prediction = FALSE, min.dist = 0.01) {
  stopifnot((!prediction && is.null(predLocs)) || (prediction && !is.null(predLocs)))
  stopifnot(min.dist > 0 || prediction)
  
  if (is.null(predLocs) && !prediction) 
    predLocs = dataLocs
  
  dataGeom <- .pointGeom(dataLocs)
  predGeom <- .pointGeom(predLocs)
  
  hasData <- inherits(dataLocs, "sf") && length(var) > 0 && !all(is.na(var))
  if (hasData) {
    if (any(is.na(match(var, .attrNames(dataLocs))))) 
      stop("The variables is not part of the data.")
  }
  
  nLocs <- length(predGeom)
  size <- min(size, length(dataGeom) + prediction)
  
  knn <- .knn(dataGeom, predGeom, size - 1, min.dist)
  allLocs <- cbind(1:nLocs, knn$index)
  allDists <- knn$dists
  
  if (hasData) {
    varValues <- .attrValues(dataLocs, var)
    if (!prediction) {
      allData <- matrix(varValues[allLocs], nLocs, size)
    } else {
      allData <- cbind(rep(NA, nLocs), 
                       matrix(varValues[allLocs[, -1]], nLocs, size - 1))
    }
    colnames(allData) <- paste(paste("N", rep(0:(size - 1), each = length(var)), sep = ""),
                               rep(var, size), sep = ".")
    
    # covariates of the central location: from the data or, for prediction, from the target
    if (length(coVar) > 0) {
      coVarSource <- if (prediction) predLocs else dataLocs
      if (!inherits(coVarSource, "sf") || !all(coVar %in% .attrNames(coVarSource)))
        stop("The covariate(s) need to be attributes of the ", 
             ifelse(prediction, "target", "data"), ".")
      coVarData <- sapply(coVar, function(cv) .attrValues(coVarSource, cv))
      coVarData <- matrix(coVarData, nrow = nLocs)
      colnames(coVarData) <- paste("N0", coVar, sep = ".")
      allData <- cbind(allData, coVarData)
    }
  } else {
    allData <- as.data.frame(matrix(NA, nLocs, size + length(coVar)))
    var <- character()
  }
  
  dimnames(allLocs) <- NULL
  return(neighbourhood(data=allData, distances=allDists, 
                       index=allLocs, var=var, coVar=coVar,
                       prediction=prediction))
}

#############
## BINNING ##
#############

# calculates lag indicies for spatial points and stores the respective separating distances
# 
# boundaries  -> are the right-side limits of the distance classes
# data --------> an sf or sfc object with POINT geometries
calcSpLagInd <- function(data, boundaries) {
  lags <- vector("list",length(boundaries))
  
  geom <- .pointGeom(data)
  dists <- .spDistMat(geom)
  
  pairs <- which(upper.tri(dists), arr.ind = TRUE)
  pairs <- pairs[order(pairs[, 1], pairs[, 2]), , drop = FALSE]
  d <- dists[pairs]
  
  # first boundary that is larger than the distance
  k <- findInterval(d, boundaries) + 1
  for (b in unique(k[k <= length(boundaries)])) {
    bool <- k == b
    lags[[b]] <- cbind(pairs[bool, 1], pairs[bool, 2], d[bool])
    dimnames(lags[[b]]) <- NULL
  }
  return(lags)
}

# the generic calc_bins, calculates bins for spatial and spatio-temporal data
setGeneric("calc_bins", function(data, var, nbins=15, boundaries=NA, cutoff=NA,
                                ..., cor.method="fasttau", plot=TRUE) {
                         standardGeneric("calc_bins") 
                         })

## calculating the spatial bins
################################

calcSpBins <- function(data, var, nbins=15, boundaries=NA, cutoff=NA, 
                       cor.method="fasttau", plot=TRUE) {

  if(is.na(cutoff)) {
    cutoff <- .bboxDiag(.pointGeom(data))/3
  }
  if(any(is.na(boundaries))) {
    boundaries <- ((1:nbins) * cutoff/nbins)
  }
    
  nbins <- length(boundaries)-1
  
  lags <- calcSpLagInd(data, boundaries)
    
  mDists <- sapply(lags, function(x) mean(x[,3]))
  np <- sapply(lags, function(x) length(x[,3]))
  varValues <- .attrValues(data, var)
  lagData <- lapply(lags, function(x) {
    lagPairs <- cbind(varValues[x[,1]], varValues[x[,2]])
    colnames(lagPairs) <- c(var, var)
    lagPairs
  })
  
  if(cor.method == "fasttau")
    lagCor <- sapply(lagData, function(x) TauMatrix(x)[1,2])
  if(cor.method %in% c("kendall","spearman","pearson"))
    lagCor <- sapply(lagData, function(x) cor(x,method=cor.method)[1,2])
  if(cor.method == "normVariogram")  
    lagCor <- sapply(lagData, function(x) 1-cor(x,method="pearson")[1,2])
  if(cor.method == "variogram")  
    lagCor <- sapply(lagData, function(x) 0.5*mean((x[,1]-x[,2])^2,na.rm=T))
    
  if(plot) { 
    plot(mDists, lagCor, xlab="distance",ylab=paste("correlation [",cor.method,"]",sep=""), 
         ylim=1.05*c(-abs(min(lagCor)), max(lagCor)), xlim=c(0,max(mDists)))
    abline(h=c(-min(lagCor),0,min(lagCor)),col="grey")
  }
  
  res <- list(np=np, meanDists = mDists, lagCor=lagCor, lags=lags)
  attr(res,"cor.method") <- cor.method
  attr(res,"variable") <- var
  return(res)
}

setMethod(calc_bins, signature("sf"), calcSpBins)

# calc bins from a (conditional) neighbourhood

calcNeighBins <- function(data, var=data@var, nbins=9, boundaries=NA, 
                          cutoff=NA, cor.method="kendall", plot=TRUE) {
  if (.isStNeighbourhood(data))
    stop("calc_bins is only available for spatial neighbourhoods.")
  dists <- data@distances
  
  corFun <- switch(cor.method,
                   fasttau=function(x) TauMatrix(x)[1,2],
                   function(x) cor(x,method=cor.method)[1,2])
  
  if (any(is.na(boundaries))) 
    boundaries <- quantile(as.vector(dists), probs=c(1:nbins/nbins))
  if(!is.na(cutoff)) {
    boundaries <- boundaries[boundaries < cutoff]
    boundaries <- unique(c(0,boundaries,cutoff))
  } else {
    boundaries <- unique(c(0,boundaries))
  }
  
  nbins <- length(boundaries)-1
  
  np <- numeric(nbins)
  moa <- numeric(nbins)
  meanDists <- numeric(nbins)

  data <- as.matrix(data@data)
  
  lagData <- list()
  
  for (i in 1:nbins) {
    bools <- (dists <= boundaries[i+1] & dists > boundaries[i])
    
    pairs <- NULL
    for(col in 1:(dim(bools)[2])) {
      pairs <- rbind(pairs, data[bools[,col],c(1,1+col)])
    }
    
    lagData[[i]] <- pairs
    moa[i] <- corFun(pairs)
    meanDists[i] <- mean(dists[bools])
    np[i] <- sum(bools)
  }
  
  if(plot) { 
    plot(meanDists, moa, xlab="distance", ylab=paste("correlation [",cor.method,"]",sep=""), 
         ylim=1.05*c(-abs(min(moa, na.rm=T)),max(moa, na.rm=T)), xlim=c(0,max(meanDists,na.rm=T)))
    abline(h=c(-min(moa),0,min(moa)),col="grey")
  }
  
  res <- list(np=np, meanDists = meanDists, lagCor=moa, lagData=lagData)
  attr(res,"cor.method") <- switch(cor.method, fasttau="kendall", cor.method)
  attr(res,"variable") <- var
  
  return(res)
}
  
setMethod(calc_bins, signature="neighbourhood", calcNeighBins)
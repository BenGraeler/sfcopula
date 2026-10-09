###############################################################
##                                                           ##
## functions based on stars preparing the use of copulas    ##
##                                                           ##
###############################################################

## calculate neighbourhood from a stars vector data cube

# returns an neighbourhood object
##################################

.stNeighbours <- function(stData, ST=NULL, spSize=4, tlags=-(0:2),
                          var=names(stData)[1], covar=character(),
                          time_steps=NA, prediction=FALSE, min_dist=0.01,
                          time_col="time") {
  stopifnot((!prediction && is.null(ST)) || (prediction && !is.null(ST)))
  stopifnot(min_dist>0 || prediction)
  stopifnot(length(tlags) > 1)

  timeSpan <- min(tlags)

  dataGeom <- .stGeom(stData)
  dataTime <- .stTime(stData)
  nTime <- length(dataTime)

  if(any(is.na(match(var,names(stData)))))
    stop("The variables is not part of stData.")
  if(length(covar)>0)
    if(any(is.na(match(covar,names(stData)))))
      stop("The covariate is not part of stData.")

  varValues <- .stValues(stData, var)

  if(!prediction) {
    if(is.na(time_steps)) {
      time_steps <- nTime+timeSpan
      reSample <- function() (1-timeSpan):nTime
    } else {
      reSample <- function() sort(sample((1-timeSpan):nTime, time_steps))
    }
    nghbrs <- .spNeighbours(dataLocs=dataGeom, var=character(), size=spSize,
                            min_dist=min_dist)
    nLocs <- length(dataGeom)*time_steps
    # draw once in advance to keep random streams in line with spcopula
    invisible(reSample())
    # each location is the centre of timeSteps neighbourhoods
    rowLoc <- rep(1:nrow(nghbrs@index), each=time_steps)
    rowTime <- unlist(lapply(1:nrow(nghbrs@index), function(i) reSample()))
    coVarValues <- lapply(covar, function(cv) .stValues(stData, cv)[cbind(nghbrs@index[rowLoc, 1], rowTime)])
  } else {
    targets <- .stTargets(ST, time_col)
    nghbrs <- .spNeighbours(dataLocs=dataGeom, predLocs=targets$geom,
                            size=spSize, var=character(), prediction=prediction,
                            min_dist=min_dist)
    nLocs <- length(targets$spInd)
    rowLoc <- targets$spInd
    rowTime <- match(as.numeric(targets$time4row), as.numeric(dataTime))
    if (any(is.na(rowTime)))
      stop("All target time instances need to be part of the time dimension of stData.")
    coVarValues <- lapply(covar, function(cv) .stTargetValues(ST, cv, targets))
  }

  nStNeighs <- (spSize-1)*length(tlags)

  stNeighData <- matrix(NA, nLocs, nStNeighs + 1 + length(covar))
  stDists <- array(NA,c(nLocs, nStNeighs, 2))
  stInd <- array(NA,c(nLocs, nStNeighs + 1, 2))

  # look up values at (location, time) pairs; outside the time frame -> NA
  lookUp <- function(sInd, tInd) {
    bool <- tInd >= 1 & tInd <= nTime & !is.na(sInd)
    res <- rep(NA, length(sInd))
    res[bool] <- varValues[cbind(sInd[bool], tInd[bool])]
    res
  }

  nbIndex <- nghbrs@index[rowLoc, , drop=FALSE]

  # the central location (unknown in the prediction case)
  if (!prediction)
    stNeighData[, 1] <- lookUp(nbIndex[, 1], rowTime)

  # the spatial neighbours of the central location for all temporal lags
  for (j in 1:length(tlags)) {
    cols <- (j-1)*(spSize-1)+2:(spSize)
    for (k in 2:spSize)
      stNeighData[, cols[k-1]] <- lookUp(nbIndex[, k], rowTime+tlags[j])
  }

  # add covariate(s) to the last column(s)
  if (length(covar) > 0) {
    coVarCols <- nStNeighs + 1 + (1:length(covar))
    stNeighData[, coVarCols] <- do.call(cbind, coVarValues)
  }

  # store spatial distances
  stDists[,,1] <- nghbrs@distances[rowLoc, rep(1:(spSize-1), length(tlags)), drop=FALSE]

  # store temporal distances
  stDists[,,2] <- matrix(rep(tlags, each=spSize-1), nLocs, nStNeighs, byrow=TRUE)

  # store space indices (the central location refers to the target in the prediction case)
  stInd[,,1] <- nbIndex[, c(1, rep(2:spSize, length(tlags))), drop=FALSE]

  # store time indices
  stInd[,,2] <- rowTime + matrix(c(0, rep(tlags, each=spSize-1)),
                                 nLocs, nStNeighs + 1, byrow=TRUE)

  if (prediction)
    stNeighData <- stNeighData[,-1,drop=FALSE]

  return(neighbourhood(as.data.frame(stNeighData), stDists, stInd, var, covar, prediction))
}


## reduction of a larger neigbopurhood based on correlation strengths
reduce_neighbours <- function(neigh, dep_fun, n,
                             prediction=neigh@prediction, drop_empty=!prediction) {
  stopifnot(n>0)
  if (!.isStNeighbourhood(neigh))
    stop("reduce_neighbours needs a spatio-temporal neighbourhood.")

  # transform distances into correlations to detect the strongest correlated ones
  dimStNeigh <- dim(neigh@distances)
  corMat <- matrix(NA, dimStNeigh[1], dimStNeigh[2])

  pb <- txtProgressBar(0, 2*dimStNeigh[1], style=3)
  for (i in 1:dimStNeigh[2]) {
    # whether neighbours are missing: set distance to NA
    if (prediction) # central location is not part of the data
      boolNA <- is.na(neigh@data[[i]])
    else {
      if(drop_empty) # neighbourrhoods with missing central location are not to be considered
        boolNA <- is.na(neigh@data[[1]]) | is.na(neigh@data[[1+i]])
      else # do not care about NA at the central location (e.g. cross-validation)
        boolNA <- is.na(neigh@data[[1+i]])
    }
    neigh@distances[boolNA,i,] <- c(NA,NA)
    tLag <- -1*neigh@distances[!boolNA,i,2][1]+1
    corMat[!boolNA,i] <- dep_fun(neigh@distances[!boolNA,i,1], tLag)
    setTxtProgressBar(pb, i*dimStNeigh[1]/dimStNeigh[2])
  }

  highCorMat <- t(apply(corMat, 1, function(x) order(x, na.last=TRUE, decreasing=TRUE)[1:n]))
  nrCM <- nrow(highCorMat)

  stNeighDataRed <- matrix(NA, nrow=nrCM, ncol=n+1+length(neigh@coVar))
  stNeighDistRed <- array(NA, dim=c(nrCM, n, 2))
  stNeighIndeRed <- array(NA, dim=c(nrCM, n+1, 2))
  if (length(neigh@coVar) > 0) {
    for (i in 1:nrCM) {
      if (prediction)
        selCol <- c(highCorMat[i,], ncol(neigh@data)-((length(neigh@coVar)-1):0))
      else
        selCol <- c(1,highCorMat[i,]+1, ncol(neigh@data)-((length(neigh@coVar)-1):0))
      stNeighDataRed[i,] <- as.numeric(neigh@data[i,selCol])
      stNeighDistRed[i,,] <- neigh@distances[i,highCorMat[i,],]
      stNeighIndeRed[i,,] <- neigh@index[i,c(1,highCorMat[i,]+1),]
      setTxtProgressBar(pb, dimStNeigh[1]+i)
    }
  } else {
    for (i in 1:nrCM) {
      if (prediction)
        selCol <- c(highCorMat[i,])
      else
        selCol <- c(1,highCorMat[i,]+1)
      stNeighDataRed[i,] <- as.numeric(neigh@data[i, selCol])
      stNeighDistRed[i,,] <- neigh@distances[i,highCorMat[i,],]
      stNeighIndeRed[i,,] <- neigh@index[i,c(1,highCorMat[i,]+1),]
      setTxtProgressBar(pb, dimStNeigh[1]+i)
    }
  }
  close(pb)

  # check whether neighbourhoods with missing central locations need to be dropped
  if (drop_empty) {
    boolNA <- !is.na(neigh@data[[1]])
    stNeighDataRed <- stNeighDataRed[boolNA,]
    stNeighDistRed <- stNeighDistRed[boolNA,,]
    stNeighIndeRed <- stNeighIndeRed[boolNA,,]
  }

  return(neighbourhood(stNeighDataRed, stNeighDistRed, stNeighIndeRed,
                         var=neigh@var, covar=neigh@coVar,
                         prediction=neigh@prediction))
}

## to be redone
# calcStNeighBins <- function(data, var="uniPM10", nbins=9, tlags=-(0:2),
#                             boundaries=NA, cutoff=NA, cor.method="fasttau") {
#   dists <- data@distances[,,1]
#
#   corFun <- switch(cor.method,
#                    fasttau=function(x) VineCopula:::fasttau(x[,1],x[,2]),
#                    function(x) cor(x,method=cor.method)[1,2])
#
#   if (any(is.na(boundaries)))
#     boundaries <- quantile(as.vector(dists), probs=c(1:nbins/nbins))
#   if(!is.na(cutoff)) {
#     boundaries <- boundaries[boundaries < cutoff]
#     boundaries <- unique(c(0,boundaries,cutoff))
#   } else {
#     boundaries <- unique(c(0,boundaries))
#   }
#
#   lagData <- NULL
#   for(tlag in tlags) { # tlag <- 0
#     tBool <- data@distances[,,2]==tlag
#     tmpLagData <- NULL
#     for(i in 1:nbins) { # i <- 1
#       sBool <- (dists <= boundaries[i + 1] & dists > boundaries[i])
#       bool <- tBool & sBool
#       pairs <- NULL
#       for (col in 1:(dim(tBool)[2])) { # col <- 1
#         if(!any(bool[, col]))
#           next
#         sInd <- data@index[bool[, col], c(1, 1 + col),1]
#         tInd <- data@index[bool[, col], c(1, 1 + col),2]
#         p1 <- apply(cbind(sInd[,1], tInd[,1]),1,
#                     function(x) data@locations[x[1], x[2],var])
#         p2 <- apply(cbind(sInd[,2], tInd[,2]),1,
#                     function(x) data@locations[x[1], x[2],var])
#         pairs <- rbind(pairs, cbind(p1,p2))
#       }
#       tmpLagData <- append(tmpLagData,list(pairs))
#     }
#     lagData <- append(lagData,list(tmpLagData))
#
#   }
#
#   lagData <- lapply(spIndices, retrieveData, tempIndices = tempIndices)
#   calcStats <- function(binnedData) {
#     cors <- NULL
#     for (i in 1:(ncol(binnedData)/2)) {
#       cors <- c(cors, cor(binnedData[, 2 * i - 1], binnedData[, 2 * i], method = cor.method, use = "pairwise.complete.obs"))
#     }
#     return(cors)
#   }
#   calcTau <- function(binnedData) {
#     cors <- NULL
#     for (i in 1:(ncol(binnedData)/2)) {
#       cors <- c(cors, VineCopula:::fasttau(binnedData[, 2 * i - 1], binnedData[, 2 * i]))
#     }
#     return(cors)
#   }
#   calcCor <- switch(cor.method, fasttau = calcTau, calcStats)
#   lagCor <- sapply(lagData, calcCor)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#   np <- numeric(0)
#   moa <- numeric(0)
#   lagData <- NULL
#   meanDists <- numeric(0)
#
#   data <- as.matrix(data@data)
#
#   for ( i in 1:nbins) {
#     bools <- (dists <= boundaries[i+1] & dists > boundaries[i])
#
#     pairs <- NULL
#     for(col in 1:(dim(bools)[2])) {
#       pairs <- rbind(pairs, data[bools[,col],c(1,1+col)])
#     }
#
#     lagData <- append(lagData, list(pairs))
#     moa <- c(moa, corFun(pairs))
#     meanDists <- c(meanDists, mean(dists[bools]))
#     np <- c(np, sum(bools))
#   }
#
#   if(plot) {
#     plot(meanDists, moa, xlab="distance", ylab=paste("correlation [",cor.method,"]",sep=""),
#          ylim=1.05*c(-abs(min(moa)),max(moa)), xlim=c(0,max(meanDists)))
#     abline(h=c(-min(moa),0,min(moa)),col="grey")
#   }
#
#   res <- list(np=np, meanDists = meanDists, lagCor=moa, lagData=lagData)
#   attr(res,"cor.method") <- switch(cor.method, fasttau="kendall", cor.method)
#   return(res)
# }
#
# setMethod(calc_bins, signature="spacetime_neighbourhood", calcStNeighBins)


# instances: number  -> number of randomly choosen temporal intances
#            NA      -> all observations
#            other   -> time indices or time stamps (Date/POSIXct) of the time dimension; tlags is set to 0 in this case.
# tlags:    numeric -> temporal shifts between obs
calcStBins <- function(data, var, nbins=15, boundaries=NA, cutoff=NA,
                       instances=NA, tlags=-(0:2), ...,
                       cor_method="fasttau", plot=FALSE) {
  cor_method <- .dotsArg(list(...), "cor.method", "cor_method", cor_method)
  geom <- .stGeom(data)
  time <- .stTime(data)
  varValues <- .stValues(data, var)

  if(is.na(cutoff))
    cutoff <- .bboxDiag(geom)/3
  if(any(is.na(boundaries)))
    boundaries <- (1:nbins) * cutoff / nbins
  if(length(instances) == 1 && is.na(instances))
    instances=length(time)

  spIndices <- calcSpLagInd(geom, boundaries)

  mDists <- sapply(spIndices, function(x) mean(x[,3]))

  lengthTime <- length(time)
  if (!is.numeric(instances) | !length(instances)==1) {
    # explicit time instances (indices or time stamps) without temporal lags
    if (!is.numeric(instances))
      instances <- match(as.numeric(instances), as.numeric(time))
    stopifnot(all(!is.na(instances)))
    tempIndices <- cbind(instances, instances)
    tlags <- 0
  }
  else {
    tempIndices <- NULL
    for (tlag in rev(tlags)) {
      if(is.na(instances))
        smplInd <- max(1,1-min(tlags)):min(lengthTime,lengthTime-min(tlags))
      else
        smplInd <- sort(sample(x=max(1,1-min(tlags)):min(lengthTime,lengthTime-min(tlags)),
                               size=min(instances,lengthTime-max(abs(tlags)))))

      tempIndices <- cbind(smplInd+tlag, tempIndices)
      tempIndices <- cbind(smplInd, tempIndices)
    }
  }

  # internal stat function
  calcStats <- function(binnedData) {
    return(c(sum(complete.cases(binnedData)),
             cor(binnedData[,1], binnedData[,2],
                 method=cor_method,
                 use="pairwise.complete.obs")))
  }

  # internal fast tau function
  calcTau <- function(tmpData) {
      tmpData <- tmpData[complete.cases(tmpData),]
    return(c(nrow(tmpData), TauMatrix(tmpData)[1,2]))
  }

  calc_cor <- switch(cor_method, fasttau=calcTau, calcStats)

  retrieveData <- function(spIndex, tempIndices, cor_fun) {
    binStats <- matrix(NA, nrow = ncol(tempIndices)/2, ncol = 2)
    for (i in 1:(ncol(tempIndices)/2)) {
      binStats[i,] <- cor_fun(cbind(as.vector(varValues[spIndex[,1], tempIndices[,2*i-1], drop=F]),
                                   as.vector(varValues[spIndex[,2], tempIndices[,2*i], drop=F])))
    }

  return(binStats)
  }

  lagStats <- lapply(spIndices, retrieveData, tempIndices=tempIndices, cor_fun=calc_cor)

  lagCor <- matrix(NA, length(tlags), nbins)
  lagNp  <- matrix(NA, length(tlags), nbins)
  for (i in 1:length(lagStats)) {
    lagNp[,i]  <- lagStats[[i]][,1]
    lagCor[,i] <- lagStats[[i]][,2]
  }

  if(plot) {
    plot(mDists, lagCor[1,],
         xlab="distance",
         ylab=paste("correlation [",cor_method,"]",sep=""),
         ylim=1.05*c(-abs(min(lagCor)), max(lagCor)),
         xlim=c(0,max(mDists)))
    abline(h=c(-min(lagCor), 0, min(lagCor)), col="grey")
  }

  res <- list(meanDists = mDists, lagCor = lagCor, lagNp=lagNp,
              lags=list(sp=spIndices, time=tempIndices))
  attr(res,"cor.method") <- cor_method
  attr(res, "variable") <- var
  return(res)
}

setMethod(calc_bins, signature(data="stars"), calcStBins)
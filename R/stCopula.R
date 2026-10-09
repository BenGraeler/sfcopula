######################################
## Spatio-Temporal Bivariate Copula ##
######################################

## constructor ##
#################

spacetime_copula <- function(components, tlags, distances=NA, dep_fun, unit="m", tres="day") {
  if(all(sapply(components, function(x) class(x)=="spatial_copula"))) {
    if(length(unique(sapply(components, function(x) x@unit))) >1 )
      stop("All spatial copulas need to have the same distance unit.")
    stopifnot(length(tlags) == length(components))
    spCopList <- components
  } else {
    spCopList <- list()
    
    if(!missing(dep_fun)) {
      getSpCop <- function(comp,dist,time) spatial_copula(comp, dist,
                                                    dep_fun=function(h) dep_fun(h, time, 1:length(tlags)), unit)
      for(i in 1:length(tlags)){
        spCopList <- append(spCopList, getSpCop(components[[i]], distances[[i]], i))
      }
    } else {
      for(i in 1:length(tlags)){
        spCopList <- append(spCopList, spatial_copula(components[[i]], distances[[i]], unit=unit))
      }
    }
  }
  
  param       <- unlist(lapply(spCopList, function(x) x@parameters))
  param.names <- unlist(lapply(spCopList, function(x) x@param.names))
  param.low   <- unlist(lapply(spCopList, function(x) x@param.lowbnd))
  param.up    <- unlist(lapply(spCopList, function(x) x@param.upbnd))
  
  new("spacetime_copula", dimension=as.integer(2), parameters=param, param.names=param.names,
      param.lowbnd=param.low, param.upbnd=param.up,
      fullname="Spatio-Temporal Copula: distance and time dependent convex combination of bivariate copulas",
      spCopList=spCopList, tlags=tlags, tres=tres)
}

## show method ##
#################

showStCopula <- function(object) {
  cat(object@fullname, "\n")
  cat("Dimension: ", object@dimension, "\n")
  cat("Copulas:\n")
  for (i in 1:length(object@spCopList)) {
    cmpCop <- object@spCopList[[i]]
    cat("  ", describeCop(cmpCop, "very short"), "at", object@tlags[i], 
        paste("[",object@tres,"]",sep=""), "\n")
    show(cmpCop)
  }
}

setMethod("show", signature("spacetime_copula"), showStCopula)

## spatial copula cdf ##
########################

pStCopula <- function (u, copula, h) {
  stopifnot(ncol(h)==2)
  stopifnot(nrow(h)==1 || nrow(h)==nrow(u))
  
  n <- nrow(u)
  tDist <- unique(h[,2])
  
  if(any(is.na(match(tDist,copula@tlags)))) 
    stop("Prediction time(s) do(es) not math the modelled time slices.")
  
  if (length(tDist)==1) {
    res <- pSpCopula(u, copula@spCopList[[match(tDist, copula@tlags)]], h[,1])
  } else {
    res <- numeric(n)
    for(t in tDist) {
      tmpInd <- h[,2]==t
      tmpCop <- copula@spCopList[[match(t, copula@tlags)]]
      res[tmpInd] <- pSpCopula(u[tmpInd,,drop=F], tmpCop, h[tmpInd,1])
    }
  }
  res
}

setMethod(pCopula, signature("numeric","spacetime_copula"), 
          function(u, copula, log, ...) pStCopula(matrix(u,ncol=2), copula, ...))
setMethod(pCopula, signature("matrix","spacetime_copula"), pStCopula)

## spatial Copula density ##
############################

dStCopula <- function (u, copula, log, h) {
  stopifnot(ncol(h)==2)
  stopifnot(nrow(h)==1 || nrow(h)==nrow(u))
  
  n <- nrow(u)
  tDist <- unique(h[,2])
  
  if(any(is.na(match(tDist,copula@tlags)))) 
    stop("Prediction time(s) do(es) not math the modelled time slices.")
  
  if (length(tDist)==1) {
    res <- dSpCopula(u, copula@spCopList[[match(tDist, copula@tlags)]], log, h[,1])
  } else {
    res <- numeric(n)
    for(t in tDist) {
      tmpInd <- h[,2]==t
      tmpCop <- copula@spCopList[[match(t, copula@tlags)]]
      res[tmpInd] <- dSpCopula(u[tmpInd,,drop=F], tmpCop, log, h[tmpInd,1])
    }
  }
  res
}

setMethod(dCopula, signature("numeric","spacetime_copula"), 
          function(u, copula, log, ...) dStCopula(matrix(u,ncol=2), copula, log=log, ...))
setMethod(dCopula, signature("matrix","spacetime_copula"), dStCopula)


## partial derivatives ##

## dduSpCopula ##
#################

dduStCopula <- function (u, copula, h) {
  stopifnot(ncol(h)==2)
  stopifnot(nrow(h)==1 || nrow(h)==nrow(u))
  
  n <- nrow(u)
  tDist <- unique(h[,2])
  
  if(any(is.na(match(tDist,copula@tlags)))) 
    stop("Prediction time(s) do(es) not match the modelled time slices.")
  
  if (length(tDist)==1) {
    res <- dduSpCopula(u, copula@spCopList[[match(tDist, copula@tlags)]], h[,1])
  } else {
    res <- numeric(n)
    for(t in tDist) {
      tmpInd <- h[,2]==t
      tmpCop <- copula@spCopList[[match(t, copula@tlags)]]
      res[tmpInd] <- dduSpCopula(u[tmpInd,,drop=F], tmpCop, h[tmpInd,1])
    }
  }
  res
}

setMethod("dduCopula", signature("numeric","spacetime_copula"), 
          function(u, copula, ...) dduStCopula(matrix(u,ncol=2), copula, ...))
setMethod("dduCopula", signature("matrix","spacetime_copula"), dduStCopula)

invdduStCopula <- function(u, copula, y, h, tol=1e-12) {
  if(!is.matrix(h)) 
    h <- matrix(h,ncol=2)
  nElem <- length(u)
  stopifnot(nElem == length(y))
  stopifnot(nrow(h) == 1 | nrow(h)==nElem)
  
  .bisect(function(v) dduStCopula(cbind(u, v), copula, h), y, tol=tol)
}

setMethod("invdduCopula", signature("numeric", "spacetime_copula"), invdduStCopula)

## ddvSpCopula ##
#################

ddvStCopula <- function (u, copula, h) {
  stopifnot(ncol(h)==2)
  stopifnot(nrow(h)==1 || nrow(h)==nrow(u))
  
  n <- nrow(u)
  tDist <- unique(h[,2])
  
  if(any(is.na(match(tDist,copula@tlags)))) 
    stop("Prediction time(s) do(es) not match the modelled time slices.")
  
  if (length(tDist)==1) {
    res <- ddvSpCopula(u, copula@spCopList[[match(tDist,copula@tlags)]], h[,1])
  } else {
    res <- numeric(n)
    for(t in tDist) {
      tmpInd <- h[,2]==t
      tmpCop <- copula@spCopList[[match(t, copula@tlags)]]
      res[tmpInd] <- ddvSpCopula(u[tmpInd,,drop=F], tmpCop, h[tmpInd,1])
    }
  }
  res
}

setMethod("ddvCopula", signature("numeric","spacetime_copula"), 
          function(u, copula, ...) ddvStCopula(matrix(u,ncol=2), copula, ...))
setMethod("ddvCopula", signature("matrix","spacetime_copula"), ddvStCopula)

invddvStCopula <- function(v, copula, y, h, tol=1e-12) {
  if(!is.matrix(h)) 
    h <- matrix(h,ncol=2)
  nElem <- length(v)
  stopifnot(nElem == length(y))
  stopifnot(nrow(h) == 1 | nrow(h)==nElem)
  
  .bisect(function(u) ddvStCopula(cbind(u, v), copula, h), y, tol=tol)
}

setMethod("invddvCopula", signature("numeric", "spacetime_copula"), invddvStCopula)

# log-likelihood by copula for all spatio-temporal lags


.loglikByStLags <- function(stBins, data, families = c(normalCopula(),
                                                             tCopula(),
                                                             claytonCopula(),
                                                             frankCopula(),
                                                             gumbelCopula()),
                                  calc_cor, lag_sub=1:length(stBins$meanDists)) {
  nTimeLags <- dim(stBins$lagCor)[1]
  if(is.null(nTimeLags))
    nTimeLags <- 1
  var <- attr(stBins, "variable")
  
  varValues <- .stValues(data, var)
  
  retrieveData <- function(spIndex, tempIndices) {
    binnedData <- NULL
    for (i in 1:(ncol(tempIndices)/2)) {
      binnedData <- cbind(binnedData, 
                          as.vector(varValues[spIndex[,1], tempIndices[,2*i-1], drop=FALSE]),
                          as.vector(varValues[spIndex[,2], tempIndices[,2*i], drop=FALSE]))
    }
    return(binnedData)
  }
  
  lagData <- lapply(stBins$lags[[1]][lag_sub], retrieveData, tempIndices=stBins$lags[[2]])
  
  tmpBins <- list(meanDists=stBins$meanDists[lag_sub])
  attr(tmpBins, "variable") <- var
  
  loglikTau <- list()
  for(j in 1:nTimeLags) {
    tmpLagData <- lapply(lagData, function(x) x[,c(2*j-1,2*j)])
    tmpLagData <- lapply(tmpLagData, function(pairs) {
      bool <- !is.na(pairs[,1]) & !is.na(pairs[,2])
      pairs[bool,]
    })
    
    if(missing(calc_cor))
      res <- loglik_by_lags.static(tmpLagData, families)
    else
      res <- loglik_by_lags.dyn(tmpBins, tmpLagData, families, 
                                     function(h) calc_cor(h, j, 1:nTimeLags))
    loglikTau[[paste("loglik",j,sep="")]] <- res
  }
  
  return(loglikTau)
}


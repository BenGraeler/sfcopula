# constructor
# dimension = "integer"     set to 2
# parameters = "numeric"    set of parameters
# param.names = "character" appropriate names
# param.lowbnd = "numeric"  appropriate lower bounds
# param.upbnd = "numeric"   appropriate upper bounds
# components="list"         list of copulas (will be automatically supplemented 
#			      by the independent copula)
# distances="numeric"       the linking distances + the range (will be assigned
#			      to the independent copula)
# unit="character"          measurement unit of distance
# depFun="function"         a optional spatial dependence function providing 
#                             Kendalls tau or Spearman's rho to calib* or exact 
#                             parameters
# combination="character"   "convex": distance dependent convex combination,
#                           "geometric": distance dependent geometric mean


spatial_copula <- function(components, distances, spDepFun, unit="m",
                           combination=c("convex", "geometric")) {
  combination <- match.arg(combination)
  if (combination == "geometric" && !missing(spDepFun))
    stop("A spatial dependence function is only supported for the convex combination.")
  
  indepCopInComponents <- sapply(components, function(x) class(x) == "indepCopula")
  if (any(indepCopInComponents) && combination == "convex") {
    components[which(indepCopInComponents)] <- list(normalCopula(0.01))
  }
  
  if (missing(spDepFun)) { 
    calibMoa <- function(copula, h) return(NULL)
  } else {
    if (is.na(match(spDepFun(NULL), c("kendall","spearman","id")))) 
      stop("spDepFun(NULL) must return 'spearman', 'kendall' or 'id'.")
    cat("The parameters of the components will be recalculated according to the provided spDepFun where possible. \nIn case no 1-1 relation is known, the copula as in components is used. \n")
    calibMoa <- switch(spDepFun(NULL), 
                       kendall=function(copula, h) iTau(copula, spDepFun(h)),
                       spearman=function(copula, h) iRho(copula, spDepFun(h)),
                       id=function(copula, h) return(h))
    
    for (i in 1:length(components)) {
      if(class(components[[i]]) != "indepCopula")
        param <- try(calibMoa(components[[i]], distances[i]),T)
      if (class(param) == "numeric")
        components[[i]]@parameters[1:length(param)] <- param
    }
  }

  param       <- unlist(lapply(components, 
                               function(x) {
                                 if(class(x)=="indepCopula") 
                                   return(NA)
                                 x@parameters}))
  param.names <- unlist(lapply(components, 
                               function(x) {
                                 if(class(x)=="indepCopula") 
                                   return(NA)
                                 x@param.names}))
  param.low   <- unlist(lapply(components, 
                               function(x) {
                                 if(class(x)=="indepCopula") 
                                   return(NA)
                                 x@param.lowbnd}))
  param.up    <- unlist(lapply(components, 
                               function(x) {
                                 if(class(x)=="indepCopula") 
                                   return(NA)
                                 x@param.upbnd}))
  
  new("spatial_copula", dimension=as.integer(2), parameters=param, param.names=param.names,
      param.lowbnd=param.low, param.upbnd=param.up,
      fullname=switch(combination,
                      convex="Spatial Copula: distance dependent convex combination of bivariate copulas",
                      geometric="Spatial Copula: distance dependent geometric mean of bivariate copulas"),
      components=components, distances=distances, calibMoa=calibMoa, unit=unit,
      combination=combination)
}

## show method
showCopula <- function(object) {
  cat(object@fullname, "\n")
  cat("Dimension: ", object@dimension, "\n")
  cat("Copulas:\n")
  for (i in 1:length(object@components)) {
    cmpCop <- object@components[[i]]
    cat("  ", describeCop(cmpCop, "very short"), "at", object@distances[i], 
        paste("[",object@unit,"]",sep=""), "\n")
  }
  if(!is.null(object@calibMoa(normalCopula(0),0))) 
    cat("A spatial dependence function is used. \n")
}

setMethod("show", signature("spatial_copula"), showCopula)

## spatial copula CDF
######################

pSpCopula <- function (u, copula, h, block=1) {
  if (missing(h)) 
    stop("Point pairs need to be provided with their separating distance \"h\".")
  if(length(h)>1 && length(h)!=nrow(u))
    stop("The distance vector must either be of the same length as rows in the data pairs or a single value.")
  
  if (copula@combination == "geometric")
    return(pCop.spGeomCop(u, copula, rep(h, length.out=nrow(u))))

  .spCopEval("p", u, copula, h)
}

setMethod(pCopula, signature("numeric","spatial_copula"), 
          function(u, copula, ...) pSpCopula(matrix(u,ncol=2),copula, ...))
setMethod(pCopula, signature("matrix","spatial_copula"), pSpCopula)

## spatial Copula density 
##########################

dSpCopula <- function (u, copula, log=FALSE, h) {
  if (missing(h)) 
    stop("Point pairs need to be provided with their separating distance \"h\".")
  if(length(h)>1 && length(h)!=nrow(u))
    stop("The distance vector must either be of the same length as rows in the data pairs or a single value.")
  
  if (copula@combination == "geometric")
    return(dCop.spGeomCop(u, copula, rep(h, length.out=nrow(u)), do.logs=log))

  .spCopEval("d", u, copula, h, log=log)
}

setMethod(dCopula, signature("numeric","spatial_copula"), 
          function(u, copula, log, ...) dSpCopula(matrix(u,ncol=2), copula, log=log, ...))
setMethod(dCopula, signature("matrix","spatial_copula"), dSpCopula)

## partial derivatives ##
## dduSpCopula
###############

dduSpCopula <- function (u, copula, h) {
  if (missing(h)) 
    stop("Point pairs need to be provided with their separating distance h.")
  if(length(h)>1 && length(h)!=nrow(u))
    stop("The distance vector must either be of the same length as rows in the data pairs or a single value.")
  
  if (copula@combination == "geometric")
    return(dduCop.spGeomCop(u, copula, rep(h, length.out=nrow(u))))

  .spCopEval("ddu", u, copula, h)
}

setMethod("dduCopula", signature("matrix","spatial_copula"), dduSpCopula)
setMethod("dduCopula", signature("numeric","spatial_copula"), 
          function(u, copula, ...) dduSpCopula(matrix(u,ncol=copula@dimension),copula, ...) )

invdduSpCopula <- function(u, copula, y, h, tol=1e-12) {
  nElem <- length(u)
  stopifnot(nElem == length(y))
  stopifnot(length(h) == 1 | length(h)==nElem)
  
  # d/du C(u, v) is increasing in v
  .bisect(function(v) dduSpCopula(cbind(u, v), copula, h), y, tol=tol)
}

setMethod("invdduCopula", signature("numeric", "spatial_copula"), invdduSpCopula)

## ddvSpCopula
###############

ddvSpCopula <- function (u, copula, h) {
  if (missing(h)) 
    stop("Point pairs need to be provided with their separating distance h.")
  if(length(h)>1 && length(h)!=nrow(u))
    stop("The distance vector must either be of the same length as rows in the data pairs or a single value.")
  
  if (copula@combination == "geometric")
    return(ddvCop.spGeomCop(u, copula, rep(h, length.out=nrow(u))))

  .spCopEval("ddv", u, copula, h)
}

setMethod("ddvCopula", signature("matrix","spatial_copula"), ddvSpCopula)
setMethod("ddvCopula", signature("numeric","spatial_copula"), 
          function(u, copula, ...) ddvSpCopula(matrix(u,ncol=copula@dimension),copula, ...) )

invddvSpCopula <- function(v, copula, y, h, tol=1e-12) {
  nElem <- length(v)
  stopifnot(nElem == length(y))
  stopifnot(length(h) == 1 | length(h) == nElem)
  
  # d/dv C(u, v) is increasing in u
  .bisect(function(u) ddvSpCopula(cbind(u, v), copula, h), y, tol=tol)
}

setMethod("invddvCopula", signature("numeric", "spatial_copula"), invddvSpCopula)

## simulation

spCop.rCop <- function(n, copula, h) {
  u <- runif(n)
  v <- invdduCopula(u, copula, y=runif(n), h=h)
  
  return(cbind(u, v))
}

setMethod("rCopula", signature("numeric", "spatial_copula"), spCop.rCop)

#############
##         ##
## FITTING ##
##         ##
#############

# two models: 
# 1) Kendall's tau driven:
#    fit curve through emp. Kendall's tau values, identify validity ranges for
#    copula families deriving parameters from the fit, fade from one family to 
#    another at borders
# 2) convex-linear combination of copulas: 
#    fit one per lag, fade from one to another

# towards the first model:

# INPUT: the stBinning
# steps
# a) fit a curve
# b) estimate bivariate copulas per lag (limited to those with some 1-1-relation 
#    to Kendall's tau')
# INTERMEDIATE RESULT
# c) select best fits based on ... e.g. log-likelihood, visual inspection
# d) compose bivariate copulas to one spatial copula
# OUTPUT: a spatial copula parametrised by distance through Kendall's tau

# towards a)
# bins   -> typically output from calc_bins
# degree -> the degree of the polynominal
# cutoff -> maximal distance that should be considered for fitting
# bounds -> the bounds of the correlation function (typically c(0,1))
# method -> the measure of association, either "kendall" or "spearman"
fitCorFunSng <- function(bins, degree, cutoff, bounds, cor.method, weighted) {
  if (weighted) {
    bins <- as.data.frame(bins[c("np","meanDists","lagCor")])
    if(!is.na(cutoff)) 
      bins <- bins[bins$meanDists <= cutoff,]
    fitCor <- lm(lagCor ~ poly(meanDists, degree), data = bins, weights=bins$np)
  } else {
    bins <- as.data.frame(bins[c("meanDists","lagCor")])
    if(!is.na(cutoff)) 
      bins <- bins[bins$meanDists <= cutoff,]
    fitCor <- lm(lagCor ~ poly(meanDists, degree), data = bins)
  }
  
  print(fitCor)
  cat("Sum of squared residuals:",sum(fitCor$residuals^2),"\n")
  
  if(cor.method=="fasttau") 
    cor.method <- "kendall"
  
  function(x) {
    if (is.null(x)) return(cor.method)
    return(pmin(bounds[2], pmax(bounds[1], 
                                eval(predict(fitCor, data.frame(meanDists=x))))))
  }
}

fit_cor_fun <- function(bins, degree=3, cutoff=NA, tlags, bounds=c(0,1), 
                      cor.method=NULL, weighted=FALSE){
  if(is.null(cor.method)) {
    if(is.null(attr(bins,"cor.method")))
      stop("Neither the bins arguments has an attribute cor.method nor is the parameter cor.method provided.") 
    else 
      cor.method <- attr(bins,"cor.method")
  } else {
    if(!is.null(attr(bins,"cor.method")) && cor.method != attr(bins,"cor.method"))
      stop("The cor.method attribute of the bins argument and the argument cor.method do not match.")
  }
  
  if(is.null(nrow(bins$lagCor))) # the spatial case
    return(fitCorFunSng(bins, degree, cutoff, bounds, cor.method, weighted))
    
  # the spatio-temporal case
  degree <- rep(degree, length.out = nrow(bins$lagCor))
  calcKTau <- list()
  for (j in 1:nrow(bins$lagCor)) {
    calcKTau[[paste("fun",j,sep="")]] <- fitCorFunSng(data.frame(np=bins$lagNp[j,],
                                                                 meanDists=bins$meanDists, 
                                                                 lagCor=bins$lagCor[j,]),
                                                      degree[j], cutoff, bounds, 
                                                      cor.method, weighted)
  }
  
  tlsort <- sort(tlags,decreasing=TRUE)
  
  corFun <- function(h, time, tlags=tlsort) {
    t <- which(tlags==time)
    calcKTau[[time]](h)
  }
  
  attr(corFun, "tlags") <- sort(tlags, decreasing=TRUE)
  return(corFun)
}


# towards b)
  
## loglikelihoods for a dynamic spatial copula
loglik_by_lags.dyn <- function(bins, lagData, families, calcCor) {
  moa <- switch(calcCor(NULL),
                kendall=function(copula, h) iTau(copula, calcCor(h)),
                spearman=function(copula, h) iRho(copula, calcCor(h)),
                id=function(copula, h) calcCor(h))
  
  loglik <- NULL
  copulas <- list()
  for (cop in families) {
    cat(describeCop(cop, "very short"),"\n")
    tmploglik <- NULL
    tmpCop <- list()
    
    pb <- txtProgressBar(0, length(bins$meanDists), style=3)
    for(i in 1:length(bins$meanDists)) {
      if(class(cop)!="indepCopula") {
        if(class(cop) == "asCopula") {
          cop <- switch(calcCor(NULL),
                        kendall=fitASC2.itau(cop, lagData[[i]], 
                                              tau=calcCor(bins$meanDists[i]))@copula,
                        spearman=fitASC2.irho(cop, lagData[[i]],
                                              rho=calcCor(bins$meanDists[i]))@copula,
                        stop(paste(calcCor(NULL), "is not yet supported.")))
          param <- cop@parameters
        } else {
          if(class(cop) == "cqsCopula") {
            cop <- switch(calcCor(NULL),
                          kendall=fitCQSec.itau(cop, lagData[[i]], 
                                                tau=calcCor(bins$meanDists[i]))@copula,
                          spearman=fitCQSec.irho(cop, lagData[[i]],
                                                rho=calcCor(bins$meanDists[i]))@copula,
                          stop(paste(calcCor(NULL), "is not yet supported.")))
            param <- cop@parameters
          } else {
            param <- moa(cop, bins$meanDists[i])
            if(!is.na(param))
              cop@parameters[1:length(param)] <- param
          }
        }
      }
      
      if(any(is.na(param)))
        tmploglik <- c(tmploglik, NA)
      else 
        tmploglik <- c(tmploglik, sum(log(.compEval("d", lagData[[i]], cop))))
      tmpCop <- append(tmpCop, cop)
      setTxtProgressBar(pb, i)
    }
    close(pb)
    loglik <- cbind(loglik, tmploglik)
    copulas[[class(cop)]] <- tmpCop
  }

  colnames(loglik) <- sapply(families, function(x) class(x)[1])

  return(list(loglik=loglik, copulas=copulas))
}

## loglikelihoods for a static spatial copula
loglik_by_lags.static <- function(lagData, families) {
  
  fits <-lapply(families, 
                function(cop) {
                  cat(describeCop(cop, "very short"), "\n")
                  lapply(lagData,
                         function(x) {
                           tryCatch(fitCopula(cop, x, estimate.variance = FALSE),
                                    error=function(e) return(NA))
                         })
                })
  
  loglik <- lapply(fits, function(x) sapply(x, function(fit) {
    if(class(fit)=="fitCopula")
      return(fit@loglik)
    else
      return(NA)
  }))
  
  loglik <- matrix(unlist(loglik),ncol=length(loglik),byrow=F)
  colnames(loglik) <- sapply(families, function(x) class(x)[1])
  
  copulas <- lapply(fits, function(x) sapply(x, function(fit) {
    if(class(fit)=="fitCopula")
      return(fit@copula)
    else
      return(NULL)
  }))

  names(copulas) <- colnames(loglik)
  
  return(list(loglik=loglik, copulas=copulas))
}

##

loglik_by_lags <- function(bins, data, families=c(normalCopula(), 
                                                       tCopula(),
                                                       claytonCopula(), frankCopula(), 
                                                       gumbelCopula()),
                                calcCor, lagSub=1:length(bins$meanDists)) {
  # spatio-temporal bins (from calc_bins on a stars cube) hold spatial and temporal lags
  if (is.list(bins$lags) && all(c("sp", "time") %in% names(bins$lags)))
    return(.loglikByStLags(bins, data, families, calcCor, lagSub))
  
  var <- attr(bins, "variable")
  
  if(missing(data)) {
    lagData <- bins$lagData
  }
  else {
    lagData <- lapply(bins$lags[lagSub], 
                      function(x) {
                        varValues <- .attrValues(data, var)
                        cbind(varValues[x[, 1]], varValues[x[, 2]])
    })
  }
  
  lagData <- lapply(lagData, 
                    function(pairs) {
                      bool <- !is.na(pairs[,1]) & !is.na(pairs[,2])
                      pairs[bool,]
                    })
  
  if(missing(calcCor))
    return(loglik_by_lags.static(lagData, families))
  else
    return(loglik_by_lags.dyn(lapply(bins, function(x) x[lagSub]),
                                          lagData, families, calcCor))
}



# towards d)
compose_spatial_copula <- function(bestFit, families, bins, calcCor, range=max(bins$meanDists)) {
  nFits <- length(bestFit)
  if(nFits > length(bins$meanDists))
    stop("There may not be less bins than best fits.\n")
  rangeIndex <- min(nFits, max(which(bins$meanDists <= range)))
  
  if (missing(calcCor)) {
    return(spatial_copula(components = as.list(families[bestFit[1:rangeIndex]]),
                    distances = bins$meanDists[1:rangeIndex], 
                    unit = "m"))
  }
  
  else {
    rangeIndex <- min(rangeIndex, which(calcCor(bins$meanDists) <= 0))
    
    return(spatial_copula(components = as.list(families[bestFit[1:rangeIndex]]),
                    distances = bins$meanDists[1:rangeIndex], 
                    unit = "m", spDepFun = calcCor))
  }
}

# in once

# bins   -> typically output from calc_bins
# cutoff -> maximal distance that should be considered for fitting
# families -> a vector of dummy copula objects of each family to be considered
#             DEFAULT: c(normal, t_df=4, clayton, frank, gumbel
# ...
# type   -> the type of curve (by now only polynominals are supported)
# degree -> the degree of the polynominal
# bounds -> the bounds of the correlation function (typically c(0,1))
# method -> the measure of association, either "kendall" or "spearman"
fit_spatial_copula <- function(bins, data, cutoff=NA, 
                        families=c(normalCopula(), tCopula(),
                                   claytonCopula(), frankCopula(),
                                   gumbelCopula()), ...) {
  calcCor <- fit_cor_fun(bins, cutoff=cutoff, ...)
  loglik <- loglik_by_lags(bins, data, families, calcCor)
  
  bestFit <- apply(apply(loglik$loglik, 1, rank),2, 
                   function(x) which(x==length(families)))
  
  return(compose_spatial_copula(bestFit, families, bins, calcCor, range=cutoff))
}


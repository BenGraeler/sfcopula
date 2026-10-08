## spatial Gaussian Copula

# "density" evaluation
spatial_gauss_loglik <- function(corFun, neigh, dataLocs, log=TRUE) {
  if (is(corFun, "spatial_gauss_copula"))
    corFun <- corFun@corFun
  neighDim <- ncol(neigh@data) - length(neigh@coVar)
  
  allDataDists <- .spDistMat(.pointGeom(dataLocs))
  
  pb <- txtProgressBar(0, nrow(neigh@data), 0, width = getOption("width") - 10, style = 3)
  
  loglik <- 0
  
  for(i in 1:nrow(neigh@data)) { # i <- 2
    setTxtProgressBar(pb, i)
    tmpDists <- allDataDists[neigh@index[i,], neigh@index[i,]]
    
    tmpCor <- corFun(tmpDists)
    
    tmpGaussCop <- normalCopula(tmpCor[lower.tri(tmpCor)], neighDim, dispstr="un")
    
    loglik <- loglik + dCopula(as.numeric(neigh@data[i, 1:neighDim]), tmpGaussCop, log=T)
  }
  close(pb)
  
  if(log)
    return(loglik)
  else
    return(exp(loglik))
}

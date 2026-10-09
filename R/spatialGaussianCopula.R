## spatial Gaussian Copula

# "density" evaluation
spatial_gauss_loglik <- function(cor_fun, neigh, data, log=TRUE) {
  if (is(cor_fun, "spatial_gauss_copula"))
    cor_fun <- cor_fun@corFun
  neighDim <- ncol(neigh@data) - length(neigh@coVar)
  
  allDataDists <- .spDistMat(.pointGeom(data))
  
  pb <- txtProgressBar(0, nrow(neigh@data), 0, width = getOption("width") - 10, style = 3)
  
  loglik <- 0
  
  for(i in 1:nrow(neigh@data)) { # i <- 2
    setTxtProgressBar(pb, i)
    tmpDists <- allDataDists[neigh@index[i,], neigh@index[i,]]
    
    tmpCor <- cor_fun(tmpDists)
    
    tmpGaussCop <- normalCopula(tmpCor[lower.tri(tmpCor)], neighDim, dispstr="un")
    
    loglik <- loglik + dCopula(as.numeric(neigh@data[i, 1:neighDim]), tmpGaussCop, log=T)
  }
  close(pb)
  
  if(log)
    return(loglik)
  else
    return(exp(loglik))
}

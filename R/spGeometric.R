# geometric means of copulas, used by spatial_copula(..., combination = "geometric")

## dists is expected to contain 0 and components is expected to contain a copula 
## for 0-distance, either perfect dependence (nugget free), or some strongly 
## correlated copula; components needs to contain one copula that is used beyond
## the range, typically the product copula

pCop.spGeomCop <- function(u, copula, h, ...) {
  dists <- copula@distances

  resLow <- numeric(nrow(u))
  resHigh <- numeric(nrow(u))
  p <- numeric(nrow(u))
  
  pairsInt <- findInterval(h, dists, left.open = TRUE, rightmost.closed = TRUE)
  
  for (int in unique(pairsInt)) {
    sel <- which(pairsInt == int)
    
    if (length(sel) == 0) 
      next;
    if (int >= length(copula@components)) {
      resLow[sel] <- resHigh[sel] <- .compEval("p", u[sel,,drop=FALSE], copula@components[[int]])
      
      p[sel] <- 0.5
    } else {
      resLow[sel] <- .compEval("p", u[sel,,drop=FALSE], copula@components[[int]])
      resHigh[sel] <- .compEval("p", u[sel,,drop=FALSE], copula@components[[int+1]])
    
      p[sel] <- (dists[int+1] - h[sel])/diff(dists[int+c(0,1)])
    }
  }
  
  resLow^p * resHigh^(1-p)
}

dCop.spGeomCop <- function(u, copula, h, do.logs=F, ...) {
  dists <- copula@distances

  dCopLow <- matrix(NA, nrow(u), 4)
  dCopHigh <- matrix(NA, nrow(u), 4)
  p <- numeric(nrow(u))
  
pairsInt <- findInterval(h, dists, left.open = TRUE, rightmost.closed = TRUE)
  
  for (int in unique(pairsInt)) {
    sel <- which(pairsInt == int)
    
    if (length(sel) == 0) 
      next;
    if (int >= length(copula@components)) {
      tmpPairs <- u[sel,,drop=FALSE]
      copLow <- copula@components[[int]]
      dCopHigh[sel,] <- dCopLow[sel,] <- cbind(.compEval("p", tmpPairs, copLow),
                                               .compEval("ddu", tmpPairs, copLow),
                                               .compEval("ddv", tmpPairs, copLow),
                                               .compEval("d", tmpPairs, copLow))
      
      p[sel] <- 0.5
    } else {
      tmpPairs <- u[sel,,drop=FALSE]
      copLow <- copula@components[[int]]
      copHigh <- copula@components[[int+1]]
      dCopLow[sel,] <- cbind(.compEval("p", tmpPairs, copLow),
                             .compEval("ddu", tmpPairs, copLow),
                             .compEval("ddv", tmpPairs, copLow),
                             .compEval("d", tmpPairs, copLow))
      dCopHigh[sel,] <- cbind(.compEval("p", tmpPairs, copHigh),
                             .compEval("ddu", tmpPairs, copHigh),
                             .compEval("ddv", tmpPairs, copHigh),
                             .compEval("d", tmpPairs, copHigh))

      p[sel] <- (dists[int+1] - h[sel])/diff(dists[int+c(0,1)])
    }
  }
  
  # the formulas below are written for L^(1-w) * H^w, the cdf is L^p * H^(1-p)
  p <- 1 - p

# f[u,v]^(-1-p)        g[u,v]^(-2+p)         ((-1+p)   (p   (g[u,v]       f^(0,1)[u,v] -f[u,v]      g^(0,1)[u,v])   (g[u,v]       f^(1,0)[u,v]- f[u,v]      g^(1,0)[u,v]) - f[u,v]        g[u,v]^2       f^(1,1)[u,v])+ p   f[u,v]^2      g[u,v]       g^(1,1)[u,v])
  res <- dCopLow[,1]^(-1-p) * dCopHigh[,1]^(-2+p) * ((-1+p) * (p * (dCopHigh[,1]*dCopLow[,3] - dCopLow[,1]*dCopHigh[,3]) * (dCopHigh[,1]*dCopLow[,2] - dCopLow[,1]*dCopHigh[,2]) - dCopLow[,1] * dCopHigh[,1]^2*dCopLow[,4]) + p * dCopLow[,1]^2*dCopHigh[,1]*dCopHigh[,4])
  
  if (do.logs)
    return(log(res))
  else 
    return(res)
}

dduCop.spGeomCop <- function(u, copula, h, do.logs=F, ...) {
  dists <- copula@distances
  
  dCopLow <- matrix(NA, nrow(u), 4)
  dCopHigh <- matrix(NA, nrow(u), 4)
  p <- numeric(nrow(u))
  
  pairsInt <- findInterval(h, dists, left.open = TRUE, rightmost.closed = TRUE)
  
  for (int in unique(pairsInt)) {
    sel <- which(pairsInt == int)
    
    if (length(sel) == 0) next;
    
    if (int >= length(copula@components)) {
      tmpPairs <- u[sel,,drop=FALSE]
      copLow <- copula@components[[int]]
      dCopHigh[sel, ] <- dCopLow[sel,] <- cbind(.compEval("p", tmpPairs, copLow),
                                                .compEval("ddu", tmpPairs, copLow))
      
      p[sel] <- 0.5
    } else {
      tmpPairs <- u[sel,,drop=FALSE]
      copLow <- copula@components[[int]]
      copHigh <- copula@components[[int+1]] # check for max 
      dCopLow[sel,] <- cbind(.compEval("p", tmpPairs, copLow),
                             .compEval("ddu", tmpPairs, copLow))
      dCopHigh[sel,] <- cbind(.compEval("p", tmpPairs, copHigh),
                              .compEval("ddu", tmpPairs, copHigh))
      
      p[sel] <- (dists[int+1] - h[sel])/diff(dists[int+c(0,1)])
    }
  }
  
  # the formulas below are written for L^(1-w) * H^w, the cdf is L^p * H^(1-p)
  p <- 1 - p

# (1-c)   f[u,v]^-c        g[u,v]^c         f^(1,0)[u,v]+ c   f[u,v]^(1-c)        g[u,v]^(-1+c)        g^(1,0)[u,v]
  (1-p) * dCopLow[,1]^-p * dCopHigh[,1]^p * dCopLow[,2] + p * dCopLow[,1]^(1-p) * dCopHigh[,1]^(p-1) * dCopHigh[,2]
}

ddvCop.spGeomCop <- function(u, copula, h) {
  dists <- copula@distances
  
  dCopLow <- matrix(NA, nrow(u), 2)
  dCopHigh <- matrix(NA, nrow(u), 2)
  p <- numeric(nrow(u))
  
  pairsInt <- findInterval(h, dists, left.open = TRUE, rightmost.closed = TRUE)
  
  for (int in unique(pairsInt)) {
    sel <- which(pairsInt == int)
    
    if (length(sel) == 0) next;
    
    if (int >= length(copula@components)) {
      tmpPairs <- u[sel,,drop=FALSE]
      copLow <- copula@components[[int]]
      dCopHigh[sel, ] <- dCopLow[sel,] <- cbind(.compEval("p", tmpPairs, copLow),
                                                .compEval("ddv", tmpPairs, copLow))
      
      p[sel] <- 0.5
    } else {
      tmpPairs <- u[sel,,drop=FALSE]
      copLow <- copula@components[[int]]
      copHigh <- copula@components[[int+1]]
      dCopLow[sel,] <- cbind(.compEval("p", tmpPairs, copLow),
                             .compEval("ddv", tmpPairs, copLow))
      dCopHigh[sel,] <- cbind(.compEval("p", tmpPairs, copHigh),
                              .compEval("ddv", tmpPairs, copHigh))
      
      p[sel] <- (dists[int+1] - h[sel])/diff(dists[int+c(0,1)])
    }
  }
  
  # the formulas below are written for L^(1-w) * H^w, the cdf is L^p * H^(1-p)
  p <- 1 - p

  # (1-c) f[u,v]^-c          g[u,v]^c         f^(0,1)[u,v]+ c f[u,v]^(1-c)      g[u,v]^(-1+c)        (g^(0,1))[u,v]
  (1-p) * dCopLow[,1]^(-p) * dCopHigh[,1]^p * dCopLow[,2] + p*dCopLow[,1]^(1-p)*dCopHigh[,1]^(p-1) * dCopHigh[,2]
  
}

##########################################################
## distance vine copula: vine copulas whose lower trees ##
## are spatial or spatio-temporal bivariate copulas     ##
##########################################################

# constructor
distance_vine_copula <- function(trees, top_cop=NULL) {
  if(!is.list(trees))
    trees <- list(trees)

  isSt <- all(sapply(trees, is, "spacetime_copula"))
  kind <- if (isSt) "Spatio-temporal" else "Spatial"

  if(is.null(top_cop)) {
    dim <- length(trees) + 1
    fullname <- paste(kind, "vine copula family with only",
                      ifelse(isSt, "spatio-temporal", "spatial"), "tree(s).")
  } else {
    dim <- top_cop@dimension + length(trees)
    fullname <- paste(kind, "vine copula family with", length(trees),
                      ifelse(isSt, "spatio-temporal", "spatial"), "tree(s).")
  }

  new("distance_vine_copula", dimension = as.integer(dim),
      parameters=numeric(), param.names = character(), param.lowbnd = numeric(),
      param.upbnd = numeric(), fullname = fullname, trees=trees, topCop=top_cop)
}

# show
showDistanceVineCopula <- function(object) {
  cat(object@fullname, "\n")
  cat("Dimension: ", object@dimension, "\n")
}

setMethod("show", signature("distance_vine_copula"), showDistanceVineCopula)

## helpers for the distances of a single tree
##############################################

.isStTree <- function(tree) is(tree, "spacetime_copula")

# bring the distances of one tree into a standard form:
# spatial: matrix [n, k]; spatio-temporal: array [n, k, 2]
.treeDists <- function(h, tree) {
  if (.isStTree(tree)) {
    if (length(dim(h)) == 2 && ncol(h) == 2)  # a single neighbour pair per row
      h <- array(h, c(nrow(h), 1, 2))
    if (length(dim(h)) != 3)
      stop("Spatio-temporal trees need distances as array [n, neighbours, 2].")
    return(h)
  }
  if(!is.matrix(h))
    h <- matrix(h, ncol=length(h))
  h
}

# distances of the i-th pair of a tree in the form the bivariate copula expects
.pairH <- function(h, i, tree) {
  if (.isStTree(tree))
    return(matrix(h[,i,], ncol=2))
  h[,i]
}

# drop one tree: log-likelihood contributions and conditioned data
.condTree <- function(u0, h, tree, loglik=TRUE) {
  h <- .treeDists(h, tree)
  l0 <- rep(0, nrow(u0))
  u1 <- matrix(NA, nrow(u0), dim(h)[2])
  for(i in 1:dim(h)[2]) {
    res <- .treeEval(u0[,c(1,i+1), drop=FALSE], tree, .pairH(h, i, tree), loglik)
    if (loglik)
      l0 <- l0 + res$d
    u1[,i] <- res$ddu
  }
  list(loglik=l0, u=u1)
}

## density
###########

dDistanceVine <- function(u, copula, log=FALSE, h) {
  if (!is.list(h))
    h <- list(h)
  stopifnot(length(copula@trees)==length(h))

  l0 <- rep(0, nrow(u)) # density contributions of the distance trees
  u0 <- u # previous level's conditional data
  for(tree in seq_along(copula@trees)) {
    res <- .condTree(u0, h[[tree]], copula@trees[[tree]])
    l0 <- l0 + res$loglik
    u0 <- res$u
  }

  if(!is.null(copula@topCop))
    l1 <- .topLogDens(u0, copula@topCop)
  else
    l1 <- 0

  if(log)
    return(l0+l1)
  else
    return(exp(l0+l1))
}

setMethod("dCopula", signature=signature("matrix","distance_vine_copula"),
          function(u, copula, log=FALSE, ...) dDistanceVine(u, copula, log=log, ...))
setMethod("dCopula", signature=signature("numeric","distance_vine_copula"),
          function(u, copula, log=FALSE, ...)
            dDistanceVine(matrix(u, ncol=copula@dimension), copula, log=log, ...))
setMethod("dCopula", signature=signature("data.frame","distance_vine_copula"),
          function(u, copula, log=FALSE, ...)
            dDistanceVine(as.matrix(u), copula, log=log, ...))

## distances of all trees
##########################

# returns a list with the distances of each tree; the first tree uses the
# distances stored in the neighbourhood, higher trees the distances between
# the neighbours of the previous tree
tree_dists <- function(neigh, data, n_trees) {
  stopifnot(is(neigh, "neighbourhood"))
  condDists <- list(neigh@distances)
  if(n_trees==1)
    return(condDists)

  geom <- .pointGeom(data)
  nNeighs <- dim(neigh@distances)[2]
  for (tree in 1:(n_trees-1)) {
    condDists[[tree+1]] <- .nextTreeDists(neigh@index, geom, tree, nNeighs - tree)
  }
  return(condDists)
}

# distances between the conditioning neighbour of tree 'tree' (column tree+1
# of the index) and the remaining neighbours of the neighbourhood
.nextTreeDists <- function(index, geom, tree, nPairs) {
  if (length(dim(index)) == 3) {
    h1 <- array(NA, c(dim(index)[1], nPairs, 2))
    for(i in seq_len(nPairs)) {
      h1[,i,1] <- .pairDists(geom, index[,tree+1,1], index[,tree+i+1,1])
      h1[,i,2] <- index[,tree+i+1,2] - index[,tree+1,2]
    }
    return(h1)
  }
  h1 <- matrix(NA, nrow(index), nPairs)
  for(i in seq_len(nPairs))
    h1[,i] <- .pairDists(geom, index[,tree+1], index[,tree+i+1])
  h1
}

## dropping a tree, returning a conditional neighbourhood
##########################################################

drop_tree <- function(neigh, data, copula) {
  stopifnot(is(neigh, "neighbourhood"))
  isSt <- .isStNeighbourhood(neigh)
  if (isSt != .isStTree(copula))
    stop("Spatio-temporal neighbourhoods need a spacetime_copula, spatial ones a spatial_copula.")

  nCoVar <- length(neigh@coVar)
  u0 <- as.matrix(neigh@data)
  nNeighs <- dim(neigh@distances)[2]
  res <- .condTree(u0[, 1:(nNeighs+1), drop=FALSE], neigh@distances, copula, loglik=FALSE)
  u1 <- res$u

  geom <- .pointGeom(data)
  h1 <- .nextTreeDists(neigh@index, geom, 1, nNeighs - 1)

  # name the conditioned variable
  varSplit <- strsplit(neigh@var,"|",fixed=TRUE)[[1]]
  cond <- suppressWarnings(as.numeric(varSplit[length(varSplit)]))

  if (isSt) {
    cond <- if (is.na(cond)) paste(neigh@var, "|0", sep = "") else paste(neigh@var, cond + 1, sep = "")
    return(neighbourhood(data = u1, distances = h1, index = neigh@index[, -1, , drop=FALSE],
                         var = cond, prediction = neigh@prediction))
  }

  # as in the spatio-temporal case, covariates are not carried over; they are
  # conditioned on the central location with cond_covariate()
  if(is.na(cond)) {
    var <- paste(neigh@var,"|0",sep="")
    colnames(u1) <- paste(paste("N", rep(1:(ncol(u1)), each=length(var)), sep=""),
                          rep(var,ncol(u1)),sep=".")
  } else {
    var <- paste(neigh@var,cond+1,sep="")
    colnames(u1) <- paste(paste("N", rep(cond:(ncol(u1)+cond-1)+2,
                                         each=length(var)), sep=""),
                          rep(var,ncol(u1)),sep=".")
  }
  return(neighbourhood(data=u1, distances=h1, index=neigh@index[,-1,drop=FALSE],
                       var=var, prediction=neigh@prediction))
}

## fitting the vine for given distance trees
##############################################

fitDistanceVine <- function(copula, data,
                            method = list(StructureSelect = FALSE,
                                          indeptest = FALSE,
                                          familyset = NA),
                            estimate.variance=FALSE) {
  if (is(data, "neighbourhood")) {
    neigh <- data
    dataLocs <- NULL
  } else {
    stopifnot(is.list(data), length(data)==2)
    neigh <- data[[1]]
    dataLocs <- data[[2]]
  }
  stopifnot(is(neigh, "neighbourhood"))
  # covariates (last columns) are not part of a distance vine
  nVars <- ncol(neigh@data) - length(neigh@coVar)
  stopifnot(copula@dimension == nVars)

  nTrees <- length(copula@trees)
  if (nTrees > 1 && is.null(dataLocs))
    stop("Fitting more than one distance tree needs the data locations: data = list(neigh, data locations).")

  dists <- if (nTrees > 1) tree_dists(neigh, dataLocs, nTrees) else list(neigh@distances)

  u0 <- as.matrix(neigh@data)[, 1:nVars, drop=FALSE] # previous level's (conditional) data
  l0 <- rep(0,nrow(u0)) # density of the distance trees
  for(tree in 1:nTrees) {
    cat("[Dropping ", tree, ". distance tree.]\n",sep="")
    res <- .condTree(u0, dists[[tree]], copula@trees[[tree]])
    l0 <- l0 + res$loglik
    u0 <- res$u
  }

  if (ncol(u0)==1) {
    cat("[No copula to be estimated at the top.]\n")
    top_cop <- NULL
    loglik <- 0
  } else if (ncol(u0)==2) {
    cat("[Estimating a single bivariate copula at the top.]\n")
    bivCop <- BiCopSelect(u0[,1],u0[,2])
    top_cop <- copulaFromFamilyIndex(bivCop$family, bivCop$par, bivCop$par2)
    loglik <- sum(dCopula(u0, top_cop, log=TRUE))
  } else {
    cat("[Estimating a",ncol(u0),"dimensional copula at the top.]\n")
    top_cop <- copula@topCop
    vineCopFit <- fitCopula(top_cop, u0, method)
    top_cop <- vineCopFit@copula
    loglik <- vineCopFit@loglik
  }

  vineCop <- distance_vine_copula(copula@trees, top_cop)

  return(new("fitCopula", estimate = vineCop@parameters, var.est = matrix(NA),
             method = paste(sapply(method, paste, collapse=", "), collapse="; "),
             loglik = sum(l0)+loglik,
             fitting.stats=list(convergence = as.integer(NA)),
             nsample = nrow(neigh@data), copula=vineCop))
}

setMethod("fitCopula", signature=signature("distance_vine_copula"), fitDistanceVine)

## conditional density of the central location
################################################

setGeneric("cond_vine", function(cond_var, dists, vine, n = 1000, ...) standardGeneric("cond_vine"))

# evaluation grid with some points in the tails
.condGrid <- function(n) {
  rat <- 50:1%x%c(1e-6,1e-5,1e-4,1e-3)
  unique(sort(c(rat, 1 - rat, 1:(n - 1)/n)))
}

# turns density values on the grid into a normalised density function on [0,1]
.condDensityFun <- function(xVals, density) {
  nx <- length(xVals)
  # the 1-e6 corners linearily to [0,1], but keep non-negative
  density <- c(max(0,2*density[1]-density[2]),
               density, max(0,2*density[nx]-density[nx-1]))
  linAppr <- approxfun(c(0, xVals, 1), density)

  # sum up the denstiy to rescale
  int <- sum(diff(c(0,xVals,1))*(0.5*diff(density)+density[-(nx+2)]))
  condVineFun <- function(u) linAppr(u)/int
  attr(condVineFun,"xVals") <- c(0,xVals,1)
  return(condVineFun)
}

condDistanceVine <- function (cond_var, dists, vine, n = 1000, ...) {
  if (missing(cond_var))
    cond_var <- .dotsArg(list(...), "condVar", "cond_var", NULL)
  if (!is.list(dists))
    dists <- list(dists)
  stopifnot(length(vine@trees)==length(dists))

  xVals <- .condGrid(n)
  repCondVar <- matrix(cond_var, ncol = length(cond_var), nrow = length(xVals), byrow = T)
  density <- dDistanceVine(cbind(xVals, repCondVar), vine, h = dists)

  .condDensityFun(xVals, density)
}

setMethod("cond_vine", signature(vine = "distance_vine_copula"), condDistanceVine)

## simulation
## Algorithm 1 from Aas et al. (2006): Pair-copula constructions of multiple dependence
## h: list with the distances of each tree; spatial trees: vectors,
##    spatio-temporal trees: matrices with columns space and time

r.distanceVine <- function(n, copula, h) {
  if (!is.null(copula@topCop))
    stop("Simulation is only implemented for vine copulas without top copula.")
  vineDim <- copula@dimension
  trees <- copula@trees

  pairH <- function(tree, j) {
    hk <- h[[tree]]
    if (.isStTree(trees[[tree]])) {
      if (length(dim(hk)) == 3)
        hk <- matrix(hk[1,,], ncol=2)
      return(matrix(hk[j,], ncol=2))
    }
    hk[j]
  }

  # all n draws at once; v[, i, j] corresponds to v[i, j] of Algorithm 1
  init <- matrix(runif(n*vineDim), n, vineDim, byrow=TRUE)
  v <- array(NA, c(n, vineDim, vineDim))
  v[,1,1] <- init[,1]
  sims <- matrix(NA, n, vineDim)
  sims[,1] <- init[,1]
  for (i in 2:vineDim) {
    v[,i,1] <- init[,i]
    for (k in (i-1):1) {
      v[,i,1] <- invddvCopula(v[,k,k], trees[[k]], y=v[,i,1], h=pairH(k, i-k))
    }
    sims[,i] <- v[,i,1]
    if(i==vineDim)
      break()
    for(j in 1:(i-1)) {
      v[,i,j+1] <- ddvCopula(cbind(v[,i,j], v[,j,j]), trees[[j]], h=pairH(j, i-j))
    }
  }

  sims
}

setMethod("rCopula", signature("numeric","distance_vine_copula"),
          function(n, copula, ...) r.distanceVine(n, copula, ...))

##################################################################
## covariate vine copula: a distance tree (spatial or spatio-   ##
## temporal) plus a covariate copula at the central location    ##
##################################################################

setClassUnion("distance_copula", c("spatial_copula", "spacetime_copula"))

validCovariateVineCopula <- function(object) {
  if (!is.function(object@coVarCop))
    return("coVarCop needs to be a function returning a bivariate copula.")
  return(validObject(object@topCop))
}

setClass("covariate_vine_copula",
         representation("copula", coVarCop="function", tree="distance_copula", topCop="copula"),
         validity = validCovariateVineCopula, contains=list("copula"))

## constructor ##
#################

# coVarCop: function of the index of the central location (location for spatial,
#           (location, time) for spatio-temporal trees) returning a bivariate copula
# tree:     spatial_copula or spacetime_copula coupling the central location with its neighbours
# topCop:   copula joining covariate|centre and neighbours|centre
covariate_vine_copula <- function(coVarCop, tree, topCop) {
  stopifnot(is(tree, "distance_copula"))
  kind <- ifelse(.isStTree(tree), "spatio-temporal", "spatial")

  new("covariate_vine_copula", dimension = as.integer(topCop@dimension+1),
      parameters=numeric(), param.names = character(), param.lowbnd = numeric(),
      param.upbnd = numeric(),
      fullname = paste("Covariate vine copula family with 1", kind, "tree."),
      coVarCop=coVarCop, tree=tree, topCop=topCop)
}

## show ##
##########

setMethod("show", signature("covariate_vine_copula"), function(object) {
  cat(object@fullname, "\n")
  cat("Dimension: ", object@dimension, "\n")
})

## covariate copula at the central locations ##
###############################################

# centre: index of the central location, a vector (spatial) or a matrix with
# columns location and time (spatio-temporal), one row per pair or a single one
.coVarCond <- function(pairs, coVarCop, centre, isSt, loglik=TRUE) {
  centre <- matrix(centre, ncol = if (isSt) 2 else 1)
  stopifnot(nrow(centre) == 1 || nrow(centre) == nrow(pairs))

  l <- numeric(nrow(pairs))
  u <- numeric(nrow(pairs))
  key <- apply(centre, 1, paste, collapse=" ")
  for (k in unique(key)) {
    rows <- if (nrow(centre) == 1) seq_len(nrow(pairs)) else which(key == k)
    cop <- coVarCop(centre[match(k, key),])
    if (loglik)
      l[rows] <- dCopula(pairs[rows,,drop=FALSE], cop, log=TRUE)
    u[rows] <- dduCopula(pairs[rows,,drop=FALSE], cop)
  }
  list(loglik=l, u=u)
}

# conditions the covariate (last data column) on the central location (first
# column); coVarCop receives the index of the central location: (location, time)
# for spatio-temporal and the location for spatial neighbourhoods
cond_covariate <- function(neigh, coVarCop) {
  stopifnot(length(neigh@coVar) == 1)
  uv <- as.matrix(neigh@data[,c(1,ncol(neigh@data))])
  .coVarCond(uv, coVarCop, .centreIndex(neigh), .isStNeighbourhood(neigh), loglik=FALSE)$u
}

.centreIndex <- function(neigh) {
  if (.isStNeighbourhood(neigh)) neigh@index[,1,] else neigh@index[,1]
}

## density ##
#############

# u: columns as in the neighbourhood data: central location, neighbours, covariate (last)
# h: distances of the neighbours: spatial [n or 1, k], spatio-temporal [n or 1, k, 2]
# centre: index of the central location (see .coVarCond)
dCovariateVine <- function(u, copula, log=FALSE, h, centre) {
  nNeighs <- ncol(u) - 2

  tree <- .condTree(u[, 1:(nNeighs+1), drop=FALSE], h, copula@tree)
  coVar <- .coVarCond(u[, c(1, nNeighs+2), drop=FALSE], copula@coVarCop, centre,
                      .isStTree(copula@tree))
  l1 <- .topLogDens(cbind(coVar$u, tree$u), copula@topCop)

  res <- tree$loglik + coVar$loglik + l1
  if(log)
    return(res)
  exp(res)
}

setMethod("dCopula", signature=signature("matrix","covariate_vine_copula"),
          function(u, copula, log=FALSE, ...) dCovariateVine(u, copula, log=log, ...))
setMethod("dCopula", signature=signature("numeric","covariate_vine_copula"),
          function(u, copula, log=FALSE, ...)
            dCovariateVine(matrix(u, ncol=copula@dimension), copula, log=log, ...))
setMethod("dCopula", signature=signature("data.frame","covariate_vine_copula"),
          function(u, copula, log=FALSE, ...)
            dCovariateVine(as.matrix(u), copula, log=log, ...))

## fitting the top copula for given tree and covariate copula ##
#################################################################

fitCovariateVine <- function(copula, data,
                             method = list(StructureSelect = FALSE,
                                           indeptest = FALSE,
                                           familyset = NA),
                             estimate.variance=FALSE) {
  neigh <- data
  stopifnot(is(neigh, "neighbourhood"), length(neigh@coVar) == 1)
  if (.isStNeighbourhood(neigh) != .isStTree(copula@tree))
    stop("Spatio-temporal neighbourhoods need a spacetime_copula tree, spatial ones a spatial_copula tree.")
  u <- as.matrix(neigh@data)
  stopifnot(copula@dimension == ncol(u))
  nNeighs <- ncol(u) - 2

  cat("[Dropping the distance tree.]\n")
  tree <- .condTree(u[, 1:(nNeighs+1), drop=FALSE], neigh@distances, copula@tree)
  cat("[Conditioning the covariate.]\n")
  coVar <- .coVarCond(u[, c(1, nNeighs+2), drop=FALSE], copula@coVarCop, .centreIndex(neigh),
                      .isStNeighbourhood(neigh))

  u1 <- cbind(coVar$u, tree$u)
  if (ncol(u1) == 2) {
    cat("[Estimating a single bivariate copula at the top.]\n")
    bivCop <- BiCopSelect(u1[,1], u1[,2])
    topCop <- copulaFromFamilyIndex(bivCop$family, bivCop$par, bivCop$par2)
    loglik <- sum(dCopula(u1, topCop, log=TRUE))
  } else {
    cat("[Estimating a",ncol(u1),"dimensional copula at the top.]\n")
    topFit <- fitCopula(copula@topCop, u1, method)
    topCop <- topFit@copula
    loglik <- topFit@loglik
  }

  cvCop <- covariate_vine_copula(copula@coVarCop, copula@tree, topCop)
  new("fitCopula", estimate = cvCop@parameters, var.est = matrix(NA),
      method = paste(sapply(method, paste, collapse=", "), collapse="; "),
      loglik = sum(tree$loglik) + sum(coVar$loglik) + loglik,
      fitting.stats=list(convergence = as.integer(NA)),
      nsample = nrow(u), copula=cvCop)
}

setMethod("fitCopula", signature=signature("covariate_vine_copula"), fitCovariateVine)

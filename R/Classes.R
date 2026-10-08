####
## an empirical survival copula representation

validEmpSurCopula <- function(object) {
  if(ncol(object@sample) != object@dimension)
    return("Dimension of the copula and the sample do not match.")
  else
    return(TRUE)
}

setClass("empSurCopula",
         representation = representation("copula", sample="matrix"),
         validity = validEmpSurCopula,
         contains = list("copula")
)

####
## the leaf copula

validLeafCopula <- function(object) {
  if (object@dimension != 2)
    return("The leaf copula only supports two dimensions.")
  
  if (any(is.na(object@parameters)))
    return("Parameter value is \"NA\".")
  else return (TRUE)
}

setClass("leafCopula",
         representation = representation("copula"),
         validity = validLeafCopula,
         contains = list("copula")
)

#####
## Hierarchical Kendall Copulas

validHkCopula <- function(object) {
  stopifnot(all(sapply(object@clusterCops, function(x) inherits(x[[1]], "copula"))))
  if (object@dimension != sum(sapply(object@clusterCops, function(x) x[[1]]@dimension))+object@nestingCop@dimension-length(object@clusterCops))
    return("The dimensions of the hierarchical Kendall copula do not match.")
  
  if(length(object@clusterCops) != length(object@kenFuns))
    return("Each cluster copula needs to have its Kendall function in 'kenFuns'.")
  
  else return (TRUE)
}

setClass("hkCopula",
         representation = representation("copula",
                                         nestingCop = "copula",
                                         clusterCops = "list",
                                         kenFuns = "list"),
         validity = validHkCopula,
         contains = list("copula")
)




## 
## the spatial copula
##
## realized as a distance dependent convex combination of biv copulas

# dimension = "numeric"     set to 2
# parameters = "numeric"    set of parameters
# param.names = "character" appropriate names
# param.lowbnd = "numeric"  appropriate lower bounds
# param.upbnd = "numeric"   appropriate upper bounds
# fullname = "character"    name printed with "show"
# components="list"         list of copulas 
# distances="numeric"       the linking distances
# unit="character"          measurement unit of distance
# depFun="function"         an optional dependence function; depFun(NULL)
#                             has to return either "spearman" or "kendall" 
#                             dependening on the moa used. Make sure depFun
#                             assings valid parameters to the copulas involved

validSpCopula <- function(object) {
  if (length(object@components) != length(object@distances)) 
    return("Length of components does not equal length of distances. \n Note: The last distance is interpreted as the range and every pair beyond the range is modelled via the last copula in components.")
  if (is.unsorted(object@distances, strictly = T))
    return("Distances have to be strictly increasing.")
  if (length(object@combination) != 1 || !object@combination %in% c("convex", "geometric"))
    return("The combination needs to be either 'convex' or 'geometric'.")
  if (object@combination == "geometric" && min(object@distances) > 0)
    return("The distance vector of a geometric spatial copula must contain 0.")
  
  check.upper <- NULL
  check.lower <- NULL
  
  nComp <- length(object@components)
  if(!is.null(object@calibMoa(normalCopula(0),0))) {
    nonIndep <- sapply(object@components[-nComp], function(x) class(x) != "indepCopula")
    for (i in (1:(nComp-1))[nonIndep]) {
      upParam <- object@calibMoa(object@components[[i]], object@distances[i+1])
      if(any(is.na(upParam))) {
        check.upper <- c(check.upper, TRUE)
      } else {
        if (class(object@components[[i]]) == "frankCopula" && upParam == 0) {
          check.upper <- c(check.upper, TRUE)
        } else {
          check.upper <- c(check.upper, FALSE)
        }
      }
        
      check.lower <- c(check.lower, is.na(object@calibMoa(object@components[[i]], c(0,object@distances)[i])))
    }
    if(sum(check.upper>0)) return(paste("Reconsider the upper boundary conditions of the following copula(s): \n",
                                        paste(sapply(object@components[check.upper], function(x) describeCop(x, "very short")), 
                                              "at", object@distances[check.upper],collapse="\n")))
    if(sum(check.lower>0)) return(paste("Reconsider the lower boundary conditions of the following copula(s): \n",
                                        paste(sapply(object@components[check.lower], function(x) describeCop(x, "very short")), 
                                              "at", object@distances[check.lower],collapse="\n")))
  }
  
  return(TRUE)
}

setClass("spatial_copula", representation = representation("copula", 
                                                     components="list",
                                                     distances="numeric", 
                                                     calibMoa="function", 
                                                     unit="character",
                                                     combination="character"),
         validity = validSpCopula, contains = list("copula"))

############################
## Spatio-Temporal Copula ##
############################

validStCopula <- function(object) {
  if(length(object@tlags) != length(object@spCopList)) return("The length of the temporal distance vector must equal the number of spatial copulas.")
  return(TRUE) # validity of any spatial_copula in spCopList is tested by the constructor, I believe
}

setClass("spacetime_copula", representation = representation("copula", 
                                                     spCopList="list", 
                                                     tlags="numeric",
                                                     tres="character"),
         validity = validStCopula, contains = list("copula"))

########################################################
## Distance Vine Copula (spatial and spatio-temporal) ##
########################################################

setClassUnion("optionalCopula", c("copula", "NULL"))

validDistanceVineCopula <- function(object) {
  if (length(object@trees) == 0)
    return("At least one distance dependent tree is needed.")
  isSp <- sapply(object@trees, is, "spatial_copula")
  isSt <- sapply(object@trees, is, "spacetime_copula")
  if (!(all(isSp) || all(isSt)))
    return("All trees need to be either spatial_copula or spacetime_copula objects.")
  if (!is.null(object@topCop))
    return(validObject(object@topCop))
  return(TRUE)
}

setClass("distance_vine_copula", representation("copula", trees="list", 
                                                topCop="optionalCopula"),
         validity = validDistanceVineCopula, contains=list("copula"))

#####################################
## neighbourhoods providing the data ##
#####################################

## spatial: distances and index are matrices
## spatio-temporal: 3-dimensional arrays with space and time in the 3rd dimension

validNeighbourhood <- function(object) {
  if(length(object@var)>1)
    return("Only a single variable name is supported.")
  dimDists <- dim(object@distances)
  dimInd <- dim(object@index)
  if (length(dimDists) != length(dimInd))
    return("Distances and index need to have the same number of dimensions.")
  
  if (nrow(object@data) != dimDists[1]) 
    return("Data and distances have unequal number of rows.")
  if (nrow(object@data) != dimInd[1]) 
    return("Data and index have unequal number of rows.")
  
  if (length(dimDists) == 2) {
    if (ncol(object@data) != dimDists[2] + 1 + length(object@coVar))
      return("Data and distances have non matching number of columns.")
    if (ncol(object@data) != dimInd[2] + length(object@coVar)) 
      return("Data and index have non matching number of columns.")
    return(TRUE)
  }
  
  if (length(dimDists) != 3 || dimDists[3] != 2 || dimInd[3] != 2)
    return("Spatio-temporal distances and index need to be arrays with 2 layers (space, time).")
  if (ncol(object@data) + object@prediction != dimInd[2] + length(object@coVar)) 
    return("Data and index have non matching number of columns.")
  if (dimDists[2]+1 != dimInd[2]) 
    return("Data and index have non matching number of columns.")
  return(TRUE)
}

setClass("neighbourhood",
         representation = representation(data = "data.frame", 
                                         distances="array", 
                                         index="array",
                                         var="character",
                                         coVar="character",
                                         prediction="logical"),         
         validity = validNeighbourhood)

#############################
## spatial Gaussian copula ##
#############################

setClass("spatial_gauss_copula", representation = representation(corFun = "function"))

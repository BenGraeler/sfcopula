## distance_vine_copula estimation following the vignette, but using a randomly 
## selected smaller subset of the original data to reduce calculation demands.
## Thus, results are likely to differ (a little) from the original study.

library("VineCopula") # load before sfcopula to not mask VC2copula's copula classes
library("sf")
library("stars")
library("sfcopula")
data("EU_RB_2005", package = "sfcopula")

## rank transformed PM10 per day (EU_RB_2005 is a stars vector data cube)
EU_RB_2005$rtPM10 <- apply(EU_RB_2005[["PM10"]], 2, 
                           function(x) rank(x, na.last = "keep") / (sum(!is.na(x)) + 1))

## spatio-temporal copula
# binning, using only 10 randomly selected temporal instances
set.seed(2005)
stBins <- calc_bins(EU_RB_2005, "rtPM10", nbins=40, tlags=-(0:2),
                   instances=10)

stDepFun <- fit_cor_fun(stBins, c(3,3,3), tlags=-(0:2))

# parameters are re-calibrated by stDepFun; avoid boundary values that
# would turn Clayton/Gumbel into independence copulas
families <- c(normalCopula(0.5), tCopula(0.5), claytonCopula(1), 
              frankCopula(1), gumbelCopula(1.5), joeBiCopula(1.5))

loglikTau <- loglik_by_lags(stBins, EU_RB_2005, families, stDepFun)

bestFitTau <- lapply(loglikTau, 
                     function(x) apply(apply(x$loglik, 1, rank), 2, which.max))

# gather the fitted copulas and representative distances
listCops <- lapply(bestFitTau, function(x) families[x[1:35]])
listDists <- rep(list(stBins$meanDists[1:35]), 3)

# build the spatio-temporal copula
stConvCop <- spacetime_copula(components=listCops, distances=listDists,
                      tlags=-(0:2), stDepFun=stDepFun)

# get the neighbours
stNeigh <- neighbours(EU_RB_2005, var="rtPM10", size=4, 
                           tlags=-(0:2), timeSteps=10, min.dist=10)
# drop neighbourhoods with missing values
stNeigh <- stNeigh[which(complete.cases(stNeigh@data))]

stVineFit <- fitCopula(distance_vine_copula(stConvCop, vineCopula(9L)), stNeigh,
                       method=list(indeptest=TRUE))
stVine <- stVineFit@copula

# retrieve the log-likelihood
stVineFit@loglik

## predict the median for two stations on three days (cross-validation style)
target <- EU_RB_2005["rtPM10", 1:2, 100:102]
predNeigh <- neighbours(EU_RB_2005, target, var="rtPM10", size=4,
                             tlags=-(0:2), prediction=TRUE, min.dist=10)
pred <- predict(stVine, predNeigh, EU_RB_2005, target,
                     margin=list(q=function(p) p), method="quantile")
cbind(observed=as.vector(t(target[["rtPM10"]])),
      predicted=as.vector(t(pred[["quantile.0.5"]])))

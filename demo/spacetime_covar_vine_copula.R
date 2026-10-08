######################################################################
# demo related to a paper (in preparation for JSS)
# Different than the study presented in the above paper, only a tempo-
# ral subset of the European air quality data is used and the set of
# copula family candidates is limited. These chnages have been neces-
# sary to maintain the "runability" of this demo.
######################################################################
library("VineCopula") # load before sfcopula to not mask VC2copula's copula classes
library("sf")
library("stars")
library("sfcopula")
library("evd")
data("EU_RB", package = "sfcopula")

# EU_RB is a stars vector data cube: 194 stations x 61 days
PM10 <- EU_RB[["PM10"]]
EMEP <- EU_RB[["EMEP"]]
nLoc <- nrow(PM10)

# estimate a GEV at each location for PM10 and EMEP
parPM10 <- matrix(NA, nLoc, 3)
parEMEP <- matrix(NA, nLoc, 3)

marPM10 <- matrix(NA, nLoc, 61)
marEMEP <- matrix(NA, nLoc, 61)

for (loc in 1:nLoc) {
  parPM10[loc, 1:3] <- fgev(PM10[loc, ])$estimate
  parEMEP[loc, 1:3] <- fgev(EMEP[loc, ])$estimate
  
  marPM10[loc,] <- pgev(PM10[loc, ], parPM10[loc,1], parPM10[loc,2], parPM10[loc,3])
  marEMEP[loc,] <- pgev(EMEP[loc, ], parEMEP[loc,1], parEMEP[loc,2], parEMEP[loc,3])
}

EU_RB$marPM10 <- marPM10
EU_RB$marEMEP <- marEMEP

########################################
## correlation between EMEP and PM10? ##
########################################

dayCor <- numeric(61)
for(day in 1:61) {
  smpl <- cbind(marPM10[, day], marEMEP[, day])
  bool <- !apply(smpl,1,function(row) any(is.na(row)))
  smpl <- smpl[bool,]
  
  dayCor[day] <- TauMatrix(smpl)[1,2]
}

weekCor <- numeric(9)
weekCop <- NULL
for(week in 1:9) {
  days <- pmin((week-1)*7+1:7,61)
  smpl <- cbind(as.vector(marPM10[, days]), as.vector(marEMEP[, days]))
  bool <- !apply(smpl,1,function(row) any(is.na(row)))
  smpl <- smpl[bool,]
  
  weekCor[week] <- TauMatrix(smpl)[1,2]
  weekCop <- append(weekCop,list(BiCopSelect(smpl[,1], smpl[,2], familyset=1:6)))
}


par(mar=c(5.1, 4.1, 4.1,6))
plot(dayCor, type="l", col="gray", xlab="day in 2005-06-01::2005-07-31", 
     ylab="Kendall's tau", main="correlation structure of PM10 and EMEP over time")
points(rep(weekCor,each=7), type="s", col="red")
segments(0:8*7+1,sapply(weekCop, function(x) x$family)/10, 1:9*7+1, sapply(weekCop, function(x) x$family)/10)
axis(4,at=1:6/10, labels=c("Gauss", "Student", "Clayton", "Gumbel", "Frank", "Joe"),las=2)
mtext("copula family",4,4.5)

#############################
# the paper starts here ... #
#############################

# define the coVariate Copula function
coVarCop <- function(stInd) {
  week <- min(ceiling(stInd[2]/7), 9)
  copulaFromFamilyIndex(weekCop[[week]]$family, weekCop[[week]]$par, 
                        weekCop[[week]]$par2)
}

## spatio-temporal copula
# binning
stBins <- calc_bins(EU_RB, "marPM10", nbins=20, tlags=-(0:2))
stDepFun <- fit_cor_fun(stBins, rep(3, 5), tlags=-(0:4))


## 
fiveColors <- c("#fed976", "#feb24c", "#fd8d3c", "#f03b20", "#bd0026")
par(mar=c(4.1, 4.1, 2.1, 1.1))
plot(stBins$meanDists/1000, stBins$lagCor[1,], 
     ylim=c(0,0.7), xlab="distance [km]", ylab="correlation [Kendall's tau]",
     col=fiveColors[5])
points(stBins$meanDists/1000, stBins$lagCor[2,], col=fiveColors[3])
points(stBins$meanDists/1000, stBins$lagCor[3,], col=fiveColors[1])
abline(h=0)
abline(h=0.025,col="grey")

fun1 <- function(x) stDepFun(x*1000, 1, 5:1)
curve(fun1, 0, 1600, add=T, col=fiveColors[5])
fun2 <- function(x) stDepFun(x*1000, 2, 5:1)
curve(fun2, 0, 1600, add=T, col=fiveColors[3])
fun3 <- function(x) stDepFun(x*1000, 3, 5:1)
curve(fun3, 0, 1600, add=T, col=fiveColors[1])

legend("topright",c("same day", "1 day before", "2 days before"),
       lty=1, pch=1, col=fiveColors[c(5,3,1)])
title("Spatio-Temporal Dependence Structure")
##

families <- c(normalCopula(), tCopula(),
              claytonCopula(), frankCopula(), gumbelCopula(), 
              joeBiCopula())

loglikTau <- loglik_by_lags(stBins, EU_RB, families, stDepFun)

bestFitTau <- lapply(loglikTau, 
                     function(x) apply(apply(x$loglik[,1:6], 1, rank),
                                       2, which.max))

bestFitTau

# define the spatio-temporal copula components
listDists <- NULL
listDists[[1]] <- stBins$meanDists[sort(unique(c(which(diff(bestFitTau$loglik1)!=0),
                                                 which(diff(bestFitTau$loglik1)!=0)+1,1,20)))]
listDists[[2]] <- stBins$meanDists[sort(unique(c(which(diff(bestFitTau$loglik2)!=0),
                                                 which(diff(bestFitTau$loglik2)!=0)+1,1,20)))]
listDists[[3]] <- stBins$meanDists[sort(unique(c(which(diff(bestFitTau$loglik3)!=0),
                                                 which(diff(bestFitTau$loglik3)!=0)+1,1,20)))]

listCops <- NULL
listCops[[1]] <- families[bestFitTau$loglik1[sort(unique(c(which(diff(bestFitTau$loglik1)!=0),
                                                           which(diff(bestFitTau$loglik1)!=0)+1,1,20)))]]
listCops[[2]] <- families[bestFitTau$loglik2[sort(unique(c(which(diff(bestFitTau$loglik2)!=0),
                                                           which(diff(bestFitTau$loglik2)!=0)+1,1,20)))]]
listCops[[3]] <- families[bestFitTau$loglik3[sort(unique(c(which(diff(bestFitTau$loglik3)!=0),
                                                           which(diff(bestFitTau$loglik3)!=0)+1,1,20)))]]

stBiCop <- spacetime_copula(components = listCops, distances = listDists, 
                    tlags=-c(0:2), stDepFun=stDepFun)


## get the neighbours
stNeigh <- neighbours(EU_RB, size=9, var="marPM10", coVar="marEMEP",
                           tlags=-(0:2), timeSteps=20, min.dist=10)
stRedNeigh <- reduce_neighbours(stNeigh, stDepFun, 5)

# condition on the spatio-temporal tree
condData <- drop_tree(stRedNeigh, EU_RB, stBiCop)

# condition the covariate on the observed phenomenon
condCoVa <- cond_covariate(stRedNeigh, coVarCop)

secTreeData <- cbind(condCoVa, as.matrix(condData@data))

vineFit <- fitCopula(vineCopula(6L), secTreeData, method=list(familyset=1:6))

stCVVC <- covariate_vine_copula(coVarCop, stBiCop, vineFit@copula)

stCVVC
## the same in one call: fitCopula drops the spatio-temporal tree, conditions
## the covariate and fits the top vine (the log-likelihood covers all parts)
stCVVCFit <- fitCopula(covariate_vine_copula(coVarCop, stBiCop, vineCopula(6L)), stRedNeigh,
                       method=list(familyset=1:6))
stCVVCFit@loglik

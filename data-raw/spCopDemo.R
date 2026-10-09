## Generates the workspace data/spCopDemo.RData following demo("spatial_copula").
## Run from the package root with sfcopula installed.
library("sf")
library("sfcopula")

data("meuse", package = "sfcopula")

meanLog <- mean(log(meuse$zinc))
sdLog <- sd(log(meuse$zinc))
meuse$marZinc <- plnorm(meuse$zinc, meanLog, sdLog)

bins <- calc_bins(meuse, var = "marZinc", nbins = 10, cutoff = 800, plot = FALSE)
calcKTauPol <- fit_cor_fun(bins, degree = 3)
environment(calcKTauPol) <- list2env(as.list(environment(calcKTauPol)),
                                     parent = globalenv())

copCandidates <- c(normalCopula(0.2), tCopula(0.2), claytonCopula(0.2),
                   frankCopula(1.2), gumbelCopula(1.2), joeBiCopula(1.5),
                   indepCopula())

loglikTau <- loglik_by_lags(bins, meuse, copCandidates, calcKTauPol)
bestFitTau <- apply(loglikTau$loglik, 1, which.max)

spCop <- spatial_copula(components = c(copCandidates[bestFitTau[1]], copCandidates[bestFitTau]),
                  distances = c(0, bins$meanDists),
                  dep_fun = calcKTauPol, unit = "m")

meuseNeigh <- neighbours(meuse, var = "marZinc", size = 5L)
meuseSpVine <- fitCopula(distance_vine_copula(spCop, vineCopula(4L)),
                         list(meuseNeigh, meuse))@copula

save(bins, bestFitTau, calcKTauPol, spCop, meuseSpVine,
     file = "data/spCopDemo.RData", compress = "xz")

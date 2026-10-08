## Creates tests/testthat/fixtures/spcopula_objects.rds: objects built with the
## original spcopula 0.2-5 (sp/spacetime) plus reference values computed by it.
## Needs spcopula, sp and spacetime installed.
suppressMessages({library(sp); library(spacetime); library(spcopula)})

data("meuse", package = "sp"); coordinates(meuse) <- ~x+y
meuse$marZinc <- plnorm(meuse$zinc, mean(log(meuse$zinc)), sd(log(meuse$zinc)))

spCop <- spCopula(components = list(claytonCopula(4), gumbelCopula(2), indepCopula()),
                  distances = c(0, 200, 500), unit = "m")
bins <- calcBins(meuse, var = "marZinc", nbins = 10, cutoff = 800, plot = FALSE)
spCopDep <- spCopula(components = list(normalCopula(0.5), claytonCopula(1), indepCopula()),
                     distances = c(0, 300, 700), spDepFun = fitCorFun(bins, degree = 3), unit = "m")
stCop <- stCopula(components = list(spCop, spCopula(list(claytonCopula(2), indepCopula()),
                                                    distances = c(0, 400), unit = "m")),
                  tlags = c(0, -1))
neigh <- getNeighbours(meuse[1:30, ], var = "marZinc", size = 4)
fit <- fitCopula(spVineCopula(spCop, vineCopula(3L)), list(neigh, meuse[1:30, ]))
spVine <- fit@copula
pureVine <- spVineCopula(list(spCop, spCop))
stVine <- stVineCopula(stCop, vineCopula(3L))
cqs <- cqsCopula(c(0.2, 0.1))
geom <- spGeomCopula(list(claytonCopula(4), gumbelCopula(2), indepCopula()), distances = c(0, 200, 500), unit = "m")

coVarCop <- function(stInd) normalCopula(0.3 + 0.1 * (stInd[2] %% 3))
covTop <- vineCopula(VineCopula::RVineMatrix(Matrix = matrix(c(3,1,2, 0,2,1, 0,0,1), 3),
                                             family = matrix(c(0,1,3, 0,0,4, 0,0,0), 3),
                                             par = matrix(c(0,0.4,1.2, 0,0,1.5, 0,0,0), 3)))
stCoVar <- stCoVarVineCopula(coVarCop, stCopula(list(spCopula(list(claytonCopula(2), indepCopula()), c(0, 400), unit = "m"),
                                                     spCopula(list(claytonCopula(1), indepCopula()), c(0, 400), unit = "m")),
                                                tlags = -(0:1)), covTop)
coVarDists <- array(c(100, 250, 0, -1), c(1, 2, 2))
coVarX <- c(0.05, 0.3, 0.5, 0.8, 0.97)

u <- cbind(c(0.2, 0.5, 0.9), c(0.3, 0.6, 0.8))
h <- c(50, 250, 600)
u3 <- matrix(c(0.2, 0.5, 0.9, 0.3, 0.6, 0.8, 0.4, 0.5, 0.7), 3)
expected <- list(
  spCop = dCopula(u, spCop, h = h),
  spCopDep = dCopula(u, spCopDep, h = h),
  stCop = dCopula(u, stCop, h = cbind(h, c(0, -1, 0))),
  spVine = dCopula(cbind(u, u), spVine, h = list(matrix(c(50, 100, 150), 1)[rep(1, 3), ])),
  pureVine = dCopula(u3, pureVine, h = list(matrix(c(50, 100), 1)[rep(1, 3), ], 120)),
  cqs = dCopula(u, cqs),
  geom = dCopula(u, geom, h = h),
  # covariate first in spcopula's conditional density
  stCoVar = condStCoVarVine(c(0.7, 0.4, 0.6), coVarDists, stCoVar, stInd = c(2, 5))(coVarX))

objects <- list(spCop = spCop, spCopDep = spCopDep, stCop = stCop, neigh = neigh, fit = fit,
                spVine = spVine, pureVine = pureVine, stVine = stVine, cqs = cqs, geom = geom,
                stCoVar = stCoVar,
                nested = list(a = spCop, b = list(neigh)))
saveRDS(list(objects = objects, expected = expected, u = u, h = h, u3 = u3,
             coVarDists = coVarDists, coVarX = coVarX),
        "tests/testthat/fixtures/spcopula_objects.rds")

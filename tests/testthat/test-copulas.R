## ported from spcopula's tests: the copula core is unaffected by sf/stars

test_that("inverse partial derivatives of the spatial copula invert", {
  data("spCopDemo", package = "sfcopula", envir = environment())
  spCop <- spatial_copula(components = list(normalCopula(), tCopula(),
                                      frankCopula(), normalCopula(), claytonCopula(),
                                      claytonCopula(), claytonCopula(), claytonCopula(),
                                      claytonCopula(), indepCopula()),
                    distances = c(0, bins$meanDists[1:9]),
                    spDepFun = calcKTauPol, unit = "m")

  y <- dduCopula(c(0.3, 0.7), spCop, h = 300)
  expect_equal(invdduCopula(0.3, spCop, y, h = 300), 0.7, tolerance = 1e-4)

  y <- ddvCopula(c(0.3, 0.7), spCop, h = 300)
  expect_equal(invddvCopula(0.7, spCop, y, h = 300), 0.3, tolerance = 1e-4)
})

test_that("inverse partial derivatives of the spatio-temporal copula invert", {
  spCopT0 <- spatial_copula(components = list(claytonCopula(8), claytonCopula(4),
                                        claytonCopula(2), claytonCopula(1),
                                        claytonCopula(0.5), indepCopula()),
                      distances = c(0, 100, 200, 300, 400, 500), unit = "km")
  spCopT1 <- spatial_copula(components = list(claytonCopula(4), claytonCopula(2),
                                        claytonCopula(1), claytonCopula(0.5),
                                        indepCopula()),
                      distances = c(0, 100, 200, 300, 400), unit = "km")
  spCopT2 <- spatial_copula(components = list(claytonCopula(2), claytonCopula(1),
                                        claytonCopula(0.5), indepCopula()),
                      distances = c(0, 100, 200, 300), unit = "km")
  stCop <- spacetime_copula(components = list(spCopT0, spCopT1, spCopT2), tlags = -(0:2))
  h <- matrix(c(150, -1), ncol = 2)

  y <- dduCopula(c(0.3, 0.7), stCop, h = h)
  expect_equal(invdduCopula(0.3, stCop, y, h = h), 0.7, tolerance = 1e-4)

  y <- ddvCopula(c(0.3, 0.7), stCop, h = h)
  expect_equal(invddvCopula(0.7, stCop, y, h = h), 0.3, tolerance = 1e-4)
})

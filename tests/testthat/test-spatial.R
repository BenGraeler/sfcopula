## reference values have been computed with spcopula 0.2-5 (sp based)

test_that("neighbourhoods match spcopula", {
  meuse <- meuseZinc()
  neigh <- neighbours(meuse, var = "marZinc", size = 5)

  expect_s4_class(neigh, "neighbourhood")
  expect_equal(dim(neigh@data), c(155, 5))
  expect_equal(neigh@index[1:3, ], matrix(c(1, 2, 3, 2, 1, 1, 3, 3, 2, 8, 8, 4, 7, 7, 7), 3))
  expect_equal(neigh@distances[1:2, 1:2],
               matrix(c(70.8378429936994, 70.8378429936994,
                        118.848643240047, 141.566238913097), 2))
  expect_equal(neigh@data[, 1], meuse$marZinc, ignore_attr = TRUE)
})

test_that("prediction neighbourhoods exclude the target and accept sfc", {
  meuse <- meuseZinc()
  pn <- neighbours(meuse[1:15, ], st_geometry(meuse[16:20, ]), var = "marZinc",
                      size = 3, prediction = TRUE, min.dist = 10)
  expect_true(all(is.na(pn@data[, 1])))
  expect_equal(pn@index[, 1], 1:5)
  expect_true(all(pn@index[, -1] <= 15))
})

test_that("spatial bins match spcopula", {
  meuse <- meuseZinc()
  bins <- calc_bins(meuse, var = "marZinc", nbins = 10, cutoff = 800, plot = FALSE)
  expect_equal(bins$np, c(25, 170, 264, 318, 349, 394, 411, 389, 426, 448))
  expect_equal(bins$lagCor[1:3], c(0.56807444771205895, 0.49449012777111551,
                                   0.31777333447200512))
})

test_that("spatial vine fit and prediction match spcopula", {
  meuse <- meuseZinc()
  qMar <- function(p) qlnorm(p, mean(log(meuse$zinc)), sd(log(meuse$zinc)))
  bins <- calc_bins(meuse, var = "marZinc", nbins = 10, cutoff = 800, plot = FALSE)
  calcKTauPol <- fit_cor_fun(bins, degree = 3)
  fams <- c(normalCopula(0.2), tCopula(0.2), claytonCopula(0.2), frankCopula(1.2),
            gumbelCopula(1.2), joeBiCopula(1.5), indepCopula())
  llTau <- quiet(loglik_by_lags(bins, meuse, fams, calcKTauPol))
  best <- apply(llTau$loglik, 1, which.max)
  spCop <- quiet(spatial_copula(components = c(fams[best[1]], fams[best]),
                          distances = c(0, bins$meanDists),
                          spDepFun = calcKTauPol, unit = "m"))

  neigh <- neighbours(meuse, var = "marZinc", size = 5)
  fit <- quiet(fitCopula(distance_vine_copula(spCop, vineCopula(4L)), list(neigh, meuse)))
  expect_equal(fit@loglik, 208.017820595265, tolerance = 1e-6)

  spVine <- distance_vine_copula(list(spCop, spCop))
  pn <- neighbours(meuse[1:15, ], meuse[16:20, ], var = "marZinc",
                      size = 3, prediction = TRUE, min.dist = 10)
  pred <- quiet(predict(spVine, pn, meuse[1:15, ], meuse[16:20, ],
                             list(q = qMar), "quantile"))
  expect_s3_class(pred, "sf")
  expect_equal(pred$quantile.0.5[1:2], c(627.993751330458, 403.470809122627),
               tolerance = 1e-6)

  predGeom <- quiet(predict(spVine, pn, meuse[1:15, ], st_geometry(meuse[16:20, ]), list(q = qMar), "quantile"))
  expect_s3_class(predGeom, "sf")
  expect_equal(predGeom$quantile.0.5, pred$quantile.0.5)
})

test_that("spatial Gaussian copula matches spcopula", {
  meuse <- meuseZinc()
  neigh <- neighbours(meuse, var = "marZinc", size = 5)
  ll <- quiet(spatial_gauss_loglik(function(h) exp(-h / 400), neigh[1:20], meuse))
  expect_equal(ll, 33.8726773626673, tolerance = 1e-8)
})

test_that("geographic coordinates use great circle distances in metres", {
  meuse <- meuseZinc()
  ll <- st_transform(meuse, 4326)
  n1 <- neighbours(meuse, var = "marZinc", size = 3)
  n2 <- neighbours(ll, var = "marZinc", size = 3)
  expect_equal(n2@index, n1@index)
  expect_equal(n2@distances, n1@distances, tolerance = 0.01)
})

test_that("non-point geometries are rejected", {
  poly <- st_sf(a = 1, geometry = st_sfc(st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 0))))))
  expect_error(neighbours(poly, var = "a", size = 2), "POINT")
})

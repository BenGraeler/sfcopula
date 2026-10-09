## the unified classes and functions

clayGumb <- function(...) list(claytonCopula(4), gumbelCopula(2), indepCopula())

test_that("spatial_copula covers convex and geometric combinations", {
  convex <- spatial_copula(clayGumb(), distances = c(0, 200, 500))
  geom <- spatial_copula(clayGumb(), distances = c(0, 200, 500), combination = "geometric")
  expect_equal(convex@unit, "m")
  expect_equal(geom@combination, "geometric")
  expect_output(show(convex), "convex combination")
  expect_output(show(geom), "geometric mean")

  # spcopula's spGeomCopula used swapped weights in its density; they agree only
  # where both weights are equal (beyond the range)
  old <- suppressWarnings(readRDS(test_path("fixtures", "spcopula_objects.rds")))
  expect_equal(dCopula(old$u, geom, h = old$h)[3], old$expected$geom[3])
  expect_false(isTRUE(all.equal(dCopula(old$u, geom, h = old$h)[1:2], old$expected$geom[1:2])))
  v <- quiet(invdduCopula(0.3, geom, y = dduCopula(c(0.3, 0.7), geom, h = 300), h = 300))
  expect_equal(v, 0.7, tolerance = 1e-4)

  expect_error(spatial_copula(clayGumb(), distances = c(10, 200, 500), combination = "geometric"),
               "must contain 0")
  expect_error(spatial_copula(clayGumb(), distances = c(0, 200, 500), dep_fun = function(h) "kendall",
                              combination = "geometric"), "convex")
})

test_that("neighbours dispatches on sf and stars", {
  meuse <- meuseZinc()
  sp <- neighbours(meuse, var = "marZinc", size = 3)
  expect_s4_class(sp, "neighbourhood")
  expect_false(sp@prediction)
  expect_equal(dim(sp@distances), c(155, 2))

  pn <- neighbours(meuse[1:10, ], meuse[11:12, ], var = "marZinc", size = 3)
  expect_true(pn@prediction) # implied by the target

  st <- neighbours(toyCube(6), var = "u", size = 3, tlags = -(0:1))
  expect_s4_class(st, "neighbourhood")
  expect_equal(dim(st@distances), c(25, 4, 2))
  expect_output(show(st), "spatio-temporal")
  expect_equal(dim(st[1:3]@distances), c(3, 4, 2))

  expect_error(calc_bins(st, plot = FALSE), "spatial neighbourhoods")
  expect_error(reduce_neighbours(sp, function(h, t) 1, 1), "spatio-temporal")
})

test_that("drop_tree passes distances between neighbours to the next tree", {
  cube <- toyCube(6)
  st <- neighbours(cube, var = "u", size = 4, tlags = -(0:1))
  stCop <- spacetime_copula(list(quiet(spatial_copula(list(claytonCopula(2), indepCopula()), c(0, 300))),
                                 quiet(spatial_copula(list(claytonCopula(1), indepCopula()), c(0, 300)))),
                            tlags = -(0:1))
  dropped <- drop_tree(st, cube, stCop)
  geom <- st_get_dimension_values(cube, "geometry")
  ind <- st@index
  for (i in 1:dim(dropped@distances)[2]) {
    expect_equal(dropped@distances[, i, 1],
                 as.numeric(st_distance(geom[ind[, 2, 1]], geom[ind[, 2 + i, 1]], by_element = TRUE)))
    expect_equal(dropped@distances[, i, 2], ind[, 2 + i, 2] - ind[, 2, 2])
  }
  expect_equal(dropped@index, ind[, -1, ])

  meuse <- meuseZinc()
  sp <- neighbours(meuse, var = "marZinc", size = 4)
  spDropped <- drop_tree(sp, meuse, spatial_copula(clayGumb(), c(0, 200, 500)))
  expect_equal(spDropped@distances, tree_dists(sp, meuse, 2)[[2]])
  expect_error(drop_tree(sp, meuse, stCop), "spatial_copula")
})

test_that("distance vines fit and simulate with several trees", {
  meuse <- meuseZinc()
  tree1 <- spatial_copula(list(claytonCopula(3), indepCopula()), c(0, 800))
  tree2 <- spatial_copula(list(gumbelCopula(1.5), indepCopula()), c(0, 800))

  # pure vine: nothing left for the top
  pure <- distance_vine_copula(list(tree1, tree2))
  expect_null(pure@topCop)
  expect_equal(pure@dimension, 3L)
  fit <- quiet(fitCopula(pure, list(neighbours(meuse, var = "marZinc", size = 3), meuse)))
  expect_true(is.finite(fit@loglik))

  # a single bivariate copula on top
  fit2 <- quiet(fitCopula(distance_vine_copula(tree1, normalCopula(0)),
                          neighbours(meuse, var = "marZinc", size = 3)))
  expect_s4_class(fit2@copula@topCop, "copula")
  expect_equal(fit2@copula@dimension, 3L)

  # simulation uses the copula of each tree (tau of Clayton(3) at 50 m: 0.6 * 750/800)
  set.seed(11)
  sims <- rCopula(300, pure, h = list(c(50, 50), 50))
  expect_equal(dim(sims), c(300, 3))
  expect_equal(cor(sims[, 1], sims[, 2], method = "kendall"),
               tau(claytonCopula(3)) * 750 / 800, tolerance = 0.1)
})

test_that("spatio-temporal vines with two trees fit and predict", {
  cube <- toyCube(8)
  set.seed(3)
  cube$u <- matrix(runif(40), 5, 8)
  spC <- function(p) quiet(spatial_copula(list(claytonCopula(p), indepCopula()), c(0, 400)))
  tree1 <- spacetime_copula(list(spC(2), spC(1)), tlags = -(0:1))
  tree2 <- spacetime_copula(list(spC(1), spC(1), spC(1)), tlags = -1:1)
  vine <- distance_vine_copula(list(tree1, tree2), vineCopula(3L))
  expect_equal(vine@dimension, 5L)

  nb <- neighbours(cube, var = "u", size = 3, tlags = -(0:1))
  fit <- quiet(fitCopula(vine, list(nb, cube), method = list(familyset = 1)))
  expect_true(is.finite(fit@loglik))

  target <- cube["u", c(1, 5), 4:5]
  pn <- neighbours(cube, target, var = "u", size = 3, tlags = -(0:1), min_dist = 10)
  pred <- quiet(predict(fit@copula, pn, cube, target, list(q = identity)))
  q <- pred[["quantile.0.5"]]
  expect_equal(dim(q), c(2, 2))
  expect_true(all(q > 0 & q < 1))
})

test_that("Gaussian copula predictions go through predict()", {
  meuse <- meuseZinc()
  qMar <- function(p) qlnorm(p, mean(log(meuse$zinc)), sd(log(meuse$zinc)))
  gauss <- spatial_gauss_copula(function(h) exp(-h / 400))
  expect_output(show(gauss), "Gaussian")
  pn <- neighbours(meuse[1:15, ], meuse[16:20, ], var = "marZinc", size = 3, min_dist = 10)
  pred <- quiet(predict(gauss, pn, meuse[1:15, ], meuse[16:20, ], list(q = qMar)))
  # spcopula 0.2-5 (differs slightly: the shared code integrates over [0, 1])
  expect_equal(pred$quantile.0.5, c(736.835271732940, 386.727876071744, 318.511183245556,
                                    325.491620605963, 333.786307886115), tolerance = 1e-5)
  expE <- quiet(predict(gauss, pn, meuse[1:15, ], meuse[16:20, ], list(q = qMar), "expectation"))
  expect_true(all(expE$expect > 0))
  expect_equal(spatial_gauss_loglik(gauss, neighbours(meuse, var = "marZinc", size = 3)[1:5], meuse),
               quiet(spatial_gauss_loglik(gauss@corFun, neighbours(meuse, var = "marZinc", size = 3)[1:5], meuse)))
})

test_that("spatio-temporal covariate vines predict", {
  cube <- toyCube(6)
  set.seed(5)
  cube$u <- matrix(runif(30), 5, 6)
  cube$cov <- matrix(runif(30), 5, 6)
  spC <- function(p) quiet(spatial_copula(list(claytonCopula(p), indepCopula()), c(0, 400)))
  stCop <- spacetime_copula(list(spC(2), spC(1)), tlags = -(0:1))
  covar_cop <- function(stInd) normalCopula(0.5)
  cvvc <- covariate_vine_copula(covar_cop, stCop, vineCopula(5L))
  target <- cube[c("cov"), c(1, 5), 3]
  pn <- neighbours(cube, target, var = "u", covar = "cov", size = 3, tlags = -(0:1), min_dist = 10)
  pred <- quiet(predict(cvvc, pn, cube, target, list(q = identity)))
  expect_true(all(pred[["quantile.0.5"]] > 0 & pred[["quantile.0.5"]] < 1))
})

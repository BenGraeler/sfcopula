test_that("spatial neighbourhoods carry covariates from the data or the target", {
  meuse <- meuseZinc()
  meuse$marCopper <- rank(meuse$copper) / (nrow(meuse) + 1)
  nb <- neighbours(meuse, var = "marZinc", covar = "marCopper", size = 4)
  expect_equal(dim(nb@data), c(155, 5))
  expect_equal(nb@data[[5]], meuse$marCopper)
  expect_equal(colnames(nb@data)[5], "N0.marCopper")

  # prediction: the covariate is taken from the target locations
  target <- meuse[150:155, ]
  target$marCopper <- seq(0.1, 0.6, by = 0.1)
  pn <- neighbours(meuse[1:149, ], target, var = "marZinc", covar = "marCopper", size = 4)
  expect_equal(pn@data[[5]], target$marCopper)
  expect_error(neighbours(meuse[1:149, ], st_geometry(target), var = "marZinc",
                          covar = "marCopper", size = 4), "attributes of the target")

  # covariates are conditioned separately and not carried over to the next tree
  cv <- cond_covariate(nb, function(loc) normalCopula(0.5))
  expect_equal(cv, dduCopula(as.matrix(nb@data[, c(1, 5)]), normalCopula(0.5)))
  dropped <- drop_tree(nb, meuse, spatial_copula(list(claytonCopula(2), indepCopula()), c(0, 800)))
  expect_length(dropped@coVar, 0)
  expect_equal(dim(dropped@data), c(155, 3))

  # distance vines ignore the covariate column
  vine <- distance_vine_copula(spatial_copula(list(claytonCopula(2), indepCopula()), c(0, 800)),
                               vineCopula(3L))
  fit <- quiet(fitCopula(vine, nb, method = list(familyset = 1)))
  expect_true(is.finite(fit@loglik))
})

## reference values computed with spcopula 0.2-5's condStCoVarVine
covarModel <- function() {
  spC <- function(p) quiet(spatial_copula(list(claytonCopula(p), indepCopula()), c(0, 400)))
  stCop <- spacetime_copula(list(spC(2), spC(1)), tlags = -(0:1))
  top <- vineCopula(VineCopula::RVineMatrix(Matrix = matrix(c(3,1,2, 0,2,1, 0,0,1), 3),
                                            family = matrix(c(0,1,3, 0,0,4, 0,0,0), 3),
                                            par = matrix(c(0,0.4,1.2, 0,0,1.5, 0,0,0), 3)))
  covariate_vine_copula(function(stInd) normalCopula(0.3 + 0.1 * (stInd[2] %% 3)), stCop, top)
}

test_that("the covariate vine density reproduces spcopula's conditional density", {
  cvvc <- covarModel()
  dists <- array(c(100, 250, 0, -1), c(1, 2, 2))
  x <- c(0.05, 0.3, 0.5, 0.8, 0.97)
  ref <- c(0.223821123348400, 1.357105292500937, 1.355851761119450,
           0.928792640378287, 0.656651976483051)
  # neighbours first, covariate last
  f <- cond_vine(c(0.4, 0.6, 0.7), dists, cvvc, centre = c(2, 5))
  expect_equal(f(x), ref)
  # the deprecated spcopula function expects the covariate first
  g <- suppressWarnings(condStCoVarVine(c(0.7, 0.4, 0.6), dists, cvvc, stInd = c(2, 5)))
  expect_equal(g(x), ref)
  expect_equal(integrate(f, 0, 1)$value, 1, tolerance = 1e-3)

  # the density is vectorised over rows with individual distances and indices
  u <- cbind(x, 0.4, 0.6, 0.7)
  h <- array(c(100, 250, 0, -1), c(1, 2, 2))[rep(1, 5), , , drop = FALSE]
  d1 <- dCopula(u, cvvc, h = h, centre = cbind(2, c(5, 5, 6, 6, 7)))
  d2 <- sapply(1:5, function(i) dCopula(u[i, , drop = FALSE], cvvc, h = h[i, , , drop = FALSE],
                                        centre = c(2, c(5, 5, 6, 6, 7)[i])))
  expect_equal(d1, d2)
  expect_equal(dCopula(u, cvvc, h = h, centre = c(2, 5), log = TRUE),
               log(dCopula(u, cvvc, h = h, centre = c(2, 5))))
})

test_that("covariate vines work with a spatial tree", {
  meuse <- meuseZinc()
  meuse$marCopper <- rank(meuse$copper) / (nrow(meuse) + 1)
  tree <- spatial_copula(list(claytonCopula(2), indepCopula()), c(0, 800))
  cvc <- covariate_vine_copula(function(loc) normalCopula(0.8), tree, vineCopula(3L))
  expect_s4_class(cvc, "covariate_vine_copula")
  expect_equal(cvc@dimension, 4L)
  expect_output(show(cvc), "spatial tree")

  nb <- neighbours(meuse, var = "marZinc", covar = "marCopper", size = 3)
  u <- as.matrix(nb@data)

  # the density is the product of tree, covariate and top copula densities
  top <- cvc@topCop
  manual <- dCopula(u[, 1:2], tree, h = nb@distances[, 1], log = TRUE) +
    dCopula(u[, 1:3, drop = FALSE][, c(1, 3)], tree, h = nb@distances[, 2], log = TRUE) +
    dCopula(u[, c(1, 4)], normalCopula(0.8), log = TRUE) +
    dCopula(cbind(dduCopula(u[, c(1, 4)], normalCopula(0.8)),
                  dduCopula(u[, 1:2], tree, h = nb@distances[, 1]),
                  dduCopula(u[, c(1, 3)], tree, h = nb@distances[, 2])), top, log = TRUE)
  expect_equal(dCopula(u, cvc, h = nb@distances, centre = nb@index[, 1], log = TRUE), manual)

  # fitting in one call
  fit <- quiet(fitCopula(cvc, nb, method = list(familyset = 1:4)))
  expect_s4_class(fit@copula, "covariate_vine_copula")
  expect_equal(fit@loglik,
               sum(dCopula(u, fit@copula, h = nb@distances, centre = nb@index[, 1], log = TRUE)),
               tolerance = 1e-6)

  # conditional density and prediction with the covariate of the target
  f <- cond_vine(c(0.3, 0.5, 0.9), c(60, 120), fit@copula, centre = 1)
  expect_equal(integrate(f, 0, 1)$value, 1, tolerance = 1e-3)

  qMar <- function(p) qlnorm(p, mean(log(meuse$zinc)), sd(log(meuse$zinc)))
  target <- meuse[150:152, ]
  pn <- neighbours(meuse[1:149, ], target, var = "marZinc", covar = "marCopper", size = 3)
  pred <- quiet(predict(fit@copula, pn, meuse[1:149, ], target, list(q = qMar)))
  expect_s3_class(pred, "sf")
  expect_true(all(pred$quantile.0.5 > 0))
  # a higher covariate shifts the prediction up (positive covariate dependence)
  high <- target
  high$marCopper <- 0.99
  pnHigh <- neighbours(meuse[1:149, ], high, var = "marZinc", covar = "marCopper", size = 3)
  predHigh <- quiet(predict(fit@copula, pnHigh, meuse[1:149, ], high, list(q = qMar)))
  expect_true(all(predHigh$quantile.0.5 > pred$quantile.0.5))

  expect_error(fitCopula(covariate_vine_copula(function(loc) normalCopula(0.8),
                                               spacetime_copula(list(tree, tree), tlags = -(0:1)),
                                               vineCopula(3L)), nb), "spacetime_copula tree")
})

test_that("old spatio-temporal covariate vines are upgraded", {
  old <- suppressWarnings(readRDS(test_path("fixtures", "spcopula_objects.rds")))
  new <- upgrade_spcopula(old$objects$stCoVar)
  expect_s4_class(new, "covariate_vine_copula")
  expect_s4_class(new@tree, "spacetime_copula")
  f <- cond_vine(c(0.4, 0.6, 0.7), old$coVarDists, new, centre = c(2, 5))
  expect_equal(f(old$coVarX), old$expected$stCoVar)
})

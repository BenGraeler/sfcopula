## objects in the fixture were created with spcopula 0.2-5, see data-raw/spcopula_fixture.R

old <- suppressWarnings(readRDS(test_path("fixtures", "spcopula_objects.rds")))

test_that("old spcopula classes are translated into sfcopula classes", {
  new <- upgrade_spcopula(old$objects)
  expect_s4_class(new$spCop, "spatial_copula")
  expect_s4_class(new$stCop, "spacetime_copula")
  expect_s4_class(new$stCop@spCopList[[1]], "spatial_copula")
  expect_s4_class(new$spVine, "distance_vine_copula")
  expect_s4_class(new$spVine@topCop, "vineCopula")
  expect_null(new$pureVine@topCop)
  expect_length(new$stVine@trees, 1)
  expect_equal(new$geom@combination, "geometric")
  expect_s4_class(new$spVine@topCop, "vineCopula")
  expect_s4_class(new$pureVine, "distance_vine_copula")
  expect_s4_class(new$stVine, "distance_vine_copula")
  expect_s4_class(new$neigh, "neighbourhood")
  expect_s4_class(new$geom, "spatial_copula")
  expect_s4_class(new$stCoVar, "covariate_vine_copula")
  expect_s4_class(new$cqs, "cqsCopula")
  expect_identical(attr(class(new$cqs), "package"), "sfcopula")
  # nested in other S4 objects and lists
  expect_s4_class(new$fit, "fitCopula")
  expect_s4_class(new$fit@copula, "distance_vine_copula")
  expect_s4_class(new$nested$a, "spatial_copula")
  expect_s4_class(new$nested$b[[1]], "neighbourhood")
  expect_named(new, names(old$objects))
})

test_that("upgraded objects reproduce the results of spcopula", {
  new <- upgrade_spcopula(old$objects)
  u <- old$u; h <- old$h; u3 <- old$u3
  exp <- old$expected
  expect_equal(dCopula(u, new$spCop, h = h), exp$spCop)
  expect_equal(dCopula(u, new$spCopDep, h = h), exp$spCopDep)
  expect_equal(dCopula(u, new$stCop, h = cbind(h, c(0, -1, 0))), exp$stCop)
  expect_equal(dCopula(cbind(u, u), new$spVine,
                       h = list(matrix(c(50, 100, 150), 1)[rep(1, 3), ])), exp$spVine)
  expect_equal(dCopula(u3, new$pureVine,
                       h = list(matrix(c(50, 100), 1)[rep(1, 3), ], 120)), exp$pureVine)
  expect_equal(dCopula(u, new$cqs), exp$cqs)
  expect_equal(dCopula(u, new$geom, h = h), exp$geom)
})

test_that("upgraded neighbourhoods work with the new functions", {
  new <- upgrade_spcopula(old$objects$neigh)
  meuse <- meuseZinc()[1:30, ]
  fresh <- neighbours(meuse, var = "marZinc", size = 4)
  expect_equal(new@distances, fresh@distances)
  expect_equal(new@index, fresh@index)
  bins <- calc_bins(new, var = "marZinc", plot = FALSE)
  expect_true(is.list(bins))
})

test_that("load_spcopula upgrades all objects of an .RData file", {
  f <- tempfile(fileext = ".RData")
  spCop <- old$objects$spCop
  neigh <- old$objects$neigh
  save(spCop, neigh, file = f)
  e <- new.env()
  expect_setequal(load_spcopula(f, envir = e), c("spCop", "neigh"))
  expect_s4_class(e$spCop, "spatial_copula")
  expect_s4_class(e$neigh, "neighbourhood")
})

test_that("current objects pass unchanged", {
  cop <- spatial_copula(list(claytonCopula(2), indepCopula()), distances = c(0, 100), unit = "m")
  expect_identical(upgrade_spcopula(cop), cop)
  expect_identical(upgrade_spcopula(1:3), 1:3)
})

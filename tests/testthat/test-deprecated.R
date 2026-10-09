test_that("deprecated spcopula-style names forward to the new functions", {
  meuse <- meuseZinc()
  expect_warning(old <- getNeighbours(meuse, var = "marZinc", size = 3), "neighbours")
  new <- neighbours(meuse, var = "marZinc", size = 3)
  expect_identical(old, new)
  expect_s4_class(new, "neighbourhood")

  expect_warning(cop <- spCopula(list(claytonCopula(2), indepCopula()),
                                 distances = c(0, 100), unit = "m"), "spatial_copula")
  expect_s4_class(cop, "spatial_copula")
})

test_that("camelCase argument names of spcopula are accepted with a warning", {
  meuse <- meuseZinc()
  new <- neighbours(meuse, var = "marZinc", size = 3, min_dist = 10)
  expect_warning(old <- neighbours(meuse, var = "marZinc", size = 3, min.dist = 10),
                 "'min.dist' is deprecated; use 'min_dist'")
  expect_identical(old, new)

  expect_warning(cop <- spatial_copula(list(claytonCopula(2), indepCopula()), c(0, 100),
                                       spDepFun = function(h) if (is.null(h)) "kendall" else 0.3),
                 "spDepFun")
  bins <- calc_bins(meuse, var = "marZinc", nbins = 5, cutoff = 800, plot = FALSE)
  expect_warning(b2 <- calc_bins(meuse, var = "marZinc", nbins = 5, cutoff = 800, plot = FALSE,
                                 cor.method = "kendall"), "cor.method")
  expect_equal(attr(b2, "cor.method"), "kendall")

  # spcopula function names with spcopula arguments: only the function is reported
  w <- character()
  withCallingHandlers(getNeighbours(meuse, var = "marZinc", size = 3, min.dist = 10),
                      warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  expect_length(w, 1)
  expect_match(w, "getNeighbours")
  expect_warning(r <- rankTransform(meuse$zinc, meuse$copper), "rank_transform")
  expect_identical(r, rank_transform(meuse$zinc, meuse$copper))
  expect_warning(f <- genInvKenFun(function(t) t), "gen_inv_ken_fun")
  expect_true(is.function(f))
  expect_warning(f2 <- gen_inv_ken_fun(kenFun = function(t) t), "kenFun")
  expect_equal(f2(0.5), f(0.5))
})

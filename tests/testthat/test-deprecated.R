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

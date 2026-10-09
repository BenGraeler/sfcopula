test_that("as_spacetime_cube builds a vector data cube from long sf tables", {
  cube <- toyCube(3)
  long <- st_as_sf(cube, long = TRUE)
  # locations are ordered by first appearance, time instances are sorted
  long <- long[order(long$time, decreasing = TRUE), ]
  rebuilt <- as_spacetime_cube(long, time = "time")
  expect_s3_class(rebuilt, "stars")
  expect_equal(rebuilt[["u"]], cube[["u"]], ignore_attr = TRUE)
  expect_equal(st_get_dimension_values(rebuilt, "time"),
               st_get_dimension_values(cube, "time"))

  # missing combinations become NA
  partial <- as_spacetime_cube(long[-1, ], time = "time")
  expect_equal(sum(is.na(partial[["u"]])), 1)
})

test_that("spatio-temporal neighbourhoods are built from stars cubes", {
  cube <- toyCube(6)
  nb <- neighbours(cube, var = "u", size= 3, tlags = -(0:1), covar = "cov")
  # 5 locations x (6 - 1) time steps, 1 + 2*2 neighbours + 1 covariate
  expect_equal(dim(nb@data), c(25, 6))
  expect_equal(dim(nb@distances), c(25, 4, 2))
  u <- cube[["u"]]
  ind <- nb@index
  for (k in 1:5)
    expect_equal(nb@data[[k]], u[cbind(ind[, k, 1], ind[, k, 2])])
  expect_equal(nb@data[[6]], 1 - nb@data[[1]])
  expect_equal(unique(nb@distances[, , 2]), matrix(c(0, 0, -1, -1), 1))
})

test_that("spatio-temporal bins and the vine pipeline match spcopula", {
  data("EU_RB", package = "sfcopula", envir = environment())
  EU_RB$rtPM10 <- apply(EU_RB[["PM10"]], 2,
                        function(x) rank(x, na.last = "keep") / (sum(!is.na(x)) + 1))
  set.seed(42)
  stBins <- calc_bins(EU_RB, "rtPM10", nbins = 20, tlags = -(0:2))
  expect_equal(stBins$lagCor[, 1], c(0.332637009502096, 0.223765302746121,
                                     0.168141296575996))
})

test_that("spatio-temporal prediction works for stars and sf targets", {
  cube <- toyCube(6)
  spCops <- lapply(list(c(4, 2), c(2, 1)), function(p)
    quiet(spatial_copula(components = c(lapply(p, claytonCopula), indepCopula()),
                   distances = c(0, 100, 200), unit = "m")))
  stCop <- spacetime_copula(components = spCops, tlags = -(0:1))
  stVine <- distance_vine_copula(stCop, vineCopula(4L))

  target <- cube["cov", c(1, 5), 3:4]
  nb <- neighbours(cube, target, size= 3, tlags = -(0:1), var = "u",
                        prediction = TRUE, min_dist = 10)
  expect_equal(nrow(nb@data), 4)
  pred <- quiet(predict(stVine, nb, cube, target, list(q = identity), "quantile"))
  expect_s3_class(pred, "stars")
  expect_true("quantile.0.5" %in% names(pred))
  q <- pred[["quantile.0.5"]]
  expect_true(all(q > 0 & q < 1))

  tsf <- st_sf(time = as.Date("2020-01-03") + c(0, 1),
               geometry = st_get_dimension_values(target, "geometry"))
  nb2 <- neighbours(cube, tsf, size= 3, tlags = -(0:1), var = "u",
                         prediction = TRUE, min_dist = 10)
  pred2 <- quiet(predict(stVine, nb2, cube, tsf, list(q = identity), "quantile"))
  expect_equal(pred2$quantile.0.5, q[cbind(1:2, 1:2)], ignore_attr = TRUE)

  bad <- tsf
  bad$time <- as.Date("2021-01-01") + 0:1
  expect_error(neighbours(cube, bad, size= 3, tlags = -(0:1), var = "u",
                               prediction = TRUE), "time")
})

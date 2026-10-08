suppressMessages({
  library(sf)
  library(stars)
})

meuseZinc <- function() {
  data("meuse", package = "sfcopula", envir = environment())
  meuse$marZinc <- plnorm(meuse$zinc, mean(log(meuse$zinc)), sd(log(meuse$zinc)))
  meuse
}

toyCube <- function(nTime = 6) {
  pts <- st_sfc(st_point(c(0, 0)), st_point(c(100, 0)), st_point(c(0, 100)),
                st_point(c(100, 100)), st_point(c(250, 50)), crs = 28992)
  time <- as.Date("2020-01-01") + 0:(nTime - 1)
  values <- matrix(seq_len(5 * nTime) / (5 * nTime + 1), 5, nTime)
  st_as_stars(list(u = values, cov = 1 - values),
              dimensions = st_dimensions(geometry = pts, time = time))
}

quiet <- function(expr) {
  out <- NULL
  invisible(capture.output(out <- suppressMessages(expr)))
  out
}

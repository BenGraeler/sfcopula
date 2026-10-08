## Converts the sp/spacetime based data sets of spcopula into sf/stars.
## Requires the packages sp and spacetime (only for this conversion).
## Run from the package root: Rscript data-raw/convert_spcopula_data.R <path to spcopula>

library(sf)
library(stars)

args <- commandArgs(trailingOnly = TRUE)
spcopulaDir <- if (length(args) > 0) args[1] else "../spcopula"

STFDF2stars <- function(x) {
  geom <- st_geometry(st_as_sf(x@sp))
  time <- as.Date(zoo::index(x@time))
  nSp <- length(geom)
  vars <- lapply(x@data, function(v) matrix(v, nrow = nSp))
  # time invariant station attributes are replicated along time
  if ("data" %in% methods::slotNames(x@sp))
    vars <- c(vars, lapply(x@sp@data, function(v) matrix(v, nrow = nSp, ncol = length(time))))
  st_as_stars(vars, dimensions = st_dimensions(geometry = geom, time = time))
}

suppressMessages({library(sp); library(spacetime)})

load(file.path(spcopulaDir, "data", "EU_RB.RData"))
EU_RB <- STFDF2stars(EU_RB)
save(EU_RB, file = "data/EU_RB.RData", compress = "xz")

load(file.path(spcopulaDir, "data", "EU_RB_2005.RData"))
EU_RB_2005 <- STFDF2stars(EU_RB_2005)
save(EU_RB_2005, file = "data/EU_RB_2005.RData", compress = "xz")

load(file.path(spcopulaDir, "data", "simulatedTriples.RData"))
save(triples, file = "data/simulatedTriples.RData", compress = "xz")

# the meuse data set from sp as sf (EPSG:28992, Dutch RD New)
data("meuse", package = "sp", envir = environment())
meuse <- st_as_sf(meuse, coords = c("x", "y"), crs = 28992)
save(meuse, file = "data/meuse.RData", compress = "xz")

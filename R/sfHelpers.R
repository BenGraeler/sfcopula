##################################################################
##                                                              ##
## internal helpers wrapping sf and stars for the copula code   ##
##                                                              ##
##################################################################

## spatial data (sf/sfc)
#########################

# returns the point geometries of an sf or sfc object
.pointGeom <- function(x) {
  if (inherits(x, "stars"))
    return(.stGeom(x))
  if (!inherits(x, c("sf", "sfc")))
    stop("Spatial data needs to be provided as 'sf' or 'sfc' object.")
  geom <- st_geometry(x)
  if (!all(st_geometry_type(geom) == "POINT"))
    stop("Only POINT geometries are supported; consider sf::st_centroid() or sf::st_point_on_surface().")
  geom
}

# the attribute columns of an sf object (or none for sfc)
.attrNames <- function(x) {
  if (inherits(x, "sf"))
    return(setdiff(names(x), attr(x, "sf_column")))
  if (inherits(x, "stars"))
    return(names(x))
  character()
}

# a single attribute as plain vector
.attrValues <- function(x, var) {
  if (!var %in% .attrNames(x))
    stop("The variable '", var, "' is not part of the data.")
  x[[var]]
}

.isLonLat <- function(geom) {
  isTRUE(st_is_longlat(geom))
}

# distance matrix between two sets of points (plain numeric, CRS units;
# metres for geographic coordinates)
.spDistMat <- function(x, y = NULL) {
  if (.isLonLat(x)) {
    d <- if (is.null(y)) st_distance(x) else st_distance(x, y)
    return(matrix(as.numeric(d), nrow(d), ncol(d)))
  }
  cx <- st_coordinates(x)[, 1:2, drop = FALSE]
  cy <- if (is.null(y)) cx else st_coordinates(y)[, 1:2, drop = FALSE]
  d2 <- outer(cx[, 1], cy[, 1], "-")^2 + outer(cx[, 2], cy[, 2], "-")^2
  sqrt(d2)
}

# element-wise distances between the points geom[i] and geom[j]
.pairDists <- function(geom, i, j) {
  if (length(i) == 0)
    return(numeric(0))
  if (.isLonLat(geom))
    return(as.numeric(st_distance(geom[i], geom[j], by_element = TRUE)))
  cc <- st_coordinates(geom)[, 1:2, drop = FALSE]
  sqrt((cc[i, 1] - cc[j, 1])^2 + (cc[i, 2] - cc[j, 2])^2)
}

# length of the diagonal of the bounding box
.bboxDiag <- function(geom) {
  bb <- st_bbox(geom)
  corners <- st_sfc(st_point(bb[c("xmin", "ymin")]),
                    st_point(bb[c("xmax", "ymax")]), crs = st_crs(geom))
  .pairDists(corners, 1, 2)
}

# k nearest neighbours of each point in predGeom among dataGeom, ignoring
# all pairs closer than min.dist (e.g. the location itself); processed in
# chunks to limit the memory footprint of the distance matrix
.knn <- function(dataGeom, predGeom, k, min.dist) {
  nData <- length(dataGeom)
  nPred <- length(predGeom)
  index <- matrix(NA_integer_, nPred, k)
  dists <- matrix(NA_real_, nPred, k)
  if (k < 1)
    return(list(index = index, dists = dists))

  chunk <- max(1, floor(1e7 / nData))
  for (start in seq(1, nPred, by = chunk)) {
    rows <- start:min(nPred, start + chunk - 1)
    d <- .spDistMat(predGeom[rows], dataGeom)
    d[d < min.dist] <- Inf
    for (r in seq_along(rows)) {
      nn <- order(d[r, ])[1:min(k, nData)]
      nn <- nn[is.finite(d[r, nn])]
      if (length(nn) > 0) {
        index[rows[r], seq_along(nn)] <- nn
        dists[rows[r], seq_along(nn)] <- d[r, nn]
      }
    }
  }
  list(index = index, dists = dists)
}

# add prediction results to the target geometry
.addPrediction <- function(target, values, name) {
  if (inherits(target, "sf")) {
    target[[name]] <- values
    return(target)
  }
  res <- st_sf(data.frame(values), geometry = st_geometry(target))
  names(res)[1] <- name
  res
}

## spatio-temporal data (stars vector data cubes)
##################################################

# names of the geometry and time dimension of a vector data cube
.stDimNames <- function(x) {
  if (!inherits(x, "stars"))
    stop("Spatio-temporal data needs to be provided as 'stars' vector data cube.")
  dims <- st_dimensions(x)
  if (length(dims) != 2)
    stop("The vector data cube needs exactly two dimensions: geometry and time.")
  isGeom <- vapply(names(dims), function(d)
    inherits(st_get_dimension_values(x, d), "sfc"), logical(1))
  if (sum(isGeom) != 1)
    stop("The vector data cube needs exactly one geometry (sfc) dimension.")
  list(geom = names(dims)[isGeom], time = names(dims)[!isGeom])
}

.stGeom <- function(x) {
  geom <- st_get_dimension_values(x, .stDimNames(x)$geom)
  if (!all(st_geometry_type(geom) == "POINT"))
    stop("Only POINT geometries are supported.")
  geom
}

.stTime <- function(x) st_get_dimension_values(x, .stDimNames(x)$time)

# attribute as matrix with locations in rows and time instances in columns
.stValues <- function(x, var) {
  if (!var %in% names(x))
    stop("The variable '", var, "' is not part of the spatio-temporal data.")
  dn <- .stDimNames(x)
  arr <- unclass(x[[var]])
  if (match(dn$geom, names(st_dimensions(x))) == 2)
    arr <- t(arr)
  dim(arr) <- c(length(.stGeom(x)), length(.stTime(x)))
  arr
}

# spatio-temporal target locations as pairs of (geometry, time instance):
# either a full vector data cube (all combinations, locations outer, time
# inner) or an sf object with one time stamp per feature
.stTargets <- function(x, timeCol) {
  if (inherits(x, "stars")) {
    geom <- .stGeom(x)
    time <- .stTime(x)
    return(list(geom = geom, time = time,
                spInd = rep(seq_along(geom), each = length(time)),
                time4row = rep(time, length(geom)),
                tInd = rep(seq_along(time), length(geom))))
  }
  if (inherits(x, "sf")) {
    if (inherits(x, "sftime"))
      timeCol <- attr(x, "time_column")
    if (!timeCol %in% names(x))
      stop("The spatio-temporal target needs a time column '", timeCol, "'.")
    geom <- .pointGeom(x)
    return(list(geom = geom, time = x[[timeCol]], spInd = seq_along(geom),
                time4row = x[[timeCol]], tInd = seq_len(nrow(x))))
  }
  stop("Spatio-temporal targets need to be a 'stars' vector data cube or an 'sf' object with a time column.")
}

# value of a variable at the target pairs (e.g. a covariate)
.stTargetValues <- function(x, var, targets) {
  if (inherits(x, "stars"))
    return(.stValues(x, var)[cbind(targets$spInd, targets$tInd)])
  .attrValues(x, var)
}

# add prediction results to a spatio-temporal target
.addStPrediction <- function(target, values, name) {
  if (inherits(target, "stars")) {
    dn <- .stDimNames(target)
    nGeom <- length(.stGeom(target))
    arr <- matrix(values, nrow = nGeom, byrow = TRUE)
    if (match(dn$geom, names(st_dimensions(target))) == 2)
      arr <- t(arr)
    target[[name]] <- arr
    return(target)
  }
  target[[name]] <- values
  target
}

## conversion
##############

# builds a vector data cube from a long sf table holding one row per
# location and time instance (missing combinations become NA)
as_spacetime_cube <- function(x, time = "time", vars = NULL) {
  stopifnot(inherits(x, "sf"))
  if (inherits(x, "sftime"))
    time <- attr(x, "time_column")
  stopifnot(time %in% names(x))
  geom <- .pointGeom(x)

  if (is.null(vars))
    vars <- setdiff(.attrNames(x), time)

  cc <- st_coordinates(geom)[, 1:2, drop = FALSE]
  key <- paste(cc[, 1], cc[, 2])
  uKey <- unique(key)
  sInd <- match(key, uKey)
  uGeom <- geom[match(uKey, key)]

  uTime <- sort(unique(x[[time]]))
  tInd <- match(x[[time]], uTime)
  if (anyDuplicated(cbind(sInd, tInd)))
    stop("Locations and time instances need to be unique combinations.")

  values <- lapply(vars, function(v) {
    arr <- matrix(x[[v]][NA_integer_][1], length(uGeom), length(uTime))
    arr[cbind(sInd, tInd)] <- x[[v]]
    arr
  })
  names(values) <- vars

  st_as_stars(values, dimensions = st_dimensions(geometry = uGeom, time = uTime))
}

##############################################################
##                                                          ##
## recovery of objects created with spcopula (sp/spacetime) ##
##                                                          ##
##############################################################

# classes of spcopula and their sfcopula counterparts
.spcopulaClassMap <- c(
  spCopula = "spatial_copula", spGeomCopula = "spatial_copula",
  stCopula = "spacetime_copula",
  mixedSpVineCopula = "distance_vine_copula",
  pureSpVineCopula = "distance_vine_copula",
  stVineCopula = "distance_vine_copula",
  stCoVarVineCopula = "covariate_vine_copula",
  neighbourhood = "neighbourhood",
  stNeighbourhood = "neighbourhood",
  # classes that kept their names
  asCopula = "asCopula", cqsCopula = "cqsCopula", tawn3pCopula = "tawn3pCopula",
  empiricalCopula = "empiricalCopula", empSurCopula = "empSurCopula",
  leafCopula = "leafCopula", mixtureCopula = "mixtureCopula",
  hkCopula = "hkCopula", trunCopula = "trunCopula")

# slot translations for classes that have been merged
.spcopulaSlotMap <- list(
  spCopula = function(slots) c(slots, list(combination = "convex")),
  spGeomCopula = function(slots) 
    c(slots, list(combination = "geometric", calibMoa = function(copula, h) return(NULL),
                  fullname = "Spatial Copula: distance dependent geometric mean of bivariate copulas")),
  mixedSpVineCopula = function(slots) {
    slots$trees <- slots$spCop; slots$spCop <- NULL; slots },
  pureSpVineCopula = function(slots) {
    slots$trees <- slots$spCop; slots$spCop <- NULL; slots$topCop <- NULL
    c(slots, list(topCop = NULL)) },
  stVineCopula = function(slots) {
    slots$trees <- list(slots$stCop); slots$stCop <- NULL; slots },
  stCoVarVineCopula = function(slots) {
    slots$tree <- slots$stCop; slots$stCop <- NULL; slots })

# the class an S4 object needs to be rebuilt as, or NA if it is fine as it is
.upgradeTarget <- function(x) {
  cls <- class(x)
  pkg <- attr(cls, "package")
  if (is.null(pkg))
    return(NA_character_)
  if (pkg == "spcopula" && cls %in% names(.spcopulaClassMap))
    return(list(class = .spcopulaClassMap[[cls]], package = "sfcopula", 
                slotMap = .spcopulaSlotMap[[cls]]))
  # copula wrappers of VineCopula moved to VC2copula
  if (pkg == "VineCopula" &&
      !is.null(getClassDef(cls, where = asNamespace("VC2copula"))))
    return(list(class = as.character(cls), package = "VC2copula"))
  NA_character_
}

.upgradeSlots <- function(x, target) {
  slots <- attributes(x)
  slots$class <- NULL
  if (is.null(target)) {
    obj <- x
  } else {
    classDef <- getClassDef(target$class, where = asNamespace(target$package))
    obj <- new(classDef)
    if (!is.null(target$slotMap))
      slots <- target$slotMap(slots)
    dropped <- setdiff(names(slots), names(classDef@slots))
    if (length(dropped) > 0)
      warning("Dropping slot(s) ", paste(dropped, collapse = ", "),
              " not present in class '", target$class, "'.", call. = FALSE)
    slots <- slots[intersect(names(slots), names(classDef@slots))]
  }
  for (s in names(slots))
    slot(obj, s, check = FALSE) <- if (is.null(slots[[s]])) NULL else upgrade_spcopula(slots[[s]])
  obj
}

# translates objects of spcopula classes (also nested in lists or slots of
# other S4 objects, e.g. fitCopula results) into the classes of sfcopula
upgrade_spcopula <- function(x) {
  if (isS4(x)) {
    target <- .upgradeTarget(x)
    if (identical(target, NA_character_)) {
      if (!is.null(attr(class(x), "package")) &&
          attr(class(x), "package") == "spcopula")
        warning("No sfcopula counterpart for class '", class(x), "'; ",
                "the object is returned unchanged.", call. = FALSE)
      # other S4 classes (e.g. fitCopula) may hold spcopula objects in their slots
      if (length(slotNames(x)) == 0 || inherits(x, "function"))
        return(x)
      return(.upgradeSlots(x, NULL))
    }
    obj <- .upgradeSlots(x, target)
    msg <- tryCatch(validObject(obj, complete = TRUE), error = function(e) conditionMessage(e))
    if (is.character(msg))
      warning("The upgraded '", target$class, "' object is not valid: ", msg, call. = FALSE)
    return(obj)
  }
  if (is.list(x) && !is.data.frame(x)) {
    attrs <- attributes(x)
    x <- lapply(x, upgrade_spcopula)
    attributes(x) <- attrs
  }
  x
}

# loads an .RData/.rda file written with spcopula and upgrades all its objects
load_spcopula <- function(file, envir = parent.frame()) {
  tmp <- new.env()
  objs <- suppressWarnings(load(file, envir = tmp))
  for (o in objs)
    assign(o, upgrade_spcopula(get(o, envir = tmp)), envir = envir)
  invisible(objs)
}

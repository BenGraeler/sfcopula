#####################################################################
## fast, vectorised evaluation of bivariate copulas                ##
##                                                                 ##
## families known to VineCopula are evaluated in compiled code     ##
## with one parameter per pair; all other families fall back to    ##
## the methods of the copula package                               ##
#####################################################################

# VineCopula family code and second parameter of a copula object (NULL if not supported)
.vcFamily <- function(cop) {
  if ("family" %in% slotNames(cop))  # copula wrappers of VC2copula
    return(list(fam = cop@family,
                par2 = if (length(cop@parameters) > 1) cop@parameters[2] else 0))
  switch(class(cop)[1],
         normalCopula = list(fam = 1, par2 = 0),
         tCopula = list(fam = 2, par2 = cop@parameters[2]),
         claytonCopula = list(fam = 3, par2 = 0),
         gumbelCopula = list(fam = 4, par2 = 0),
         frankCopula = list(fam = 5, par2 = 0),
         NULL)
}

# parameters at which the one-parameter families equal the independence copula
.vcIndep <- function(fam, par) {
  switch(as.character(fam %% 10),
         "1" = par == 0, "3" = par == 0, "4" = par == 1, "5" = par == 0, "6" = par == 1,
         rep(FALSE, length(par)))
}

# parameter ranges of VineCopula for the families 1-6 (other families: trusted)
.vcValid <- function(fam, par, par2) {
  switch(as.character(fam),
         "1" = abs(par) < 1,
         "2" = abs(par) < 1 & par2 > 2,
         "3" = par > 0,
         "4" = par >= 1,
         "5" = par != 0 & is.finite(par),
         "6" = par >= 1,
         is.finite(par))
}

.vcCall <- function(type, u, fam, par, par2) {
  switch(type,
         p = VineCopula::BiCopCDF(u[,1], u[,2], fam, par, par2, check.pars = FALSE),
         d = VineCopula::BiCopPDF(u[,1], u[,2], fam, par, par2, check.pars = FALSE),
         ddu = VineCopula::BiCopHfunc1(u[,1], u[,2], fam, par, par2, check.pars = FALSE),
         ddv = VineCopula::BiCopHfunc2(u[,1], u[,2], fam, par, par2, check.pars = FALSE))
}

# evaluates cdf ("p"), density ("d") or a partial derivative ("ddu", "ddv") of a
# bivariate copula for all rows of u; par: NULL (parameters of cop) or a matrix
# with the leading parameters of each row
.compEval <- function(type, u, cop, par = NULL) {
  n <- nrow(u)
  if (n == 0)
    return(numeric(0))
  if (is(cop, "indepCopula"))
    return(switch(type, p = u[,1]*u[,2], d = rep(1, n), ddu = u[,2], ddv = u[,1]))

  res <- numeric(n)
  rows <- seq_len(n)
  vc <- .vcFamily(cop)
  if (!is.null(vc) && (is.null(par) || ncol(par) == 1)) {
    par1 <- if (is.null(par)) rep(cop@parameters[1], n) else par[,1]
    fam <- rep(vc$fam, n)
    indep <- .vcIndep(vc$fam, par1)
    ok <- indep | .vcValid(vc$fam, par1, vc$par2)
    fam[indep] <- 0
    if (any(ok))
      res[ok] <- .vcCall(type, u[ok,,drop=FALSE], fam[ok], par1[ok], rep(vc$par2, sum(ok)))
    if (type == "d")
      res[.onBoundary(u)] <- 0
    rows <- which(!ok)
    if (length(rows) == 0)
      return(res)
  }

  # fallback: methods of the copula package, grouped by parameter
  fun <- switch(type, p = pCopula, d = dCopula, ddu = dduCopula, ddv = ddvCopula)
  if (is.null(par)) {
    res[rows] <- fun(u[rows,,drop=FALSE], cop)
    return(res)
  }
  key <- apply(par[rows,,drop=FALSE], 1, paste, collapse=" ")
  for (k in unique(key)) {
    sel <- rows[key == k]
    tmpCop <- cop
    tmpCop@parameters[1:ncol(par)] <- par[sel[1],]
    res[sel] <- fun(u[sel,,drop=FALSE], tmpCop)
  }
  res
}

# parameters of a component calibrated by the spatial dependence function for
# each distance (once per distinct distance); a matrix with one row per distance
.calibrate <- function(cop, h, calibMoa) {
  uh <- unique(h)
  par <- tryCatch(calibMoa(cop, uh), error = function(e) NULL)
  if (is.numeric(par) && length(par) == length(uh) && length(cop@parameters) >= 1 &&
      !inherits(cop, c("asCopula", "cqsCopula"))) {
    parMat <- matrix(par, ncol = 1)
  } else {
    # families with several calibrated parameters: one distance at a time
    parList <- lapply(uh, function(x) calibMoa(cop, x))
    parMat <- do.call(rbind, parList)
  }
  parMat[match(h, uh), , drop = FALSE]
}

# evaluation of a spatial copula with convex combination for all pairs;
# several types (e.g. c("d", "ddu")) share one calibration and return a list
.spCopEval <- function(type, u, copula, h, log = FALSE) {
  n <- nrow(u)
  h <- rep(h, length.out = n)
  dists <- copula@distances
  nDists <- length(dists)
  depFun <- !is.null(copula@calibMoa(normalCopula(0), 0))
  comps <- copula@components

  # a matrix with one column per type
  evalComp <- function(cop, rows, calibrate) {
    par <- NULL
    if (calibrate && !is(cop, "indepCopula"))
      par <- .calibrate(cop, h[rows], copula@calibMoa)
    vapply(type, function(tp) .compEval(tp, u[rows,,drop=FALSE], cop, par), numeric(length(rows)))
  }

  res <- matrix(0, n, length(type))
  int <- findInterval(h, dists)
  for (k in unique(int)) {
    rows <- which(int == k)
    if (k == 0) {          # below the first distance
      res[rows,] <- evalComp(comps[[1]], rows, depFun)
    } else if (k >= nDists) { # beyond the range
      res[rows,] <- evalComp(comps[[nDists]], rows, FALSE)
    } else {
      lowerCop <- comps[[k]]
      upperCop <- comps[[k+1]]
      if (depFun && class(lowerCop) == class(upperCop)) {
        res[rows,] <- evalComp(lowerCop, rows, TRUE)
      } else {
        w <- (dists[k+1] - h[rows]) / (dists[k+1] - dists[k])
        res[rows,] <- w * evalComp(lowerCop, rows, depFun) +
          (1 - w) * evalComp(upperCop, rows, depFun)
      }
    }
  }

  if (log)
    res[, type == "d"] <- log(res[, type == "d"])
  if (length(type) == 1)
    return(res[,1])
  setNames(lapply(seq_along(type), function(i) res[,i]), type)
}

# log-density and first partial derivative of a distance tree's copula (spatial
# or spatio-temporal) for all pairs; h: vector (spatial) or matrix [n or 1, 2]
.treeEval <- function(u, tree, h, loglik = TRUE) {
  if (is(tree, "spacetime_copula")) {
    h <- matrix(h, ncol = 2)
    res <- list(d = numeric(nrow(u)), ddu = numeric(nrow(u)))
    for (t in unique(h[,2])) {
      if (is.na(match(t, tree@tlags)))
        stop("Prediction time(s) do(es) not match the modelled time slices.")
      sel <- if (nrow(h) == 1) rep(TRUE, nrow(u)) else h[,2] == t
      hSp <- if (nrow(h) == 1) h[1,1] else h[sel,1]
      tmp <- .treeEval(u[sel,,drop=FALSE], tree@spCopList[[match(t, tree@tlags)]], hSp, loglik)
      res$ddu[sel] <- tmp$ddu
      if (loglik) res$d[sel] <- tmp$d
    }
    return(res)
  }
  if (is(tree, "spatial_copula")) {
    if (tree@combination == "geometric") {
      h <- rep(h, length.out = nrow(u))
      return(list(d = if (loglik) dCop.spGeomCop(u, tree, h, do.logs = TRUE) else NULL,
                  ddu = dduCop.spGeomCop(u, tree, h)))
    }
    if (!loglik)
      return(list(ddu = .spCopEval("ddu", u, tree, h)))
    return(.spCopEval(c("d", "ddu"), u, tree, h, log = TRUE))
  }
  list(d = if (loglik) dCopula(u, tree, h = h, log = TRUE) else NULL,
       ddu = dduCopula(u, tree, h = h))
}

# rows on or outside the boundary of the unit hypercube; like the copula
# package, densities are 0 there
.onBoundary <- function(u) {
  out <- rowSums(u <= 0 | u >= 1) > 0
  out & !is.na(out)
}

# log-density of a top copula (vine copulas directly through VineCopula)
.topLogDens <- function(u, cop) {
  if (is(cop, "vineCopula")) {
    res <- rep(-Inf, nrow(u))
    inside <- !.onBoundary(u)
    if (any(inside))
      res[inside] <- VineCopula::RVineLogLik(u[inside,,drop=FALSE], cop@RVM, separate = TRUE)$loglik
    return(res)
  }
  dCopula(u, cop, log = TRUE)
}

# vectorised bisection solving fn(x) = y on (0, 1) for an increasing fn
.bisect <- function(fn, y, tol = 1e-12, maxit = 60) {
  lo <- rep(0, length(y))
  hi <- rep(1, length(y))
  for (it in 1:maxit) {
    mid <- (lo + hi) / 2
    below <- fn(mid) < y
    lo[below] <- mid[below]
    hi[!below] <- mid[!below]
    if (max(hi - lo) < tol)
      break
  }
  (lo + hi) / 2
}

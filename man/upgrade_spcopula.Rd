\name{upgrade_spcopula}
\alias{upgrade_spcopula}
\alias{load_spcopula}
\title{Translate objects created with spcopula into sfcopula classes}
\description{
Objects saved with the \pkg{spcopula} package (e.g. \code{spCopula}, \code{stCopula}, \code{mixedSpVineCopula}, \code{neighbourhood}) cannot be used by \pkg{sfcopula}, because their classes are renamed and their class definitions refer to \pkg{spcopula}. \code{upgrade_spcopula} rebuilds such objects as their \pkg{sfcopula} counterparts. \code{load_spcopula} loads an \code{.RData} file and upgrades all objects in it.
}
\usage{
upgrade_spcopula(x)
load_spcopula(file, envir = parent.frame())
}
\arguments{
  \item{x}{any R object. Objects of \pkg{spcopula} classes are translated, including when nested in lists or in slots of other S4 objects (e.g. the \code{copula} slot of a \code{fitCopula} result). All other objects are returned unchanged.}
  \item{file}{path to an \code{.RData}/\code{.rda} file written with \pkg{spcopula}.}
  \item{envir}{the environment to assign the upgraded objects to.}
}
\details{
The classes are translated as follows: \code{spCopula} to \code{\linkS4class{spatial_copula}} (with \code{combination = "convex"}), \code{spGeomCopula} to \code{\linkS4class{spatial_copula}} (with \code{combination = "geometric"}), \code{stCopula} to \code{\linkS4class{spacetime_copula}}, \code{mixedSpVineCopula}, \code{pureSpVineCopula} and \code{stVineCopula} to \code{\linkS4class{distance_vine_copula}} (the spatial copula(s) or the spatio-temporal copula become the \code{trees}; pure spatial vines get a \code{NULL} \code{top_cop}), \code{stCoVarVineCopula} to \code{\linkS4class{covariate_vine_copula}}, and \code{neighbourhood} and \code{stNeighbourhood} to \code{\linkS4class{neighbourhood}}. Copula families that kept their names (e.g. \code{asCopula}, \code{cqsCopula}, \code{empiricalCopula}) are re-registered with \pkg{sfcopula}, and copula wrappers created with old versions of \pkg{VineCopula} are re-registered with \pkg{VC2copula}.

Slots are copied (and renamed where classes have been merged); slots that the new class does not define are dropped with a warning, and a warning is issued if the upgraded object does not pass the validity check.

Neighbourhoods only hold indices into the data. They need to be used together with the same data, converted to \code{sf} or \code{stars} in the same order (see \code{\link{as_spacetime_cube}}). Note that distances of geographic coordinates are in kilometres in \pkg{spcopula} but in metres in \pkg{sfcopula}.

Functions stored in objects (e.g. the spatial dependence function of a \code{spatial_copula}) keep working when \pkg{sfcopula} is attached, which also attaches \pkg{copula}.
}
\value{
\code{upgrade_spcopula}: the upgraded object. \code{load_spcopula}: invisibly, the names of the loaded objects.
}
\seealso{\code{\link{sfcopula-deprecated}} for the renamed functions.}
\examples{
\dontrun{
# a workspace saved with spcopula
load_spcopula("my_spcopula_workspace.RData")

# or single objects
spCop <- upgrade_spcopula(readRDS("my_spCopula.rds"))
}
cop <- spatial_copula(list(claytonCopula(2), indepCopula()), distances = c(0, 100))
identical(upgrade_spcopula(cop), cop) # current objects pass unchanged
}
\keyword{spatial}

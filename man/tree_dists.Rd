\name{tree_dists}
\alias{tree_dists}

\title{
Calculating the conditional (spatio-temporal) distances
}
\description{
All distances of all involved neighbours are calculated according to multiple spatial or spatio-temporal trees forming a C-Vine.
}
\usage{tree_dists(neigh, data, n.trees)}

\arguments{
  \item{neigh}{the initial \code{\linkS4class{neighbourhood}}}
  \item{data}{the initial data locations (the same as used to generate \code{neigh}): an \code{\link[sf]{sf}} or \code{sfc} object with POINT geometries in the spatial case, the \code{stars} vector data cube or its \code{sfc} geometries in the spatio-temporal case.}
  \item{n.trees}{the number of distance dependent trees}
}

\value{
A \code{list} with one entry per tree holding the distances of the pairs involved in the conditional data pairs. The first entry are the distances of \code{neigh}. For tree \eqn{i+1}, the distances between the \eqn{i}-th neighbour and the remaining neighbours are provided: a matrix \code{[n, k-i]} in the spatial case and an array \code{[n, k-i, 2]} holding the spatial and temporal distances in the spatio-temporal case.
}

\author{
Benedikt Graeler
}

\seealso{
\code{\linkS4class{neighbourhood}}, \code{\link{drop_tree}}, \code{\linkS4class{distance_vine_copula}}
}

\examples{
library("sf")
data("meuse", package = "sfcopula")

meuseNeigh <- neighbours(meuse, var = "zinc", size = 4)
dists <- tree_dists(meuseNeigh, meuse, 3)
lapply(dists, head, n = 2)
}

\keyword{ spatial }
\keyword{ multivariate }
\keyword{distribution}

\name{loglik_by_lags}
\alias{loglik_by_lags}

\title{Log-likelihoods by copula family and spatial or spatio-temporal lag class}

\description{
This function calculates the log-likelihood for a set of provided copula families per lag class. The copulas' parameters are either fitted by a provided distance dependent function, or through the function \code{\link{fitCopula}} per lag and copula family. Spatial as well as spatio-temporal bins (as derived by \code{\link{calc_bins}} from a \code{stars} vector data cube) are supported; for the latter, the spatial lags are investigated for each temporal lag.
}
\usage{
loglik_by_lags(bins, data, families = c(normalCopula(),
                                             tCopula(),
                                             claytonCopula(),
                                             frankCopula(),
                                             gumbelCopula()),
                    calc_cor, lag_sub=1:length(bins$meanDists), ...)
}

\arguments{
  \item{bins}{Typically, the output of the \code{\link{calc_bins}} function providing columns \code{meanDists} to tune the copula's parameter and \code{lags} (or \code{lagData}) to evaluate the log-likelihood. Spatio-temporal bins are identified by their entries \code{lags$sp} and \code{lags$time}.}
  \item{data}{the spatial data set (\code{\link[sf]{sf}}) or the \code{stars} vector data cube used to derive the bins from for data retrieval. If missing (spatial case only), the entry \code{lagData} of \code{bins} is used.}
  \item{families}{A vector of representatives of the bivariate copula families that should be investigated.}
  \item{calc_cor}{The tuning function for the copulas' parameters. Typically, the output of the \code{\link{fit_cor_fun}} function. Its either based on Kendall's tau, Spearman's rho or the identity. This function needs to return one of "kendall", "spearman" or "id" when passed "NULL". This parameter might be "missing" for static convex spatial copulas. In the spatio-temporal case, a function of the spatial distance and the temporal lag index, \code{calcCor(h, t, tlags)}, as returned by \code{\link{fit_cor_fun}} with \code{tlags}.}
  \item{lag_sub}{A set of indices indicating which of the bins to use.}
  \item{\dots}{deprecated argument names of \pkg{spcopula} (\code{calcCor}, \code{lagSub}, \code{stBins}) are accepted with a warning.}
}

\value{A list containing a matrix (\code{loglik}) of spatial lags (rows) and copula family names (columns) holding the calculated log-likelihood value and a list of the corresponding copula fits. In the spatio-temporal case, a list holding such a result per temporal lag (named \code{loglik1}, \code{loglik2}, ...).}

\author{Benedikt Graeler}

\seealso{
See also \code{\link{calc_bins}}, \code{\link{fit_cor_fun}} or \code{\link{spatial_copula}}.
}

\examples{
# load spatial data
library("sf")
data("meuse", package = "sfcopula")

# drop margins
meuse$marZinc <- plnorm(meuse$zinc, mean(log(meuse$zinc)), sd(log(meuse$zinc)))

# load parts of the workspace from demo("spatial_copula", package = "sfcopula")
data("spCopDemo", package = "sfcopula") 

calcKTauPol <- fit_cor_fun(bins, degree=3)

loglikTau <- loglik_by_lags(bins, meuse, calc_cor=calcKTauPol)
loglikTau$loglik

# the spatio-temporal case
library("stars")
data("EU_RB", package = "sfcopula")

# rank transformed margins
EU_RB$rtPM10 <- array(rank(EU_RB$PM10, na.last = "keep") /
                        (sum(!is.na(EU_RB$PM10)) + 1),
                      dim(EU_RB$PM10))

stBins <- calc_bins(EU_RB, var="rtPM10", nbins=2, cutoff=2e5,
                   instances=10, tlags=-(0:1))

loglikTauSt <- loglik_by_lags(stBins, EU_RB,
                              families=c(normalCopula(), frankCopula()))
loglikTauSt$loglik1$loglik
}

\keyword{sfcopula}
\keyword{log-likelihood}
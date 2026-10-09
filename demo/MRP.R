library("sfcopula")

## get the data
data("simulatedTriples")

## rank order transformation
peakVol <- rank_transform(triples[,1],triples[,3])
colnames(peakVol) <- c("Qp","Vp")
plot(peakVol, asp=1)

# Kendall's tau correlation
cor(triples,method="kendall")

# estimate a Clayton copula by means of maximum likelihood. It captures the
# strong (lower tail) dependence of peak and volume (Kendall's tau of 0.85).
# Note: a BB7 copula cannot be fitted to this sample, as VineCopula restricts
# its second parameter to at most 6, which limits Kendall's tau to about 0.75.
copQV <- fitCopula(claytonCopula(), peakVol, method="ml")@copula
copQV
tau(copQV)

# we use a design return period of 10 years
# the MAR-case: given the 0.9 quantile of one marginal, what is the 
# corresponding quantile of the second variable with the given joint return 
# period of 1/(1-0.9)
v_MAR <- c(0.9, invdduCopula(0.9, copQV, 0.9))
v_MAR
pCopula(v_MAR, copQV)

## the analytical Kendall distribution
kendallFunQV <- get_kendall_distr(copQV)

# the Kendall distribution value (fraction of pairs having a smaller copula value than "t")
kendallFunQV(t=0.9)

curve(kendallFunQV, from=0, to=1, asp=1)

# the critical level of the KEN2-RP for 10 years
t_KEN2 <- critical_level(kendallFunQV, 10, mu=1)
t_KEN2

# the corresponding KEN2-RP for the OR-RP
kendall_rp(kendallFunQV, cl=0.9, mu=1, copula=copQV)

# the empirical counterparts: fraction of events beyond the OR and KEN2 levels
mean(pCopula(peakVol, copQV) > 0.9)
mean(pCopula(peakVol, copQV) > t_KEN2)

# illustrating the critical lines (empirically)
contour(copQV,pCopula,levels=c(0.9,t_KEN2),
        xlim=c(0.8,1), ylim=c(0.8,1), n=100, asp=1, col="blue")


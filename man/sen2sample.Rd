\name{sen2sample}
\alias{sen2sample}
\concept{Sensitivity analysis}
\concept{Observational study}
\concept{Causal inference}
\concept{Two-sample sensitivity analysis}
\concept{rbounds}
\title{
Treatment versus control sensitivity analysis without strata.
}
\description{
Performs a two-sample, treated-versus-control sensitivity analysis without strata or matching.  The method is described in Rosenbaum and Krieger (1990) and Rosenbaum (2002, Section 4.6).  The example in those references is used below to illustrate use of the sen2sample() function.
}
\usage{
sen2sample(sc, z, gamma = 1, alternative = "greater", method="BU")
}
\arguments{
  \item{sc}{
A vector of scored outcomes for one stratum.  For instance, for Wilcoxon's rank sum test,
these would be the ranks of the outcomes.
}
  \item{z}{
Treatment indicators, with length(z)=length(sc).  Here, z[i]=1 if i is treated and z[i]=0 if i is control.
}
  \item{gamma}{
The sensitivity parameter \eqn{\Gamma}, where \eqn{\Gamma \ge 1}.
}
  \item{alternative}{
If alternative="greater", then the test rejects for large values of the test statistic.  If alternative="less" then the test rejects for small values of the test statistic.  In a sensitivity analysis, it is safe but somewhat conservative to perform a two-sided test at level \eqn{\alpha} by doing two one-sided tests each at level \eqn{\alpha}/2.
}
  \item{method}{
If method="RK" or if method="BU", exact expectations and variances are used in a large sample approximation.  Methods "RK" and "BU" should give the same answer, but "RK" uses formulas from Rosenbaum and Krieger (1990), while "BU" obtains exact moments for the extended hypergeometric distribution using the BiasedUrn package and then applies Proposition 20, page 155, section 4.7.4 of Rosenbaum (2002).  In contrast, method="LS" does not use exact expectations and variances, but rather uses the large sample approximations in section 4.6.4 of Rosenbaum (2002).  Finally, method="AD" uses method="LS" for large strata and method="BU" for smaller strata.
}
}

\value{
  \item{sc}{
A vector of scored outcomes for one stratum.  For instance, for Wilcoxon's rank sum test,
these would be the ranks of the outcomes in the current stratum.
}
  \item{z}{
Treatment indicators, with length(z)=length(sc).  Here, z[i]=1 if i is treated and z[i]=0 if i is control.
}
  \item{g}{
The sensitivity parameter \eqn{\Gamma}, where \eqn{\Gamma \ge 1}.
}
}
\references{
Rosenbaum, P. R. and Krieger, A. M. (1990). Sensitivity of two-sample permutation inferences in observational studies.  Journal of the American Statistical Association, 85, 493-498.

Rosenbaum, P. R. (2002). Observational Studies (2nd edition). New York: Springer.  Section 4.6.

Skerfving, S., Hansson, K., Mangs, C., Lindsten, J., & Ryman, N. (1974). Methylmercury-induced chromosome damage in man. Environmental Research, 7, 83-98.
}
\author{
Paul R. Rosenbaum
}
\note{
The example is from Table 1, page 497, of Rosenbaum and Krieger (1990).  The example is also Table 4.15, page 146, in Rosenbaum (2002).  The data are orignally from Skerfving et al. (1974).
}

\seealso{
For binary responses, use the sensitivity2x2xk package.  For matched responses, use one of the following packages: sensitivitymult, sensitivitymv, sensitivitymw, sensitivityfull.
}
\examples{

mercury<-c(5.3, 15, 11, 5.8, 17, 7, 8.5, 9.4, 7.8, 12, 8.7, 4, 3, 12.2, 6.1, 10.2,
           100, 70, 196, 69, 370, 270, 150, 60, 330, 1100, 40, 100, 70,
           150, 200, 304, 236, 178, 41, 120, 330, 62, 12.8)
z<-c(rep(0,16),rep(1,23))
CuCells<-c(2.7, .5, 0, 0, 5, 0, 0, 1.3, 0, 1.8, 0, 0, 1.0, 1.8,
           0, 3.1, .7, 4.6, 0, 1.7, 5.2, 0, 5, 9.5, 2, 3, 1, 3.5,
           2, 5, 5.5, 2, 3, 4, 0, 2, 2.2, 0, 2)

#Reproduces Rosenbaum and Krieger (1990), page 497
sen2sample(rank(mercury),z,gamma=5)
#Reproduces Rosenbaum and Krieger (1990), page 497
sen2sample(rank(CuCells),z,gamma=2)
(551.500000-492.334479)/sqrt(1153.775252) #Computation of the deviate
#Intermediate calculations: expectation and variance are in row 21.
evall(rank(CuCells),z,2,method="RK")

#The following three examples, if run, reproduce the
#calculations in the final paragraph of page 145,
#Section 4.6.6 of Rosenbaum (2002) Observational Studies, 2nd Ed.
#The first calculation uses large sample approximations
#to expectations and variances.
sen2sample(rank(mercury),z,gamma=2,method="LS")
#The next two calculations use exact expectations and variances
sen2sample(rank(mercury),z,gamma=2,method="RK")
sen2sample(rank(mercury),z,gamma=2,method="BU")
}

\keyword{ htest }
\keyword{ robust }

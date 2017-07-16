\name{evall}
\alias{evall}
\title{
Compute expectations and variances for one stratum.
}
\description{
Of limited interest to most users, the evall() function plays an internal role in 2-sample and stratified sensitivity analyses.  The expectation and variance returned by the evall() function are defined in the third paragraph of section 4, page 495, of Rosenbaum and Krieger (1990).  The function evall() calls the function ev() to determine the expectation and variance of the test statistic for an unobserved covariate u with length(z)-m 0's followed by m 1's, doing this for m=1,...,length(z)-1.
}
\usage{
evall(sc, z, g, method)
}
\arguments{
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
  \item{method}{
If method="RK" or if method="BU", exact expectations and variances are used in a large sample approximation.  Methods "RK" and "BU" should give the same answer, but "RK" uses formulas from Rosenbaum and Krieger (1990), while "BU" obtains exact moments for the extended hypergeometric distribution using the BiasedUrn package and then applies Proposition 20, page 155, section 4.7.4 of Rosenbaum (2002).  In contrast, method="LS" does not use exact expectations and variances, but rather uses the large sample approximations in section 4.6.4 of Rosenbaum (2002).  Finally, method="AD" uses method="LS" for large strata and method="BU" for smaller strata.
}
}
\details{
The evall() function is called by the sen2sample() function and the senstrat() function.
}
\value{
A data.frame with length(z)-1 rows and three columns.  The first column, m, gives the number of 1's in the unobserved covariate vector, u.  The second column, expect, and the third column, var, give the expectation and variance of the test statistic for this u.
}
\references{
Rosenbaum, P. R. and Krieger, A. M. (1990). Sensitivity of two-sample permutation inferences in observational studies.  Journal of the American Statistical Association, 85, 493-498.

Rosenbaum, P. R. (2002). Observational Studies (2nd edition). New York: Springer.  Section 4.6.
}
\author{
Paul R. Rosenbaum
}

\note{
The example is from Table 1, page 497, of Rosenbaum and Krieger (1990).  The example is also Table 4.15, page 146, in Rosenbaum (2002).  The example refers to Cu cells.  The data are orignally from Skerfving et al. (1974).
}

\examples{
z<-c(rep(0,16),rep(1,23))
CuCells<-c(2.7, .5, 0, 0, 5, 0, 0, 1.3, 0, 1.8, 0, 0, 1.0, 1.8,
           0, 3.1, .7, 4.6, 0, 1.7, 5.2, 0, 5, 9.5, 2, 3, 1, 3.5,
           2, 5, 5.5, 2, 3, 4, 0, 2, 2.2, 0, 2)
evall(rank(CuCells),z,2,"RK")
}

\keyword{ htest }
\keyword{ robust }

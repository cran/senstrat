\name{ev}
\alias{ev}

\title{
Computes the null expectation and variance for one stratum.
}
\description{
Of limited interest to most users, the ev function plays an internal role in 2-sample and stratified sensitivity analyses.  The expectation and variance returned by the ev function are defined in the third paragraph of section 4, page 495, of Rosenbaum and Krieger (1990).
}
\usage{
ev(sc, z, m, g, method)
}
\arguments{
  \item{sc}{
A vector of scored outcomes for one stratum.  For instance, for Wilcoxon's rank sum test,
these would be the ranks of the outcomes in the current stratum.
}
  \item{z}{
Treatment indicators, with length(z)=length(sc).  Here, z[i]=1 if i is treated and z[i]=0 if i is control.
}
  \item{m}{
The unobserved covariate u has length(z)-m 0's followed by m 1's.
}
  \item{g}{
The sensitivity parameter \eqn{\Gamma}, where \eqn{\Gamma \ge 1}.
}
  \item{method}{
If method="RK" or if method="BU", exact expectations and variances are used in a large sample approximation.  Methods "RK" and "BU" should give the same answer, but "RK" uses formulas from Rosenbaum and Krieger (1990), while "BU" obtains exact moments for the extended hypergeometric distribution using the BiasedUrn package and then applies Proposition 20, page 155, section 4.7.4 of Rosenbaum (2002).  In contrast, method="LS" does not use exact expectations and variances, but rather uses the large sample approximations in section 4.6.4 of Rosenbaum (2002).  Finally, method="AD" uses method="LS" for large strata and method="BU" for smaller strata.
}
}
\details{
The function ev() is called by the function evall().
}
\value{
\item{expect }{Null expectation of the test statistic.}
\item{vari }{Null variance of the test statistic.}
}
\references{

Rosenbaum, P. R. and Krieger, A. M. (1990). Sensitivity of two-sample permutation inferences in observational studies.  Journal of the American Statistical Association, 85, 493-498.

Rosenbaum, P. R. (2002). Observational Studies (2nd edition). New York: Springer.  Section 4.6.

}
\author{
Paul R. Rosenbaum
}



\examples{
ev(1:5,c(0,1,0,1,0),3,2,"RK")
ev(1:5,c(0,1,0,1,0),3,2,"BU")
}
\keyword{ htest }
\keyword{ robust }

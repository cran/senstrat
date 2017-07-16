\name{zeta}
\alias{zeta}
\title{
zeta function in sensitivity analysis
}
\description{
Of limited interest to most users, the zeta function plays an internal role in 2-sample and stratified sensitivity analyses.  The zeta function is equation (8), page 495, in Rosenbaum and Krieger (1990).
}
\usage{
zeta(bigN, n, m, g)
}
\arguments{
  \item{bigN}{
Total sample size in this stratum.
}
  \item{n}{
Treated sample size in this stratum.
}
  \item{m}{
The number of 1's in the vector u of unobserved covariates.  Here, u
has bigN-m 0's followed by m 1's.
}
  \item{g}{
The sensitivity parameter \eqn{\Gamma}, where \eqn{\Gamma \ge 1}.
}
}

\value{
The value of the zeta function.
}
\references{

Rosenbaum, P. R. and Krieger, A. M. (1990). Sensitivity of two-sampler permutation inferences in observational studies.  Journal of the American Statistical Association, 85, 493-498.

Rosenbaum, P. R. (2002). Observational Studies (2nd edition). New York: Springer.  Section 4.6.

}
\author{
Paul R. Rosenbaum
}
\note{
The zeta function is called by computep.
}


\examples{
zeta(10,5,6,2)
}


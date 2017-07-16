\name{computep}
\alias{computep}
\title{
Computes individual and pairwise treatment assignment probabilities.
}
\description{
Of limited interest to most users, the computep function plays an internal role in 2-sample and stratified sensitivity analyses.  The computep function is equations (9) and (10), page 496, in Rosenbaum and Krieger (1990).
}
\usage{
computep(bigN, n, m, g)
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
\item{p1 }{Equation (9), page 496, in Rosenbaum and Krieger (1990) evaluated with u[i]=1.}
\item{p0 }{Equation (9), page 496, in Rosenbaum and Krieger (1990) evaluated with u[i]=0.}
\item{p11 }{Equation (10), page 496, in Rosenbaum and Krieger (1990) evaluated with u[i]=1, u[j]=1.}
\item{p10 }{Equation (10), page 496, in Rosenbaum and Krieger (1990) evaluated with u[i]=1, u[j]=0.}
\item{p00 }{Equation (10), page 496, in Rosenbaum and Krieger (1990) evaluated with u[i]=0, u[j]=0.}
}
\references{

Rosenbaum, P. R. and Krieger, A. M. (1990). Sensitivity of two-sample permutation inferences in observational studies.  Journal of the American Statistical Association, 85, 493-498.

Rosenbaum, P. R. (2002). Observational Studies (2nd edition). New York: Springer.  Section 4.6.

}
\author{
Paul R. Rosenbaum
}
\note{
The function computep is called by the function ev.
}


\examples{
computep(10,5,6,2)
}


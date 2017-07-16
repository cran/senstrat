\name{mscores}
\alias{mscores}
\concept{M-estimate}
\concept{Stratification}
\concept{M-test}
\title{
Computes M-Scores for Two-Sample or Stratified Permutation Inference.
}
\description{
Computes M-scores for two-sample or stratified permutation inference and sensitivity analyses.  For instance, the scores may be used in functions sen2sample() and senstrat().
}
\usage{
mscores(y, z, st = NULL, inner = 0, trim = 3, lambda = 0.5,
     tau = 0)
}
\arguments{
  \item{y}{
A vector outcome.  An error will result if y contains NAs.
}
  \item{z}{
Treatment indicators, with length(z)=length(y).  Here, z[i]=1 if i is treated and z[i]=0 if i is control.
}
  \item{st}{
For an unstratified two-sample comparison, st = NULL.  For a stratified comparison, st is a vector of stratum indicators with length(st)=length(y).  The vector st may be numeric, such as 1,2,..., or it may be a factor.
}
  \item{inner}{
See trim below.
}
  \item{trim}{
The two parameters, inner>=0 and trim>inner, determine the odd \eqn{\psi}-function of the M-statistic.  Huber favored an M-statistic similar to a trimmed mean, with inner=0, trim>0, so that \eqn{\psi(y)} is max(-trim,min(trim,y)).  Setting inner>0 ignores y near zero and is most useful for matched pairs; see Rosenbaum (2013).  In general, \eqn{\psi(y)} is sign(y)*trim*min(trim-inner,max(0,abs(y)-inner))/(trim-inner), so it is zero for y in the interval [-inner,inner], and rises linearly from 0 to trim on the interval y=inner to y=trim, and equals trim for y>trim.  If trim=Inf, then no trimming is done, and \eqn{\psi(y)} = y.
}
  \item{lambda}{
The \eqn{\psi}-function is applied to a treated-minus-control difference in responses after scaling by the lambda quantile of the within-strata absolute differences.  Typically, lambda=1/2 for the median.
}
  \item{tau}{
If tau=0, then the null hypothesis is Fisher's sharp null hypothesis of no treatment effect.  If tau is nonzero, then tau is subtracted from the y's for treated responses before scoring the y's.  If tau is nonzero, the null hypothesis is that the treatment has a constant additive effect tau, the same for all strata.
}
}
\value{
A vector of length(y) containing the M-scores.  The M-scores may be used in senstrat() or sen2Sample().
}
\references{

Huber, P. (1981) Robust Statistics. New York: Wiley, 1981.  M-statistics were developed by Huber.

Maritz, J. S. (1979) Exact robust confidence intervals for location. Biometrika 1979, 66, 163-166.  Maritz proposed small adjustments to M-statistics for matched pairs that permit them to be used in exact permutation tests and confidence intervals.

Rosenbaum, P. R. (2007) Sensitivity analysis for m-estimates, tests and confidence intervals in matched observational studies.  Biometrics, 2007, 63, 456-464. <doi:10.1111/j.1541-0420.2006.00717.x>  Extends the method of Maritz (1979) to matching with multiple controls and to sensitivity analysis in observational studies.

Rosenbaum, P. R. (2013) Impact of multiple matched controls on design sensitivity in observational studies. Biometrics, 2013, 69, 118-127.  Computes the design sensitivity of M-statistics and proposes the trim parameter for use with matched pairs to increase design sensitivity.

Rosenbaum, P. R. (2014) Weighted M-statistics with superior design sensitivity in matched observational studies with multiple controls. Journal of the American Statistical Association 109, 1145-1158.  Discusses weights for matched sets to increase design sensitivity.

Rosenbaum, P. R. (2017) Sensitivity analysis for stratified comparisons in an observational study
of the effect of smoking on homocysteine levels.  Manuscript.  Among other things, extends the method
of Maritz (1979) to stratified comparisons.
}
\author{
Paul R. Rosenbaum
}

\note{
The ith response in stratum s, R[si], is compared to the jth response in stratum s, R[sj], as
\eqn{\psi}((R[si]-R[sj])/\eqn{\sigma}), and for each fixed i these values are averaged over the n[s]-1 choices of j in stratum s, where n[s] is the size of stratum s, thereby producing the M-score for R[si].  Here, \eqn{\sigma} is the lambda quantile, usually the median, of |R[si]-R[sj]|, taken over all within stratum comparisons.  This extension to stratified comparisons of the method of Maritz (1979) for matched pairs is described in Rosenbaum (2007, 2017).
}

\seealso{
The packages senitivitymv, sensitivitymw and sensitivityfull use M-scores in matched sets.  The M-scores from those packages are similar, but are weighted differently, particularly when matched sets have varying sizes.
}

\examples{
data("homocyst")
attach(homocyst)
sc<-mscores(log2(homocysteine),z,st=stf)
par(mfrow=c(1,2))
boxplot(log2(homocysteine)~z,main="Data")
boxplot(sc~z,main="Mscores")
detach(homocyst)
}

\keyword{ htest }
\keyword{ robust }

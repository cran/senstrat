\name{hodgeslehmann}
\alias{hodgeslehmann}

\title{
Computes Hodges-Lehmann Aligned Ranks.
}
\description{
Computes Hodges-Lehmann (1962) aligned ranks for use in stratified sensitivity analyses.  For instance, the scores, sc, used in function senstrat() might be Hodges-Lehmann aligned ranks.
}
\usage{
hodgeslehmann(y,z,st,align="median",tau=0)
}
\arguments{
  \item{y}{
A vector outcome.  An error will result if y contains NAs.
}
  \item{z}{
Treatment indicators, with length(z)=length(y).  Here, z[i]=1 if i is treated and z[i]=0 if i is control.
}
  \item{st}{
A vector of stratum indicators with length(st)=length(y).  The vector st may be numeric, such as 1,2,..., or it may be a factor.
}
  \item{align}{
Within each stratum, the observations are aligned or centered by subtracting a location measure computed from observations in that stratum.  If align="median", then the median is subtracted.  If align="mean", then the mean is subtracted.   If align="huber", then the Huber's m-estimate is subtracted, as computed by the function huber() in the MASS package.  If align="hl", then the one-sample Hodges-Lehmann estimate is subtracted, as suggested by Mehrotra et al. (2010), and as computed by the wilcox.test() function in the stats package.  The wilcox.test() command may generate numerous warnings about situations that are not hazardous, such as ties.  When align="hl", warnings produced by wilcox.test() are suppressed.
}
  \item{tau}{
If tau=0, then the null hypothesis is Fisher's sharp null hypothesis of no treatment effect.  If tau is nonzero, then tau is subtracted from the y's for treated responses before aligning and ranking.  If tau is nonzero, the null hypothesis is that the treatment has a constant additive effect tau, the same for all strata.
}
}
\value{
A vector of length(y) containing the aligned ranks.
}
\references{

Hodges, J. L. and Lehmann, E. L. (1962) Rank methods for combination of independent experiments in analysis of variance.  Annals of Mathematical Statistics, 33, 482-497.

Lehmann, E. L. (1975) Nonparametrics.  San Francisco: Holden-Day.

Mehrotra, D. V., Lu, X., and Xiaoming, L. (2010).  Rank-based analysis of stratified experiments: alternatives to the van Elteren test. American Statistician, 64, 121-130.

}
\author{
Paul R. Rosenbaum
}

\examples{
data("homocyst")
attach(homocyst)
sc<-hodgeslehmann(log2(homocysteine),z,stf,align="hl")
summary(sc)
length(sc)
detach(homocyst)
}

\keyword{ htest }
\keyword{ robust }

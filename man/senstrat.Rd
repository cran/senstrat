\name{senstrat}
\alias{senstrat}
\concept{Sensitivity analysis}
\concept{Observational study}
\concept{Causal inference}
\concept{Stratification}
\concept{M-estimate}
\concept{Hodges-Lehmann aligned ranks}
\concept{rbounds}
\title{
Sensitivity Analysis With Strata.
}
\description{
Performs a sensitivity analysis with strata.   The underlying method is described in Rosenbaum (2017),
and the main example illustrates calculations from that paper.  Use sen2sample() if there are no strata.
}
\usage{
senstrat(sc, z, st, gamma = 1, alternative = "greater",
       level = 0.05, method="BU", detail = FALSE)
}
\arguments{
  \item{sc}{
A vector of scored outcomes.  For instance, these scored outcomes might be produced by the hodgeslehmann() function or the mscores() function.  An error will result if sc contains NAs.
}
  \item{z}{
Treatment indicators, with length(z)=length(sc).  Here, z[i]=1 if i is treated and z[i]=0 if i is control.
}
  \item{st}{
Vector of stratum indicators, with length(st)=length(sc).  The vector st may be numeric, say 1, 2, ..., or it may be a factor.  A factor will be converted to integers for computations.  If there is only one stratum with data, then it is better to use sen2sample, not senstrat, and a warning will be given.
}
  \item{gamma}{
The sensitivity parameter \eqn{\Gamma}, where \eqn{\Gamma \ge 1}.
}
  \item{alternative}{
If alternative="greater", then the test rejects for large values of the test statistic.  If alternative="less" then the test rejects for small values of the test statistic.  In a sensitivity analysis, it is safe but somewhat conservative to perform a two-sided test at level \eqn{\alpha} by doing two one-sided tests each at level \eqn{\alpha}/2.
}
  \item{level}{
The sensitivity analysis is for a test of the null hypothesis of no treatment effect performed with the stated level, conventionally level=0.05.  If there is no treatment effect, so the null hypothesis is true, and if the bias in treatment assignment is at most \eqn{\Gamma}, then the chance that the sensitivity analysis will falsely reject the null hypothesis is at most 0.05 if level=0.05.  It would be common to report that rejection of the null hypothesis at conventional level=0.05 is insensitive to a bias of \eqn{\Gamma} if this \eqn{\Gamma} is the largest \eqn{\Gamma} leading to rejection.  Determining this largest \eqn{\Gamma} would entail running senstrat several times with different values of \eqn{\Gamma}.  It is perfectly reasonable, if less conventional, to conduct a sensitivity analysis with some other level, say level=0.01, and this might be necessary if one were testing several hypotheses, correcting for multiple testing using the Bonferroni or Holm procedures.
}
  \item{method}{
If method="RK" or if method="BU", exact expectations and variances are used in a large sample approximation.  Methods "RK" and "BU" should give the same answer, but "RK" uses formulas from Rosenbaum and Krieger (1990), while "BU" obtains exact moments for the extended hypergeometric distribution using the BiasedUrn package and then applies Proposition 20, page 155, section 4.7.4 of Rosenbaum (2002).  In contrast, method="LS" does not use exact expectations and variances, but rather uses the large sample approximations in section 4.6.4 of Rosenbaum (2002).  Finally, method="AD" uses method="LS" for large strata and method="BU" for smaller strata.
}
  \item{detail}{
If detail=FALSE, concise practical output is produced.  The option detail=TRUE provides additional details about the computations, which may be useful in understanding the computations or in trouble shooting, but the additional details are not useful in data analysis.
}
}
\details{
The method uses a Normal approximation to the distribution of the test statistic.  If method is not "LS", then this approximation is suitable for either a few strata containing many people or for many strata each stratum containing only a few people.  In contrast, method="LS" is useful only if every single stratum contains a large sample.
}
\value{
\item{Conclusion }{An English sentence stating the conclusion of the sensitivity analysis.  This sentence says whether the null hypothesis has been rejected at the stated level, defaulting to level=0.05, in the presence of a bias of at most \eqn{\Gamma}.  This statement is correct, if perhaps ever so slightly conservative.  This statement should be understood as the conclusion produced by senstrat, with the remaining output being descriptive or approximate.}
\item{Result }{Numeric results, including an approximate P-value, the deviate that
was compared with the Normal distribution to produce the P-value, the test statistic formed as the sum of the scores for treated individuals, its null expectation and variance, and the value of \eqn{\Gamma}.}
\item{Description }{The number of nondegenerate strata used in computations and the total number of treated individuals and controls in those strata.}
\item{StrataUse }{An English sentence stating whether all strata were used or alternatively describing degenerate strata that were not used.  See the Note.}
\item{LinearBound}{If detail=TRUE, then the Result above is labeled as the LinearBound; it is a safe but perhaps slightly conservative P-value.}
\item{Separable}{If detail=TRUE, then separable approximation of Gastwirth et al. (2000) is reported; it is a slightly liberal P-value.  In many if not most examples, LinearBound and Separable are in close agreement, so the issue of liberal versus conservative does not arise.  The stated Conclusion is based on the conservative LinearBound.}
\item{Remark}{If detail=TRUE, then Remark contains an English sentence commenting upon the agreement of the LinearBound and the Separable approximation.}
\item{lambda}{If detail=TRUE, then the values of the separable \eqn{\lambda(b)} and the linear Taylor bound  \eqn{\lambda(b)+\sum\eta_{s}} from Rosenbaum (2017).  The Remark above says there is agreement if these two quantities have the same sign.}
}
\references{

Gastwirth, J. L., Krieger, A. M., and Rosenbaum, P. R. (2000) Asymptotic separability in sensitivity analysis.  Journal of the Royal Statistical Society B 2000, 62, 545-556. <doi:10.1111/1467-9868.00249>

Hansen, B. B. (2004) Full matching in an observational study of coaching for the SAT.  Journal of the American Statistical Association, 99, 609-618.  Application of full matching as an alternative to conventional stratification.  See also Hansen's R package optmatch.

Rosenbaum, P. R. and Krieger, A. M. (1990). Sensitivity of two-sample permutation inferences in observational studies.  Journal of the American Statistical Association, 85, 493-498.

Rosenbaum, P. R. (1991) A characterization of optimal designs for observational studies. Journal of the Royal Statistical Society B, 53, 597-610. Introduces full matching as an alternative to conventional stratification.

Rosenbaum, P. R. (2002). Observational Studies (2nd edition). New York: Springer.  Section 4.6.

Rosenbaum, P. R. (2007) Sensitivity analysis for m-estimates, tests and confidence intervals in matched observational studies.  Biometrics, 2007, 63, 456-464. <doi:10.1111/j.1541-0420.2006.00717.x> See the erpcp example below.

Rosenbaum, P. R. (2014) Weighted M-statistics with superior design sensitivity in matched observational studies with multiple controls. Journal of the American Statistical Association, 109, 1145-1158.
<doi:10.1080/01621459.2013.879261>  Contains the mercury example below, 397 matched triples.

Rosenbaum, P. R. and Small, D. S. (2017) An adaptive Mantel-Haenszel test for sensitivity analysis in observational studies.  Biometrics, 73, 422–430.  <doi:10.1111/biom.12591>  The 2x2x2 BRCA example from Satagopan et al. (2001) in this paper can be used to compare the senstrat() function and the mhLS() function in the sensitivity2x2xk package for the Mantel-Haenszel test.  The 2x2x2 table is in the documentation for mhLS(), but must be reformated as individual data for use by senstrat.  With binary outcomes, the extreme unobserved covariate is known from theory.  Not knowing this theory, senstrat() computes the same answer as mhLS() for gamma=7.

Rosenbaum, P. R. (2017) Sensitivity analysis for stratified comparisons in an observational study
of the effect of smoking on homocysteine levels.  Manuscript.

Satagopan, J. M., Offit, K., Foulkes, W., Robson, M. E. Wacholder, S., Eng, C. M., Karp, S. E. and Begg, C. B. (2001). The lifetime risks of breast cancer in Ashkenazi Jewish carriers of brca1 and brca2 mutations. Cancer Epidemology, Biomarkers and Prevention, 10, 467-473.

Werfel, U., Langen, V., Eickhoff, I. et al. Elevated DNA strand breakage frequencies in lymphocytes of welders exposed to chromium and nickel. Carcinogenesis, 1998, 19, 413-418.  Used in the erpcp example below.

}
\author{
Paul R. Rosenbaum
}
\note{
Strata that contain only treated subjects or only controls do not affect permutation inferences.  These strata are removed before computations begin.  Inclusion of these strata would not alter the permutation inference.  A message will indicate whether any strata have been removed; see StrataUse in the value section above.  You can avoid strata that do not contribute by using full matching in place of conventional stratification; see Rosenbaum (1991) and Hansen (2004) and R packages optmatch and sensitivityfull.

The output produces a rejection or acceptance of the null hypothesis at a stated level=\eqn{\alpha} in the presence of a bias of at most \eqn{\Gamma}.  This statement is entirely safe, in the sense that it is at worst a tad conservative, falsely rejecting a true null hypothesis with probability at most \eqn{\alpha} in the presence of a bias of at most \eqn{\Gamma}.  To produce a true P-value, you would need to run the program several times to find the smallest level=\eqn{\alpha} that leads to rejection, and the P-value produced in this standard way would share the property of the test in being, at worst, slightly conservative.  To save time, the output contains an approximate P-value that agrees with the accept/reject decision, but if this P-value is much smaller than the level -- say, rejection at 0.05 with a P-value of 0.00048, then unlike the reject/accept decison, the 0.00048 P-value may not be conservative.  I have never found the approximate P-value to be misleading.  However, having seen an approximate P-value of 0.00048, it is easy to check whether you are formally entitled to reject at \eqn{\alpha}=0.0005 by rerunning the program with level=0.0005 and basing the conclusion on the reject/accept decision at level \eqn{\alpha}=0.0005.

When there are many small strata, Gastwirth, Krieger and Rosenbaum (2000, GKR) proposed a separable approximation to the sensitivity bound.  In principle, this separable approximation is a tad liberal: it does not find the absolute worst unobserved covariate u, but rather a very bad u, such that as the number of strata increases the difference between the worst u and the very bad u becomes negligible.  The current function senstrat() improves upon the separable approximation in the following way.  This improvement is discussed in Rosenbaum (2017).  It makes a one-step linear Taylor correction to the separable approximation which is guaranteed to be slightly conservative, rather than slightly liberal, so it is always safe to use: it falsely rejects at level = 0.05 with probability at most 0.05 in the presence of a bias of at most \eqn{\Gamma}.  More precisely, unlike the method of GKR, the one-step LinearBound correction does not require many small strata: in large samples, it falsely rejects at level=0.05 with probability at most 0.05 whether there are few or many strata, even if some of the strata are much larger than others.  If detail=FALSE, conclusions are based on the LinearBound without further comments.  This is reasonable, because the LinearBound is safe to use in all cases, being at worst slightly conservative.  If detail=TRUE, the LinearBound and the Separable approximation are compared.  Usually, the LinearBound and the Separable approximation yield conclusions that are very close, providing some reassurance that the LinearBound is not very conservative and the Separable approximation is not very liberal.  The option detail=TRUE is an aid to someone who wants to understand the LinearBound, but it is not a tool required for data analysis.
}

\seealso{
If outcomes are binary, then use the sensitivity2x2xk package.  If there are no strata -- that is, if everyone is in the same stratum, so there is a single stratum -- then use the function sen2sample() in this package.  If the strata are matched pairs or matched sets, then use one of the packages sensitivitymult, sensitivityfull, sensitivitymv, sensitivitymw.
}
\examples{
data("homocyst")
attach(homocyst)
sc<-hodgeslehmann(log2(homocysteine),z,stf,align="hl")
senstrat(sc,z,stf,gamma=1.8)
# Compare this with:
senstrat(sc,z,stf,gamma=1.8,detail=TRUE)
# With detail=TRUE, it is seen that the separable and Taylor bounds
# on the maximum P-value are nearly identical.  The Taylor upper bound
# is safe -- i.e., at worst conservative -- in all cases.
detach(homocyst)
#
# The following example compares senmw in the sensitivitymw package
# to senstrat in an example with 397 matched triples, one treated,
# two controls.  We expect the separable approximation to work well
# with S=397 small strata, and indeed the results are identical.
library(sensitivitymw)
data(mercury)
senmw(mercury,gamma=15)
# Reformat mercury for use by senstrat().
z<-c(rep(1,397),rep(0,397),rep(0,397))
st<-rep(1:397,3)
y<-as.vector(as.matrix(mercury))
sc<-mscores(y,z,st=st)
senstrat(sc,z,st,gamma=15,detail=TRUE)
# The separable approximation from senmw() and senstrat() are identical,
# as they should be, and the Taylor approximation in senstrat()
# makes no adjustment to the separable approximation.
#
# The following example from the sensitivitymw package
# is for 39 matched pairs, so the separable algorithm
# and the Taylor approximation are not needed, yet
# they both provide exactly the correct answer.
library(sensitivitymw)
data(erpcp)
senmw(erpcp,gamma=3)
# Reformat erpcp for use by senstrat().
z<-c(rep(1,39),rep(0,39))
st<-rep(1:39,2)
y<-as.vector(as.matrix(erpcp))
sc<-mscores(y,z,st=st)
senstrat(sc,z,st,gamma=3,detail=TRUE)
}

\keyword{ htest }
\keyword{ robust }

\name{homocyst}
\alias{homocyst}
\docType{data}
\title{
Homocysteine levels in daily smokers and never smokers.
}
\description{
Data from NHANES 2005-2006 concerning homocysteine levels in daily smokers (z=1) and never smokers (z=0), aged 20 and older.  Daily smokers smoked every day for the last 30 days, smoking an average of at least 10 cigarettes per day.  Never smokers smoked fewer than 100 cigarettes in their lives, do not smoke now, and had no tobacco use in the previous 5 days.
}
\usage{data("homocyst")}
\format{
  A data frame with 2475 observations on the following 10 variables.
  \describe{
    \item{\code{SEQN}}{2005-2006 NHANES ID number.}
    \item{\code{homocysteine}}{Homocysteine level, umol/L.  Based on LBXHCY.}
    \item{\code{z}}{z=1 for a daily smoker, z=0 for a never smoker.  Based on SMQ020, SMQ040, SMD641, SMD650, SMQ680.}
    \item{\code{stf}}{A factor for strata indicating female, age, education, BMI and poverty.}
    \item{\code{st}}{Numeric strata indicating female, age, education, BMI and poverty.}
    \item{\code{female}}{1=female, 0=male.  Based on RIAGENDR}
    \item{\code{age3}}{Three age categories, 20-39, 40-50, >=60.  Based on RIDAGEYR.}
    \item{\code{ed3}}{Three education categories, <High School, High School, at least some College.  Based on DMDEDUC2.}
    \item{\code{bmi3}}{Three of the body-mass-index, BMI, <30, [30,35), >= 35.  Based on BMXBMI.}
    \item{\code{pov2}}{TRUE=income at least twice the poverty level, FALSE otherwise.}
  }
}
\details{
Bazzano et al. (2003) noted higher homocysteine levels in smokers than in nonsmokers.  See also Pimentel et al. (2016) for a related analysis.  The example below is from Rosenbaum (2017).
}
\source{
NHANES, the US National Health and Nutrition Examination Survey, 2005-2006.
}
\references{
Bazzano, L. A., He, J., Muntner, P., Vupputuri, S. and Whelton, P. K. (2003) Relationship between cigarette smoking and novel risk factors for cardiovascular disease in the United States.  Annals of Internal Medicine, 138, 891-897.

Pimentel, S. D., Small, D. S. and Rosenbaum, P. R. (2016) Constructed second control groups and attenuation of unmeasured biases.  Journal of the American Statistical Association, 111, 1157-1167.

Rosenbaum, P. R. (2017) Sensitivity analysis for stratified comparisons in an observational study
of the effect of smoking on homocysteine levels.  Manuscript.
}
\examples{
data(homocyst)
#Homocysteine levels for daily smokers and nonsmokers.
boxplot(log(homocyst$homocysteine)~homocyst$z)
}
\keyword{datasets}

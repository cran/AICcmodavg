\name{importance}
\Rdversion{1.1}
\alias{importance}
\alias{print.importance}
\title{
Compute Importance Values of Variable
}
\description{
This function calculates the relative importance of variables (w+) based
on the sum of Akaike weights (model probabilities) of the models that
include the variable.  Note that this measure of evidence is only
appropriate when the variable appears in the same number of models as
those that do not include the variable.
}
\usage{
importance(cand.set, parm, modnames, c.hat = 1, second.ord = TRUE,
           nobs = NULL, parm.type = NULL)  
}

\arguments{
  \item{cand.set}{
a list storing each of the models in the candidate model set.
}
  \item{parm}{
the parameter of interest for which a measure of relative importance is
  required. 
}
  \item{modnames}{
a character vector of model names to facilitate the identification of
each model in the model selection table.
}
  \item{c.hat}{
value of overdispersion parameter (i.e., variance inflation factor) such
as that obtained from 'c_hat'.  Note that values of c.hat different from
1 are only appropriate for binomial GLM's with trials > 1 (i.e.,
success/trial or cbind(success, failure) syntax), with Poisson
GLM's, or single-season occupancy models (MacKenzie et al. 2002).  If
c.hat > 1, 'AICc' will return the quasi-likelihood analogue of the
information criterion requested.  This option is not supported for
generalized linear mixed models of the 'mer' class. 
}
  \item{second.ord}{
logical. If TRUE, the function returns the second-order Akaike
  information criterion (i.e., AICc).
}
  \item{nobs}{
this argument allows to specify a numeric value other than total sample
size to compute the AICc (i.e., 'nobs' defaults to total number of
observations).  This is relevant only for mixed models or various models
of 'unmarkedFit' classes where sample size is not straightforward.  In
such cases, one might use total number of observations or number of
independent clusters (e.g., sites) as the value of 'nobs'.
}
\item{parm.type}{
  this argument specifies the parameter type on which the effect size 
  will be computed and is only relevant for models of 'unmarkedFitOccu',
  'unmarkedFitColExt', 'unmarkedFitOccuRN', 'unmarkedFitPCount',
  'unmarkedFitPCO', 'unmarkedFitDS', 'unmarkedFitGDS', and
  'unmarkedFitOccuFP' classes.  The character strings supported vary
  with the type of model fitted.  For 'unmarkedFitOccu' objects, either
  'psi' or 'detect' can be supplied to indicate whether the parameter is
  on occupancy or detectability, respectively.  For 'unmarkedFitColExt',
  possible values are 'psi', 'gamma', 'epsilon', and 'detect', for
  parameters on occupancy in the inital year, colonization, extinction,
  and detectability, respectively.  For 'unmarkedFitOccuRN' objects,
  either 'lambda or 'detect' can be entered for abundance and
  detectability parameters, respectively.  For 'unmarkedFitPCount'
  objects, 'lambda' or 'detect' denote parameters on abundance and
  detectability, respectively.  For 'unmarkedFitPCO' objects, one can
  enter 'lambda', 'gamma', 'omega', or 'detect', to specify parameters
  on abundance, recruitment, apparent survival, and detectability,
  respectively. For 'unmarkedFitDS' and 'unmarkedFitGDS' objects, only
  'lambda' is supported for the moment. For 'unmarkedFitOccuFP' objects,
  one can specify 'psi', 'detect', or 'fp', for occupancy,
  detectability, and probability of assigning false-positives,
  respectively.
}
}

\value{
'importance' returns an object of class 'importance' consisting of the
following components:
\item{parm}{the parameter for which an importance value is required.} 
\item{w.plus}{the parameter for which an importance value is required.} 
\item{w.minus}{the sum of Akaike weights for the models that exclude the
  parameter of interest} 
}
\references{
Burnham, K. P., and Anderson, D. R. (2002) \emph{Model Selection and
  Multimodel Inference: a practical information-theoretic
  approach}. Second edition. Springer: New York.

MacKenzie, D. I., Nichols, J. D., Lachman, G. B., Droege, S., Royle,
J. A., Langtimm, C. A. (2002) Estimating site occupancy rates when
detection probabilities are less than one. \emph{Ecology} \bold{83},
2248--2255.
}
\author{
Marc J. Mazerolle
}

\seealso{
\code{\link{AICc}}, \code{\link{aictab}}, \code{\link{c_hat}},
\code{\link{confset}}, \code{\link{evidence}}, \code{\link{modavg}}, 
\code{\link{modavg.shrink}}, \code{\link{modavgpred}} 
}
\examples{
##example on Orthodont data set in nlme
require(nlme)

##set up candidate model list
Cand.models <- list( )
Cand.models[[1]] <- lme(distance ~ age, data = Orthodont, method = "ML")
##random is ~ age | Subject
Cand.models[[2]] <- lme(distance ~ age + Sex, data = Orthodont,
                        random = ~ 1, method = "ML")
Cand.models[[3]] <- lme(distance ~ 1, data = Orthodont, random = ~ 1,
                        method = "ML") 
Cand.models[[4]] <- lme(distance ~ Sex, data = Orthodont, random = ~ 1,
                        method = "ML") 

##create a vector of model names
Modnames <- paste("mod", 1:length(Cand.models), sep = "")

importance(cand.set = Cand.models, parm = "age", modnames = Modnames,
           second.ord = TRUE, nobs = NULL)
##round to 4 digits after decimal point
print(importance(cand.set = Cand.models, parm = "age", modnames = Modnames,
                 second.ord = TRUE, nobs = NULL), digits = 4)



##single-season occupancy model example modified from ?occu
require(unmarked)
##single season
data(frogs)
pferUMF <- unmarkedFrameOccu(pfer.bin)
## add some fake covariates for illustration
siteCovs(pferUMF) <- data.frame(sitevar1 = rnorm(numSites(pferUMF)),
                                sitevar2 = rnorm(numSites(pferUMF))) 
     
## observation covariates are in site-major, observation-minor order
obsCovs(pferUMF) <- data.frame(obsvar1 = rnorm(numSites(pferUMF) *
                                 obsNum(pferUMF))) 

##set up candidate model set
fm1 <- occu(~ obsvar1 ~ sitevar1, pferUMF)
fm2 <- occu(~ 1 ~ sitevar1, pferUMF)
fm3 <- occu(~ obsvar1 ~ sitevar2, pferUMF)
fm4 <- occu(~ 1 ~ sitevar2, pferUMF)
Cand.mods <- list(fm1, fm2, fm3, fm4)
Modnames <- c("fm1", "fm2", "fm3", "fm4")

##compute importance value for 'sitevar1' on occupancy
importance(cand.set = Cand.mods, modnames = Modnames, parm = "sitevar1",
           parm.type = "psi")
##compute importance value for 'obsvar1' on detectability
importance(cand.set = Cand.mods, modnames = Modnames, parm = "obsvar1",
           parm.type = "detect")
detach(package:unmarked)
}
\keyword{models}

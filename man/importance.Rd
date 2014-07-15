\name{importance}
\Rdversion{1.1}
\alias{importance}
\alias{importance.default}
\alias{importance.AICaov.lm}
\alias{importance.AICsclm.clm}
\alias{importance.AICclmm}
\alias{importance.AICclogit.coxph}
\alias{importance.AICcoxme}
\alias{importance.AICcoxph}
\alias{importance.AICglm.lm}
\alias{importance.AICglmerMod}
\alias{importance.AICgls}
\alias{importance.AIClm}
\alias{importance.AIClme}
\alias{importance.AIClmekin}
\alias{importance.AICmaxlikeFit.list}
\alias{importance.AICmer}
\alias{importance.AICmultinom.nnet}
\alias{importance.AICnlmerMod}
\alias{importance.AICpolr}
\alias{importance.AICrlm.lm}
\alias{importance.AICunmarkedFitColExt}
\alias{importance.AICunmarkedFitOccu}
\alias{importance.AICunmarkedFitOccuFP}
\alias{importance.AICunmarkedFitOccuRN}
\alias{importance.AICunmarkedFitPCount}
\alias{importance.AICunmarkedFitPCO}
\alias{importance.AICunmarkedFitDS}
\alias{importance.AICunmarkedFitGDS}
\alias{importance.AICunmarkedFitMPois}
\alias{importance.AICunmarkedFitGMM}
\alias{importance.AICvglm}
\alias{importance.AICzeroinfl}
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
importance(cand.set, parm, modnames = NULL, second.ord = TRUE,
           nobs = NULL, \dots)

\method{importance}{AICaov.lm}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, ...)

\method{importance}{AICsclm.clm}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, ...)

\method{importance}{AICclmm}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, ...)

\method{importance}{AICclogit.coxph}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, ...)

\method{importance}{AICcoxme}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, ...)

\method{importance}{AICcoxph}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, ...)

\method{importance}{AICglm.lm}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, c.hat = 1, ...)

\method{importance}{AICglmerMod}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, ...)

\method{importance}{AICgls}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, ...)

\method{importance}{AIClm}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, ...)

\method{importance}{AIClme}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, ...)

\method{importance}{AIClmekin}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, ...)

\method{importance}{AICmaxlikeFit.list}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, c.hat = 1, ...)

\method{importance}{AICmer}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, ...)

\method{importance}{AICmultinom.nnet}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, c.hat = 1, ...)

\method{importance}{AICnlmerMod}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, ...)

\method{importance}{AICpolr}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, ...)

\method{importance}{AICrlm.lm}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, ...)

\method{importance}{AICunmarkedFitColExt}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, c.hat = 1, parm.type = NULL, ...)

\method{importance}{AICunmarkedFitOccu}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, c.hat = 1, parm.type = NULL, ...)

\method{importance}{AICunmarkedFitOccuFP}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, c.hat = 1, parm.type = NULL, ...)

\method{importance}{AICunmarkedFitOccuRN}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, c.hat = 1, parm.type = NULL, ...)

\method{importance}{AICunmarkedFitPCount}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, c.hat = 1, parm.type = NULL, ...)

\method{importance}{AICunmarkedFitPCO}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, c.hat = 1, parm.type = NULL, ...)

\method{importance}{AICunmarkedFitDS}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, c.hat = 1, parm.type = NULL, ...)

\method{importance}{AICunmarkedFitGDS}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, c.hat = 1, parm.type = NULL, ...)

\method{importance}{AICunmarkedFitMPois}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, c.hat = 1, parm.type = NULL, ...)

\method{importance}{AICunmarkedFitGMM}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, c.hat = 1, parm.type = NULL, ...)

\method{importance}{AICvglm}(cand.set, parm, modnames = NULL, 
           second.ord = TRUE, nobs = NULL, c.hat = 1, ...)

\method{importance}{AICzeroinfl}(cand.set, parm, modnames = NULL, 
           second.ord = TRUE, nobs = NULL, ...)
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
each model in the model selection table. If \code{NULL}, the function
uses the names in the cand.set list of candidate models.  If no names
appear in the list, generic names (e.g., \code{Mod1}, \code{Mod2}) are
supplied in the table in the same order as in the list of candidate models.
}
  \item{second.ord}{
logical. If \code{TRUE}, the function returns the second-order Akaike 
information criterion (i.e., AICc). 
}
  \item{nobs}{
this argument allows to specify a numeric value other than total sample
size to compute the AICc (i.e., \code{nobs} defaults to total number of
observations).  This is relevant only for mixed models or various models
of \code{unmarkedFit} classes where sample size is not straightforward.  In
such cases, one might use total number of observations or number of
independent clusters (e.g., sites) as the value of \code{nobs}.
}
\item{c.hat}{
  value of overdispersion parameter (i.e., variance inflation factor) such
  as that obtained from \code{c_hat}.  Note that values of \code{c.hat}
  different from 1 are only appropriate for binomial GLM's with trials >
  1 (i.e., success/trial or cbind(success, failure) syntax), with
  Poisson GLM's, or single-season occupancy models (MacKenzie et
  al. 2002).  If c.hat > 1, \code{AICc} will return the quasi-likelihood
  analogue of the information criterion requested.  This option is not
  supported for generalized linear mixed models of the \code{mer} and
  \code{merMod} classes.
}
\item{parm.type}{
  this argument specifies the parameter type on which the effect size 
  will be computed and is only relevant for models of
  \code{unmarkedFitOccu}, \code{unmarkedFitColExt},
  \code{unmarkedFitOccuFP}, \code{unmarkedFitOccuRN},
  \code{unmarkedFitMPois}, \code{unmarkedFitGPC},
  \code{unmarkedFitPCount}, \code{unmarkedFitPCO}, \code{unmarkedFitDS},
  \code{unmarkedFitGDS}, and \code{unmarkedFitGMM} classes.  The
  character strings supported vary with the type of model fitted.  For
  \code{unmarkedFitOccu} objects, either \code{psi} or \code{detect} can
  be supplied to indicate whether the parameter is on occupancy or
  detectability, respectively.  For \code{unmarkedFitColExt}, possible
  values are \code{psi}, \code{gamma}, \code{epsilon}, and
  \code{detect}, for parameters on occupancy in the inital year,
  colonization, extinction, and detectability, respectively.  For
  \code{unmarkedFitOccuFP} objects, one can specify \code{psi},
  \code{detect}, or \code{fp}, for occupancy, detectability, and
  probability of assigning false-positives, respectively.  For
  \code{unmarkedFitOccuRN} objects, either \code{lambda} or
  \code{detect} can be entered for abundance and detectability
  parameters, respectively.  For \code{unmarkedFitPCount} and
  \code{unmarkedFitMPois} objects, \code{lambda} or \code{detect} denote
  parameters on abundance and detectability, respectively.  For
  \code{unmarkedFitPCO} objects, one can enter \code{lambda}, 
  \code{gamma}, \code{omega}, or \code{detect}, to specify parameters on
  abundance, recruitment, apparent survival, and detectability,
  respectively.  For \code{unmarkedFitDS} objects, only \code{lambda} is
  supported for the moment. For \code{unmarkedFitGDS} objects, \code{lambda}
  and \code{phi} denote abundance and availability, respectively.  For
  \code{unmarkedFitGMM} and \code{unmarkedFitGPC} objects,
  \code{lambda}, \code{phi}, and \code{detect} denote abundance,
  availability, and detectability, respectively.
}
\item{\dots}{
  additional arguments passed to the function.
}
}

\value{
\code{importance} returns an object of class \code{importance}
consisting of the following components:
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
\code{\link{modavgShrink}}, \code{\link{modavgPred}} 
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
if(require(unmarked)) {
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
}
\keyword{models}

\name{modavgShrink}
\Rdversion{1.1}
\alias{modavgShrink}
\alias{modavgShrink.default}
\alias{modavgShrink.AICaov.lm}
\alias{modavgShrink.AICsclm.clm}
\alias{modavgShrink.AICclmm}
\alias{modavgShrink.AICcoxme}
\alias{modavgShrink.AICcoxph}
\alias{modavgShrink.AICglm.lm}
\alias{modavgShrink.AICgls}
\alias{modavgShrink.AIClm}
\alias{modavgShrink.AIClme}
\alias{modavgShrink.AIClmekin}
\alias{modavgShrink.AICmaxlikeFit.list}
\alias{modavgShrink.AICmer}
\alias{modavgShrink.AICglmerMod}
\alias{modavgShrink.AIClmerMod}
\alias{modavgShrink.AICmultinom.nnet}
\alias{modavgShrink.AICpolr}
\alias{modavgShrink.AICrlm.lm}
\alias{modavgShrink.AICvglm}
\alias{modavgShrink.AICzeroinfl}
\alias{modavgShrink.AICunmarkedFitOccu}
\alias{modavgShrink.AICunmarkedFitColExt}
\alias{modavgShrink.AICunmarkedFitOccuRN}
\alias{modavgShrink.AICunmarkedFitPCount}
\alias{modavgShrink.AICunmarkedFitPCO}
\alias{modavgShrink.AICunmarkedFitDS}
\alias{modavgShrink.AICunmarkedFitGDS}
\alias{modavgShrink.AICunmarkedFitOccuFP}
\alias{modavgShrink.AICunmarkedFitMPois}
\alias{modavgShrink.AICunmarkedFitFitGMM}
\alias{modavgShrink.AICunmarkedFitFitGPC}
\alias{print.modavgShrink}

\title{
Compute Model-averaged Parameter Estimate with Shrinkage (Multimodel Inference)
}
\description{
This function computes an alternative version of model-averaging
parameter estimates that consists in shrinking estimates toward 0 to
reduce model selection bias as in Burnham and Anderson (2002, p. 152),
Anderson (2008, pp. 130-132) and Lukacs et al. (2010).  Specifically,
models without the parameter of interest have an estimate and variance
of 0.  \code{modavgShrink} also returns unconditional standard errors
and unconditional confidence intervals as described in Buckland et
al. (1997) and Burnham and Anderson (2002). 
}
\usage{
modavgShrink(cand.set, parm, modnames = NULL, second.ord = TRUE,
              nobs = NULL, uncond.se = "revised", conf.level = 0.95,
              \dots)
\method{modavgShrink}{AICaov.lm}(cand.set, parm, modnames = NULL,
        second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, \dots)

\method{modavgShrink}{AICsclm.clm}(cand.set, parm, modnames = NULL,
        second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, \dots)

\method{modavgShrink}{AICclmm}(cand.set, parm, modnames = NULL,
        second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, \dots)

\method{modavgShrink}{AICcoxph}(cand.set, parm, modnames = NULL,
        second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, \dots)

\method{modavgShrink}{AICglm.lm}(cand.set, parm, modnames = NULL,
        second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, c.hat = 1, gamdisp = NULL, \dots)

\method{modavgShrink}{AICgls}(cand.set, parm, modnames = NULL,
        second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, \dots)

\method{modavgShrink}{AIClm}(cand.set, parm, modnames = NULL,
        second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, \dots)

\method{modavgShrink}{AIClme}(cand.set, parm, modnames = NULL,
        second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, \dots)

\method{modavgShrink}{AIClmekin}(cand.set, parm, modnames = NULL,
        second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, \dots)

\method{modavgShrink}{AICmer}(cand.set, parm, modnames = NULL,
        second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, \dots)

\method{modavgShrink}{AICglmerMod}(cand.set, parm, modnames = NULL,
        second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, \dots)

\method{modavgShrink}{AIClmerMod}(cand.set, parm, modnames = NULL,
        second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, \dots)

\method{modavgShrink}{AICmaxlikeFit.list}(cand.set, parm, modnames = NULL,
        second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, c.hat = 1, \dots)

\method{modavgShrink}{AICmultinom.nnet}(cand.set, parm, modnames =
        NULL, second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, c.hat = 1, \dots)

\method{modavgShrink}{AICpolr}(cand.set, parm, modnames = NULL,
        second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, \dots)

\method{modavgShrink}{AICrlm.lm}(cand.set, parm, modnames = NULL,
        second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, \dots)

\method{modavgShrink}{AICvglm}(cand.set, parm, modnames = NULL,
        second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, c.hat = 1, \dots)

\method{modavgShrink}{AICzeroinfl}(cand.set, parm, modnames = NULL,
        second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, \dots)

\method{modavgShrink}{AICunmarkedFitOccu}(cand.set, parm, modnames =
        NULL, second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, c.hat = 1, parm.type = NULL, \dots)

\method{modavgShrink}{AICunmarkedFitColExt}(cand.set, parm, modnames
           = NULL, second.ord = TRUE, nobs = NULL, uncond.se =
           "revised", conf.level = 0.95, c.hat = 1, parm.type = NULL,
           \dots)

\method{modavgShrink}{AICunmarkedFitOccuRN}(cand.set, parm, modnames
        = NULL, second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, c.hat = 1, parm.type = NULL, \dots)

\method{modavgShrink}{AICunmarkedFitPCount}(cand.set, parm, modnames
        = NULL, second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, c.hat = 1, parm.type = NULL, \dots)

\method{modavgShrink}{AICunmarkedFitPCO}(cand.set, parm, modnames =
        NULL, second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, c.hat = 1, parm.type = NULL, \dots)

\method{modavgShrink}{AICunmarkedFitDS}(cand.set, parm, modnames =
        NULL, second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, c.hat = 1, parm.type = NULL, \dots)

\method{modavgShrink}{AICunmarkedFitGDS}(cand.set, parm, modnames =
        NULL, second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, c.hat = 1, parm.type = NULL, \dots)

\method{modavgShrink}{AICunmarkedFitOccuFP}(cand.set, parm, modnames
        = NULL, second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, c.hat = 1, parm.type = NULL, \dots)

\method{modavgShrink}{AICunmarkedFitMPois}(cand.set, parm, modnames
        = NULL, second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, c.hat = 1, parm.type = NULL, \dots)

\method{modavgShrink}{AICunmarkedFitFitGMM}(cand.set, parm, modnames
        = NULL, second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, c.hat = 1, parm.type = NULL, \dots)

\method{modavgShrink}{AICunmarkedFitFitGPC}(cand.set, parm, modnames
        = NULL, second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, c.hat = 1, parm.type = NULL, \dots)
}


%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cand.set}{
a list storing each of the models in the candidate model set.
}
  \item{parm}{
the parameter of interest, enclosed between quotes, for which a
model-averaged estimate is required.  For a categorical variable, 
the label of the estimate must be included as it appears in the output
(see 'Details' below). 
}
  \item{modnames}{
a character vector of model names to facilitate the identification of
each model in the model selection table. If \code{NULL}, the function
uses the names in the cand.set list of candidate models. If no names
appear in the list, generic names (e.g., \code{Mod1}, \code{Mod2}) are
supplied in the table in the same order as in the list of candidate
models.
}
\item{second.ord}{
  logical.  If \code{TRUE}, the function returns the second-order Akaike
  information criterion (i.e., AICc).
}
\item{nobs}{
  this argument allows to specify a numeric value other than total
  sample size to compute the AICc (i.e., \code{nobs} defaults to total
  number of observations).  This is relevant only for mixed models or
  various models of \code{unmarkedFit} classes where sample size is not
  straightforward.  In such cases, one might use total number of
  observations or number of independent clusters (e.g., sites) as the
  value of \code{nobs}.
}
\item{uncond.se}{
  either, \code{"old"}, or \code{"revised"}, specifying the equation
  used to compute the unconditional standard error of a model-averaged
  estimate.  With \code{uncond.se = "old"}, computations are based on
  equation 4.9 of Burnham and Anderson (2002), which was the former way
  to compute unconditional standard errors.  With \code{uncond.se =
    "revised"}, equation 6.12 of Burnham and Anderson (2002) is used.
  Anderson (2008, p. 111) recommends use of the revised version for the
  computation of unconditional standard errors and it is now the
  default.  Note that versions of package AICcmodavg < 1.04 used the old
  method to compute unconditional standard errors.  
}
  \item{conf.level}{
the confidence level requested for the computation of unconditional
confidence intervals.
}
  \item{c.hat}{
value of overdispersion parameter (i.e., variance inflation factor) such
as that obtained from \code{c_hat}.  Note that values of c.hat different
from 1 are only appropriate for binomial GLM's with trials > 1 (i.e.,
success/trial or cbind(success, failure) syntax), with Poisson GLM's, 
single-season occupancy models (MacKenzie et al. 2002), or
\emph{N}-mixture models (Royle 2004, Dail and Madsen 2011).  If
\code{c.hat} > 1, \code{modavgShrink} will return the quasi-likelihood
analogue of the information criteria requested and multiply the
variance-covariance matrix of the estimates by this value (i.e., SE's
are multiplied by \code{sqrt(c.hat)}).  This option is not supported for
generalized linear mixed models of the \code{mer} or \code{merMod}
classes.
}
  \item{gamdisp}{
if gamma GLM is used, the dispersion parameter should be specified here
to apply the same value to each model.
}
\item{parm.type}{
   this argument specifies the parameter type on which the effect size
  will be computed and is only relevant for models of
  \code{unmarkedFitOccu}, \code{unmarkedFitColExt},
  \code{unmarkedFitOccuFP}, \code{unmarkedFitOccuRN}, 
  \code{unmarkedFitMPois}, \code{unmarkedFitPCount},
  \code{unmarkedFitPCO}, \code{unmarkedFitDS}, \code{unmarkedFitGDS},
  \code{unmarkedFitGMM}, and \code{unmarkedFitGPC} classes.  The
  character strings supported vary with the type of model fitted.  For
  \code{unmarkedFitOccu} objects, either \code{psi} or \code{detect} can
  be supplied to indicate whether the parameter is on occupancy or
  detectability, respectively.  For \code{unmarkedFitColExt}, possible
  values  are \code{psi}, \code{gamma}, \code{epsilon}, and
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
  supported for the moment.  For \code{unmarkedFitGDS}, \code{lambda}
  and \code{phi} denote abundance and availability, respectively.  For
  \code{unmarkedFitGMM} and \code{unmarkedFitGPC} objects,
  \code{lambda}, \code{phi}, and \code{detect} denote abundance,
  availability, and detectability, respectively.
}
\item{\dots}{
  additional arguments passed to the function.
}
}

\details{
The parameter for which a model-averaged estimate is requested must be
specified with the \code{parm} argument and must be identical to its
label in the model output (e.g., from \code{summary}).  For factors, one
must specify the name of the variable and the level of interest.  The
shrinkage version of model averaging is only appropriate for cases where
each parameter is given an equal weighting in the model (i.e., each
parameter must appear the same number of times in the models) and has
the same interpretation across all models.  As a result, models with
interaction terms or polynomial terms are not supported by
\code{modavgShrink}. 

\code{modavgShrink} is implemented for a list containing objects of
\code{clm}, \code{clmm}, \code{clogit}, \code{coxme}, \code{coxph},
\code{glm}, \code{gls}, \code{lm}, \code{lme}, \code{lmekin},
\code{maxlikeFit}, \code{mer}, \code{glmerMod}, \code{lmerMod},
\code{multinom}, \code{polr}, \code{rlm}, \code{vglm}, \code{zeroinfl}
classes as well as various models of \code{unmarkedFit} classes.
}
\value{
\code{modavgShrink} creates an object of class \code{modavgShrink}
with the following components:

\item{Parameter}{the parameter for which a model-averaged estimate with
  shrinkage was obtained} 
\item{Mod.avg.table}{the model selection table based on models including
  the parameter of interest} 
\item{Mod.avg.beta}{the model-averaged estimate based on all models}
\item{Uncond.SE}{the unconditional standard error for the model-averaged
  estimate (as opposed to the conditional SE based on a single model)}
\item{Conf.level}{the confidence level used to compute the confidence
  interval} 
\item{Lower.CL}{the lower confidence limit}
\item{Upper.CL}{the upper confidence limit}

}
\references{
Anderson, D. R. (2008) \emph{Model-based Inference in the Life Sciences:
a primer on evidence}. Springer: New York.

Buckland, S. T., Burnham, K. P., Augustin, N. H. (1997) Model selection:
an integral part of inference. \emph{Biometrics} \bold{53}, 603--618.  

Burnham, K. P., Anderson, D. R. (2002) \emph{Model Selection and
Multimodel Inference: a practical information-theoretic
approach}. Second edition. Springer: New York.

Burnham, K. P., Anderson, D. R. (2004) Multimodel inference:
understanding AIC and BIC in model selection. \emph{Sociological
Methods and Research} \bold{33}, 261--304.

Dail, D., Madsen, L. (2011) Models for estimating abundance from 
repeated counts of an open population. \emph{Biometrics} \bold{67},
577--587. 

Lukacs, P. M., Burnham, K. P., Anderson, D. R. (2010) Model selection
bias and Freedman's paradox. \emph{Annals of the Institute of
Statistical Mathematics} \bold{62}, 117--125. 

MacKenzie, D. I., Nichols, J. D., Lachman, G. B., Droege, S., Royle,
J. A., Langtimm, C. A. (2002) Estimating site occupancy rates when
detection probabilities are less than one. \emph{Ecology} \bold{83},
2248--2255.

Mazerolle, M. J. (2006) Improving data analysis in herpetology: using
Akaike's Information Criterion (AIC) to assess the strength of
biological hypotheses. \emph{Amphibia-Reptilia} \bold{27}, 169--180. 

Royle, J. A. (2004) \emph{N}-mixture models for estimating population
size from spatially replicated counts. \emph{Biometrics} \bold{60},
108--115.
}
\author{
Marc J. Mazerolle
}

\seealso{
\code{\link{AICc}}, \code{\link{aictab}}, \code{\link{c_hat}},
\code{\link{importance}}, \code{\link{confset}}, \code{\link{evidence}},
\code{\link{modavg}}, \code{\link{modavgPred}}
}
\examples{
##cement example in Burnham and Anderson 2002
data(cement)
##setup same model set as in Table 3.2, p. 102         
Cand.models <- list( )
Cand.models[[1]] <- lm(y ~ x1 + x2, data = cement)
Cand.models[[2]] <- lm(y ~ x1 + x2 + x4, data = cement)          
Cand.models[[3]] <- lm(y ~ x1 + x2 + x3, data = cement)
Cand.models[[4]] <- lm(y ~ x1 + x4, data = cement)
Cand.models[[5]] <- lm(y ~ x1 + x3 + x4, data = cement)
Cand.models[[6]] <- lm(y ~ x2 + x3 + x4, data = cement)
Cand.models[[7]] <- lm(y ~ x1 + x2 + x3 + x4, data = cement)
Cand.models[[8]] <- lm(y ~ x3 + x4, data = cement)
Cand.models[[9]] <- lm(y ~ x2 + x3, data = cement)
Cand.models[[10]] <- lm(y ~ x4, data = cement)
Cand.models[[11]] <- lm(y ~ x2, data = cement)
Cand.models[[12]] <- lm(y ~ x2 + x4, data = cement)
Cand.models[[13]] <- lm(y ~ x1, data = cement)
Cand.models[[14]] <- lm(y ~ x1 + x3, data = cement)
Cand.models[[15]] <- lm(y ~ x3, data = cement)

##vector of model names
Modnames <- paste("mod", 1:15, sep="")

##AICc          
aictab(cand.set = Cand.models, modnames = Modnames)

##compute model-averaged estimate with shrinkage - each parameter
##appears 8 times in the models 
modavgShrink(cand.set = Cand.models, modnames = Modnames, parm = "x1")

##compare against classic model-averaging
modavg(cand.set = Cand.models, modnames = Modnames, parm = "x1")
##note that model-averaged estimate with shrinkage is closer to 0 than
##with the classic version

##remove a few models from the set and run again
Cand.unbalanced <- Cand.models[-c(3, 14, 15)]

##set up model names
Modnames <- paste("mod", 1:length(Cand.unbalanced), sep="")

##issues an error because some parameters appear more often than others
\dontrun{modavgShrink(cand.set = Cand.unbalanced,
                       modnames = Modnames, parm = "x1")}



##example on Orthodont data set in nlme
require(nlme)

##set up candidate model list
##age and sex parameters appear in the same number of models
##same number of models with and without these parameters
Cand.models <- list( )
Cand.models[[1]] <- lme(distance ~ age, data = Orthodont, method = "ML") 
##random is ~ age | Subject as it is a grouped data frame
Cand.models[[2]] <- lme(distance ~ age + Sex, data = Orthodont,
                        random = ~ 1, method = "ML")
Cand.models[[3]] <- lme(distance ~ 1, data = Orthodont, random = ~ 1, 
                        method = "ML") 
Cand.models[[4]] <- lme(distance ~ Sex, data = Orthodont, random = ~ 1,
                        method = "ML")  

##create a vector of model names
Modnames <- paste("mod", 1:length(Cand.models), sep = "")

##compute importance values for age
imp.age <- importance(cand.set = Cand.models, parm = "age",
                      modnames = Modnames, second.ord = TRUE,
                      nobs = NULL)

##compute shrinkage version of model averaging on age
mod.avg.age.shrink <- modavgShrink(cand.set = Cand.models,
                                    parm = "age", modnames = Modnames,
                                    second.ord = TRUE, nobs = NULL)

##compute classic version of model averaging on age
mod.avg.age.classic <- modavg(cand.set = Cand.models, parm = "age",
                              modnames = Modnames, second.ord = TRUE,
                              nobs = NULL)

##correspondence between shrinkage version and classic version of
##model averaging 
mod.avg.age.shrink$Mod.avg.beta/imp.age$w.plus
mod.avg.age.classic$Mod.avg.beta
detach(package:nlme)



\dontrun{
##example of N-mixture model modified from ?pcount
if(require(unmarked)) {
data(mallard)
mallardUMF <- unmarkedFramePCount(mallard.y, siteCovs = mallard.site,
                                  obsCovs = mallard.obs)
##set up models so that each variable on abundance appears twice
fm.mall.one <- pcount(~ ivel + date  ~ length + forest, mallardUMF,
                      K = 30)
fm.mall.two <- pcount(~ ivel + date  ~ elev + forest, mallardUMF,
                      K = 30)
fm.mall.three <- pcount(~ ivel + date  ~ length + elev, mallardUMF,
                        K = 30)

##model list and names
Cands <- list(fm.mall.one, fm.mall.two, fm.mall.three)
Modnames <- c("length + forest", "elev + forest", "length + elev")

##compute model-averaged estimate with shrinkage for elev on abundance
modavgShrink(cand.set = Cands, modnames = Modnames, parm = "elev",
              parm.type = "lambda")
detach(package:unmarked)
}
}

}
\keyword{models}

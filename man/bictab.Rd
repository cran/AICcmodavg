\name{bictab}
\Rdversion{1.1}
\alias{bictab}
\alias{bictab.default}
\alias{bictab.AICaov.lm}
\alias{bictab.AICbetareg}
\alias{bictab.AICsclm.clm}
\alias{bictab.AICclm}
\alias{bictab.AICclmm}
\alias{bictab.AICcoxme}
\alias{bictab.AICcoxph}
\alias{bictab.AICfitdist}
\alias{bictab.AICfitdistr}
\alias{bictab.AICglm.lm}
\alias{bictab.AICglmmTMB}
\alias{bictab.AICgls}
\alias{bictab.AICgnls.gls}
\alias{bictab.AIChurdle}
\alias{bictab.AIClavaan}
\alias{bictab.AIClm}
\alias{bictab.AIClme}
\alias{bictab.AIClmekin}
\alias{bictab.AICmaxlikeFit.list}
\alias{bictab.AICmer}
\alias{bictab.AIClmerMod}
\alias{bictab.AIClmerModLmerTest}
\alias{bictab.AICglmerMod}
\alias{bictab.AICnlmerMod}
\alias{bictab.AICmultinom.nnet}
\alias{bictab.AICnlme.lme}
\alias{bictab.AICnls}
\alias{bictab.AICpolr}
\alias{bictab.AICrlm.lm}
\alias{bictab.AICsurvreg}
\alias{bictab.AICunmarkedFitOccu}
\alias{bictab.AICunmarkedFitColExt}
\alias{bictab.AICunmarkedFitOccuRN}
\alias{bictab.AICunmarkedFitPCount}
\alias{bictab.AICunmarkedFitPCO}
\alias{bictab.AICunmarkedFitDS}
\alias{bictab.AICunmarkedFitGDS}
\alias{bictab.AICunmarkedFitOccuFP}
\alias{bictab.AICunmarkedFitMPois}
\alias{bictab.AICunmarkedFitGMM}
\alias{bictab.AICunmarkedFitGPC}
\alias{bictab.AICunmarkedFitOccuMulti}
\alias{bictab.AICunmarkedFitOccuMS}
\alias{bictab.AICunmarkedFitOccuTTD}
\alias{bictab.AICunmarkedFitMMO}
\alias{bictab.AICunmarkedFitDSO}
\alias{bictab.AICunmarkedFitGOccu}
\alias{bictab.AICunmarkedFitOccuComm}
\alias{bictab.AICvglm}
\alias{bictab.AICzeroinfl}
\alias{print.bictab}
\title{
Create Model Selection Tables Based on BIC
}
\description{
This function creates a model selection table based on the Bayesian
information criterion (Schwarz 1978, Burnham and Anderson 2002).  The
table ranks the models based on the BIC and also provides delta BIC and
BIC model weights.  The function adjusts for overdispersion in model
selection  by using the QBIC when \code{c.hat > 1}.  \code{bictab}
selects the appropriate function to create the model selection table
based on the object class.  The current version works with lists
containing objects of \code{aov}, \code{betareg}, \code{clm},
\code{clmm}, \code{clogit}, \code{coxme}, \code{coxph}, \code{fitdist},
\code{fitdistr}, \code{glm}, \code{glmmTMB}, \code{gls}, \code{gnls},
\code{hurdle}, \code{lavaan}, \code{lm}, \code{lme}, \code{lmekin},
\code{maxlikeFit}, \code{mer}, \code{merMod}, \code{lmerModLmerTest},
\code{multinom}, \code{nlme}, \code{nls}, \code{polr}, \code{rlm},
\code{survreg}, \code{vglm}, and \code{zeroinfl} classes as well as
various models of \code{unmarkedFit} classes but does not yet allow
mixing of different classes.
}

\usage{
bictab(cand.set, modnames = NULL, nobs = NULL,
       sort = TRUE, \dots)

\method{bictab}{AICaov.lm}(cand.set, modnames = NULL,
        nobs = NULL, sort = TRUE, \dots)

\method{bictab}{AICbetareg}(cand.set, modnames = NULL,
        nobs = NULL, sort = TRUE, \dots)

\method{bictab}{AICsclm.clm}(cand.set, modnames = NULL,
        nobs = NULL, sort = TRUE, \dots)

\method{bictab}{AICclmm}(cand.set, modnames = NULL,
        nobs = NULL, sort = TRUE, \dots)

\method{bictab}{AICclm}(cand.set, modnames = NULL,
        nobs = NULL, sort = TRUE, \dots)

\method{bictab}{AICcoxme}(cand.set, modnames = NULL,
        nobs = NULL, sort = TRUE, \dots)

\method{bictab}{AICcoxph}(cand.set, modnames = NULL,
        nobs = NULL, sort = TRUE, \dots)

\method{bictab}{AICfitdist}(cand.set, modnames = NULL,
        nobs = NULL, sort = TRUE, \dots)

\method{bictab}{AICfitdistr}(cand.set, modnames = NULL,
        nobs = NULL, sort = TRUE, \dots)

\method{bictab}{AICglm.lm}(cand.set, modnames = NULL,
        nobs = NULL, sort = TRUE, c.hat = 1, \dots)  

\method{bictab}{AICglmmTMB}(cand.set, modnames = NULL,
        nobs = NULL, sort = TRUE, c.hat = 1, \dots)  

\method{bictab}{AICgls}(cand.set, modnames = NULL,
        nobs = NULL, sort = TRUE, \dots)

\method{bictab}{AICgnls.gls}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, \dots)

\method{bictab}{AIChurdle}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, \dots)

\method{bictab}{AIClavaan}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, \dots)

\method{bictab}{AIClm}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, \dots)

\method{bictab}{AIClme}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, \dots)

\method{bictab}{AIClmekin}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, \dots)

\method{bictab}{AICmaxlikeFit.list}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, c.hat = 1, \dots) 

\method{bictab}{AICmer}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, \dots)

\method{bictab}{AIClmerMod}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, \dots)

\method{bictab}{AIClmerModLmerTest}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, \dots)

\method{bictab}{AICglmerMod}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, \dots) 

\method{bictab}{AICnlmerMod}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, \dots)

\method{bictab}{AICmultinom.nnet}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, c.hat = 1, \dots)

\method{bictab}{AICnlme.lme}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, \dots) 

\method{bictab}{AICnls}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, \dots)

\method{bictab}{AICpolr}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, \dots)

\method{bictab}{AICrlm.lm}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, \dots)

\method{bictab}{AICsurvreg}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, \dots)

\method{bictab}{AICunmarkedFitOccu}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, c.hat = 1, \dots)

\method{bictab}{AICunmarkedFitColExt}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, c.hat = 1, \dots)

\method{bictab}{AICunmarkedFitOccuRN}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, c.hat = 1, \dots)

\method{bictab}{AICunmarkedFitPCount}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, c.hat = 1, \dots)

\method{bictab}{AICunmarkedFitPCO}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, c.hat = 1, \dots)

\method{bictab}{AICunmarkedFitDS}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, c.hat = 1, \dots)

\method{bictab}{AICunmarkedFitGDS}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, c.hat = 1, \dots)

\method{bictab}{AICunmarkedFitOccuFP}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, c.hat = 1, \dots)

\method{bictab}{AICunmarkedFitMPois}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, c.hat = 1, \dots)

\method{bictab}{AICunmarkedFitGMM}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, c.hat = 1, \dots)

\method{bictab}{AICunmarkedFitGPC}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, c.hat = 1, \dots)

\method{bictab}{AICunmarkedFitOccuMulti}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, c.hat = 1, \dots)

\method{bictab}{AICunmarkedFitOccuMS}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, c.hat = 1, \dots)

\method{bictab}{AICunmarkedFitOccuTTD}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, c.hat = 1, \dots)

\method{bictab}{AICunmarkedFitMMO}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, c.hat = 1, \dots)

\method{bictab}{AICunmarkedFitDSO}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, c.hat = 1, \dots)

\method{bictab}{AICunmarkedFitGOccu}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, c.hat = 1, \dots)

\method{bictab}{AICunmarkedFitOccuComm}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, c.hat = 1, \dots)

\method{bictab}{AICvglm}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, c.hat = 1, \dots)

\method{bictab}{AICzeroinfl}(cand.set, modnames = NULL,
         nobs = NULL, sort = TRUE, \dots)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cand.set}{
a list storing each of the models in the candidate model set.
}
  \item{modnames}{
a character vector of model names to facilitate the identification of
each model in the model selection table. If \code{NULL}, the function
uses the names in the cand.set list of candidate models (i.e., a named
list).  If no names appear in the list and no character vector is
provided, generic names (e.g., \code{Mod1}, \code{Mod2}) are supplied in
the table in the same order as in the list of candidate models.
}
\item{nobs}{
this argument allows to specify a numeric value other than total sample
size to compute the BIC (i.e., \code{nobs} defaults to total number of
observations).  This is relevant only for mixed models or various models
of \code{unmarkedFit} classes where sample size is not straightforward.
In such cases, one might use total number of observations or number of
independent clusters (e.g., sites) as the value of \code{nobs}.
}
\item{sort}{
logical.  If \code{TRUE}, the model selection table is ranked according
to the BIC values.
}

\item{c.hat}{
value of overdispersion parameter (i.e., variance inflation factor) such
as that obtained from \code{c_hat}.  Note that values of c.hat different
from 1 are only appropriate for binomial GLM's with trials > 1 (i.e.,
success/trial or cbind(success, failure) syntax), with Poisson GLM's, 
single-season occupancy models (MacKenzie et al. 2002), dynamic
occupancy models (MacKenzie et al. 2003), or \emph{N}-mixture models
(Royle 2004, Dail and Madsen 2011).  If \code{c.hat} > 1, \code{bictab}
will return the quasi-likelihood analogue of the BIC (QBIC) and multiply
the variance-covariance matrix of the estimates by this value (i.e.,
SE's are multiplied by \code{sqrt(c.hat)}).  This option is not
supported for generalized linear mixed models of the \code{mer} or
\code{merMod} classes.
}
\item{\dots}{
  additional arguments passed to the function.
}
}
\details{
BIC tends to favor simpler models than AIC whenever \emph{n} > 8
(Schwarz 1978, Link and Barker 2006, Anderson 2008).  BIC assigns
uniform prior probabilities across all models (i.e., equal 1/R), whereas
in AIC and AICc, prior probabilities increase with sample size (Burnham
and Anderson 2004, Link and Barker 2010).  Some authors argue that BIC
requires the true model to be included in the model set, whereas AIC or
AICc does not (Burnham and Anderson 2002).  However, Link and Barker
(2006, 2010) consider both as assuming that a model in the model set
approximates truth. 

\code{bictab} internally creates a new class for the \code{cand.set}
list of candidate models, according to the contents of the list.  The
current function is implemented for \code{clogit}, \code{coxme},
\code{coxph}, \code{fitdist}, \code{fitdistr}, \code{glm},
\code{glmmTMB}, \code{gls}, \code{gnls}, \code{hurdle}, \code{lavaan},
\code{lm}, \code{lme}, \code{lmekin}, \code{maxlikeFit}, \code{mer},
\code{merMod}, \code{lmerModLmerTest}, \code{multinom}, \code{nlme},
\code{nls}, \code{polr}, \code{rlm}, \code{survreg}, \code{vglm}, and
\code{zeroinfl} classes as well as various \code{unmarkedFit} classes.
The function constructs a model selection table based on BIC.
}
\value{
\code{bictab} creates an object of class \code{bictab} with the
following components:

\item{Modname}{the name of each model of the candidate model set.}
\item{K }{the number of estimated parameters for each model.}
\item{(Q)BIC}{the Bayesian information criterion for each model.}
\item{Delta_(Q)BIC}{the delta BIC component.}
\item{ModelLik }{the relative likelihood of the model given the
  data (exp(-0.5*delta[i])).  This is not to be confused with the
  likelihood of the parameters given the data.  The relative likelihood
  can then be normalized across all models to get the model probabilities.} 
\item{(Q)BICWt}{the BIC model weights, also termed "model probabilities"
  (Burnham and Anderson 2002, Link and Barker 2006, Anderson 2008).
  These measures indicate the level of support (i.e., weight of
  evidence) in favor of any given model being the most parsimonious
  among the candidate model set.} 
\item{Cum.Wt }{the cumulative BIC weights.  These are only meaningful
if results in table are sorted in decreasing order of BIC weights (i.e.,
\code{sort = TRUE}).}
\item{c.hat}{if c.hat was specified as an argument, it is included in
  the table.}
\item{LL}{the log-likelihood of each model.}
\item{Res.LL}{if parameters are estimated by restricted
  maximum-likelihood (REML), the restricted log-likelihood of each
  model.} 
}
\references{
Anderson, D. R. (2008) \emph{Model-based Inference in the Life Sciences:
a primer on evidence}. Springer: New York.

Burnham, K. P., Anderson, D. R. (2002) \emph{Model Selection and
Multimodel Inference: a practical information-theoretic
approach}. Second edition. Springer: New York.

Burnham, K. P., Anderson, D. R. (2004) Multimodel inference:
understanding AIC and BIC in model selection. \emph{Sociological Methods
and Research} \bold{33}, 261--304.

Dail, D., Madsen, L. (2011) Models for estimating abundance from 
repeated counts of an open population. \emph{Biometrics} \bold{67},
577--587. 

Link, W. A., Barker, R. J. (2006) Model weights and the foundations of
multimodel inference. \emph{Ecology} \bold{87}, 2626--2635.

Link, W. A., Barker, R. J. (2010) Bayesian Inference with Ecological
Applications. Academic Press: Boston.  

MacKenzie, D. I., Nichols, J. D., Lachman, G. B., Droege, S., Royle,
J. A., Langtimm, C. A. (2002) Estimating site occupancy rates when
detection probabilities are less than one. \emph{Ecology} \bold{83},
2248--2255.

MacKenzie, D. I., Nichols, J. D., Hines, J. E., Knutson, M. G.,
Franklin, A. B. (2003) Estimating site occupancy, colonization, and
local extinction when a species is detected imperfectly. \emph{Ecology}
\bold{84}, 2200--2207.

Royle, J. A. (2004) \emph{N}-mixture models for estimating population
size from spatially replicated counts. \emph{Biometrics} \bold{60},
108--115.

Schwarz, G. (1978) Estimating the dimension of a model. \emph{Annals of
  Statistics} \bold{6}, 461--464.

}
\author{
Marc J. Mazerolle
}
\seealso{
\code{\link{aictab}}, \code{\link{bictabCustom}}, \code{\link{confset}},
\code{\link{evidence}}, \code{\link{importance}}, \code{\link{useBIC}},
}
\examples{
##Mazerolle (2006) frog water loss example
data(dry.frog)

##setup a subset of models of Table 1
Cand.models <- list( )
Cand.models[[1]] <- lm(log_Mass_lost ~ Shade + Substrate +
                       cent_Initial_mass + Initial_mass2,
                       data = dry.frog)
Cand.models[[2]] <- lm(log_Mass_lost ~ Shade + Substrate +
                       cent_Initial_mass + Initial_mass2 +
                       Shade:Substrate, data = dry.frog)
Cand.models[[3]] <- lm(log_Mass_lost ~ cent_Initial_mass +
                       Initial_mass2, data = dry.frog)
Cand.models[[4]] <- lm(log_Mass_lost ~ Shade + cent_Initial_mass +
                       Initial_mass2, data = dry.frog)
Cand.models[[5]] <- lm(log_Mass_lost ~ Substrate + cent_Initial_mass +
                       Initial_mass2, data = dry.frog)

##create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(Cand.models), sep = " ")

##generate BIC table
bictab(cand.set = Cand.models, modnames = Modnames, sort = TRUE)
##round to 4 digits after decimal point and give log-likelihood
print(bictab(cand.set = Cand.models, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)



\dontrun{
##Burnham and Anderson (2002) flour beetle data
data(beetle)
##models as suggested by Burnham and Anderson p. 198          
Cand.set <- list( )
Cand.set[[1]] <- glm(Mortality_rate ~ Dose, family =
                     binomial(link = "logit"), weights = Number_tested,
                     data = beetle)
Cand.set[[2]] <- glm(Mortality_rate ~ Dose, family =
                     binomial(link = "probit"), weights = Number_tested,
                     data = beetle)
Cand.set[[3]] <- glm(Mortality_rate ~ Dose, family =
                     binomial(link ="cloglog"), weights = Number_tested,
                     data = beetle)

##set up named list
names(Cand.set) <- c("logit", "probit", "cloglog")

##compare models
##model names will be taken from the list if modnames is not specified
bictab(cand.set = Cand.set)
}


##two-way ANOVA with interaction
data(iron)
##full model
m1 <- lm(Iron ~ Pot + Food + Pot:Food, data = iron)
##additive model
m2 <- lm(Iron ~ Pot + Food, data = iron)
##null model
m3 <- lm(Iron ~ 1, data = iron)

##candidate models
Cand.aov <- list(m1, m2, m3)
Cand.names <- c("full", "additive", "null")
bictab(Cand.aov, Cand.names)



##single-season occupancy model example modified from ?occu
\dontrun{
require(unmarked)
##single season example modified from ?occu
data(frogs)
pferUMF <- unmarkedFrameOccu(pfer.bin)
##add fake covariates
siteCovs(pferUMF) <- data.frame(sitevar1 = rnorm(numSites(pferUMF)),
                                sitevar2 = runif(numSites(pferUMF))) 
     
##observation covariates
obsCovs(pferUMF) <- data.frame(obsvar1 = rnorm(numSites(pferUMF) *
                                 obsNum(pferUMF))) 

##set up candidate model set
fm1 <- occu(~ obsvar1 ~ sitevar1, pferUMF)
fm2 <- occu(~ 1 ~ sitevar1, pferUMF)
fm3 <- occu(~ obsvar1 ~ sitevar2, pferUMF)
fm4 <- occu(~ 1 ~ sitevar2, pferUMF)

##assemble models in named list (alternative to using 'modnames' argument)
Cand.mods <- list("fm1" = fm1, "fm2" = fm2, "fm3" = fm3, "fm4" = fm4)

##compute table based on QBIC that accounts for c.hat
bictab(cand.set = Cand.mods, c.hat = 3.9)

detach(package:unmarked)
}
}

\keyword{models}

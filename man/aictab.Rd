\name{aictab}
\Rdversion{1.1}
\alias{aictab}
\alias{aictab.default}
\alias{aictab.AICaov.lm}
\alias{aictab.AICsclm.clm}
\alias{aictab.AICclmm}
\alias{aictab.AICcoxme}
\alias{aictab.AICcoxph}
\alias{aictab.AICglm.lm}
\alias{aictab.AICgls}
\alias{aictab.AIClm}
\alias{aictab.AIClme}
\alias{aictab.AIClmekin}
\alias{aictab.AICmaxlikeFit.list}
\alias{aictab.AICmer}
\alias{aictab.AIClmerMod}
\alias{aictab.AICglmerMod}
\alias{aictab.AICnlmerMod}
\alias{aictab.AICmultinom.nnet}
\alias{aictab.AICnlme}
\alias{aictab.AICnls}
\alias{aictab.AICpolr}
\alias{aictab.AICrlm.lm}
\alias{aictab.AICunmarkedFitOccu}
\alias{aictab.AICunmarkedFitColExt}
\alias{aictab.AICunmarkedFitOccuRN}
\alias{aictab.AICunmarkedFitPCount}
\alias{aictab.AICunmarkedFitPCO}
\alias{aictab.AICunmarkedFitDS}
\alias{aictab.AICunmarkedFitGDS}
\alias{aictab.AICunmarkedFitOccuFP}
\alias{aictab.AICunmarkedFitMPois}
\alias{aictab.AICunmarkedFitGMM}
\alias{aictab.AICunmarkedFitGPC}
\alias{aictab.AICvglm}
\alias{aictab.AICzeroinfl}
\alias{print.aictab}
\title{
Create Model Selection Tables
}
\description{
This function creates a model selection table based on one of the
following information criteria:  AIC, AICc, QAIC, QAICc.  The table
ranks the models based on the selected information criteria and also
provides delta AIC and Akaike weights.  \code{aictab} selects the
appropriate function to create the model selection table based on the 
object class.  The current version works with lists containing objects
of \code{clm}, \code{clmm}, \code{clogit}, \code{coxme}, \code{coxph},
\code{glm}, \code{gls}, \code{lm}, \code{lme}, \code{lmekin},
\code{maxlikeFit}, \code{mer}, \code{merMod}, \code{multinom},
\code{nlme}, \code{nls}, \code{polr}, \code{rlm}, \code{vglm}, and
\code{zeroinfl} classes as well as various models of \code{unmarkedFit}
classes but does not yet allow mixing of different classes.
}
\usage{
aictab(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL,
       sort = TRUE, \dots)

\method{aictab}{AICaov.lm}(cand.set, modnames = NULL,
        second.ord = TRUE, nobs = NULL, sort = TRUE, \dots)

\method{aictab}{AICsclm.clm}(cand.set, modnames = NULL,
       second.ord = TRUE, nobs = NULL, sort = TRUE, \dots)

\method{aictab}{AICclmm}(cand.set, modnames = NULL,
       second.ord = TRUE, nobs = NULL, sort = TRUE, \dots)

\method{aictab}{AICcoxme}(cand.set, modnames = NULL,
       second.ord = TRUE, nobs = NULL, sort = TRUE, \dots)

\method{aictab}{AICcoxph}(cand.set, modnames = NULL,
        second.ord = TRUE, nobs = NULL, sort = TRUE, \dots)

\method{aictab}{AICglm.lm}(cand.set, modnames = NULL,
        second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1, \dots)  

\method{aictab}{AICgls}(cand.set, modnames = NULL,
        second.ord = TRUE, nobs = NULL, sort = TRUE, \dots)

\method{aictab}{AIClm}(cand.set, modnames = NULL,
        second.ord = TRUE, nobs = NULL, sort = TRUE, \dots)

\method{aictab}{AIClme}(cand.set, modnames = NULL,
        second.ord = TRUE, nobs = NULL, sort = TRUE, \dots)

\method{aictab}{AIClmekin}(cand.set, modnames = NULL,
        second.ord = TRUE, nobs = NULL, sort = TRUE, \dots)

\method{aictab}{AICmaxlikeFit.list}(cand.set, modnames = NULL,
        second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1, \dots) 

\method{aictab}{AICmer}(cand.set, modnames = NULL,
        second.ord = TRUE, nobs = NULL, sort = TRUE, \dots)

\method{aictab}{AIClmerMod}(cand.set, modnames = NULL,
        second.ord = TRUE, nobs = NULL, sort = TRUE, \dots)

\method{aictab}{AICglmerMod}(cand.set, modnames = NULL,
        second.ord = TRUE, nobs = NULL, sort = TRUE, \dots) 

\method{aictab}{AICnlmerMod}(cand.set, modnames = NULL,
        second.ord = TRUE, nobs = NULL, sort = TRUE, \dots)

\method{aictab}{AICmultinom.nnet}(cand.set, modnames = NULL,
        second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1, \dots)

\method{aictab}{AICnlme}(cand.set, modnames = NULL,
        second.ord = TRUE, nobs = NULL, sort = TRUE, \dots) 

\method{aictab}{AICnls}(cand.set, modnames = NULL,
        second.ord = TRUE, nobs = NULL, sort = TRUE, \dots)

\method{aictab}{AICpolr}(cand.set, modnames = NULL,
        second.ord = TRUE, nobs = NULL, sort = TRUE, \dots)

\method{aictab}{AICrlm.lm}(cand.set, modnames = NULL,
        second.ord = TRUE, nobs = NULL, sort = TRUE, \dots)

\method{aictab}{AICunmarkedFitOccu}(cand.set, modnames = NULL,
        second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1, \dots)

\method{aictab}{AICunmarkedFitColExt}(cand.set, modnames = NULL,
        second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1, \dots)

\method{aictab}{AICunmarkedFitOccuRN}(cand.set, modnames = NULL,
        second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1, \dots)

\method{aictab}{AICunmarkedFitPCount}(cand.set, modnames = NULL,
        second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1, \dots)

\method{aictab}{AICunmarkedFitPCO}(cand.set, modnames = NULL,
        second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1, \dots)

\method{aictab}{AICunmarkedFitDS}(cand.set, modnames = NULL,
        second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1, \dots)

\method{aictab}{AICunmarkedFitGDS}(cand.set, modnames = NULL,
        second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1, \dots)

\method{aictab}{AICunmarkedFitOccuFP}(cand.set, modnames = NULL,
        second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1, \dots)

\method{aictab}{AICunmarkedFitMPois}(cand.set, modnames = NULL,
        second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1, \dots)

\method{aictab}{AICunmarkedFitGMM}(cand.set, modnames = NULL,
        second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1, \dots)

\method{aictab}{AICunmarkedFitGPC}(cand.set, modnames = NULL,
        second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1, \dots)

\method{aictab}{AICvglm}(cand.set, modnames = NULL,
        second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1, \dots)

\method{aictab}{AICzeroinfl}(cand.set, modnames = NULL,
        second.ord = TRUE, nobs = NULL, sort = TRUE, \dots)
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
  \item{second.ord}{
logical.  If \code{TRUE}, the function returns the second-order
Akaike information criterion (i.e., AICc).
}
  \item{nobs}{
this argument allows to specify a numeric value other than total sample
size to compute the AICc (i.e., \code{nobs} defaults to total number of
observations).  This is relevant only for mixed models or various models
of \code{unmarkedFit} classes where sample size is not straightforward.
In such cases, one might use total number of observations or number of
independent clusters (e.g., sites) as the value of \code{nobs}.
}
\item{sort}{
logical.  If \code{TRUE}, the model selection table is ranked according
to the (Q)AIC(c) values.
}
\item{c.hat}{
value of overdispersion parameter (i.e., variance inflation factor) such
as that obtained from \code{c_hat}.  Note that values of c.hat different
from 1 are only appropriate for binomial GLM's with trials > 1 (i.e.,
success/trial or cbind(success, failure) syntax), with Poisson
GLM's, single-season occupancy models (MacKenzie et al. 2002), or
\emph{N}-mixture models (Royle 2004, Dail and Madsen 2011).  If
c.hat > 1, \code{AICc} will return the quasi-likelihood analogue of the
information criterion requested.  This option is not supported for
generalized linear mixed models of the \code{mer} and \code{merMod}
classes.
}
\item{\dots}{
  additional arguments passed to the function.
}
}
\details{
\code{aictab} internally creates a new class for the \code{cand.set}
list of candidate models, according to the contents of the list. The
current function is implemented for \code{clogit}, \code{coxme},
\code{coxph}, \code{glm}, \code{gls}, \code{lm}, \code{lme},
\code{lmekin}, \code{maxlikeFit}, \code{mer}, \code{merMod},
\code{multinom}, \code{nlme}, \code{nls}, \code{polr}, \code{rlm},
\code{vglm}, and \code{zeroinfl} classes as well as various
\code{unmarkedFit} classes. The function constructs a model selection
table based on one of the four information criteria: AIC, AICc, QAIC,
and QAICc. 

Ten guidelines for model selection:

1) Carefully construct your candidate model set.  Each model
should represent a specific (interesting) hypothesis to test.

2) Keep your candidate model set short.  It is ill-advised to consider
as many models as there are data.

3) Check model fit.  Use your global model (most complex model) or
subglobal models to determine if the assumptions are valid.  If none of
your models fit the data well, information criteria will only indicate
the most parsimonious of the poor models.  

4) Avoid data dredging (i.e., looking for patterns after an initial round
of analysis). 

5) Avoid overfitting models.  You should not estimate too many
parameters for the number of observations available in the sample.

6) Be careful of missing values.  Remember that values that are missing
only for certain variables change the data set and sample size,
depending on which variable is included in any given model.  I suggest
to remove missing cases before starting model selection. 

7) Use the same response variable for all models of the candidate model
set.  It is inappropriate to run some models with a transformed response
variable and others with the untransformed variable.  A workaround is to
use a different link function for some models (i.e., identity vs log link).

8) When dealing with models with overdispersion, use the same value of
c-hat for all models in the candidate model set.  For binomial models
with trials > 1 (i.e., success/trial or cbind(success, failure) syntax)
or with Poisson GLM's, you should estimate the c-hat from the most
complex model (global model).  If c-hat > 1, you should use the same
value for each model of the candidate model set (where appropriate) and
include it in the count of parameters (K).  Similarly, for negative
binomial models, you should estimate the dispersion parameter from the
global model and use the same value across all models.  

9) Burnham and Anderson (2002) recommend to avoid mixing the
information-theoretic approach and notions of significance (i.e., P
values).  It is best to provide estimates and a measure of their
precision (standard error, confidence intervals).

10) Determining the ranking of the models is just the first step.
Akaike weights sum to 1 for the entire model set and can be interpreted
as the weight of evidence in favor of a given model being the best one
given the candidate model set considered and the data at hand. Models
with large Akaike weights have strong support.  Evidence ratios,
importance values, and confidence sets for the best model are all
measures that assist in interpretation.  In cases where the top ranking
model has an  Akaike weight > 0.9, one can base inference on this single
most parsimonious model.  When many models rank highly (i.e., delta
(Q)AIC(c) < 4), one should model-average the parameters of interest
appearing in the top models.  Model averaging consists in making
inference based on the whole set of candidate models, instead of basing
conclusions on a single 'best' model.  It is an elegant way of making
inference based on the information contained in the entire model set.
}
\value{
\code{aictab} creates an object of class \code{aictab} with the
following components:

\item{Modname}{the names of each model of the candidate model set.}
\item{K }{the number of estimated parameters for each model.}
\item{(Q)AIC(c) }{the information criteria requested for each model
  (AICc, AICc, QAIC, QAICc).}
\item{Delta_(Q)AIC(c) }{the appropriate delta AIC component depending on
  the information criteria selected.}
\item{ModelLik }{the relative likelihood of the model given the
  data (exp(-0.5*delta[i])).  This is not to be confused with the
  likelihood of the parameters given the data.  The relative likelihood
  can then be normalized across all models to get the model probabilities.} 
\item{(Q)AIC(c)Wt}{the Akaike weights, also termed "model
  probabilities" sensu Burnham and Anderson (2002) and Anderson (2008).
  These measures indicate the level of support (i.e., weight of
  evidence) in favor of any given model being the most parsimonious
  among the candidate model set.}
\item{Cum.Wt }{the cumulative Akaike weights.  These are only meaningful
if results in table are sorted in decreasing order of Akaike weights
(i.e., sort = TRUE).}
\item{c.hat}{if c.hat was specified as an argument, it is included in
  the table.}
\item{LL}{if c.hat = 1 and parameters estimated by maximum likelihood,
  the log-likelihood of each model.}
\item{Quasi.LL}{if c.hat > 1, the quasi log-likelihood of each model.}
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
\code{\link{AICc}}, \code{\link{confset}}, \code{\link{c_hat}},
\code{\link{evidence}}, \code{\link{importance}}, \code{\link{modavg}}, 
\code{\link{modavgShrink}}, \code{\link{modavgPred}} 
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

##generate AICc table
aictab(cand.set = Cand.models, modnames = Modnames, sort = TRUE)
##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = Cand.models, modnames = Modnames, sort = TRUE),
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

##check c-hat
c_hat(Cand.set[[1]])
c_hat(Cand.set[[2]])
c_hat(Cand.set[[3]])
##lowest value of c-hat < 1 for these non-nested models, thus use
##c.hat = 1 
       
##set up named list
names(Cand.set) <- c("logit", "probit", "cloglog")

##compare models
##model names will be taken from the list if modnames is not specified
res.table <- aictab(cand.set = Cand.set, second.ord = FALSE)
##note that delta AIC and Akaike weights are identical to Table 4.7
print(res.table, digits = 2, LL = TRUE) #print table with 2 digits and
##print log-likelihood in table
print(res.table, digits = 4, LL = FALSE) #print table with 4 digits and
##do not print log-likelihood
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
aictab(Cand.aov, Cand.names)



\dontrun{
##modified example of Cox regression from ?coxph
require(survival)
##Create a simple data set for a time-dependent model
test2 <- list(start=c(1, 2, 5, 2, 1, 7, 3, 4, 8, 8),
              stop =c(2, 3, 6, 7, 8, 9, 9, 9,14,17),
              event=c(1, 1, 1, 1, 1, 1, 1, 0, 0, 0),
              x    =c(1, 0, 0, 1, 0, 1, 1, 1, 0, 0),
              height = c(12.3, 10.5, 9.2, 5.6, 8.9,
                11.0, 16.1, 10.2, 9.9, 14.8))
m.cox <- coxph(Surv(start, stop, event) ~ x, test2)
m.cox2 <- coxph(Surv(start, stop, event) ~ x + height, test2)
Cands <- list(m.cox, m.cox2)
Mods <- c("x", "additive")
aictab(Cands, Mods)
}


##mer class example modified from ?glmer
if(require(lme4)) {
##create proportion of incidence
cbpp$prop <- cbpp$incidence/cbpp$size
gm1 <- glmer(prop ~ period + (1 | herd), family = binomial, weights =
             size, data = cbpp)
gm2 <- glmer(prop ~ 1 + (1 | herd), family = binomial, weights =
             size, data = cbpp)

Cands <- list(gm1, gm2)
Model.names <- c("effect of period", "no effect of period")
aictab(cand.set = Cands, modnames = Model.names, sort = TRUE)
detach(package:lme4)
}


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

##compute table
aictab(cand.set = Cand.mods, modnames = Modnames,
       second.ord = TRUE)

detach(package:unmarked)
}
}

\keyword{models}

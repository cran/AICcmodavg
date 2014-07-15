\name{modavgEffect}
\Rdversion{1.1}
\alias{modavgEffect}
\alias{modavgEffect.default}
\alias{modavgEffect.AICaov.lm}
\alias{modavgEffect.AICglm.lm}
\alias{modavgEffect.AICgls}
\alias{modavgEffect.AIClm}
\alias{modavgEffect.AIClme}
\alias{modavgEffect.AICmer}
\alias{modavgEffect.AICglmerMod}
\alias{modavgEffect.AIClmerMod}
\alias{modavgEffect.AICrlm.lm}
\alias{modavgEffect.AICunmarkedFitOccu}
\alias{modavgEffect.AICunmarkedFitColExt}
\alias{modavgEffect.AICunmarkedFitOccuRN}
\alias{modavgEffect.AICunmarkedFitPCount}
\alias{modavgEffect.AICunmarkedFitPCO}
\alias{modavgEffect.AICunmarkedFitDS}
\alias{modavgEffect.AICunmarkedFitGDS}
\alias{modavgEffect.AICunmarkedFitOccuFP}
\alias{modavgEffect.AICunmarkedFitMPois}
\alias{modavgEffect.AICunmarkedFitGMM}
\alias{modavgEffect.AICunmarkedFitGPC}
\alias{print.modavgEffect}

\title{
Compute Model-averaged Effect Sizes (Multimodel Inference on Group Differences)
}
\description{
  This function model-averages the effect size between two groups
  defined by a categorical variable based on the entire model set and
  computes the unconditional standard error and unconditional confidence
  intervals as described in Buckland et al. (1997) and Burnham and
  Anderson (2002). This can be particularly useful when dealing with
  data from an experiment (e.g., ANOVA) and when the focus is to
  determine the effect of a given factor. This is an
  information-theoretic alternative to multiple comparisons (e.g.,
  Burnham et al. 2011).
}
\usage{
modavgEffect(cand.set, modnames = NULL, newdata, second.ord = TRUE,
              nobs = NULL, uncond.se = "revised", conf.level = 0.95,
              \dots)

\method{modavgEffect}{AICaov.lm}(cand.set, modnames = NULL, newdata,
        second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, \dots)

\method{modavgEffect}{AICglm.lm}(cand.set, modnames = NULL, newdata, 
        second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, type = "response", c.hat = 1, gamdisp = NULL,
        \dots)

\method{modavgEffect}{AICgls}(cand.set, modnames = NULL, newdata,
        second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, \dots)

\method{modavgEffect}{AIClm}(cand.set, modnames = NULL, newdata,
        second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, \dots)

\method{modavgEffect}{AIClme}(cand.set, modnames = NULL, newdata,
        second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, \dots) 

\method{modavgEffect}{AICmer}(cand.set, modnames = NULL, newdata,
        second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, type = "response", \dots)

\method{modavgEffect}{AICglmerMod}(cand.set, modnames = NULL,
        newdata, second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, type = "response", \dots)

\method{modavgEffect}{AIClmerMod}(cand.set, modnames = NULL,
        newdata, second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, \dots)

\method{modavgEffect}{AICrlm.lm}(cand.set, modnames = NULL, newdata,
        second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, \dots)

\method{modavgEffect}{AICunmarkedFitOccu}(cand.set, modnames = NULL,
        newdata, second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, type = "response", c.hat = 1, parm.type =
        NULL, \dots)

\method{modavgEffect}{AICunmarkedFitColExt}(cand.set, modnames =
        NULL, newdata, second.ord = TRUE, nobs = NULL, uncond.se =
        "revised", conf.level = 0.95, type = "response", c.hat = 1,
        parm.type = NULL, \dots)

\method{modavgEffect}{AICunmarkedFitOccuRN}(cand.set, modnames =
        NULL, newdata, second.ord = TRUE, nobs = NULL, uncond.se =
        "revised", conf.level = 0.95, type = "response", c.hat = 1,
        parm.type = NULL, \dots)

\method{modavgEffect}{AICunmarkedFitPCount}(cand.set, modnames =
        NULL, newdata, second.ord = TRUE, nobs = NULL, uncond.se =
        "revised", conf.level = 0.95, type = "response", c.hat = 1,
        parm.type = NULL, \dots)

\method{modavgEffect}{AICunmarkedFitPCO}(cand.set, modnames = NULL,
        newdata, second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, type = "response", c.hat = 1, parm.type =
        NULL, \dots)

\method{modavgEffect}{AICunmarkedFitDS}(cand.set, modnames = NULL,
        newdata, second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, type = "response", c.hat = 1, parm.type =
        NULL, \dots)

\method{modavgEffect}{AICunmarkedFitGDS}(cand.set, modnames = NULL,
        newdata, second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, type = "response", c.hat = 1, parm.type =
        NULL, \dots)

\method{modavgEffect}{AICunmarkedFitOccuFP}(cand.set, modnames =
        NULL, newdata, second.ord = TRUE, nobs = NULL, uncond.se =
        "revised", conf.level = 0.95, type = "response", c.hat = 1,
        parm.type = NULL, \dots)

\method{modavgEffect}{AICunmarkedFitMPois}(cand.set, modnames =
        NULL, newdata, second.ord = TRUE, nobs = NULL, uncond.se =
        "revised", conf.level = 0.95, type = "response", c.hat = 1,
        parm.type = NULL, \dots)

\method{modavgEffect}{AICunmarkedFitGMM}(cand.set, modnames =
       NULL, newdata, second.ord = TRUE, nobs = NULL, uncond.se =
       "revised", conf.level = 0.95, type = "response", c.hat = 1,
       parm.type = NULL, \dots)

\method{modavgEffect}{AICunmarkedFitGPC}(cand.set, modnames =
        NULL, newdata, second.ord = TRUE, nobs = NULL, uncond.se =
        "revised", conf.level = 0.95, type = "response", c.hat = 1,
        parm.type = NULL, \dots)
}


%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cand.set}{
    a list storing each of the models in the candidate model set.
  }
  \item{modnames}{
a character vector of model names to facilitate the identification of
each model in the model selection table. If \code{NULL}, the function
uses the names in the cand.set list of candidate models. If no names
appear in the list, generic names (e.g., \code{Mod1}, \code{Mod2}) are
supplied in the table in the same order as in the list of candidate
models.
}
  \item{newdata}{
a data frame with two rows and where the columns correspond to the
explanatory variables specified in the candidate models. Note that this
data set must have the same structure as that of the original data frame
for which we want to make predictions, specifically, the same variable
type and names that appear in the original data set.  Each row of the
data set defines one of the two groups compared.  The first row in
\code{newdata} defines the first group, whereas the second row defines
the second group.  The effect size is computed as the prediction in the
first row minus the prediction in the second row (first row - second 
row). Only the column relating to the grouping variable can change value
and all others must be held constant for the comparison (see 'Details').
}
  \item{second.ord}{
logical.  If \code{TRUE}, the function returns the second-order Akaike
information criterion (i.e., AICc).
}
  \item{nobs}{
this argument allows the specification of a numeric value other than
total sample size to compute the AICc (i.e., \code{nobs} defaults to
total number of observations).  This is relevant only for mixed models
or various models of \code{unmarkedFit} classes where sample size is not
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
the confidence level (1 - alpha) requested for the computation of
unconditional confidence intervals.  To obtain confidence intervals
corrected for multiple comparisons between pairs of treatments, it is
possible to adjust the alpha level according to various strategies such
as the Bonferroni correction (Dunn 1961).
}
  \item{type}{
the scale of prediction requested, one of \code{"response"} or
\code{"link"} (only relevant for \code{glm}, \code{mer}, and
\code{unmarkedFit} classes).  Note that the value \code{"terms"} is not
defined for \code{modavgEffect}).  
}
  \item{c.hat}{
value of overdispersion parameter (i.e., variance inflation factor) such  
as that obtained from \code{c_hat}.  Note that values of \code{c.hat}
different from 1 are only appropriate for binomial GLM's with trials > 1
(i.e., success/trial or cbind(success, failure) syntax), with Poisson
GLM's, single-season occupancy models (MacKenzie et al. 2002), or 
\emph{N}-mixture models (Royle 2004, Dail and Madsen 2011).  If
\code{c.hat > 1}, \code{modavgEffect} will return the quasi-likelihood
analogue of the information criteria requested and multiply the
variance-covariance matrix of the estimates by this value (i.e., SE's
are multiplied by \code{sqrt(c.hat)}.  This option is not supported for
generalized linear mixed models of the \code{mer} class.
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
  character strings supported vary with the type of model fitted. 
  For \code{unmarkedFitOccu} objects, either \code{psi} or \code{detect}
  can be supplied to indicate whether the parameter is on occupancy or
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
The strategy used here to compute effect sizes is to work from the
\code{newdata} object to create two predictions from a given model and
compute the differences and standard errors between both values. This
step is executed for each model in the candidate model set, to obtain a
model-averaged estimate of the effect size and unconditional standard
error. As a result, the \code{newdata} argument is restricted to two
rows, each for a given prediction. To specify each group, the values
entered in the column for each explanatory variable must be identical,
except for the grouping variable. A sensible choice of value for the
explanatory variables is the average of the variable. 

Model-averaging effect sizes is most useful in true experiments (e.g.,
ANOVA-type designs), where one wants to obtain the best estimate of
effect size given the support of each candidate model. This can be
considered as a information-theoretic analog of traditional multiple
comparisons, except that the information contained in the entire model
set is used instead of being restricted to a single model.  See
'Examples' below for applications.

\code{modavgEffect} calls the appropriate method depending on the class
of objects in the list.  The current classes supported include
\code{glm}, \code{gls}, \code{lm}, \code{lme}, \code{mer},
\code{glmerMod}, \code{lmerMod}, \code{rlm}, as well as
\code{unmarkedFitOccu}, \code{unmarkedFitColExt},
\code{unmarkedFitOccuFP}, \code{unmarkedFitOccuRN},
\code{unmarkedFitPCount}, \code{unmarkedFitPCO}, \code{unmarkedFitDS},
\code{unmarkedFitGDS}, \code{unmarkedFitMPois}, \code{unmarkedFitGMM},
and \code{unmarkedFitGPC} classes.
}
\value{
The result is an object of class \code{modavgEffect} with the following
components: 

\item{Group.variable}{the grouping variable defining the two groups
  compared} 
\item{Group1}{the first group considered in the comparison}
\item{Group2}{the second group considered in the comparison}
\item{Type}{the scale on which the model-averaged effect size was
  computed (e.g., response or link)}
\item{Mod.avg.table}{the full model selection table including the
entire set of candidate models}
\item{Mod.avg.eff}{the model-averaged effect size based on the entire
  candidate model set}
\item{Uncond.SE}{the unconditional standard error for the model-averaged
  effect size}
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

Burnham, K. P., Anderson, D. R., Huyvaert, K. P. (2011) AIC model
selection and multimodel inference in behaviorial ecology: some
background, observations and comparisons. \emph{Behavioral Ecology and
Sociobiology} \bold{65}, 23--25. 

Dail, D., Madsen, L. (2011) Models for estimating abundance from 
repeated counts of an open population. \emph{Biometrics} \bold{67},
577--587.

Dunn, O. J. (1961) Multiple comparisons among means. \emph{Journal of the
American Statistical Association} \bold{56}, 52--64.

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
\code{\link{confset}}, \code{\link{evidence}}, \code{\link{importance}},
\code{\link{modavgShrink}}, \code{\link{modavgPred}}
}
\examples{
##heights (cm) of plants grown under two fertilizers, Ex. 9.5 from
##Zar (1984)
heights <- data.frame(Height = c(48.2, 54.6, 58.3, 47.8, 51.4, 52.0,
                        55.2, 49.1, 49.9, 52.6, 52.3, 57.4, 55.6, 53.2,
                        61.3, 58.0, 59.8, 54.8),
                      Fertilizer = c(rep("old", 10), rep("new", 8)))

##run linear model hypothesizing an effect of fertilizer
m1 <- lm(Height ~ Fertilizer, data = heights)

##run null model (no effect of fertilizer)
m0 <- lm(Height ~ 1, data = heights)

##assemble models in list
Cands <- list(m1, m0)
Modnames <- c("Fert", "null")

##compute model selection table to compare
##both hypotheses
aictab(cand.set = Cands, modnames = Modnames)
##note that model with fertilizer effect is much better supported
##than the null

##compute model-averaged effect sizes: one model hypothesizes a
##difference of 0, whereas the other assumes a difference

##prepare newdata object from which differences between groups
##will be computed
##the first row of the newdata data.frame relates to the first group,
##whereas the second row corresponds to the second group
pred.data <- data.frame(Fertilizer = c("new", "old"))

##compute best estimate of effect size accounting for model selection
##uncertainty
modavgEffect(cand.set = Cands, modnames = Modnames,
              newdata = pred.data)


##classical one-way ANOVA type-design
##generate data for two groups and control
set.seed(seed = 15)
y <- round(c(rnorm(n = 15, mean = 10, sd = 5),
       rnorm(n = 15, mean = 15, sd = 5),
       rnorm(n = 15, mean = 12, sd = 5)), digits = 2)
##groups
group <- c(rep("cont", 15), rep("trt1", 15), rep("trt2", 15))

##combine in data set
aov.data <- data.frame(Y = y, Group = group)
rm(y, group)

##run model with group effect
lm.eff <- lm(Y ~ Group, data = aov.data)
##null model
lm.0 <- lm(Y ~ 1, data = aov.data)

##compare both models
Cands <- list(lm.eff, lm.0)
Mods <- c("group effect", "no group effect")
aictab(cand.set = Cands, modnames = Mods)
##model with group effect has most of the weight

##compute model-averaged effect sizes
##trt1 - control
modavgEffect(cand.set = Cands, modnames = Modnames,
              newdata = data.frame(Group = c("trt1", "cont")))
##trt1 differs from cont

##trt2 - control
modavgEffect(cand.set = Cands, modnames = Modnames,
              newdata = data.frame(Group = c("trt2", "cont")))
##trt2 does not differ from cont



##two-way ANOVA type design, Ex. 13.1 (Zar 1984) of plasma calcium
##concentration (mg/100 ml) in birds as a function of sex and hormone
##treatment
birds <- data.frame(Ca = c(16.87, 16.18, 17.12, 16.83, 17.19, 15.86,
                      14.92, 15.63, 15.24, 14.8, 19.07, 18.77, 17.63,
                      16.99, 18.04, 17.2, 17.64, 17.89, 16.78, 16.92,
                      32.45, 28.71, 34.65, 28.79, 24.46, 30.54, 32.41,
                      28.97, 28.46, 29.65),
                    Sex = c("M", "M", "M", "M", "M", "F", "F", "F", "F",
                      "F", "M", "M", "M", "M", "M", "F", "F", "F", "F",
                      "F", "M", "M", "M", "M", "M", "F", "F", "F", "F",
                      "F"),
                    Hormone = as.factor(c(1, 1, 1, 1, 1, 1, 1, 1, 1,
                      1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3,
                      3, 3, 3, 3, 3)))


##candidate models
##interactive effects
m.inter <- lm(Ca ~ Sex + Hormone + Sex:Hormone, data = birds)

##additive effects
m.add <- lm(Ca ~ Sex + Hormone, data = birds)

##Sex only
m.sex <- lm(Ca ~ Sex, data = birds)

##Hormone only
m.horm <- lm(Ca ~ Hormone, data = birds)

##null
m.0 <- lm(Ca ~ 1, data = birds)

##model selection
Cands <- list(m.inter, m.add, m.sex, m.horm, m.0)
Mods <- c("interaction", "additive", "sex only", "horm only", "null")
aictab(Cands, Mods)
##there is some support for a hormone only treatment, but also for
##additive effects

##compute model-averaged effects of sex, and set the other variable
##to a constant value
##M - F
sex.data <- data.frame(Sex = c("M", "F"), Hormone = c("1", "1"))
modavgEffect(Cands, Mods, newdata = sex.data)
##no support for a sex main effect

##hormone 1 - 3, but set Sex to a constant value
horm1.data <- data.frame(Sex = c("M", "M"), Hormone = c("1", "3"))
modavgEffect(Cands, Mods, newdata = horm1.data)

##hormone 2 - 3, but set Sex to a constant value
horm2.data <- data.frame(Sex = c("M", "M"), Hormone = c("2", "3"))
modavgEffect(Cands, Mods, newdata = horm2.data)



##Poisson regression with anuran larvae example from Mazerolle (2006)
data(min.trap)
##assign "UPLAND" as the reference level as in Mazerolle (2006)          
min.trap$Type <- relevel(min.trap$Type, ref = "UPLAND") 

##set up candidate models          
Cand.mod <- list( )
##global model          
Cand.mod[[1]] <- glm(Num_anura ~ Type + log.Perimeter,
                     family = poisson, offset = log(Effort),
                     data = min.trap) 
Cand.mod[[2]] <- glm(Num_anura ~ log.Perimeter, family = poisson,
                     offset = log(Effort), data = min.trap) 
Cand.mod[[3]] <- glm(Num_anura ~ Type, family = poisson,
                     offset = log(Effort), data = min.trap) 
Cand.mod[[4]] <- glm(Num_anura ~ 1, family = poisson,
                     offset = log(Effort), data = min.trap) 
          
##check c-hat for global model
vif.hat <- c_hat(Cand.mod[[1]]) #uses Pearson's chi-square/df

##assign names to each model
Modnames <- c("type + logperim", "type", "logperim", "intercept only") 

##compute model-averaged estimate of difference between abundance at bog
##pond and upland pond
##create newdata object to make predictions
pred.data <- data.frame(Type = c("BOG", "UPLAND"),
                        log.Perimeter = mean(min.trap$log.Perimeter),
                        Effort = mean(min.trap$Effort))
modavgEffect(Cand.mod, Modnames, newdata = pred.data, c.hat = vif.hat, type = "response")
##little suport for an pond type effect


##mixed linear model example from ?nlme
library(nlme)
Cand.models <- list( )
Cand.models[[1]] <- lme(distance ~ age, data = Orthodont, method="ML") # random is ~ age
Cand.models[[2]] <- lme(distance ~ age + Sex, data = Orthodont, random = ~ 1, method="ML")
Cand.models[[3]] <-lme(distance ~ 1, data = Orthodont, random = ~ 1, method="ML")
Cand.models[[4]] <-lme(distance ~ Sex, data = Orthodont, random = ~ 1, method="ML")

Modnames <- c("age", "age + sex", "null", "sex")

data.other <- data.frame(age = mean(Orthodont$age), Sex = factor(c("Male", "Female")))
modavgEffect(cand.set = Cand.models, modnames = Modnames, newdata = data.other,
              conf.level = 0.95, second.ord = TRUE, nobs = NULL, uncond.se = "revised")


\dontrun{
##site occupancy analysis example
library(unmarked)
##single season model
data(frogs)
pferUMF <- unmarkedFrameOccu(pfer.bin)
##create a bogus site group
site.group <- c(rep(1, times = nrow(pfer.bin)/2), rep(0, nrow(pfer.bin)/2))

## add some fake covariates for illustration
siteCovs(pferUMF) <- data.frame(site.group, sitevar1 =
                                rnorm(numSites(pferUMF)),
                                sitevar2 = runif(numSites(pferUMF)))
     
## observation covariates are in site-major, observation-minor order
obsCovs(pferUMF) <- data.frame(obsvar1 = rnorm(numSites(pferUMF) * obsNum(pferUMF)))
     
fm1 <- occu(~ obsvar1 ~ site.group, pferUMF)
fm2 <- occu(~ obsvar1 ~ 1, pferUMF)

Cand.mods <- list(fm1, fm2)
Modnames <- c("fm1", "fm2")

##model selection table
aictab(cand.set = Cand.mods, modnames = Modnames, second.ord = TRUE)

##model-averaged effect sizes comparing site.group 1 - site.group 0
newer.dat <- data.frame(site.group = c(0, 1))

modavgEffect(cand.set = Cand.mods, modnames = Modnames, type = "response",
              second.ord = TRUE, newdata = newer.dat, parm.type = "psi")
##no support for an effect of site group (not surprising given that we
##created this bogus variable just for this example


##single season N-mixture models
data(mallard)
##this variable was created to illustrate the use of modavgEffect( )
##with detection variables
mallard.site$site.group <- c(rep(1, 119), rep(0, 120))
mallardUMF <- unmarkedFramePCount(mallard.y, siteCovs = mallard.site,
                                  obsCovs = mallard.obs)
siteCovs(mallardUMF)
tmp.covs <- obsCovs(mallardUMF)
obsCovs(mallardUMF)$date2 <- tmp.covs$date^2
(fm.mall <- pcount(~ site.group ~ length + elev + forest, mallardUMF, K=30))
(fm.mallb <- pcount(~ 1 ~ length + elev + forest, mallardUMF, K=30))
     
Cands <- list(fm.mall, fm.mallb)
Modnames <- c("one", "null")

##model averaged effect size of site.group 1 - site.group 0 on response
##scale (point estimate)
modavgEffect(Cands, Modnames, newdata = data.frame(site.group = c(0, 1)),
              parm.type = "detect", type = "response")

##model averaged effect size of site.group 1 - site.group 0 on link
##scale (here, logit link)
modavgEffect(Cands, Modnames, newdata = data.frame(site.group = c(0, 1)),
              parm.type = "detect", type = "link")

detach(package:unmarked)
}
}
\keyword{models}

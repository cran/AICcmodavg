\name{aictab}
\Rdversion{1.1}
\alias{aictab}
\alias{aictab.glm}
\alias{aictab.gls}
\alias{aictab.lme}
\alias{aictab.mer}
\alias{aictab.mult}
\alias{aictab.polr}
\alias{aictab.unmarked}
\alias{print.aictab}
\title{
Create Model Selection Tables
}
\description{
This function creates a model selection table based on one of the
following information criteria:  AIC, AICc, QAIC, QAICc.  The table
ranks the models based on the selected information criteria and also
provides delta AIC and Akaike weights.  'aictab' selects the
appropriate function to create the model selection table based on the 
object class.  The current version works with objects of 'lm', 'glm',
'gls', 'lme', 'mer', 'multinom', 'polr' classes as well as various
models of 'unmarkedFit' classes but does not yet allow mixing of
different classes.  
}
\usage{
aictab(cand.set, modnames, sort = TRUE, c.hat = 1, second.ord = TRUE,
       nobs = NULL)

aictab.glm(cand.set, modnames, sort = TRUE, c.hat = 1,
           second.ord = TRUE, nobs = NULL)

aictab.gls(cand.set, modnames, sort = TRUE, second.ord = TRUE,
           nobs = NULL) 

aictab.lme(cand.set, modnames, sort = TRUE, second.ord = TRUE,
           nobs = NULL) 

aictab.mer(cand.set, modnames, sort = TRUE, second.ord = TRUE,
           nobs = NULL) 

aictab.mult(cand.set, modnames, sort = TRUE, c.hat = 1,
            second.ord = TRUE, nobs = NULL)

aictab.polr(cand.set, modnames, sort = TRUE, second.ord = TRUE,
            nobs = NULL) 

aictab.unmarked(cand.set, modnames, sort = TRUE, c.hat = 1,
                second.ord = TRUE, nobs = NULL) 

}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cand.set}{
a list storing each of the models in the candidate model set.
}
  \item{modnames}{
a character vector of model names to facilitate the identification of
each model in the model selection table.
}
  \item{sort}{
logical.  If true, the model selection table is ranked according to the
(Q)AIC(c) values.
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
logical.  If TRUE, the function returns the second-order
Akaike information criterion (i.e., AICc).
}
  \item{nobs}{
this argument allows to specify a numeric value other than total sample
size to compute the AICc (i.e., 'nobs' defaults to total number of
observations).  This is relevant only for mixed models or various models
of 'unmarkedFit' classes where sample size is not straightforward.  In
such cases, one might use total number of observations or number of
independent clusters (e.g., sites) as the value of 'nobs'.
}
}
\details{
'aictab' is a function that calls either 'aictab.glm', 'aictab.gls',
'aictab.lme', 'aictab.mer', 'aictab.mult', 'aictab.polr', or
'aictab.unmarked' depending on the class of the object.  The current
function is implemented for 'lm', 'glm', 'gls', 'lme', 'mer', 'multinom',
and 'polr' classes as well as 'unmarkedFitOccu', 'unmarkedFitColExt',
'unmarkedFitOccuRN', 'unmarkedFitPCount', 'unmarkedFitPCO' classes.  The
function constructs a model selection table based on one of the four
information criteria:  AIC, AICc, QAIC, and QAICc.

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
'aictab', 'aictab.glm', 'aictab.lme', 'aictab.mer', 'aictab.mult', 
'aictab.polr', and 'aictab.unmarked' create an object of class 'aictab'
with the following components: 

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
\item{(Q)AIC(c)wt}{the Akaike weights, also termed "model
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

MacKenzie, D. I., Nichols, J. D., Lachman, G. B., Droege, S., Royle,
J. A., Langtimm, C. A. (2002) Estimating site occupancy rates when
detection probabilities are less than one. \emph{Ecology} \bold{83},
2248--2255.

Mazerolle, M. J. (2006) Improving data analysis in herpetology: using
Akaike's Information Criterion (AIC) to assess the strength of
biological hypotheses. \emph{Amphibia-Reptilia} \bold{27}, 169--180.
}
\author{
Marc J. Mazerolle
}
\seealso{
\code{\link{AICc}}, \code{\link{confset}}, \code{\link{c_hat}},
\code{\link{evidence}}, \code{\link{importance}}, \code{\link{modavg}}, 
\code{\link{modavg.shrink}}, \code{\link{modavgpred}} 
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
       
Modnames <- paste("Mod", 1:length(Cand.set), sep = " ")
res.table <- aictab(cand.set = Cand.set, modnames = Modnames,
                    second.ord = FALSE)
##note that delta AIC and Akaike weights are identical to Table 4.7
print(res.table, digits = 2, LL = TRUE) #print table with 2 digits and
##print log-likelihood in table
print(res.table, digits = 4, LL = FALSE) #print table with 4 digits and
##do not print log-likelihood



##mer class example modified from ?glmer
require(lme4)
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

##compute table
aictab(cand.set = Cand.mods, modnames = Modnames,
       second.ord = TRUE)

detach(package:unmarked)
}

\keyword{models}

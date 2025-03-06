\name{AICcmodavg-package}
\Rdversion{1.1}
\alias{AICcmodavg-package}
\alias{AICcmodavg}
\docType{package}
\title{
Model Selection and Multimodel Inference Based on (Q)AIC(c)
}
\description{
Description:  This package includes functions to create model selection
tables based on Akaike's information criterion (AIC) and the
second-order AIC (AICc), as well as their quasi-likelihood counterparts
(QAIC, QAICc).  The package also features functions to conduct classic
model averaging (multimodel inference) for a given parameter of interest
or predicted values, as well as a shrinkage version of model averaging
parameter estimates.  Other handy functions enable the computation of
relative variable importance, evidence ratios, and confidence sets for
the best model.  The present version supports Cox proportional hazards
models and conditional logistic regression (\code{coxph} and
\code{coxme} classes), linear models (\code{lm} class), generalized
linear models (\code{glm}, \code{glm.nb}, \code{vglm}, \code{hurdle},
and \code{zeroinfl} classes), linear models fit by generalized least
squares (\code{gls} class), linear mixed models (\code{lme} class),
generalized linear mixed models (\code{mer}, \code{merMod}, and
\code{glmmTMB} classes), multinomial and ordinal logistic regressions
(\code{multinom}, \code{polr}, \code{clm}, and \code{clmm} classes),
robust regression models (\code{rlm} class), beta regression models
(\code{betareg} class), parametric survival models (\code{survreg}
class), nonlinear models (\code{nls} and \code{gnls} classes), nonlinear
mixed models (\code{nlme} and \code{nlmerMod} classes), univariate models
(\code{fitdist} and \code{fitdistr} classes), and certain types of
latent variable models (\code{lavaan} class).  The package also supports
various models of \code{unmarkedFit} and \code{maxLikeFit} classes
estimating demographic parameters after accounting for imperfect
detection probabilities.  Some functions also allow the creation of
model selection tables for Bayesian models of the \code{bugs} and
\code{rjags} classes.  Objects following model selection and multimodel
inference can be formatted to LaTeX using \code{xtable} methods included
in the package. 
}

\details{
\tabular{ll}{
Package: \tab AICcmodavg\cr
Type: \tab Package\cr
Version: \tab 2.3-4\cr
Date: \tab 2025-03-03\cr
License: \tab GPL (>=2 )\cr
LazyLoad: \tab yes\cr
}

Many functions of the package require a list of models as the input to
conduct model selection and multimodel inference.  Thus, you should
start by organizing the output of the models in a list (See 'Examples'
below).

This package contains several useful functions for model selection and
multimodel inference for several model classes:
\describe{
\item{\code{\link{AICc}}}{Computes AIC, AICc, and their quasi-likelihood counterparts (QAIC, QAICc).}
\item{\code{\link{aictab}}}{Constructs model selection tables with
  number of parameters, AIC, delta AIC, Akaike weights or variants based on AICc, QAIC, and QAICc for a set of candidate models.}
\item{\code{\link{bictab}}}{Constructs model selection tables with number of parameters, BIC, delta BIC, BIC weights for a set of
  candidate models.}
\item{\code{\link{boot.wt}}}{Computes summary statistics from
  detection histories.}
\item{\code{\link{confset}}}{Determines the confidence set for the
  best model based on one of three criteria.}
\item{\code{\link{DIC}}}{Extracts DIC.}
\item{\code{\link{dictab}}}{Constructs model selection tables with
  number of parameters, DIC, delta DIC, DIC weights for a set of
  candidate models.}
\item{\code{\link{evidence}}}{Computes the evidence ratio between the
  highest-ranked model based on the information criteria selected and a
  lower-ranked model.}
\item{\code{\link{importance}}}{Computes importance values (w+) for
  the support of a given parameter among set of candidate models.}
\item{\code{\link{modavg}}}{Computes model-averaged estimate,
  unconditional standard error, and unconditional confidence interval of
  a parameter of interest among a set of candidate models.}
\item{\code{\link{modavgEffect}}}{Computes model-averaged effect
  sizes between groups based on the entire candidate model set.}
\item{\code{\link{modavgShrink}}}{Computes shrinkage version of
  model-averaged estimate, unconditional standard error, and
  unconditional confidence interval of a parameter of interest among
  entire set of candidate models.} 
\item{\code{\link{modavgPred}}}{Computes model-average predictions, 
  unconditional SE's, and confidence intervals among entire set of
  candidate models.}
\item{\code{\link{multComp}}}{Performs multiple comparisons across
  levels of a factor in a model selection framework.}
\item{\code{\link{useBIC}}}{Computes BIC or a quasi-likelihood
  counterparts (QBIC).}
}

For models not yet supported by the functions above, the following
can be useful for model selection and multimodel inference conducted
from input values supplied by the user:
\describe{
  \item{\code{\link{AICcCustom}}}{Computes AIC, AICc, QAIC, and QAICc from
    user-supplied input values of log-likelihood and number of
    parameters.}
  \item{\code{\link{aictabCustom}}}{Creates model selection tables based on
    (Q)AIC(c) from user-supplied input values of log-likelihood and
    number of parameters.}
  \item{\code{\link{bictabCustom}}}{Creates model selection tables based on
    (Q)BIC from user-supplied input values of log-likelihood and number
    of parameters.}
  \item{\code{\link{ictab}}}{Creates model selection tables from user-supplied
    values of an information criterion.}
  \item{\code{\link{modavgCustom}}}{Computes model-averaged parameter estimate
    based on (Q)AIC(c) from user-supplied input values of
    log-likelihood, number of parameters, parameter estimates, and
    standard errors.}
  \item{\code{\link{modavgIC}}}{Computes model-averaged parameter estimate
    from user-supplied values of information criterion, parameter
    estimates, and standard errors.}
  \item{\code{\link{useBICCustom}}}{Computes BIC and QBIC from user-supplied
    input values of log-likelihoods and number of parameters.}
}


A number of functions for model diagnostics are available:
\describe{
\item{\code{\link{c_hat}}}{Estimates variance inflation factor for
  binomial or Poisson GLM's based on various estimators.}
\item{\code{\link{checkConv}}}{Checks the convergence information of
  the algorithm for the model.}
\item{\code{\link{checkParms}}}{Checks the occurrence of parameter
  estimates with high standard errors in a model.}
\item{\code{\link{countDist}}}{Computes summary statistics from
  distance sampling data.}
\item{\code{\link{countHist}}}{Computes summary statistics from
  count history data.} 
\item{\code{\link{covDiag}}}{Computes covariance diagnostics for
  lambda in \emph{N}-mixture models.}
\item{\code{\link{detHist}}}{Computes summary statistics from
  detection histories.}
\item{\code{\link{detTime}}}{Computes summary statistics from
  time-to-detection data.}
\item{\code{\link{extractCN}}}{Extracts condition number from models of
  certain classes.}
\item{\code{\link{mb.gof.test}}}{Computes the MacKenzie and Bailey
  goodness-of-fit test for single season and dynamic occupancy models
  using the Pearson chi-square statistic.}
\item{\code{\link{Nmix.gof.test}}}{Computes goodness-of-fit test for
  \emph{N}-mixture models based on the Pearson chi-square statistic.}
}

Other utility functions include:
\describe{
\item{\code{\link{anovaOD}}}{Computes likelihood-ratio test statistic
  corrected for overdispersion between two models.}
\item{\code{\link{extractLL}}}{Extracts log-likelihood from
  models of certain classes.}
\item{\code{\link{extractSE}}}{Extracts standard errors from models of
  certain classes and adds the labels.}
\item{\code{\link{extractX}}}{Extracts the predictors and associated
  information on variables from a list of candidate models.}
\item{\code{\link{fam.link.mer}}}{Extracts the distribution family and 
  link function from a generalized linear mixed model of classes \code{mer}
  and \code{merMod}.}
\item{\code{\link{predictSE}}}{Computes predictions and associated 
  standard errors models of certain classes.}
\item{\code{\link{summaryOD}}}{Displays summary of model output
  adjusted for overdispersion.}
\item{\link[xtable]{xtable}}{Extends methods from \pkg{xtable} to
  format various objects resulting from model selection and multimodel
  inference to LaTeX or HTML tables.}
}
}

\author{
  Marc J. Mazerolle <marc.mazerolle@uqat.ca>.
}

\references{
  Anderson, D. R. (2008) \emph{Model-based inference in the life
  sciences: a primer on evidence}. Springer: New York. 

  Burnham, K. P., and Anderson, D. R. (2002) \emph{Model selection and
  multimodel inference: a practical information-theoretic
  approach}. Second edition. Springer: New York.

  Burnham, K. P., Anderson, D. R. (2004) Multimodel inference:
  understanding AIC and BIC in model selection. \emph{Sociological
  Methods and Research} \bold{33}, 261--304.  

  Mazerolle, M. J. (2006) Improving data analysis in herpetology: using
  Akaike's Information Criterion (AIC) to assess the strength of
  biological hypotheses. \emph{Amphibia-Reptilia} \bold{27}, 169--180. 
}

\examples{
##Example 1:  Poisson GLM with offset
##anuran larvae example from Mazerolle (2006) 
data(min.trap)
##assign "UPLAND" as the reference level as in Mazerolle (2006)          
min.trap$Type <- relevel(min.trap$Type, ref = "UPLAND") 

##set up candidate models in a list
Cand.mod <- list()
##global model          
Cand.mod[[1]] <- glm(Num_anura ~ Type + log.Perimeter + Num_ranatra,
                     family = poisson, offset = log(Effort),
                     data = min.trap) 
Cand.mod[[2]] <- glm(Num_anura ~ Type + log.Perimeter, family = poisson,
                     offset = log(Effort), data = min.trap) 
Cand.mod[[3]] <- glm(Num_anura ~ Type + Num_ranatra, family = poisson,
                     offset = log(Effort), data = min.trap) 
Cand.mod[[4]] <- glm(Num_anura ~ Type, family = poisson,
                     offset = log(Effort), data = min.trap) 
Cand.mod[[5]] <- glm(Num_anura ~ log.Perimeter + Num_ranatra,
                     family = poisson, offset = log(Effort),
                     data = min.trap) 
Cand.mod[[6]] <- glm(Num_anura ~ log.Perimeter, family = poisson,
                     offset = log(Effort), data = min.trap) 
Cand.mod[[7]] <- glm(Num_anura ~ Num_ranatra, family = poisson,
                     offset = log(Effort), data = min.trap) 
Cand.mod[[8]] <- glm(Num_anura ~ 1, family = poisson,
                     offset = log(Effort), data = min.trap) 
          
##check c-hat for global model
c_hat(Cand.mod[[1]], method = "pearson") #uses Pearson's chi-square/df
##note the very low overdispersion: in this case, the analysis could be
##conducted without correcting for c-hat as its value is reasonably close
##to 1  

##output of model corrected for overdispersion
summaryOD(Cand.mod[[1]], c.hat = 1.04)

##assign names to each model
Modnames <- c("type + logperim + invertpred", "type + logperim",
              "type + invertpred", "type", "logperim + invertpred",
              "logperim", "invertpred", "intercept only") 

##model selection table based on AICc
aictab(cand.set = Cand.mod, modnames = Modnames)

##compute evidence ratio
evidence(aictab(cand.set = Cand.mod, modnames = Modnames))

##compute confidence set based on 'raw' method
confset(cand.set = Cand.mod, modnames = Modnames, second.ord = TRUE,
        method = "raw")  

##compute importance value for "TypeBOG" - same number of models
##with vs without variable
importance(cand.set = Cand.mod, modnames = Modnames, parm = "TypeBOG") 

##compute model-averaged estimate of "TypeBOG" using the natural average
modavg(cand.set = Cand.mod, modnames = Modnames, parm = "TypeBOG")

##compute model-averaged estimate of "TypeBOG" using shrinkage estimator
##same number of models with vs without variable
modavgShrink(cand.set = Cand.mod, modnames = Modnames,
             parm = "TypeBOG")

##compute model-averaged predictions for two types of ponds
##create a data set for predictions
dat.pred <- data.frame(Type = factor(c("BOG", "UPLAND")),
                       log.Perimeter = mean(min.trap$log.Perimeter),
                       Num_ranatra = mean(min.trap$Num_ranatra),
                       Effort = mean(min.trap$Effort))

##model-averaged predictions across entire model set
modavgPred(cand.set = Cand.mod, modnames = Modnames,
           newdata = dat.pred, type = "response")

##compute model-averaged effect size between two groups
##'newdata' data frame must be limited to two rows
modavgEffect(cand.set = Cand.mod, modnames = Modnames,
             newdata = dat.pred, type = "link")


\dontrun{
##Example 2:  single-season occupancy model example modified from ?occu
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

##check detection history data from data object
detHist(pferUMF)

##set up candidate model set
fm1 <- occu(~ obsvar1 ~ sitevar1, pferUMF)
##check detection history data from model object
detHist(fm1)

fm2 <- occu(~ 1 ~ sitevar1, pferUMF)
fm3 <- occu(~ obsvar1 ~ sitevar2, pferUMF)
fm4 <- occu(~ 1 ~ sitevar2, pferUMF)
Cand.models <- list(fm1, fm2, fm3, fm4)

##assign names to elements in list
##alternative to using 'modnames' argument
names(Cand.models) <- c("fm1", "fm2", "fm3", "fm4")

##check GOF of global model and estimate c-hat
mb.gof.test(fm4, nsim = 100) #nsim should be > 1000

##check for high SE's in models
lapply(Cand.models, checkParms, simplify = FALSE)

##compute table
print(aictab(cand.set = Cand.models,
             second.ord = TRUE), digits = 4)

##export as LaTeX table
if(require(xtable)) {
xtable(aictab(cand.set = Cand.models,
              second.ord = TRUE))
}

##compute evidence ratio
evidence(aictab(cand.set = Cand.models))
##evidence ratio between top model vs lowest-ranked model
evidence(aictab(cand.set = Cand.models), model.high = "fm2", model.low = "fm3")

##compute confidence set based on 'raw' method
confset(cand.set = Cand.models, second.ord = TRUE,
        method = "raw")  

##compute importance value for "sitevar1" on occupancy
##same number of models with vs without variable
importance(cand.set = Cand.models, parm = "sitevar1",
           parm.type = "psi") 

##compute model-averaged estimate of "sitevar1" on occupancy
##(natural average)
modavg(cand.set = Cand.models, parm = "sitevar1",
       parm.type = "psi")

##compute model-averaged estimate of "sitevar1"
##(shrinkage estimator)
##same number of models with vs without variable
modavgShrink(cand.set = Cand.models,
             parm = "sitevar1", parm.type = "psi")

##compute model-average predictions
##check explanatory variables appearing in models
extractX(Cand.models, parm.type = "psi")

##create a data set for predictions
dat.pred <- data.frame(sitevar1 = seq(from = min(siteCovs(pferUMF)$sitevar1),
                         to = max(siteCovs(pferUMF)$sitevar1), by = 0.5),
                       sitevar2 = mean(siteCovs(pferUMF)$sitevar2))

##model-averaged predictions of psi across range of values
##of sitevar1 and entire model set
modavgPred(cand.set = Cand.models, newdata = dat.pred,
           parm.type = "psi")
detach(package:unmarked)
}

\dontrun{
##Example 3:  example with user-supplied values of log-likelihoods and
##number of parameters

##vector with model LL's
LL <- c(-38.8876, -35.1783, -64.8970)
     
##vector with number of parameters
Ks <- c(7, 9, 4)
     
##create a vector of names to trace back models in set
Modnames <- c("Cm1", "Cm2", "Cm3")

##generate AICc table
aictabCustom(logL = LL, K = Ks, modnames = Modnames, nobs = 121,
             sort = TRUE) 
##generate AIC table
aictabCustom(logL = LL, K = Ks, modnames = Modnames,
             second.ord = FALSE, nobs = 121, sort = TRUE)

##model averaging parameter estimate
##vector of beta estimates for a parameter of interest
model.ests <- c(0.0478, 0.0480, 0.0478)
     
##vector of SE's of beta estimates for a parameter of interest
model.se.ests <- c(0.0028, 0.0028, 0.0034)
     
##compute model-averaged estimate and unconditional SE based on AICc
modavgCustom(logL = LL, K = Ks, modnames = Modnames, 
             estimate = model.ests, se = model.se.ests, nobs = 121)
##compute model-averaged estimate and unconditional SE based on BIC
modavgCustom(logL = LL, K = Ks, modnames = Modnames, 
             estimate = model.ests, se = model.se.ests, nobs = 121,
             useBIC = TRUE)
}

\dontrun{
##Example 4:  example with user-supplied values of information criterion
##model selection based on WAIC

##WAIC values
waic <- c(105.74, 107.36, 108.24, 100.57)
##number of effective parameters
effK <- c(7.45, 5.61, 6.14, 6.05)

##create a vector of names to trace back models in set
Modnames <- c("global model", "interactive model",
              "additive model", "invertpred model")

##generate WAIC model selection table
ictab(ic = waic, K = effK, modnames = Modnames,
      sort = TRUE, ic.name = "WAIC")

##compute model-averaged estimate
##vector of predictions
Preds <- c(0.106, 0.137, 0.067, 0.050)
##vector of SE's for prediction
Ses <- c(0.128, 0.159, 0.054, 0.039)

##compute model-averaged estimate and unconditional SE based on WAIC
modavgIC(ic = waic, K = effK, modnames = Modnames, 
         estimate = Preds, se = Ses,
         ic.name = "WAIC")

##export as LaTeX table
if(require(xtable)) {
##model-averaged estimate and confidence interval
xtable(modavgIC(ic = waic, K = effK, modnames = Modnames, 
       estimate = Preds, se = Ses,
       ic.name = "WAIC"))
##model selection table with estimate and SE's from each model
xtable(modavgIC(ic = waic, K = effK, modnames = Modnames, 
       estimate = Preds, se = Ses,
       ic.name = "WAIC"), print.table = TRUE)
}
}
}

\keyword{ models }

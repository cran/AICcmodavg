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
the best model.  The present version works with Cox proportional hazards
models and conditional logistic regression (\code{coxph} and
\code{coxme} classes), linear models (\code{lm} class), generalized
linear models (\code{glm}, \code{vglm}, and \code{zeroinfl} classes),
linear models fit by generalized least squares (\code{gls} class),
linear mixed models (\code{lme} class), generalized linear mixed models
(\code{mer} and \code{merMod} classes), multinomial and ordinal logistic
regressions (\code{multinom}, \code{polr}, \code{clm}, and \code{clmm}
classes), robust regression models (\code{rlm} class), nonlinear models
(\code{nls} class), and nonlinear mixed models (\code{nlme} and
\code{nlmerMod} classes).  The package also supports various models of
\code{unmarkedFit} and \code{maxLikeFit} classes estimating demographic
parameters after accounting for imperfect detection probabilities.  Some
functions also allow the creation of model selection tables for Bayesian
models of the \code{bugs} and \code{rjags} classes.
}
\details{
\tabular{ll}{
Package: \tab AICcmodavg\cr
Type: \tab Package\cr
Version: \tab 2.00\cr
Date: \tab 2014-07-14\cr
License: \tab GPL (>=2 )\cr
LazyLoad: \tab yes\cr
}
This package contains several useful functions for model selection and
multimodel inference:
\itemize{
\item{\code{\link{AICc}}}  {Computes AIC, AICc, and their
  quasi-likelihood counterparts (QAIC, QAICc).}
\item{\code{\link{aictab}}}  {Constructs model selection tables with
  number of parameters, AIC, delta AIC, Akaike weights or variants based
  on other AICc, QAIC, and QAICc for a set of candidate models.}
\item{\code{\link{boot.wt}}}  {Computes model selection relative
  frequencies based on the bootstrap.}
\item{\code{\link{c_hat}}}  {Computes an estimate of variance inflation
  factor for binomial or Poisson GLM's based on Pearson's chi-square.}
\item{\code{\link{confset}}}  {Determines the confidence set for the
  best model based on one of three criteria.}
\item{\code{\link{DIC}}}  {Extracts DIC.}
\item{\code{\link{dictab}}}  {Constructs model selection tables with
  number of parameters, DIC, delta DIC, DIC weights for a set of
  candidate models.}
\item{\code{\link{evidence}}}  {Computes the evidence ratio between the
  highest-ranked model based on the information criteria selected and a
  lower-ranked model.}
\item{\code{\link{extractLL}}}  {Extracts log-likelihood from
  models of certain classes.}
\item{\code{\link{extractSE}}}  {Extracts standard errors from models of
  certain classes and adds the labels.}
\item{\code{\link{fam.link.mer}}}  {Extracts the distribution family and 
  link function from a generalized linear mixed model of classes \code{mer}
  and \code{merMod}.} 
\item{\code{\link{importance}}}  {Computes importance values (w+) for
  the support of a given parameter among set of candidate models.}
\item{\code{\link{mb.gof.test}}}  {Computes the MacKenzie and Bailey
  goodness-of-fit test for single season occupancy models using the
  Pearson chi-square statistic.}
\item{\code{\link{modavg}}}  {Computes model-averaged estimate,
  unconditional standard error, and unconditional confidence interval of
  a parameter of interest among a set of candidate models.}
\item{\code{\link{modavgEffect}}}  {Computes model-averaged effect
  sizes between groups based on the entire candidate model set.}
\item{\code{\link{modavgShrink}}}  {Computes shrinkage version of
  model-averaged estimate, unconditional standard error, and
  unconditional confidence interval of a parameter of interest among a
  set of candidate models.} 
\item{\code{\link{modavgPred}}}  {Computes model-average predictions and 
  unconditional SE's among entire set of candidate models.}
\item{\code{\link{multComp}}}  {Performs multiple comparisons in a
  model selection framework.}
\item{\code{\link{Nmix.gof.test}}}  {Computes goodness-of-fit test for
  \emph{N}-mixture models based on the Pearson chi-square statistic.}
\item{\code{\link{predictSE}}}  {Computes predictions and associated 
  standard errors models of certain classes.}
}
}
\author{
  Marc J. Mazerolle <marc.mazerolle@uqat.ca>.
}
\references{
  Anderson, D. R. (2008) \emph{Model-based inference in the life sciences:
    a primer on evidence}. Springer: New York.

  Burnham, K. P., and Anderson, D. R. (2002) \emph{Model selection and
    multimodel inference: a practical information-theoretic approach}. Second
  edition. Springer: New York. 

  Burnham, K. P., Anderson, D. R. (2004) Multimodel inference:
  understanding AIC and BIC in model selection. \emph{Sociological
    Methods and Research} \bold{33}, 261--304.

  Mazerolle, M. J. (2006) Improving data analysis in herpetology: using
  Akaike's Information Criterion (AIC) to assess the strength of
  biological hypotheses. \emph{Amphibia-Reptilia} \bold{27}, 169--180. 
}

\examples{
##anuran larvae example from Mazerolle (2006) - Poisson GLM with offset
data(min.trap)
##assign "UPLAND" as the reference level as in Mazerolle (2006)          
min.trap$Type <- relevel(min.trap$Type, ref = "UPLAND") 

##set up candidate models          
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
c_hat(Cand.mod[[1]]) #uses Pearson's chi-square/df
##note the very low overdispersion: in this case, the analysis could be
##conducted without correcting for c-hat as its value is reasonably close
##to 1  

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

##compute model-averaged estimate of "TypeBOG"
modavg(cand.set = Cand.mod, modnames = Modnames, parm = "TypeBOG")

##compute model-averaged estimate of "TypeBOG" with shrinkage
##same number of models with vs without variable
modavgShrink(cand.set = Cand.mod, modnames = Modnames,
             parm = "TypeBOG")

##compute model-average predictions for two types of ponds
##create a data set for predictions
dat.pred <- data.frame(Type = factor(c("BOG", "UPLAND")),
                       log.Perimeter = mean(min.trap$log.Perimeter),
                       Num_ranatra = mean(min.trap$Num_ranatra),
                       Effort = mean(min.trap$Effort))

##model-averaged predictions across entire model set
modavgPred(cand.set = Cand.mod, modnames = Modnames,
           newdata = dat.pred)


\dontrun{
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
Cand.models <- list(fm1, fm2, fm3, fm4)
Modnames <- c("fm1", "fm2", "fm3", "fm4")

##compute table
print(aictab(cand.set = Cand.models, modnames = Modnames,
       second.ord = TRUE), digits = 4)

##compute evidence ratio
evidence(aictab(cand.set = Cand.models, modnames = Modnames))
##evidence ratio between top model vs lowest-ranked model
evidence(aictab(cand.set = Cand.models, modnames = Modnames), model.high = "fm2", model.low = "fm3")

##compute confidence set based on 'raw' method
confset(cand.set = Cand.models, modnames = Modnames, second.ord = TRUE,
        method = "raw")  

##compute importance value for "sitevar1" on occupancy
##same number of models with vs without variable
importance(cand.set = Cand.models, modnames = Modnames, parm = "sitevar1",
           parm.type = "psi") 

##compute model-averaged estimate of "sitevar1" on occupancy
modavg(cand.set = Cand.models, modnames = Modnames, parm = "sitevar1",
       parm.type = "psi")

##compute model-averaged estimate of "sitevar1" with shrinkage
##same number of models with vs without variable
modavgShrink(cand.set = Cand.models, modnames = Modnames,
             parm = "sitevar1", parm.type = "psi")

##compute model-average predictions for two types of ponds
##create a data set for predictions
dat.pred <- data.frame(sitevar1 = seq(from = min(siteCovs(pferUMF)$sitevar1),
                         to = max(siteCovs(pferUMF)$sitevar1), by = 0.5),
                       sitevar2 = mean(siteCovs(pferUMF)$sitevar2))

##model-averaged predictions of psi across range of values
##of sitevar1 and entire model set
modavgPred(cand.set = Cand.models, modnames = Modnames,
           newdata = dat.pred, parm.type = "psi")
detach(package:unmarked)
}
}
}

\keyword{ models }

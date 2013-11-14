\name{modavgpred}
\Rdversion{1.1}
\alias{modavgpred}
\alias{modavgpred.glm}
\alias{modavgpred.gls}
\alias{modavgpred.lme}
\alias{modavgpred.mer}
\alias{modavgpred.merMod}
\alias{modavgpred.rlm}
\alias{modavgpred.unmarked}
\alias{print.modavgpred}
\title{
Compute Model-averaged Predictions
}
\description{
This function computes the model-averaged predictions and unconditional
standard errors based on the entire candidate model set.  The function
is currently implemented for 'glm', 'gls', 'lm', 'lme', 'mer', and
'merMod' object classes that are stored in a list as well as various
models of 'unmarkedFit' classes.
}
\usage{
modavgpred(cand.set, modnames, newdata, type = "response", c.hat = 1,
           gamdisp = NULL, second.ord = TRUE, nobs = NULL,
           uncond.se = "revised", parm.type = NULL) 

modavgpred.glm(cand.set, modnames, newdata, type = "response",
               c.hat = 1, gamdisp = NULL, second.ord = TRUE,
               nobs = NULL, uncond.se = "revised") 

modavgpred.gls(cand.set, modnames, newdata, second.ord = TRUE,
               nobs = NULL, uncond.se = "revised")

modavgpred.lme(cand.set, modnames, newdata, second.ord = TRUE,
               nobs = NULL, uncond.se = "revised")

modavgpred.mer(cand.set, modnames, newdata, type = "response",
               c.hat = 1, second.ord = TRUE, nobs = NULL,
               uncond.se = "revised")

modavgpred.merMod(cand.set, modnames, newdata, type = "response",
                  c.hat = 1, second.ord = TRUE, nobs = NULL,
                  uncond.se = "revised")

modavgpred.rlm(cand.set, modnames, newdata, second.ord = TRUE,
               nobs = NULL, uncond.se = "revised")

modavgpred.unmarked(cand.set, modnames, newdata, second.ord = TRUE,
                    type = "response", c.hat = 1, nobs = NULL,
                    uncond.se = "revised", parm.type = NULL)

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
  \item{newdata}{
a data frame with the same structure as that of the original data frame
for which we want to make predictions.
}
  \item{type}{
the scale of prediction requested, one of "response" or "link".  The
latter is only relevant for 'glm', 'mer', and 'unmarkedFit' classes.
Note that the value "terms" is not defined for 'modavgpred'.
}
  \item{c.hat}{
value of overdispersion parameter (i.e., variance inflation factor) such
as that obtained from 'c_hat'.  Note that values of c.hat different from
1 are only appropriate for binomial GLM's with trials > 1 (i.e.,
success/trial or cbind(success, failure) syntax), with Poisson GLM's, 
single-season occupancy models (MacKenzie et al. 2002), or
\emph{N}-mixture models (Royle 2004, Dail and Madsen 2011).  If
c.hat > 1, 'modavgpred' will return the quasi-likelihood analogue of the
information criteria requested and multiply the variance-covariance
matrix of the estimates by this value (i.e., SE's are multiplied by
sqrt(c.hat)).  This option is not supported for generalized linear mixed
models of the 'mer' class.
}
\item{gamdisp}{
  the value of the gamma dispersion parameter.
}
\item{second.ord}{
  logical.  If TRUE, the function returns the second-order
  Akaike information criterion (i.e., AICc).
}
  \item{nobs}{
  this argument allows to specify a numeric value other than total
  sample size to compute the AICc (i.e., 'nobs' defaults to total number
  of observations).  This is relevant only for mixed models or various
  models of 'unmarkedFit' classes where sample size is not
  straightforward.  In such cases, one might use total number of
  observations or number of independent clusters (e.g., sites) as the
  value of 'nobs'.
}
  \item{uncond.se}{either, "old", or "revised", specifying the equation
used to compute the unconditional standard error of a model-averaged
estimate.  With uncond.se = "old", computations are based on equation
4.9 of Burnham and Anderson (2002), which was the former way to compute
unconditional standard errors.  With uncond.se = "revised", equation
6.12 of Burnham and Anderson (2002) is used.  Anderson (2008, p. 111)
recommends use of the revised version for the computation of
unconditional standard errors and it is now the default.  Note that
versions of package AICcmodavg < 1.04 used the old method to compute
unconditional standard errors. 
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

\details{
  'modavgpred' is a function that calls 'modavgpred.gls',
  'modavgpred.glm', 'modavgpred.lme', 'modavgpred.mer',
  'modavgpred.merMod', 'modavgpred.rlm', or 'modavgpred.unmarked'
  depending on the class of the object.  The candidate models must be
  stored in a list. Note that a data frame from which to make
  predictions must be supplied with the 'newdata' argument and that all
  variables appearing in the model set must appear in this data frame.
  Variables must be of the same type as in the original analysis (e.g.,
  factor, numeric).
  
  One can compute unconditional confidence intervals around the
  predictions from the elements returned by 'modavgpred'.  The classic
  computation based on asymptotic normality of the estimator is
  appropriate to estimate confidence intervals of beta estimates (i.e.,
  estimates on the linear predictor scale).  For predictions of some types
  of response variables (e.g., discrete values such as counts, or binary
  variables), the normal approximation may be inappropriate.  In such
  cases, it is often better to compute the confidence intervals on the
  linear predictor scale and then back-transform the limits to the scale
  of the response variable.  Burnham and Anderson (2002, p. 164) suggest
  alternative methods of computing unconditional confidence intervals for
  small degrees of freedom with profile likelihood intervals or
  bootstrapping.
}
\value{
'modavgpred' returns an object of class 'modavgpred' with the
following components:
\item{type}{
  the scale of predicted values (response or link) for 'glm', 'mer',
  'merMod', or 'unmarkedFit' classes.
}
\item{mod.avg.pred}{
  the model-averaged prediction over the entire candidate model set.
}
\item{uncond.se}{
  the unconditional standard error of each model-averaged prediction.
}

}
\references{
Anderson, D. R. (2008) \emph{Model-based Inference in the Life Sciences:
  a primer on evidence}. Springer: New York.

Burnham, K. P., Anderson, D. R. (2002) \emph{Model Selection and
  Multimodel Inference: a practical information-theoretic
  approach}. Second edition. Springer: New York.  

Dail, D., Madsen, L. (2011) Models for estimating abundance from 
repeated counts of an open population. \emph{Biometrics} \bold{67},
577--587. 

MacKenzie, D. I., Nichols, J. D., Lachman, G. B., Droege, S., Royle,
J. A., Langtimm, C. A. (2002) Estimating site occupancy rates when
detection probabilities are less than one. \emph{Ecology} \bold{83},
2248--2255.

Royle, J. A. (2004) \emph{N}-mixture models for estimating population
size from spatially replicated counts. \emph{Biometrics} \bold{60},
108--115.
}
\author{
Marc J. Mazerolle
}

\seealso{
\code{\link{AICc}}, \code{\link{aictab}}, \code{\link{importance}},
\code{\link{c_hat}}, \code{\link{confset}}, \code{\link{evidence}}, 
\code{\link{modavg}}, \code{\link{modavg.effect}},
\code{\link{modavg.shrink}}, \code{\link{predict}},
\code{\link{predict.glm}}, \code{\link{predictSE.gls}},
\code{\link{predictSE.lme}}, \code{\link{predictSE.mer}} 
}

\examples{
##example from subset of models in Table 1 in Mazerolle (2006)
data(dry.frog)

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
Cand.models[[4]] <- lm(log_Mass_lost ~ Shade + cent_Initial_mass +
                       Initial_mass2, data = dry.frog)
Cand.models[[5]] <- lm(log_Mass_lost ~ Substrate + cent_Initial_mass +
                       Initial_mass2, data = dry.frog)

##setup model names
Modnames <- paste("mod", 1:length(Cand.models), sep = "")

##compute model-averaged value and unconditional SE of predicted log of
##mass lost for frogs of average mass in shade for each substrate type

##first create data set to use for predictions
new.dat <- data.frame(Shade = c(1, 1, 1),
                      cent_Initial_mass = c(0, 0, 0),
                      Initial_mass2 = c(0, 0, 0),
                      Substrate = c("SOIL", "SPHAGNUM", "PEAT")) 

##compare unconditional SE's using both methods
modavgpred(cand.set = Cand.models, modnames = Modnames,
           newdata = new.dat, type = "response", uncond.se = "old")
modavgpred(cand.set = Cand.models, modnames = Modnames,
           newdata = new.dat, type = "response", uncond.se = "revised")
##round to 4 digits after decimal point
print(modavgpred(cand.set = Cand.models, modnames = Modnames,
                 newdata = new.dat, type = "response",
                 uncond.se = "revised"), digits = 4)



##Gamma glm
##clotting data example from 'gamma.shape' in MASS package of
##Venables and Ripley (2002, Modern applied statistics with
##S. Springer-Verlag: New York.)
clotting <- data.frame(u = c(5, 10, 15, 20, 30, 40, 60, 80, 100),
                       lot1 = c(118, 58, 42, 35, 27, 25, 21, 19, 18),
                       lot2 = c(69, 35, 26, 21, 18, 16, 13, 12, 12))
clot1 <- glm(lot1 ~ log(u), data = clotting, family = Gamma)

library(MASS)
gamma.dispersion(clot1) #dispersion parameter
gamma.shape(clot1) #reciprocal of dispersion parameter ==
##shape parameter 
summary(clot1, dispersion = gamma.dispersion(clot1))  #better

##create list with models
Cand <- list( )
Cand[[1]] <- glm(lot1 ~ log(u), data = clotting, family = Gamma)
Cand[[2]] <- glm(lot1 ~ 1, data = clotting, family = Gamma)

##create vector of model names
Modnames <- paste("mod", 1:length(Cand), sep = "")

##compute model-averaged predictions on scale of response variable for
##all observations
modavgpred(cand.set = Cand, modnames = Modnames, newdata = clotting,
           gamdisp = gamma.dispersion(clot1), type = "response") 

##compute model-averaged predictions on scale of linear predictor
modavgpred(cand.set = Cand, modnames = Modnames, newdata = clotting,
           gamdisp = gamma.dispersion(clot1), type = "link")

##compute model-averaged predictions on scale of linear predictor
\dontrun{
modavgpred(cand.set = Cand, modnames = Modnames, newdata = clotting,
           gamdisp = gamma.dispersion(clot1), type = "terms") #returns an error
##because type = "terms" is not defined for 'modavgpred'
}
\dontrun{
modavgpred(cand.set = Cand, modnames = Modnames, newdata = clotting,
           type = "terms") #returns an error because
##no gamma dispersion parameter was specified (i.e., 'gamdisp' missing) 
}


##example of model-averaged predictions from N-mixture model
##each variable appears twice in the models - this is a bit longer
\dontrun{
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
fm.mall.four <- pcount(~ ivel + date  ~ 1, mallardUMF, K = 30)

##model list
Cands <- list(fm.mall.one, fm.mall.two, fm.mall.three, fm.mall.four)
Modnames <- c("length + forest", "elev + forest", "length + elev",
              "null")

##compute model-averaged predictions of abundance for values of elev
modavgpred(cand.set = Cands, modnames = Modnames, newdata =
           data.frame(elev = seq(from = -1.4, to = 2.4, by = 0.1),
                      length = 0, forest = 0), parm.type = "lambda",
           type = "response")

##compute model-averaged predictions of detection for values of ivel
modavgpred(cand.set = Cands, modnames = Modnames, newdata =
           data.frame(ivel = seq(from = -1.75, to = 5.9, by = 0.5),
                      date = 0), parm.type = "detect",
           type = "response")
}



##example of model-averaged abundance from distance model
##this is a bit longer
\dontrun{
data(linetran) #example from ?distsamp
     
ltUMF <- with(linetran, {
  unmarkedFrameDS(y = cbind(dc1, dc2, dc3, dc4),
                  siteCovs = data.frame(Length, area, habitat),
                  dist.breaks = c(0, 5, 10, 15, 20),
                  tlength = linetran$Length * 1000, survey = "line", unitsIn = "m")
})
     
## Half-normal detection function. Density output (log scale). No covariates.
fm1 <- distsamp(~ 1 ~ 1, ltUMF)
     
## Halfnormal. Covariates affecting both density and and detection.
fm2 <- distsamp(~area + habitat ~ habitat, ltUMF)

## Hazard function. Covariates affecting both density and and detection.
fm3 <- distsamp(~area + habitat ~ habitat, ltUMF, keyfun="hazard")

##assemble model list
Cands <- list(fm1, fm2, fm3)
Modnames <- paste("mod", 1:length(Cands), sep = "")

##model-average predictions on abundance
modavgpred(cand.set = Cands, modnames = Modnames, parm.type = "lambda", type = "link",
           newdata = data.frame(area = mean(linetran$area), habitat = c("A", "B")))
detach(package:unmarked)
}
}


##example using Orthodont data set from Pinheiro and Bates (2000)
require(nlme)

##set up candidate models
m1 <- gls(distance ~ age, correlation = corCompSymm(value = 0.5, form = ~ 1 | Subject),
          data = Orthodont, method = "ML")

m2 <- gls(distance ~ 1, correlation = corCompSymm(value = 0.5, form = ~ 1 | Subject),
          data = Orthodont, method = "ML")

##assemble in list
Cand.models <- list(m1, m2)
##model names
Modnames <- c("age effect", "null model")

##model selection table
aictab(cand.set = Cand.models, modnames = Modnames)

##model-averaged predictions
modavgpred(cand.set = Cand.models, modnames = Modnames, newdata =
data.frame(age = c(8, 10, 12, 14)))
detach(package:nlme)

}

\keyword{models}

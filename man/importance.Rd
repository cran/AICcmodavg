\name{importance}
\Rdversion{1.1}
\alias{importance}
\alias{importance.default}
\alias{importance.AICaov.lm}
\alias{importance.AICbetareg}
\alias{importance.AICsclm.clm}
\alias{importance.AICclmm}
\alias{importance.AICclm}
\alias{importance.AICclogit.coxph}
\alias{importance.AICcoxme}
\alias{importance.AICcoxph}
\alias{importance.AICglm.lm}
\alias{importance.AICglmmTMB}
\alias{importance.AICglmerMod}
\alias{importance.AIClmerModLmerTest}
\alias{importance.AICgls}
\alias{importance.AIChurdle}
\alias{importance.AIClm}
\alias{importance.AIClme}
\alias{importance.AIClmekin}
\alias{importance.AICmaxlikeFit.list}
\alias{importance.AICmer}
\alias{importance.AICmultinom.nnet}
\alias{importance.AICnegbin.glm.lm}
\alias{importance.AICnlmerMod}
\alias{importance.AICpolr}
\alias{importance.AICrlm.lm}
\alias{importance.AICsurvreg}
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
\alias{importance.AICunmarkedFitGPC}
\alias{importance.AICunmarkedFitOccuMulti}
\alias{importance.AICunmarkedFitOccuMS}
\alias{importance.AICunmarkedFitOccuTTD}
\alias{importance.AICunmarkedFitMMO}
\alias{importance.AICunmarkedFitDSO}
\alias{importance.AICunmarkedFitGOccu}
\alias{importance.AICunmarkedFitOccuComm}
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

\method{importance}{AICbetareg}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, ...)

\method{importance}{AICsclm.clm}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, ...)

\method{importance}{AICclm}(cand.set, parm, modnames = NULL,
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

\method{importance}{AIClmerModLmerTest}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, ...)

\method{importance}{AICglmmTMB}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, c.hat = 1, ...)

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

\method{importance}{AICnegbin.glm.lm}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, ...)

\method{importance}{AICnlmerMod}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, ...)

\method{importance}{AICpolr}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, ...)

\method{importance}{AICrlm.lm}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, ...)

\method{importance}{AICsurvreg}(cand.set, parm, modnames = NULL,
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

\method{importance}{AICunmarkedFitGPC}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, c.hat = 1, parm.type = NULL, ...)

\method{importance}{AICunmarkedFitOccuMulti}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, c.hat = 1, parm.type = NULL, ...)

\method{importance}{AICunmarkedFitOccuMS}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, c.hat = 1, parm.type = NULL,
       ...)

\method{importance}{AICunmarkedFitOccuTTD}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, c.hat = 1, parm.type = NULL,
       ...)

\method{importance}{AICunmarkedFitMMO}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, c.hat = 1, parm.type = NULL,
       ...)

\method{importance}{AICunmarkedFitDSO}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, c.hat = 1, parm.type = NULL,
       ...)

\method{importance}{AICunmarkedFitGOccu}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, c.hat = 1, parm.type = NULL,
       ...)

\method{importance}{AICunmarkedFitOccuComm}(cand.set, parm, modnames = NULL,
           second.ord = TRUE, nobs = NULL, c.hat = 1, parm.type = NULL,
       ...)

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
as that obtained from \code{c_hat}.  Note that values of c.hat different
from 1 are only appropriate for binomial GLM's with trials > 1 (i.e.,
success/trial or cbind(success, failure) syntax), with Poisson GLM's, 
single-season occupancy models (MacKenzie et al. 2002), dynamic
occupancy models (MacKenzie et al. 2003), or \emph{N}-mixture models
(Royle 2004, Dail and Madsen 2011).  If \code{c.hat} > 1,
\code{importance} will return the quasi-likelihood analogue of the
information criteria requested and multiply the variance-covariance
matrix of the estimates by this value (i.e., SE's are multiplied by
\code{sqrt(c.hat)}).  This option is not supported for generalized
linear mixed models of the \code{mer} or \code{merMod} classes.
}
\item{parm.type}{this argument specifies the parameter type on which
  the variable of interest will be extracted and is only relevant for
  models of \code{unmarkedFit} classes.  The character strings supported
  vary with the type of model fitted.  For \code{unmarkedFitOccu},
  \code{unmarkedFitOccuMulti}, \code{unmarkedFitOccuComm} objects,
  either \code{psi} or \code{detect} can be supplied to indicate whether
  the parameter is on occupancy or detectability, respectively.  For
  \code{unmarkedFitColExt} objects, possible values are \code{psi}, 
  \code{gamma}, \code{epsilon}, and \code{detect}, for parameters on
  occupancy in the inital year, colonization, extinction, and
  detectability, respectively.  For \code{unmarkedFitOccuTTD} objects,
  possible values are \code{psi}, \code{gamma}, \code{epsilon}, and
  \code{detect}, for parameters on occupancy in the inital year,
  colonization, extinction, and time-to-dection (lambda rate parameter),
  respectively.  For \code{unmarkedFitOccuFP} objects, one can specify
  \code{psi}, \code{detect}, \code{falsepos}, and \code{certain}, for
  occupancy, detectability, probability of assigning false-positives,
  and probability detections are certain, respectively.  For
  \code{unmarkedFitOccuMS} objects, possible values are \code{psi},
  \code{phi}, or \code{detect}, denoting occupancy, transition, and
  detection probabilities, respectively.  For \code{unmarkedFitOccuRN}
  objects, either \code{lambda} or \code{detect} can be entered for
  abundance and detectability parameters, respectively.  For
  \code{unmarkedFitPCount} and \code{unmarkedFitMPois} objects,
  \code{lambda} or \code{detect} denote parameters on abundance and
  detectability, respectively.  For \code{unmarkedFitPCO},
  \code{unmarkedFitMMO}, and \code{unmarkedFitDSO} objects, one can
  enter \code{lambda}, \code{gamma}, \code{omega}, \code{iota}, or
  \code{detect}, to specify parameters on abundance, recruitment,
  apparent survival, immigration, and detectability, respectively.  For
  \code{unmarkedFitDS} objects, \code{lambda} and \code{detect} are
  supported.  For \code{unmarkedFitGDS}, \code{lambda}, \code{phi}, and
  \code{detect} denote abundance, availability, and detection
  probability, respectively.  For \code{unmarkedFitGMM} and
  \code{unmarkedFitGPC} objects, \code{lambda}, \code{phi}, and
  \code{detect} denote abundance, availability, and detectability,
  respectively.  For \code{unmarkedFitGOccu} objects, possible values
  are \code{psi}, \code{phi}, or \code{detect}, denoting occupancy,
  availability, and detection probabilities, respectively. 
}
\item{\dots}{
  additional arguments passed to the function.
}
}

\value{
\code{importance} returns an object of class \code{importance}
consisting of the following components:
\item{parm}{the parameter for which an importance value is required.} 
\item{w.plus}{the sum of Akaike weights for the models that include the
  parameter of interest.}
\item{w.minus}{the sum of Akaike weights for the models that exclude the
  parameter of interest.} 
}
\references{
Burnham, K. P., and Anderson, D. R. (2002) \emph{Model Selection and
  Multimodel Inference: a practical information-theoretic
  approach}. Second edition. Springer: New York.

Dail, D., Madsen, L. (2011) Models for estimating abundance from
  repeated counts of an open population. \emph{Biometrics} \bold{67},
  577--587. 
  
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
\dontrun{
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
detach(package:nlme)
}


##single-season occupancy model example modified from ?occu
\dontrun{
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


##example with multispecies occupancy modify from ?occuMulti
##Simulate 3 species data
N <- 80
nspecies <- 3
J <- 4
     
occ_covs <- as.data.frame(matrix(rnorm(N * 10),ncol=10))
names(occ_covs) <- paste('par',1:10,sep='')

det_covs <- list()
for (i in 1:nspecies){
    det_covs[[i]] <- matrix(rnorm(N*J),nrow=N)
}
names(det_covs) <- paste('par',1:nspecies,sep='')

##True vals
beta <- c(0.5,0.2,0.4,0.5,-0.1,-0.3,0.2,0.1,-1,0.1)
f1 <- beta[1] + beta[2]*occ_covs$par1
f2 <- beta[3] + beta[4]*occ_covs$par2
f3 <- beta[5] + beta[6]*occ_covs$par3
f4 <- beta[7]
f5 <- beta[8]
f6 <- beta[9]
f7 <- beta[10]
f <- cbind(f1,f2,f3,f4,f5,f6,f7)
z <- expand.grid(rep(list(1:0),nspecies))[,nspecies:1]
colnames(z) <- paste('sp',1:nspecies,sep='')
dm <- model.matrix(as.formula(paste0("~.^",nspecies,"-1")),z)

psi <- exp(f %*% t(dm))
psi <- psi/rowSums(psi)

##True state
ztruth <- matrix(NA,nrow=N,ncol=nspecies)
for (i in 1:N){
    ztruth[i,] <- as.matrix(z[sample(8,1,prob=psi[i,]),])
}
     
p_true <- c(0.6,0.7,0.5)

## fake y data
y <- list()
     
for (i in 1:nspecies){
    y[[i]] <- matrix(NA,N,J)
    for (j in 1:N){
        for (k in 1:J){
            y[[i]][j,k] <- rbinom(1,1,ztruth[j,i]*p_true[i])
        }
    }
}
names(y) <- c('coyote','tiger','bear')
     
##Create the unmarked data object
data <- unmarkedFrameOccuMulti(y=y,siteCovs=occ_covs,obsCovs=det_covs)

## Formulas for state and detection processes
## Length should match number/order of columns in fDesign
occFormulas <- c('~par1 + par2','~par2','~par3','~1','~1','~1','~1')
occFormulas2 <- c('~par1 + par3','~par1 + par2','~par1 + par2 + par3',
                  "~ 1", "~1", "~ 1", "~1")

##Length should match number/order of species in data@ylist
detFormulas <- c('~1','~1','~1')
     
fit <- occuMulti(detFormulas,occFormulas,data)
fit2 <- occuMulti(detFormulas,occFormulas2,data)

##importance
importance(cand.set = list(fit, fit2), parm = "[coyote] par2",
           parm.type = "psi")


detach(package:unmarked)
}
}

\keyword{models}

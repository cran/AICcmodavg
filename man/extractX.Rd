\name{extractX}
\Rdversion{1.1}
\alias{extractX}
\alias{extractX.default}
\alias{extractX.AICaov.lm}
\alias{extractX.AICglm.lm}
\alias{extractX.AICglmmTMB}
\alias{extractX.AICgls}
\alias{extractX.AIClm}
\alias{extractX.AIClme}
\alias{extractX.AICglmerMod}
\alias{extractX.AIClmerMod}
\alias{extractX.AIClmerModLmerTest}
\alias{extractX.AICrlm.lm}
\alias{extractX.AICsurvreg}
\alias{extractX.AICunmarkedFitOccu}
\alias{extractX.AICunmarkedFitColExt}
\alias{extractX.AICunmarkedFitOccuRN}
\alias{extractX.AICunmarkedFitPCount}
\alias{extractX.AICunmarkedFitPCO}
\alias{extractX.AICunmarkedFitDS}
\alias{extractX.AICunmarkedFitGDS}
\alias{extractX.AICunmarkedFitOccuFP}
\alias{extractX.AICunmarkedFitMPois}
\alias{extractX.AICunmarkedFitGMM}
\alias{extractX.AICunmarkedFitGPC}
\alias{extractX.AICunmarkedFitOccuMulti}
\alias{extractX.AICunmarkedFitOccuMS}
\alias{extractX.AICunmarkedFitOccuTTD}
\alias{extractX.AICunmarkedFitMMO}
\alias{extractX.AICunmarkedFitDSO}
\alias{extractX.AICunmarkedFitGOccu}
\alias{extractX.AICunmarkedFitOccuComm}
\alias{print.extractX}
\title{
Extract Predictors from Candidate Model List
}
\description{
This function extracts the predictors used in candidate models.  The
function is currently implemented for \code{glm}, \code{glmmTMB},
\code{gls}, \code{lm}, \code{lme}, \code{merMod},
\code{lmerModLmerTest}, \code{rlm}, \code{survreg} object classes that
are stored in a list as well as various models of \code{unmarkedFit}
classes.
}
\usage{
extractX(cand.set, \dots)

\method{extractX}{AICaov.lm}(cand.set, \dots) 

\method{extractX}{AICglm.lm}(cand.set, \dots)

\method{extractX}{AICglmmTMB}(cand.set, \dots)

\method{extractX}{AIClm}(cand.set, \dots) 

\method{extractX}{AICgls}(cand.set, \dots) 

\method{extractX}{AIClme}(cand.set, \dots)

\method{extractX}{AICglmerMod}(cand.set, \dots)

\method{extractX}{AIClmerMod}(cand.set, \dots)

\method{extractX}{AIClmerModLmerTest}(cand.set, \dots)

\method{extractX}{AICrlm.lm}(cand.set, \dots) 

\method{extractX}{AICsurvreg}(cand.set, \dots) 

\method{extractX}{AICunmarkedFitOccu}(cand.set,
        parm.type = NULL, \dots)

\method{extractX}{AICunmarkedFitColExt}(cand.set,
        parm.type = NULL, \dots)  

\method{extractX}{AICunmarkedFitOccuRN}(cand.set,
        parm.type = NULL, \dots)

\method{extractX}{AICunmarkedFitPCount}(cand.set,
        parm.type = NULL, \dots) 

\method{extractX}{AICunmarkedFitPCO}(cand.set,
        parm.type = NULL, \dots)

\method{extractX}{AICunmarkedFitDS}(cand.set,
        parm.type = NULL, \dots) 

\method{extractX}{AICunmarkedFitGDS}(cand.set,
        parm.type = NULL, \dots)

\method{extractX}{AICunmarkedFitOccuFP}(cand.set,
        parm.type = NULL, \dots)

\method{extractX}{AICunmarkedFitMPois}(cand.set,
        parm.type = NULL, \dots)

\method{extractX}{AICunmarkedFitGMM}(cand.set,
        parm.type = NULL, \dots)

\method{extractX}{AICunmarkedFitGPC}(cand.set,
        parm.type = NULL, \dots)

\method{extractX}{AICunmarkedFitOccuMulti}(cand.set,
        parm.type = NULL, \dots)

\method{extractX}{AICunmarkedFitOccuMS}(cand.set,
        parm.type = NULL, \dots)

\method{extractX}{AICunmarkedFitOccuTTD}(cand.set,
        parm.type = NULL, \dots)

\method{extractX}{AICunmarkedFitMMO}(cand.set,
        parm.type = NULL, \dots)

\method{extractX}{AICunmarkedFitDSO}(cand.set,
        parm.type = NULL, \dots)

\method{extractX}{AICunmarkedFitGOccu}(cand.set,
        parm.type = NULL, \dots)

\method{extractX}{AICunmarkedFitOccuComm}(cand.set,
        parm.type = NULL, \dots)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cand.set}{
    a list storing each of the models in the candidate model set.
  }
  \item{parm.type}{this argument specifies the parameter type on which
    the predictors will be extracted and is only relevant for models of
    \code{unmarkedFit} classes.  The character strings supported
    vary with the type of model fitted.  For \code{unmarkedFitOccu}, 
    \code{unmarkedFitOccuMulti}, and \code{unmarkedFitOccuComm} objects,
    either \code{psi} or \code{detect} can be supplied to indicate
    whether the parameter is on occupancy or detectability,
    respectively.  For \code{unmarkedFitColExt} and
    \code{unmarkedFitOccuTTD}, possible values  are \code{psi},
    \code{gamma}, \code{epsilon}, and \code{detect}, for parameters on
    occupancy in the inital year, colonization, extinction, and
    detectability, respectively.  For \code{unmarkedFitOccuFP} objects,
    one can specify \code{psi}, \code{detect}, \code{falsepos}, and
    \code{certain}, for occupancy, detectability, probability of
    assigning false-positives, and probability detections are certain,
    respectively.  For \code{unmarkedFitOccuMS} objects, possible values
    are \code{psi}, \code{phi}, or \code{detect}, denoting occupancy,
    transition, and detection probabilities, respectively.  For
    \code{unmarkedFitOccuRN} objects, either \code{lambda} or
    \code{detect} can be entered for abundance and detectability
    parameters, respectively.  For \code{unmarkedFitPCount} and
    \code{unmarkedFitMPois} objects, \code{lambda} or \code{detect}
    denote parameters on abundance and detectability, respectively.  For
    \code{unmarkedFitPCO}, \code{unmarkedFitMMO}, and
    \code{unmarkedFitDSO} objects, one can enter \code{lambda},
    \code{gamma}, \code{omega}, \code{iota}, or \code{detect}, to
    specify parameters on abundance, recruitment, apparent survival,
    immigration, and detectability, respectively.  For
    \code{unmarkedFitDS} objects, \code{lambda} and \code{detect} are
    supported.  For \code{unmarkedFitGDS}, \code{lambda}, \code{phi},
    and \code{detect} denote abundance, availability, and detection
    probability, respectively.  For \code{unmarkedFitGMM} and
    \code{unmarkedFitGPC} objects, \code{lambda}, \code{phi}, and
    \code{detect} denote abundance, availability, and detectability,
    respectively. For \code{unmarkedFitGOccu} objects, possible values
    are \code{psi}, \code{phi}, or \code{detect}, denoting occupancy,
    availability, and detection probabilities, respectively.  
}
\item{\dots}{
  additional arguments passed to the function.
}
}

\details{
  The candidate models must be stored in a list.  The results of
  \code{extractX} are useful in preparing a \code{newdata}
  data frame to use in computing model-averaged predictions with
  \code{modavgPred} or differences between groups with
  \code{modavgEffect} (Burnham and Anderson 2002, Anderson 2008, Burnham
  et al. 2011).
}
\value{
\code{extractX} returns an object of class \code{extractX} with the 
following components:
\item{predictors}{
  a character vector of the names of the predictors included in the
  model, excluding the intercept term. 
}
\item{data}{
  a data frame or, in the case of \code{unmarkedFit} objects, a list of
  data frames (e.g., obsCovs, siteCovs, yearlySiteCovs).
}
}
\references{
Anderson, D. R. (2008) \emph{Model-based Inference in the Life Sciences:
a primer on evidence}. Springer: New York.

Burnham, K. P., Anderson, D. R. (2002) \emph{Model Selection and
  Multimodel Inference: a practical information-theoretic
  approach}. Second edition. Springer: New York.  

Burnham, K. P., Anderson, D. R., Huyvaert, K. P. (2011) AIC model
selection and multimodel inference in behaviorial ecology: some
background, observations and comparisons. \emph{Behavioral Ecology and 
  Sociobiology} \bold{65}, 23--25.

Mazerolle, M. J. (2006) Improving data analysis in herpetology: using
Akaike's Information Criterion (AIC) to assess the strength of
biological hypotheses. \emph{Amphibia-Reptilia} \bold{27}, 169--180.

Pinheiro, J. C., Bates, D. M. (2000). \emph{Mixed-effects Models in S and
  S-PLUS}. Springer Verlag: New York.

Royle, J. A. (2004) \emph{N}-mixture models for estimating population
size from spatially replicated counts. \emph{Biometrics} \bold{60},
108--115.
}
\author{
Marc J. Mazerolle
}

\seealso{
  \code{\link{extractCN}}, \code{\link{extractSE}},
  \code{\link{modavgPred}}, \code{\link{modavgCustom}},
  \code{\link{modavgEffect}}, \code{\link{predict}},
  \code{\link{predictSE}}
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
##assign names
names(Cand.models) <- paste(1:length(Cand.models))

##extract predictors from candidate model set
orig.data <- extractX(cand.set = Cand.models)
orig.data
str(orig.data)

\dontrun{
##model-averaged prediction with original variables
modavgPred(Cand.models, newdata = orig.data$data)
}


##example of model-averaged predictions from N-mixture model (e.g., Royle 2004)
##modified from ?pcount
##each variable appears twice on lambda in the models
\dontrun{
require(unmarked)    
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
names(Cands) <- c("length + forest", "elev + forest", "length + elev",
                  "null")

##extract predictors on lambda
lam.dat <- extractX(cand.set = Cands, parm.type = "lambda")
lam.dat
str(lam.dat)

##extract predictors on detectability
extractX(cand.set = Cands, parm.type = "detect")

##model-averaged predictions on lambda
##extract data
siteCovs <- lam.dat$data$siteCovs
##create vector of forest values
forest <- seq(min(siteCovs$forest),
              max(siteCovs$forest),
              length.out = 40)
dframe <- data.frame(forest = forest,
                     length = mean(siteCovs$length),
                     elev = mean(siteCovs$elev))
modavgPred(Cands, parm.type = "lambda",
           newdata = dframe)
detach(package:unmarked)
}


##example of model-averaged abundance from distance model
\dontrun{
require(unmarked)
data(linetran) #example from ?distsamp
     
ltUMF <- with(linetran, {
  unmarkedFrameDS(y = cbind(dc1, dc2, dc3, dc4),
                  siteCovs = data.frame(Length, area, habitat),
                  dist.breaks = c(0, 5, 10, 15, 20),
                  tlength = linetran$Length * 1000, survey = "line",
                  unitsIn = "m")
})
     
## Half-normal detection function. Density output (log scale). No covariates.
fm1 <- distsamp(~ 1 ~ 1, ltUMF)
     
## Halfnormal. Covariates affecting both density and and detection.
fm2 <- distsamp(~area + habitat ~ habitat, ltUMF)

## Hazard function. Covariates affecting both density and and detection.
fm3 <- distsamp(~area + habitat ~ habitat, ltUMF, keyfun="hazard")

##assemble model list
Cands <- list(fm1, fm2, fm3)

##model-average predictions on abundance
extractX(cand.set = Cands, parm.type = "lambda")
detach(package:unmarked)
}



##example using Orthodont data set from Pinheiro and Bates (2000)
\dontrun{
require(nlme)

##set up candidate models
m1 <- gls(distance ~ age, correlation = corCompSymm(value = 0.5, form = ~ 1 | Subject),
          data = Orthodont, method = "ML")

m2 <- gls(distance ~ 1, correlation = corCompSymm(value = 0.5, form = ~ 1 | Subject),
          data = Orthodont, method = "ML")

##assemble in list
Cand.models <- list("age effect" = m1, "null model" = m2)

##model-averaged predictions
extractX(cand.set = Cand.models)
detach(package:nlme)
}
}
\keyword{models}

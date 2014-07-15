\name{modavg}
\Rdversion{1.1}
\alias{modavg}
\alias{modavg.default}
\alias{modavg.AICaov.lm}
\alias{modavg.AICsclm.clm}
\alias{modavg.AICclmm}
\alias{modavg.AICcoxme}
\alias{modavg.AICcoxph}
\alias{modavg.AICglm.lm}
\alias{modavg.AICgls}
\alias{modavg.AIClm}
\alias{modavg.AIClme}
\alias{modavg.AIClmekin}
\alias{modavg.AICmaxlikeFit.list}
\alias{modavg.AICmer}
\alias{modavg.AIClmerMod}
\alias{modavg.AICglmerMod}
\alias{modavg.AICmultinom.nnet}
\alias{modavg.AICpolr}
\alias{modavg.AICrlm.lm}
\alias{modavg.AICvglm}
\alias{modavg.AICzeroinfl}
\alias{modavg.AICunmarkedFitOccu}
\alias{modavg.AICunmarkedFitColExt}
\alias{modavg.AICunmarkedFitOccuRN}
\alias{modavg.AICunmarkedFitPCount}
\alias{modavg.AICunmarkedFitPCO}
\alias{modavg.AICunmarkedFitDS}
\alias{modavg.AICunmarkedFitGDS}
\alias{modavg.AICunmarkedFitOccuFP}
\alias{modavg.AICunmarkedFitMPois}
\alias{modavg.AICunmarkedFitFitGMM}
\alias{modavg.AICunmarkedFitFitGPC}
\alias{print.modavg}

\title{
Compute Model-averaged Parameter Estimate (Multimodel Inference)
}
\description{
This function model-averages the estimate of a parameter of interest
among a set of candidate models, computes the unconditional standard
error and unconditional confidence intervals as described in Buckland et
al. (1997) and Burnham and Anderson (2002). This model-averaged estimate
is also referred to as a natural average of the estimate by Burnham and
Anderson (2002, p. 152).
}
\usage{
modavg(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, 
       uncond.se = "revised", conf.level = 0.95, exclude = NULL, warn =
       TRUE, \dots) 

\method{modavg}{AICaov.lm}(cand.set, parm, modnames = NULL, second.ord =
        TRUE, nobs = NULL, uncond.se = "revised", conf.level = 0.95,
        exclude = NULL, warn = TRUE, \dots)

\method{modavg}{AICsclm.clm}(cand.set, parm, modnames = NULL,
        second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, exclude = NULL, warn = TRUE, \dots)

\method{modavg}{AICclmm}(cand.set, parm, modnames = NULL, second.ord 
        = TRUE, nobs = NULL, uncond.se = "revised", conf.level = 0.95,
        exclude = NULL, warn = TRUE, \dots)

\method{modavg}{AICcoxme}(cand.set, parm, modnames = NULL, second.ord
        = TRUE, nobs = NULL, uncond.se = "revised", conf.level = 0.95,
        exclude = NULL, warn = TRUE, \dots)

\method{modavg}{AICcoxph}(cand.set, parm, modnames = NULL, second.ord
        = TRUE, nobs = NULL, uncond.se = "revised", conf.level = 0.95,
        exclude = NULL, warn = TRUE, \dots)

\method{modavg}{AICglm.lm}(cand.set, parm, modnames = NULL,
        second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, exclude = NULL, warn = TRUE, c.hat = 1,
        gamdisp = NULL, \dots)

\method{modavg}{AICgls}(cand.set, parm, modnames = NULL, second.ord =
           TRUE, nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           exclude = NULL, warn = TRUE, \dots)

\method{modavg}{AIClm}(cand.set, parm, modnames = NULL, second.ord =
        TRUE, nobs = NULL, uncond.se = "revised", conf.level = 0.95,
        exclude = NULL, warn = TRUE, \dots)

\method{modavg}{AIClme}(cand.set, parm, modnames = NULL, second.ord =
        TRUE, nobs = NULL, uncond.se = "revised", conf.level = 0.95,
        exclude = NULL, warn = TRUE, \dots)

\method{modavg}{AIClmekin}(cand.set, parm, modnames = NULL,
         second.ord = TRUE, nobs = NULL, uncond.se = "revised",
         conf.level = 0.95, exclude = NULL, warn = TRUE, \dots)

\method{modavg}{AICmaxlikeFit.list}(cand.set, parm, modnames = NULL,
        second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, exclude = NULL, warn = TRUE, c.hat = 1,
        \dots)

\method{modavg}{AICmer}(cand.set, parm, modnames = NULL, second.ord =
        TRUE, nobs = NULL, uncond.se = "revised", conf.level = 0.95,
        exclude = NULL, warn = TRUE, \dots)

\method{modavg}{AIClmerMod}(cand.set, parm, modnames = NULL,
        second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, exclude = NULL, warn = TRUE, \dots)

\method{modavg}{AICglmerMod}(cand.set, parm, modnames = NULL,
        second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, exclude = NULL, warn = TRUE, \dots)

\method{modavg}{AICmultinom.nnet}(cand.set, parm, modnames = NULL, 
        second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, exclude = NULL, warn = TRUE, c.hat = 1,
        \dots)

\method{modavg}{AICpolr}(cand.set, parm, modnames = NULL, second.ord
        = TRUE, nobs = NULL, uncond.se = "revised", conf.level = 0.95,
        exclude = NULL, warn = TRUE, \dots)

\method{modavg}{AICrlm.lm}(cand.set, parm, modnames = NULL,
        second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, exclude = NULL, warn = TRUE, \dots)

\method{modavg}{AICvglm}(cand.set, parm, modnames = NULL, second.ord
         = TRUE, nobs = NULL, uncond.se = "revised", conf.level = 0.95,
         exclude = NULL, warn = TRUE, c.hat = 1, \dots)

\method{modavg}{AICzeroinfl}(cand.set, parm, modnames = NULL,
         second.ord = TRUE, nobs = NULL, uncond.se = "revised",
         conf.level = 0.95, exclude = NULL, warn = TRUE, \dots)

\method{modavg}{AICunmarkedFitOccu}(cand.set, parm, modnames = NULL, 
        second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, exclude = NULL, warn = TRUE, c.hat = 1,
        parm.type = NULL, \dots)

\method{modavg}{AICunmarkedFitColExt}(cand.set, parm, modnames =
        NULL, second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, exclude = NULL, warn = TRUE, c.hat = 1,
        parm.type = NULL, \dots)

\method{modavg}{AICunmarkedFitOccuRN}(cand.set, parm, modnames =
        NULL, second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, exclude = NULL, warn = TRUE, c.hat = 1,
        parm.type = NULL, \dots)

\method{modavg}{AICunmarkedFitPCount}(cand.set, parm, modnames =
        NULL, second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, exclude = NULL, warn = TRUE, c.hat = 1,
        parm.type = NULL, \dots)

\method{modavg}{AICunmarkedFitPCO}(cand.set, parm, modnames = NULL,
        second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, exclude = NULL, warn = TRUE, c.hat = 1,
        parm.type = NULL, \dots)

\method{modavg}{AICunmarkedFitDS}(cand.set, parm, modnames = NULL,
        second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, exclude = NULL, warn = TRUE, c.hat = 1,
        parm.type = NULL, \dots)

\method{modavg}{AICunmarkedFitGDS}(cand.set, parm, modnames = NULL,
        second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, exclude = NULL, warn = TRUE, c.hat = 1,
        parm.type = NULL, \dots)

\method{modavg}{AICunmarkedFitOccuFP}(cand.set, parm, modnames =
        NULL, second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, exclude = NULL, warn = TRUE, c.hat = 1,
        parm.type = NULL, \dots)

\method{modavg}{AICunmarkedFitMPois}(cand.set, parm, modnames =
        NULL, second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, exclude = NULL, warn = TRUE, c.hat = 1,
        parm.type = NULL, \dots)

\method{modavg}{AICunmarkedFitFitGMM}(cand.set, parm, modnames =
       NULL, second.ord = TRUE, nobs = NULL, uncond.se = "revised",
       conf.level = 0.95, exclude = NULL, warn = TRUE, c.hat = 1,
       parm.type = NULL, \dots)

\method{modavg}{AICunmarkedFitFitGPC}(cand.set, parm, modnames =
        NULL, second.ord = TRUE, nobs = NULL, uncond.se = "revised",
        conf.level = 0.95, exclude = NULL, warn = TRUE, c.hat = 1,
        parm.type = NULL, \dots)
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
\item{exclude}{
this argument excludes models based on the terms specified for the 
computation of a model-averaged estimate of \code{parm}.  The
\code{exclude} argument is set to \code{NULL} by default and does not
exclude any models other than those without the \code{parm}.  When
\code{parm} is a main effect but is also involved in
interactions/polynomial terms in some models, one should specify the
interaction/polynomial terms as a list to exclude models with these
terms from the computation of model-averaged estimate of the main effect
(e.g., \code{exclude = list("sex:mass", "mass2")}).  See 'Details' 
and 'Examples' below.
}
  \item{warn}{
logical.  If \code{TRUE}, \code{modavg} performs a check and isssues a
warning when the value in \code{parm} occurs more than once in any given
model.  This is a check for potential interaction/polynomial terms in
the model when such terms are constructed with the usual operators
(e.g., \code{I( )} for polynomial terms, \code{:} for interaction
terms).
}
  \item{c.hat}{
value of overdispersion parameter (i.e., variance inflation factor) such
as that obtained from \code{c_hat}.  Note that values of c.hat different
from 1 are only appropriate for binomial GLM's with trials > 1 (i.e.,
success/trial or cbind(success, failure) syntax), with Poisson GLM's, 
single-season occupancy models (MacKenzie et al. 2002), or
\emph{N}-mixture models (Royle 2004, Dail and Madsen 2011).  If
c.hat > 1, \code{modavg} will return the quasi-likelihood analogue of
the information criteria requested and multiply the variance-covariance
matrix of the estimates by this value (i.e., SE's are multiplied by
sqrt(c.hat)). This option is not supported for generalized linear mixed
models of the \code{mer} or \code{merMod} classes.
}
  \item{gamdisp}{
if gamma GLM is used, the dispersion parameter should be specified here
to apply the same value to each model.
}
\item{parm.type}{  this argument specifies the parameter type on which the effect size 
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
The parameter for which a model-averaged estimate is requested must be
specified with the \code{parm} argument and must be identical to its
label in the model output (e.g., from \code{summary}).  For factors, one
must specify the name of the variable and the level of interest.
\code{modavg} includes checks to find variations of interaction terms
specified in the \code{parm} and \code{exclude} arguments.  However, to
avoid problems, one should specify interaction terms consistently for
all models: e.g., either \code{a:b} or \code{b:a} for all models, but
not a mixture of both. 

You must exercise caution when some models include interaction or
polynomial terms, because main effect terms do not have the same
interpretation when they also appear in an interaction/polynomial term
in the same model.  In such cases, one should exclude models containing 
interaction terms where the main effect is involved with the
\code{exclude} argument of \code{modavg}.  Note that \code{modavg}
checks for potential cases of multiple instances of a variable appearing
more than once in a given model (presumably in an interaction) and
issues a warning.  To correctly compute the model-averaged estimate of a
main effect involved in interaction/polynomial terms, specify the
interaction terms(s) that should not appear in the same model with the
\code{exclude} argument.  This will effectively exclude models from the
computation of the model-averaged estimate.  

When \code{warn = TRUE}, \code{modavg} looks for matches among the
labels of the estimates with \code{identical}.  It then compares the
results to partial matches with \code{regexpr}, and issues a warning
whenever they are different.  As a result, \code{modavg} may issue a
warning when some variables or levels of categorical variables have
nested names (e.g., \code{treat}, \code{treat10}; \code{L}, \code{TL}).
When this warning is only due to the presence of similarly named
variables in the models (and NOT due to interaction terms), you can
suppress this warning by setting \code{warn = FALSE}.

\code{modavg} is implemented for a list containing objects of
\code{clm}, \code{clmm}, \code{clogit}, \code{coxme}, \code{coxph},
\code{glm}, \code{gls}, \code{lm}, \code{lme}, \code{lmekin},
\code{maxlikeFit}, \code{mer}, \code{glmerMod}, \code{lmerMod},
\code{multinom}, \code{polr}, \code{rlm}, \code{vglm}, \code{zeroinfl}
classes as well as various models of \code{unmarkedFit} classes.
}
\value{
\code{modavg} creates an object of class \code{modavg} with the following
components: 

\item{Parameter}{the parameter for which a model-averaged estimate was
  obtained} 
\item{Mod.avg.table}{the reduced model selection table based on models
  including the parameter of interest}
\item{Mod.avg.beta}{the model-averaged estimate based on all models
  including the parameter of interest (see 'Details' above regarding the
  exclusion of models where parameter of interest is involved in an
  interaction)}
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
\code{\link{modavgEffect}}, \code{\link{modavgShrink}},
\code{\link{modavgPred}} 
}
\examples{
##anuran larvae example from Mazerolle (2006)
data(min.trap)
##assign "UPLAND" as the reference level as in Mazerolle (2006)          
min.trap$Type <- relevel(min.trap$Type, ref = "UPLAND") 

##set up candidate models          
Cand.mod <- list( )
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

##compute model-averaged estimate of TypeBOG
modavg(parm = "TypeBOG", cand.set = Cand.mod, modnames = Modnames)
##round to 4 digits after decimal point
print(modavg(parm = "TypeBOG", cand.set = Cand.mod,
             modnames = Modnames), digits = 4)

##compute c-hat estimate based on residual deviance as in Mazerolle
##(2006)
Cand.mod[[1]]$deviance/Cand.mod[[1]]$df.residual

##compute model-averaged estimate of TypeBOG as in Table 4 of
##Mazerolle (2006) 
modavg(parm = "TypeBOG", cand.set = Cand.mod, modnames = Modnames,
       c.hat = 1.11) 



##example with similarly-named variables and interaction terms
set.seed(seed = 4)
resp <- rnorm(n = 40, mean = 3, sd = 1)
size <- rep(c("small", "medsmall", "high", "medhigh"), times = 10)
set.seed(seed = 4)
mass <- rnorm(n = 40, mean = 2, sd = 0.1)
mass2 <- mass^2
age <- rpois(n = 40, lambda = 3.2)
agecorr <- rpois(n = 40, lambda = 2) 
sizecat <- rep(c("a", "ab"), times = 20)
data1 <- data.frame(resp = resp, size = size, sizecat = sizecat,
                    mass = mass, mass2 = mass2, age = age,
                    agecorr = agecorr)

##set up models in list
Cand <- list( )
Cand[[1]] <- lm(resp ~ size + agecorr, data = data1)
Cand[[2]] <- lm(resp ~ size + mass + agecorr, data = data1)
Cand[[3]] <- lm(resp ~ age + mass, data = data1)
Cand[[4]] <- lm(resp ~ age + mass + mass2, data = data1)
Cand[[5]] <- lm(resp ~ mass + mass2 + size, data = data1)
Cand[[6]] <- lm(resp ~ mass + mass2 + sizecat, data = data1)
Cand[[7]] <- lm(resp ~ sizecat, data = data1)
Cand[[8]] <- lm(resp ~ sizecat + mass + sizecat:mass, data = data1)
Cand[[9]] <- lm(resp ~ agecorr + sizecat + mass + sizecat:mass,
                 data = data1) 

##create vector of model names
Modnames <- paste("mod", 1:length(Cand), sep = "")

aictab(cand.set = Cand, modnames = Modnames, sort = TRUE) #correct

##as expected, issues warning as mass occurs sometimes with "mass2" or
##"sizecatab:mass" in some of the models
\dontrun{modavg(cand.set = Cand, parm = "mass", modnames = Modnames)}

##no warning issued, because "age" and "agecorr" never appear in same model
modavg(cand.set = Cand, parm = "age", modnames = Modnames)

##as expected, issues warning because warn=FALSE, but it is a very bad
##idea in this example since "mass" occurs with "mass2" and "sizecat:mass"
##in some of the models - results are INCORRECT
\dontrun{modavg(cand.set = Cand, parm = "mass", modnames = Modnames,
                warn = FALSE)}

##correctly excludes models with quadratic term and interaction term
##results are CORRECT
modavg(cand.set = Cand, parm = "mass", modnames = Modnames,
       exclude = list("mass2", "sizecat:mass")) 

##correctly computes model-averaged estimate because no other parameter
##occurs simultaneously in any of the models
modavg(cand.set = Cand, parm = "sizesmall", modnames = Modnames) #correct

##as expected, issues a warning because "sizecatab" occurs sometimes in
##an interaction in some models
\dontrun{modavg(cand.set = Cand, parm = "sizecatab",
                modnames = Modnames) }

##exclude models with "sizecat:mass" interaction - results are CORRECT
modavg(cand.set = Cand, parm = "sizecatab", modnames = Modnames,
       exclude = list("sizecat:mass"))



##example with generalized linear mixed model
##modified example from ?glmer
require(lme4)
##create proportion of incidence
cbpp$prop <- cbpp$incidence/cbpp$size
##add bogus variable
cbpp$randx <- rnorm(n = nrow(cbpp), mean = 12, sd = 3)

##model with period
gm1 <- glmer(prop ~ period + (1 | herd), family = binomial,
             weights = size, data = cbpp)
gm2 <- glmer(prop ~ period + randx + (1 | herd), family = binomial,
             weights = size, data = cbpp)
##model without period
gm3 <- glmer(prop ~ 1 + (1 | herd), family = binomial,
             weights = size, data = cbpp)
Cands <- list(gm1, gm2, gm3)
Modnames <- c("period", "period + randx", "null")

##model selection
aictab(cand.set = Cands, modnames = Modnames)
##model average for difference between period 1 vs period 4
modavg(cand.set = Cands, modnames = Modnames, parm = "period4")
detach(package:lme4)


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
modavg(Cands, Mods, parm = "x")



##example with multiple-season occupancy model modified from ?colext
##this is a bit longer
\dontrun{
require(unmarked)
data(frogs)
umf <- formatMult(masspcru)
obsCovs(umf) <- scale(obsCovs(umf))
siteCovs(umf) <- rnorm(numSites(umf))
yearlySiteCovs(umf) <- data.frame(year = factor(rep(1:7,
                                    numSites(umf))))

##set up model with constant transition rates
fm <- colext(psiformula = ~ 1, gammaformula = ~ 1, epsilonformula = ~ 1,
             pformula = ~ JulianDate + I(JulianDate^2), data = umf,
             control = list(trace=1, maxit=1e4))

##model with with year-dependent transition rates
fm.yearly <- colext(psiformula = ~ 1, gammaformula = ~ year,
                    epsilonformula = ~ year,
                    pformula = ~ JulianDate + I(JulianDate^2),
                    data = umf)

##store in list and assign model names
Cand.mods <- list(fm, fm.yearly)
Modnames <- c("psi1(.)gam(.)eps(.)p(Date + Date2)",
              "psi1(.)gam(Year)eps(Year)p(Date + Date2)")

##compute model-averaged estimate of occupancy in the first year
modavg(cand.set = Cand.mods, modnames = Modnames, parm = "(Intercept)",
       parm.type = "psi")

##compute model-averaged estimate of Julian Day squared on detectability
modavg(cand.set = Cand.mods, modnames = Modnames,
       parm = "I(JulianDate^2)", parm.type = "detect")
}


##example of model-averaged estimate of area from distance model
##this is a bit longer
\dontrun{
data(linetran) #example modified from ?distsamp
     
ltUMF <- with(linetran, {
  unmarkedFrameDS(y = cbind(dc1, dc2, dc3, dc4),
                  siteCovs = data.frame(Length, area, habitat),
                  dist.breaks = c(0, 5, 10, 15, 20),
                  tlength = linetran$Length * 1000, survey = "line", unitsIn = "m")
})
     
## Half-normal detection function. Density output (log scale). No covariates.
fm1 <- distsamp(~ 1 ~ 1, ltUMF)
     
## Halfnormal. Covariates affecting both density and detection.
fm2 <- distsamp(~ area + habitat ~ area + habitat, ltUMF)

## Hazard function. Covariates affecting both density and detection.
fm3 <- distsamp(~ habitat ~ area + habitat, ltUMF, keyfun="hazard")

##assemble model list
Cands <- list(fm1, fm2, fm3)
Modnames <- paste("mod", 1:length(Cands), sep = "")

##model-average estimate of area on abundance
modavg(cand.set = Cands, modnames = Modnames, parm = "area", parm.type = "lambda")
detach(package:unmarked)
}
}
\keyword{models}

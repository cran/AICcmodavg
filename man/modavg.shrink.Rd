\name{modavg.shrink}
\Rdversion{1.1}
\alias{modavg.shrink}
\alias{modavg.shrink.glm}
\alias{modavg.shrink.gls}
\alias{modavg.shrink.lme}
\alias{modavg.shrink.mer}
\alias{modavg.shrink.mult}
\alias{modavg.shrink.polr}
\alias{print.modavg.shrink}

\title{
Compute Model-averaged Parameter Estimate with Shrinkage (Multimodel Inference)
}
\description{
This function computes an alternative version of model-averaging
parameter estimates that consists in shrinking estimates toward 0 as in
Burnham and Anderson (2002) and Anderson (2008).  Specifically, models
without the parameter of interest have an estimate and variance of
0. 'modavg.shrink' also returns unconditional standard error and
unconditional confidence intervals as described in Buckland et
al. (1997) and Burnham and Anderson (2002).
}
\usage{
modavg.shrink(cand.set, parm, modnames, c.hat = 1, gamdisp = NULL,
conf.level = 0.95, second.ord = TRUE, nobs = NULL,
uncond.se = "revised")

modavg.shrink.glm(cand.set, parm, modnames, c.hat = 1, gamdisp = NULL,
conf.level = 0.95, second.ord = TRUE, nobs = NULL,
uncond.se = "revised")

modavg.shrink.gls(cand.set, parm, modnames, conf.level = 0.95,
second.ord = TRUE, nobs = NULL, uncond.se = "revised")

modavg.shrink.lme(cand.set, parm, modnames, conf.level = 0.95,
second.ord = TRUE, nobs = NULL, uncond.se = "revised")

modavg.shrink.mer(cand.set, parm, modnames, conf.level = 0.95,
second.ord = TRUE, nobs = NULL, uncond.se = "revised")

modavg.shrink.mult(cand.set, parm, modnames, c.hat = 1,
conf.level = 0.95, second.ord = TRUE, nobs = NULL,
uncond.se = "revised") 

modavg.shrink.polr(cand.set, parm, modnames, conf.level = 0.95,
second.ord = TRUE, nobs = NULL, uncond.se = "revised")

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
each model in the model selection table.
}
  \item{c.hat}{
value of overdispersion parameter (i.e., variance inflation factor) such
as that obtained from 'c_hat'.  Note that values of c.hat different from
1 are only appropriate either for binomial GLM's with trials > 1 (i.e.,
success/trial or cbind(success, failure) syntax) or for Poisson GLM's.
If c.hat > 1, 'aictab' will return the quasi-likelihood analogue of the
information criteria requested and multiply the variance-covariance
matrix of the estimates by this value (i.e., SE's are multiplied by
sqrt(c.hat)).  This option is not supported for generalized linear
mixed models of the 'mer' class.
}
  \item{gamdisp}{
if gamma GLM is used, the dispersion parameter should be specified here
to apply the same value to each model.
}
  \item{conf.level}{
the confidence level requested for the computation of unconditional
confidence intervals.
}
\item{second.ord}{
  logical.  If TRUE, the function returns the second-order Akaike
  information criterion (i.e., AICc).
}
\item{nobs}{
  this argument allows to specify a numeric value other than total
  sample size to compute the AICc (i.e., 'nobs' defaults to total number
  of observations).  This is relevant only for linear mixed models where
  sample size is not straightforward. In such cases, one might use total
  number of observations or number of independent clusters as the value
  of 'nobs'.
}
\item{uncond.se}{
  either, "old", or "revised", specifying the equation used to compute
  the unconditional standard error of a model-averaged estimate.  With
  uncond.se = "old", computations are based on equation 4.9 of Burnham
  and Anderson (2002), which was the former way to compute unconditional
  standard errors.  With uncond.se = "revised", equation 6.12 of Burnham
  and Anderson (2002) is used.  Anderson (2008, p. 111) recommends use
  of the revised version for the computation of unconditional standard
  errors and it is now the default.  Note that versions of package
  AICcmodavg < 1.04 used the old method to compute unconditional
  standard errors.  
}
}

\details{
The parameter must be specified as it appears in the model output.  The
shrinkage version of model averaging is only appropriate for cases where
each parameter is given an equal weighting in the model (i.e., each
parameter must appear the same number of times in the models) and has
the same interpretation across all models.  As a result, models with
interaction terms or polynomial terms are not supported by
'modavg.shrink'.

'modavg.shrink' is a function that calls 'modavg.glm.shrink',
'modavg.gls.shrink', 'modavg.lme.shrink', 'modavg.mer.shrink',
'modavg.mult.shrink', or 'modavg.polr.shrink' depending on the class of
the object.  The current function is implemented for a list containing
objects of 'lm', 'glm', 'gls', 'lme', 'mer', 'multinom', and 'polr'
classes. 
}
\value{
'modavg.shrink', 'modavg.shrink.glm', 'modavg.shrink.gls',
'modavg.shrink.lme', 'modavg.shrink.mer', 'modavg.shrink.mult', and
'modavg.shrink.polr' create an object of class 'modavg' with the
following components:  

\item{Parameter}{the parameter for which a model-averaged estimate with
  shrinkage was obtained} 
\item{Mod.avg.table}{the model selection table based on models including
  the parameter of interest} 
\item{Mod.avg.beta}{the model-averaged estimate based on all models}
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

Mazerolle, M. J. (2006) Improving data analysis in herpetology: using
Akaike's Information Criterion (AIC) to assess the strength of
biological hypotheses. \emph{Amphibia-Reptilia} \bold{27}, 169--180. 

}
\author{
Marc J. Mazerolle
}

\seealso{
\code{\link{AICc}}, \code{\link{aictab}}, \code{\link{c_hat}},
\code{\link{importance}}, \code{\link{confset}}, \code{\link{evidence}},
\code{\link{modavg}}, \code{\link{modavgpred}}
}
\examples{
##cement example in Burnham and Anderson 2002
data(cement)
##setup same model set as in Table 3.2          
Cand.models <- list()
Cand.models[[1]] <- lm(y ~ x1 + x2, data = cement)
Cand.models[[2]] <- lm(y ~ x1 + x2 + x4, data = cement)          
Cand.models[[3]] <- lm(y ~ x1 + x2 + x3, data = cement)
Cand.models[[4]] <- lm(y ~ x1 + x4, data = cement)
Cand.models[[5]] <- lm(y ~ x1 + x3 + x4, data = cement)
Cand.models[[6]] <- lm(y ~ x2 + x3 + x4, data = cement)
Cand.models[[7]] <- lm(y ~ x1 + x2 + x3 + x4, data = cement)
Cand.models[[8]] <- lm(y ~ x3 + x4, data = cement)
Cand.models[[9]] <- lm(y ~ x2 + x3, data = cement)
Cand.models[[10]] <- lm(y ~ x4, data = cement)
Cand.models[[11]] <- lm(y ~ x2, data = cement)
Cand.models[[12]] <- lm(y ~ x2 + x4, data = cement)
Cand.models[[13]] <- lm(y ~ x1, data = cement)
Cand.models[[14]] <- lm(y ~ x1 + x3, data = cement)
Cand.models[[15]] <- lm(y ~ x3, data = cement)

Modnames <- paste("mod", 1:15, sep="")

##AICc          
aictab(cand.set = Cand.models, modnames = Modnames)

##compute model-averaged estimate with shrinkage - each parameter
##appears 8 times in the models 
modavg.shrink(cand.set = Cand.models, modnames = Modnames, parm = "x1")

##compare against classic model-averaging
modavg(cand.set = Cand.models, modnames = Modnames, parm = "x1")
##note that model-averaged estimate with shrinkage is closer to 0 than
##with the classic version

##remove a few models from the set and run again
Cand.unbalanced <- Cand.models[-c(3, 14, 15)]

##set up model names
Modnames <- paste("mod", 1:length(Cand.unbalanced), sep="")

##issues an error because some parameters appear more often than others
\dontrun{modavg.shrink(cand.set = Cand.unbalanced, modnames = Modnames,
parm = "x1")}


##example on Orthodont data set in nlme
require(nlme)

##set up candidate model list
##age and sex parameters appear in the same number of models
##same number of models with and without these parameters
Cand.models <- list()
Cand.models[[1]] <- lme(distance ~ age, data = Orthodont, method = "ML") 
##random is ~ age | Subject as it is a grouped data frame
Cand.models[[2]] <- lme(distance ~ age + Sex, data = Orthodont, random =
~ 1, method = "ML")
Cand.models[[3]] <- lme(distance ~ 1, data = Orthodont, random = ~ 1, 
method = "ML") 
Cand.models[[4]] <- lme(distance ~ Sex, data = Orthodont, random = ~ 1,
method = "ML")  

##create a vector of model names
Modnames <- NULL
for (i in 1:length(Cand.models)) {
Modnames[i] <- paste("mod", i, sep = "")
}

##compute importance values for age
imp.age <- importance(cand.set = Cand.models, parm = "age",
modnames = Modnames, second.ord = TRUE, nobs = NULL)

##compute shrinkage version of model averaging on age
mod.avg.age.shrink <- modavg.shrink(cand.set = Cand.models,
parm = "age", modnames = Modnames, second.ord = TRUE, nobs = NULL)

##compute classic version of model averaging on age
mod.avg.age.classic <- modavg(cand.set = Cand.models, parm = "age",
modnames = Modnames, second.ord = TRUE, nobs = NULL)

##correspondence between shrinkage version and classic version of
##model averaging 
mod.avg.age.shrink$Mod.avg.beta/imp.age$w.plus
mod.avg.age.classic$Mod.avg.beta

}
\keyword{models}
\name{modavgCustom}
\Rdversion{1.1}
\alias{modavgCustom}
\alias{print.modavgCustom}

\title{
Compute Model-averaged Parameter Estimate from User-supplied Input Based
on (Q)AIC(c)
}
\description{
This function model-averages the estimate of a parameter of interest
among a set of candidate models, and computes the unconditional standard
error and unconditional confidence intervals as described in Buckland et
al. (1997) and Burnham and Anderson (2002). 
}
\usage{
modavgCustom(logL, K, modnames = NULL, estimate, se, second.ord = TRUE,
             nobs = NULL, uncond.se = "revised", conf.level = 0.95,
             c.hat = 1, useBIC = FALSE) 
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{logL}{
    a vector of log-likelihood values for the models in the candidate
    model set.
  }
  \item{K}{
    a vector containing the number of estimated parameters for each
    model in the candidate model set.
  }
  \item{modnames}{
a character vector of model names to facilitate the identification of
each model in the model selection table.  If \code{NULL}, the function
uses the names in the cand.set list of candidate models. If no names
appear in the list, generic names (e.g., \code{Mod1}, \code{Mod2}) are
supplied in the table in the same order as in the list of candidate
models.
}
\item{estimate}{
  a vector of estimates for each of the models in the candidate model
  set.  Estimates can be either beta estimates for a parameter of
  interest or a single prediction from each model.
}
\item{se}{
  a vector of standard errors for each of the estimates appearing in the
  \code{estimate} vector.
}  
\item{second.ord}{
logical.  If \code{TRUE}, the function returns the second-order
Akaike information criterion (i.e., AICc).  This argument is ignored if
\code{useBIC = TRUE}.
}
\item{nobs}{
  the sample size required to compute the AICc, QAICc, BIC, or QBIC.
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
  default.
}
\item{conf.level}{
the confidence level (\eqn{1 - \alpha}) requested for the computation of
unconditional confidence intervals.
}
  \item{c.hat}{
value of overdispersion parameter (i.e., variance inflation factor) such
as that obtained from \code{c_hat}.  Note that values of \code{c.hat}
different from 1 are only appropriate for binomial GLM's with trials > 1
(i.e., success/trial or cbind(success, failure) syntax), with Poisson
GLM's, single-season and dynamic occupancy models (MacKenzie et
al. 2002, 2003), \emph{N}-mixture models (Royle 2004, Dail and Madsen
2011), or capture-mark-recapture models (e.g., Lebreton et al. 1992).
If \code{c.hat} > 1, \code{modavgCustom} will return the
quasi-likelihood analogue of the information criteria requested and
multiply the variance-covariance matrix of the estimates by this value
(i.e., SE's are multiplied by \code{sqrt(c.hat)}).
}
\item{useBIC}{
logical.  If \code{TRUE}, the function returns the Bayesian information
criterion (BIC) when \code{c.hat = 1} or the quasi-likelihood BIC (QBIC)
when \code{c.hat > 1}.
}
}

\details{
\code{modavgCustom} computes a model-averaged estimate from the vector
of parameter estimates specified in \code{estimate}.  Estimates and
their associated standard errors must be specified in the same order as
the log-likelihood, number of estimated parameters, and model names.
Estimates provided may be for a parameter of interest (i.e., beta
estimates) or predictions from each model.  This function is most useful
when model input is imported into R from other software (e.g., Program
MARK, PRESENCE) or for model classes that are not yet supported by the
other model averaging functions such as \code{modavg} or
\code{modavgPred}.
}

\value{
\code{modavgCustom} creates an object of class \code{modavgCustom} with
the following components: 

\item{Mod.avg.table}{the model selection table}
\item{Mod.avg.est}{the model-averaged estimate}
\item{Uncond.SE}{the unconditional standard error for the model-averaged
  estimate}
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

Dail, D., Madsen, L. (2011) Models for estimating abundance from 
repeated counts of an open population. \emph{Biometrics} \bold{67},
577--587.

Lebreton, J.-D., Burnham, K. P., Clobert, J., Anderson, D. R. (1992)
Modeling survival and testing biological hypotheses using marked
animals: a unified approach with case-studies. \emph{Ecological
  Monographs} \bold{62}, 67--118. 

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
\code{\link{AICcCustom}}, \code{\link{aictabCustom}},
\code{\link{bictabCustom}}, \code{\link{modavg}},
\code{\link{modavgIC}}, \code{\link{modavgShrink}},
\code{\link{modavgPred}}
}
\examples{
\dontrun{
##model averaging parameter estimate (natural average)
##vector with model LL's
LL <- c(-38.8876, -35.1783, -64.8970)

##vector with number of parameters
Ks <- c(7, 9, 4)

##create a vector of names to trace back models in set
Modnames <- c("Cm1", "Cm2", "Cm3")

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


##model-averaging with shrinkage based on AICc
##set up candidate models
data(min.trap)
Cand.mod <- list( )
##global model          
Cand.mod[[1]] <- glm(Num_anura ~ Type + log.Perimeter,
                     family = poisson, offset = log(Effort),
                     data = min.trap) 
Cand.mod[[2]] <- glm(Num_anura ~ Type + Num_ranatra, family = poisson,
                     offset = log(Effort), data = min.trap) 
Cand.mod[[3]] <- glm(Num_anura ~ log.Perimeter + Num_ranatra,
                     family = poisson, offset = log(Effort), data = min.trap)
Model.names <- c("Type + log.Perimeter", "Type + Num_ranatra",
                 "log.Perimeter + Num_ranatra")
##model-averaged estimate with shrinkage (glm model type is already supported)
modavgShrink(cand.set = Cand.mod, modnames = Model.names,
             parm = "log.Perimeter")

##equivalent manual version of model-averaging with shrinkage
##this is especially useful when model classes are not supported
##extract vector of LL
LLs <- sapply(Cand.mod, FUN = function(i) logLik(i)[1])
##extract vector of K
Ks <- sapply(Cand.mod, FUN = function(i) attr(logLik(i), "df"))
##extract betas
betas <- sapply(Cand.mod, FUN = function(i) coef(i)["log.Perimeter"])
##second model does not include log.Perimeter
betas[2] <- 0
##extract SE's
ses <- sapply(Cand.mod, FUN = function(i) sqrt(diag(vcov(i))["log.Perimeter"]))
ses[2] <- 0
##model-averaging with shrinkage based on AICc
modavgCustom(logL = LLs, K = Ks, modnames = Model.names,
             nobs = nrow(min.trap), estimate = betas, se = ses)
##model-averaging with shrinkage based on BIC
modavgCustom(logL = LLs, K = Ks, modnames = Model.names,
             nobs = nrow(min.trap), estimate = betas, se = ses,
             useBIC = TRUE)
}
}
\keyword{models}

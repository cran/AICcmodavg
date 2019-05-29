\name{useBIC}
\Rdversion{1.1}
\alias{useBIC}
\alias{useBIC.default}
\alias{useBIC.aov}
\alias{useBIC.betareg}
\alias{useBIC.clm}
\alias{useBIC.clmm}
\alias{useBIC.coxme}
\alias{useBIC.coxph}
\alias{useBIC.fitdist}
\alias{useBIC.fitdistr}
\alias{useBIC.glm}
\alias{useBIC.glmmTMB}
\alias{useBIC.gls}
\alias{useBIC.gnls}
\alias{useBIC.hurdle}
\alias{useBIC.lavaan}
\alias{useBIC.lm}
\alias{useBIC.lme}
\alias{useBIC.lmekin}
\alias{useBIC.maxlikeFit}
\alias{useBIC.mer}
\alias{useBIC.merMod}
\alias{useBIC.lmerModLmerTest}
\alias{useBIC.multinom}
\alias{useBIC.nlme}
\alias{useBIC.nls}
\alias{useBIC.polr}
\alias{useBIC.rlm}
\alias{useBIC.survreg}
\alias{useBIC.unmarkedFit}
\alias{useBIC.vglm}
\alias{useBIC.zeroinfl}
\title{
  Computing BIC or QBIC
}
\description{
  Functions to compute the Bayesian information criterion (BIC) or a
  quasi-likelihood analogue (QBIC).
}
\usage{
useBIC(mod, return.K = FALSE, nobs = NULL, \dots) 

\method{useBIC}{aov}(mod, return.K = FALSE, nobs = NULL, \dots) 

\method{useBIC}{betareg}(mod, return.K = FALSE, nobs = NULL, \dots) 

\method{useBIC}{clm}(mod, return.K = FALSE, nobs = NULL, \dots) 

\method{useBIC}{clmm}(mod, return.K = FALSE, nobs = NULL, \dots)  

\method{useBIC}{coxme}(mod, return.K = FALSE, nobs = NULL, \dots) 

\method{useBIC}{coxph}(mod, return.K = FALSE, nobs = NULL, \dots) 

\method{useBIC}{fitdist}(mod, return.K = FALSE, nobs = NULL, \dots)

\method{useBIC}{fitdistr}(mod, return.K = FALSE, nobs = NULL, \dots)

\method{useBIC}{glm}(mod, return.K = FALSE, nobs = NULL, c.hat = 1,
\dots)

\method{useBIC}{glmmTMB}(mod, return.K = FALSE, nobs = NULL, c.hat = 1,
\dots)

\method{useBIC}{gls}(mod, return.K = FALSE, nobs = NULL, \dots)

\method{useBIC}{gnls}(mod, return.K = FALSE, nobs = NULL, \dots)

\method{useBIC}{hurdle}(mod, return.K = FALSE, nobs = NULL, \dots)

\method{useBIC}{lavaan}(mod, return.K = FALSE, nobs = NULL, \dots)

\method{useBIC}{lm}(mod, return.K = FALSE, nobs = NULL, \dots)

\method{useBIC}{lme}(mod, return.K = FALSE, nobs = NULL, \dots)

\method{useBIC}{lmekin}(mod, return.K = FALSE, nobs = NULL, \dots) 

\method{useBIC}{maxlikeFit}(mod, return.K = FALSE, nobs = NULL, c.hat =
1, \dots)

\method{useBIC}{mer}(mod, return.K = FALSE, nobs = NULL, \dots)

\method{useBIC}{merMod}(mod, return.K = FALSE, nobs = NULL, \dots)  

\method{useBIC}{lmerModLmerTest}(mod, return.K = FALSE, nobs = NULL, \dots)  

\method{useBIC}{multinom}(mod, return.K = FALSE, nobs = NULL, c.hat = 1,
\dots)

\method{useBIC}{nlme}(mod, return.K = FALSE, nobs = NULL, \dots)

\method{useBIC}{nls}(mod, return.K = FALSE, nobs = NULL, \dots)

\method{useBIC}{polr}(mod, return.K = FALSE, nobs = NULL, \dots)

\method{useBIC}{rlm}(mod, return.K = FALSE, nobs = NULL, \dots)

\method{useBIC}{survreg}(mod, return.K = FALSE, nobs = NULL, \dots)

\method{useBIC}{unmarkedFit}(mod, return.K = FALSE, nobs = NULL, c.hat =
1, \dots)

\method{useBIC}{vglm}(mod, return.K = FALSE, nobs = NULL, c.hat = 1,
\dots)

\method{useBIC}{zeroinfl}(mod, return.K = FALSE, nobs = NULL, \dots)
}

\arguments{
  \item{mod}{
    an object of class \code{aov}, \code{betareg}, \code{clm},
    \code{clmm}, \code{clogit}, \code{coxme}, \code{coxph},
    \code{fitdist}, \code{fitdistr}, \code{glm}, \code{glmmTMB},
    \code{gls}, \code{gnls}, \code{hurdle}, \code{lavaan}, \code{lm},
    \code{lme}, \code{lmekin}, \code{maxlikeFit}, \code{mer},
    \code{merMod}, \code{lmerModLmerTest}, \code{multinom}, \code{nlme},
    \code{nls}, \code{polr}, \code{rlm}, \code{survreg}, \code{vglm},
    \code{zeroinfl}, and various \code{unmarkedFit} classes containing
    the output of a model. 
}
  
  \item{return.K}{
    logical.  If \code{FALSE}, the function returns the information
    criterion specified.  If \code{TRUE}, the function returns K (number
    of estimated parameters) for a given model.
}
\item{nobs}{
  this argument allows to specify a numeric value other than total
  sample size to compute the BIC (i.e., \code{nobs} defaults to total
  number of observations).  This is relevant only for mixed models or
  various models of \code{unmarkedFit} classes where sample size is not
  straightforward.  In such cases, one might use total number of
  observations or number of independent clusters (e.g., sites) as the 
  value of \code{nobs}.
}
 \item{c.hat}{
  value of overdispersion parameter (i.e., variance inflation factor)
  such as that obtained from \code{c_hat}.  Note that values of c.hat
  different from 1 are only appropriate for binomial GLM's with trials
  > 1 (i.e., success/trial or cbind(success, failure) syntax), with
  Poisson GLM's, single-season occupancy models (MacKenzie et al. 2002),
  dynamic occupancy models (MacKenzie et al. 2003), or \emph{N}-mixture
  models (Royle 2004, Dail and Madsen 2011).  If \code{c.hat} > 1,
  \code{useBIC} will return the quasi-likelihood analogue of the
  information criteria requested and multiply the variance-covariance
  matrix of the estimates by this value (i.e., SE's are multiplied by
  \code{sqrt(c.hat)}).  This option is not supported for generalized
  linear mixed models of the \code{mer} or \code{merMod} classes.
}
\item{\dots}{
  additional arguments passed to the function.
  }
}
\details{
  \code{useBIC} computes the Bayesian information criterion (BIC,
    Schwarz 1978): \deqn{BIC = -2 * log-likelihood + K * log(n),} where
    the log-likelihood is the maximum log-likelihood of the model, \emph{K}
    corresponds to the number of estimated parameters, and \emph{n}
    corresponds to the sample size of the data set.

    In the presence of overdispersion, a quasi-likelihood analogue of
    the BIC (QBIC) will be computed, as \deqn{QBIC = \frac{-2 *
    log-likelihood}{c-hat} + K * log(n),} where \emph{c-hat} is the
    overdispersion parameter specified by the user with the argument
    \code{c.hat}.  Note that BIC or QBIC values are meaningful to select
    among \code{gls} or \code{lme} models fit by maximum likelihood.
    BIC or QBIC based on REML are valid to select among different models
    that only differ in their random effects (Pinheiro and Bates 2000). 
}
\value{
  \code{useBIC} returns the BIC or the number of estimated parameters,
    depending on the values of the arguments.
}

\note{
  The actual (Q)BIC values are not really interesting in themselves, as
  they depend directly on the data, parameters estimated, and likelihood
  function.  Furthermore, a single value does not tell much about model
  fit.  Information criteria become relevant when compared to one
  another for a given data set and set of candidate models. 
%
%  Note that for robust regression models of class \code{rlm}, the
%  AIC is computed based on Tharmaratnam and Claeskens (2013).  The
%  second-order AIC is not yet implemented for the \code{rlm} class.
}

\references{
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

  MacKenzie, D. I., Nichols, J. D., Hines, J. E., Knutson, M. G.,
  Franklin, A. B. (2003) Estimating site occupancy, colonization, and
  local extinction when a species is detected imperfectly. \emph{Ecology}
  \bold{84}, 2200--2207.

  Pinheiro, J. C., Bates, D. M. (2000) \emph{Mixed-effect models in S and
  S-PLUS}. Springer Verlag: New York.  

  Royle, J. A. (2004) \emph{N}-mixture models for estimating population
  size from spatially replicated counts. \emph{Biometrics} \bold{60}, 
  108--115.

  Schwarz, G. (1978) Estimating the dimension of a model. \emph{Annals of
  Statistics} \bold{6}, 461--464.
%  Tharmaratnam, K., Claeskens, G. (2013) A comparison of robust
%  versions of the AIC based on M-, S- and MM-estimators. \emph{Statistics}
%  \bold{47}, 216--235. 

}
\author{
  Marc J. Mazerolle
}

\seealso{
\code{\link{AICc}}, \code{\link{bictab}},
\code{\link{bictabCustom}}, \code{\link{useBICCustom}}
}


\examples{
##cement data from Burnham and Anderson (2002, p. 101)
data(cement)
##run multiple regression - the global model in Table 3.2
glob.mod <- lm(y ~ x1 + x2 + x3 + x4, data = cement)

##compute BIC with full likelihood
useBIC(glob.mod, return.K = FALSE)

##compute BIC for mixed model on Orthodont data set in Pinheiro and
##Bates (2000)
\dontrun{
require(nlme)
m1 <- lme(distance ~ age, random = ~1 | Subject, data = Orthodont,
          method= "ML")
useBIC(m1, return.K = FALSE)
}

}
\keyword{models}
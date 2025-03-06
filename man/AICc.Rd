\name{AICc}
\Rdversion{1.1}
\alias{AICc}
\alias{AICc.default}
\alias{AICc.aov}
\alias{AICc.betareg}
\alias{AICc.clm}
\alias{AICc.clmm}
\alias{AICc.coxme}
\alias{AICc.coxph}
\alias{AICc.fitdist}
\alias{AICc.fitdistr}
\alias{AICc.glm}
\alias{AICc.glmmTMB}
\alias{AICc.gls}
\alias{AICc.gnls}
\alias{AICc.hurdle}
\alias{AICc.lavaan}
\alias{AICc.lm}
\alias{AICc.lme}
\alias{AICc.lmekin}
\alias{AICc.maxlikeFit}
\alias{AICc.mer}
\alias{AICc.merMod}
\alias{AICc.lmerModLmerTest}
\alias{AICc.multinom}
\alias{AICc.negbin}
\alias{AICc.nlme}
\alias{AICc.nls}
\alias{AICc.polr}
\alias{AICc.rlm}
\alias{AICc.survreg}
\alias{AICc.unmarkedFit}
\alias{AICc.vglm}
\alias{AICc.zeroinfl}
\title{
  Computing AIC, AICc, QAIC, and QAICc
}
\description{
  Functions to compute Akaike's information criterion (AIC), the
  second-order AIC (AICc), as well as their quasi-likelihood
  counterparts (QAIC, QAICc).
}
\usage{
AICc(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, \dots) 

\method{AICc}{aov}(mod, return.K = FALSE, second.ord = TRUE,
     nobs = NULL, \dots)

\method{AICc}{betareg}(mod, return.K = FALSE, second.ord = TRUE,
     nobs = NULL, \dots)

\method{AICc}{clm}(mod, return.K = FALSE, second.ord = TRUE,
     nobs = NULL, \dots) 

\method{AICc}{clmm}(mod, return.K = FALSE, second.ord = TRUE,
     nobs = NULL, \dots)  

\method{AICc}{coxme}(mod, return.K = FALSE, second.ord = TRUE,
     nobs = NULL, \dots) 

\method{AICc}{coxph}(mod, return.K = FALSE, second.ord = TRUE,
     nobs = NULL, \dots) 

\method{AICc}{fitdist}(mod, return.K = FALSE, second.ord = TRUE,
     nobs = NULL, \dots)

\method{AICc}{fitdistr}(mod, return.K = FALSE, second.ord = TRUE,
     nobs = NULL, \dots)

\method{AICc}{glm}(mod, return.K = FALSE, second.ord = TRUE,
     nobs = NULL, c.hat = 1, \dots)

\method{AICc}{glmmTMB}(mod, return.K = FALSE, second.ord = TRUE,
     nobs = NULL, c.hat = 1, \dots)

\method{AICc}{gls}(mod, return.K = FALSE, second.ord = TRUE,
     nobs = NULL, \dots)

\method{AICc}{gnls}(mod, return.K = FALSE, second.ord = TRUE,
     nobs = NULL, \dots)

\method{AICc}{hurdle}(mod, return.K = FALSE, second.ord = TRUE,
     nobs = NULL, \dots)

\method{AICc}{lavaan}(mod, return.K = FALSE, second.ord = TRUE,
     nobs = NULL, \dots)

\method{AICc}{lm}(mod, return.K = FALSE, second.ord = TRUE,
     nobs = NULL, \dots)

\method{AICc}{lme}(mod, return.K = FALSE, second.ord = TRUE,
     nobs = NULL, \dots)

\method{AICc}{lmekin}(mod, return.K = FALSE, second.ord = TRUE,
     nobs = NULL, \dots) 

\method{AICc}{maxlikeFit}(mod, return.K = FALSE, second.ord = TRUE,
     nobs = NULL, c.hat = 1, \dots)

\method{AICc}{mer}(mod, return.K = FALSE, second.ord = TRUE,
     nobs = NULL, \dots)

\method{AICc}{merMod}(mod, return.K = FALSE, second.ord = TRUE,
     nobs = NULL, \dots) 

\method{AICc}{lmerModLmerTest}(mod, return.K = FALSE, second.ord = TRUE,
     nobs = NULL, \dots) 

\method{AICc}{multinom}(mod, return.K = FALSE, second.ord = TRUE,
     nobs = NULL, c.hat = 1, \dots)

\method{AICc}{negbin}(mod, return.K = FALSE, second.ord = TRUE,
     nobs = NULL, \dots)

\method{AICc}{nlme}(mod, return.K = FALSE, second.ord = TRUE,
     nobs = NULL, \dots)

\method{AICc}{nls}(mod, return.K = FALSE, second.ord = TRUE,
     nobs = NULL, \dots)

\method{AICc}{polr}(mod, return.K = FALSE, second.ord = TRUE,
     nobs = NULL, \dots)

\method{AICc}{rlm}(mod, return.K = FALSE, second.ord = TRUE,
     nobs = NULL, \dots)

\method{AICc}{survreg}(mod, return.K = FALSE, second.ord = TRUE,
     nobs = NULL, \dots)

\method{AICc}{unmarkedFit}(mod, return.K = FALSE, second.ord = TRUE,
     nobs = NULL, c.hat = 1, \dots)

\method{AICc}{vglm}(mod, return.K = FALSE, second.ord = TRUE,
     nobs = NULL, c.hat = 1, \dots)

\method{AICc}{zeroinfl}(mod, return.K = FALSE, second.ord = TRUE,
     nobs = NULL, \dots)
}

\arguments{
  \item{mod}{
    an object of class \code{aov}, \code{betareg}, \code{clm},
    \code{clmm}, \code{clogit}, \code{coxme}, \code{coxph},
    \code{fitdist}, \code{fitdistr}, \code{glm}, \code{glmmTMB},
    \code{gls}, \code{gnls}, \code{hurdle}, \code{lavaan}, \code{lm},
    \code{lme}, \code{lmekin}, \code{maxlikeFit}, \code{mer},
    \code{merMod}, \code{lmerModLmerTest}, \code{multinom},
    \code{negbin}, \code{nlme}, \code{nls}, \code{polr}, \code{rlm},
    \code{survreg}, \code{vglm}, \code{zeroinfl}, and various
    \code{unmarkedFit} classes containing the output of a model.
  }
  
  \item{return.K}{
    logical.  If \code{FALSE}, the function returns the information
    criterion specified.  If \code{TRUE}, the function returns K (number
    of estimated parameters) for a given model.
}
\item{second.ord}{ 
  logical.  If \code{TRUE}, the function returns the second-order Akaike
  information criterion (i.e., AICc). 
}
\item{nobs}{
  this argument allows to specify a numeric value other than total
  sample size to compute the AICc (i.e., \code{nobs} defaults to total number
  of observations).  This is relevant only for mixed models or various
  models of \code{unmarkedFit} classes where sample size is not
  straightforward.  In such cases, one might use total number of
  observations or number of independent clusters (e.g., sites) as the 
  value of \code{nobs}.
}
 \item{c.hat}{
  value of overdispersion parameter (i.e., variance inflation factor) such
as that obtained from \code{c_hat}.  Note that values of c.hat different
from 1 are only appropriate for binomial GLM's with trials > 1 (i.e.,
success/trial or cbind(success, failure) syntax), with Poisson GLM's, 
single-season occupancy models (MacKenzie et al. 2002), dynamic
occupancy models (MacKenzie et al. 2003), or \emph{N}-mixture models
(Royle 2004, Dail and Madsen 2011).  If \code{c.hat} > 1,
\code{AICc} will return the quasi-likelihood analogue of the
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
  \code{AICc} computes one of the following four information criteria:
  
  Akaike's information criterion (AIC, Akaike 1973), \deqn{-2 *
  log-likelihood + 2 * K,} where the log-likelihood is the maximum
  log-likelihood of the model and \emph{K} corresponds to the number of
  estimated parameters.
  
  Second-order or small sample AIC (AICc, Sugiura   1978, Hurvich and
  Tsai 1989, 1991), \deqn{-2 * log-likelihood + 2 * K  * (n/(n - K - 1)),}
  where \emph{n} is the sample size of the data set.
  
  Quasi-likelihood AIC (QAIC, Burnham and Anderson 2002), \deqn{QAIC =
  \frac{-2 * log-likelihood}{c-hat} + 2 * K,} where \emph{c-hat} is the
  overdispersion parameter specified by the user with the argument
  \code{c.hat}.

  Quasi-likelihood AICc (QAICc, Burnham and Anderson 2002), \deqn{QAIC = 
  \frac{-2 * log-likelihood}{c-hat} + 2 * K * (n/(n - K - 1))}. 

  Note that AIC and AICc values are meaningful to select among
  \code{gls} or \code{lme} models fit by maximum likelihood. AIC and
  AICc based on REML are valid to select among different models that
  only differ in their random effects (Pinheiro and Bates 2000).
}
\value{
  \code{AICc} returns the AIC, AICc, QAIC, or QAICc, or the number of
  estimated parameters, depending on the values of the arguments.
}

\note{
  The actual (Q)AIC(c) values are not really interesting in themselves,
  as they depend directly on the data, parameters estimated, and
  likelihood function.  Furthermore, a single value does not tell much
  about model fit.  Information criteria become relevant when compared
  to one another for a given data set and set of candidate models.
%
%  Note that for robust regression models of class \code{rlm}, the
%  AIC is computed based on Tharmaratnam and Claeskens (2013).  The
%  second-order AIC is not yet implemented for the \code{rlm} class.
}

\references{
  Akaike, H. (1973) Information theory as an extension of the maximum
  likelihood principle. In:  \emph{Second International Symposium on
  Information Theory}, pp. 267--281. Petrov, B.N., Csaki, F., Eds,
  Akademiai Kiado, Budapest.

  Anderson, D. R. (2008) \emph{Model-based Inference in the Life
  Sciences: a primer on evidence}. Springer: New York.

  Burnham, K. P., Anderson, D. R. (2002) \emph{Model Selection and
  Multimodel Inference: a practical information-theoretic
  approach}. Second edition. Springer: New York.

  Burnham, K. P., Anderson, D. R. (2004) Multimodel inference:
  understanding AIC and BIC in model selection. \emph{Sociological
  Methods and Research} \bold{33}, 261--304.

  Dail, D., Madsen, L. (2011) Models for estimating abundance from
  repeated counts of an open population. \emph{Biometrics} \bold{67},
  577--587. 
  
  Hurvich, C. M., Tsai, C.-L. (1989) Regression and time series model
  selection in small samples. \emph{Biometrika} \bold{76}, 297--307. 

  Hurvich, C. M., Tsai, C.-L. (1991) Bias of the corrected AIC
  criterion for underfitted regression and time series
  models. \emph{Biometrika} \bold{78}, 499--509. 

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

  Sugiura, N. (1978) Further analysis of the data by Akaike's
  information criterion and the finite corrections. \emph{Communications
    in Statistics:  Theory and Methods} \bold{A7}, 13--26.

%  Tharmaratnam, K., Claeskens, G. (2013) A comparison of robust
%  versions of the AIC based on M-, S- and MM-estimators. \emph{Statistics}
%  \bold{47}, 216--235. 

}
\author{
  Marc J. Mazerolle
}

\seealso{
\code{\link{AICcCustom}}, \code{\link{aictab}}, \code{\link{confset}},
\code{\link{importance}}, \code{\link{evidence}}, \code{\link{c_hat}},
\code{\link{modavg}}, \code{\link{modavgShrink}},
\code{\link{modavgPred}}, \code{\link{useBIC}}, 
}


\examples{
##cement data from Burnham and Anderson (2002, p. 101)
data(cement)
##run multiple regression - the global model in Table 3.2
glob.mod <- lm(y ~ x1 + x2 + x3 + x4, data = cement)

##compute AICc with full likelihood
AICc(glob.mod, return.K = FALSE)

##compute AIC with full likelihood 
AICc(glob.mod, return.K = FALSE, second.ord = FALSE)
##note that Burnham and Anderson (2002) did not use full likelihood
##in Table 3.2 and that the MLE estimate of the variance was
##rounded to 2 digits after decimal point  


##compute AICc for mixed model on Orthodont data set in Pinheiro and
##Bates (2000)
\dontrun{
require(nlme)
m1 <- lme(distance ~ age, random = ~1 | Subject, data = Orthodont,
          method= "ML")
AICc(m1, return.K = FALSE)
}

}
\keyword{models}

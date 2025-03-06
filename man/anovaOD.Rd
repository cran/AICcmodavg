\name{anovaOD}
\Rdversion{1.1}
\alias{anovaOD}
\alias{anovaOD.glm}
\alias{anovaOD.unmarkedFitOccu}
\alias{anovaOD.unmarkedFitColExt}
\alias{anovaOD.unmarkedFitOccuRN}
\alias{anovaOD.unmarkedFitPCount}
\alias{anovaOD.unmarkedFitPCO}
\alias{anovaOD.unmarkedFitDS}
\alias{anovaOD.unmarkedFitGDS}
\alias{anovaOD.unmarkedFitOccuFP}
\alias{anovaOD.unmarkedFitMPois}
\alias{anovaOD.unmarkedFitGMM}
\alias{anovaOD.unmarkedFitGPC}
\alias{anovaOD.unmarkedFitOccuMulti}
\alias{anovaOD.unmarkedFitOccuMS}
\alias{anovaOD.unmarkedFitOccuTTD}
\alias{anovaOD.unmarkedFitMMO}
\alias{anovaOD.unmarkedFitDSO}
\alias{anovaOD.unmarkedFitGOccu}
\alias{anovaOD.unmarkedFitOccuComm}
\alias{anovaOD.glmerMod}
\alias{anovaOD.maxlikeFit}
\alias{anovaOD.multinom}
\alias{anovaOD.vglm}
\alias{print.anovaOD}

\title{
  Likelihood-Ratio Test Corrected for Overdispersion
}
\description{
Compute likelihood-ratio test between a given model and a simpler model.
}
\usage{
anovaOD(mod.simple, mod.complex, c.hat = 1, 
        nobs = NULL, \dots)

\method{anovaOD}{glm}(mod.simple, mod.complex, c.hat = 1,  
        nobs = NULL, \dots)

\method{anovaOD}{unmarkedFitOccu}(mod.simple, mod.complex, c.hat = 1,  
        nobs = NULL, \dots)

\method{anovaOD}{unmarkedFitColExt}(mod.simple, mod.complex, c.hat = 1,  
        nobs = NULL, \dots)

\method{anovaOD}{unmarkedFitOccuRN}(mod.simple, mod.complex, c.hat = 1,  
        nobs = NULL, \dots)

\method{anovaOD}{unmarkedFitPCount}(mod.simple, mod.complex, c.hat = 1,  
        nobs = NULL, \dots)

\method{anovaOD}{unmarkedFitPCO}(mod.simple, mod.complex, c.hat = 1,  
        nobs = NULL, \dots)

\method{anovaOD}{unmarkedFitDS}(mod.simple, mod.complex, c.hat = 1,  
        nobs = NULL, \dots)

\method{anovaOD}{unmarkedFitGDS}(mod.simple, mod.complex, c.hat = 1,  
        nobs = NULL, \dots)

\method{anovaOD}{unmarkedFitOccuFP}(mod.simple, mod.complex, c.hat = 1,  
        nobs = NULL, \dots)

\method{anovaOD}{unmarkedFitMPois}(mod.simple, mod.complex, c.hat = 1,  
        nobs = NULL, \dots)

\method{anovaOD}{unmarkedFitGMM}(mod.simple, mod.complex, c.hat = 1,  
        nobs = NULL, \dots)

\method{anovaOD}{unmarkedFitGPC}(mod.simple, mod.complex, c.hat = 1,  
        nobs = NULL, \dots)

\method{anovaOD}{unmarkedFitOccuMS}(mod.simple, mod.complex, c.hat = 1,  
        nobs = NULL, \dots)

\method{anovaOD}{unmarkedFitOccuTTD}(mod.simple, mod.complex, c.hat = 1,  
        nobs = NULL, \dots)

\method{anovaOD}{unmarkedFitMMO}(mod.simple, mod.complex, c.hat = 1,  
        nobs = NULL, \dots)

\method{anovaOD}{unmarkedFitDSO}(mod.simple, mod.complex, c.hat = 1,  
        nobs = NULL, \dots)

\method{anovaOD}{unmarkedFitGOccu}(mod.simple, mod.complex, c.hat = 1,  
        nobs = NULL, \dots)

\method{anovaOD}{unmarkedFitOccuComm}(mod.simple, mod.complex, c.hat = 1,  
        nobs = NULL, \dots)

\method{anovaOD}{glmerMod}(mod.simple, mod.complex, c.hat = 1,  
        nobs = NULL, \dots)

\method{anovaOD}{maxlikeFit}(mod.simple, mod.complex, c.hat = 1,  
        nobs = NULL, \dots)

\method{anovaOD}{multinom}(mod.simple, mod.complex, c.hat = 1,  
        nobs = NULL, \dots)

\method{anovaOD}{vglm}(mod.simple, mod.complex, c.hat = 1,  
        nobs = NULL, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mod.simple}{
    an object of class \code{glm}, \code{glmmTMB},
    \code{maxlikeFit}, \code{mer}, \code{merMod}, \code{multinom},  
    \code{vglm}, and various \code{unmarkedFit} classes containing 
    the output of a model.  This model should be a simpler version of
    \code{mod.complex} resulting from a deletion of certain terms (i.e.,
    nested model).
  }
  
  \item{mod.complex}{
    an object of the same class as \code{mod.simple}.
  }

  \item{c.hat}{
    value of overdispersion parameter (i.e., variance inflation factor)
    such as that obtained from \code{c_hat}, \code{mb.gof.test}, or
    \code{Nmix.gof.test}.  Typically, this value should be computed for
    the most complex model and applied to simpler models.
  }

  \item{nobs}{
    the number of observations used in the analysis.  If \code{nobs =
    NULL}, the total number of rows are used as the sample size to
    compute the residual degrees of freedom as \eqn{nobs - K}, where
    \emph{K} is the number of estimated parameters.  This is relevant
    only for mixed models or various models of \code{unmarkedFit}
    classes where sample size is not straightforward.  In such cases,
    one might use total number of observations or number of independent
    clusters (e.g., sites) as the value of \code{nobs}.
  }
  \item{\dots}{
    additional arguments passed to the function.
  }

}  
\details{
  This function applies a correction for overdispersion on the
  likelihood-ratio test between a model and its simpler counterpart.
  The simpler model must be nested within the more complex model,
  typically as the result of deleting terms.  You should supply the
  \code{c.hat} value of the most complex of the two models you are
  comparing.

  When \eqn{1 < c.hat < 4}, the likelihood-ratio test is computed
  as:
  
  \deqn{LR = \frac{-2 * (LL.simple - LL.complex)}{(K.complex - K.simple)
      * c.hat}}
  
  where \emph{LL.simple} and \emph{LL.complex} are the log-likelihoods
  of the simple and complex models, respectively, and where
  \emph{K.complex} and \emph{K.simple} are the number of estimated
  parameters in each model.  The test statistic is approximately
  distributed as \eqn{F_{K.complex - K.simple, n - K.complex}}, where
  \emph{n} is the number of observations (i.e., \code{nobs}) used in the
  analysis (Venables and Ripley 2002).

  When \code{nobs = NULL}, the number of observations is based on the
  number of rows of the data frame used in the analysis.  For mixed
  models or various models of \code{unmarkedFit}, sample size is less
  straightforward, and \code{nobs} could be based on the total number of
  observations or on the number of independent clusters (e.g., sites),
  among other choices.

  When \code{c.hat = 1}, the likelihood-ratio test simplifies to:

  \deqn{LR = -2 * (LL.simple - LL.complex)}
  
  where in this case the test statistic is distributed as a
  \eqn{\chi^2_{K.complex - K.simple}} (McCullagh and Nelder 1989).

  The function supports different model types such as Poisson GLM's and
  GLMM's, single-season and dynamic occupancy models (MacKenzie et
  al. 2002, 2003),  and various \emph{N}-mixture models (Royle 2004,
  Dail and Madsen 2011).
}
\value{
\code{anovaOD} returns an object of class \code{anovaOD} as a list with
the following components:

\item{form.simple}{a character string of the parameters estimated in \code{mod.simple}.}
\item{form.complex}{a character string of the parameters estimated in \code{mod.complex}.}
\item{c.hat}{the \code{c.hat} estimate used to adjust the
  likelihood-ratio test.} 
\item{devMat}{a matrix storing as columns the number of parameters
  estimated (\code{K}), the log-likelihood of each model \code{logLik},
  the difference in estimated parameters between the two models
  (\code{Kdiff}), minus twice the difference in log-likelihoods between
  the models (\code{-2LL}), the test statistic, and the associated
  P-value. 
}
}

\references{
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

  Mazerolle, M. J. (2006) Improving data analysis in herpetology: using
  Akaike's Information Criterion (AIC) to assess the strength of
  biological hypotheses. \emph{Amphibia-Reptilia} \bold{27}, 169--180. 
  
  McCullagh, P., Nelder, J. A. (1989) \emph{Generalized Linear
    Models}. Second edition. Chapman and Hall: New York.

  Royle, J. A. (2004) \emph{N}-mixture models for estimating population
  size from spatially replicated counts. \emph{Biometrics} \bold{60}, 
  108--115.

  Venables, W. N., Ripley, B. D. (2002) \emph{Modern Applied Statistics
  with S}. Second edition. Springer-Verlag: New York.

}
\author{
Marc J. Mazerolle
}

\seealso{
  \code{\link{c_hat}}, \code{\link{mb.gof.test}},
  \code{\link{Nmix.gof.test}}, \code{\link{summaryOD}} 
}

\examples{
##anuran larvae example from Mazerolle (2006)
data(min.trap)
##assign "UPLAND" as the reference level as in Mazerolle (2006)          
min.trap$Type <- relevel(min.trap$Type, ref = "UPLAND") 

##run model
m1 <- glm(Num_anura ~ Type + log.Perimeter + Num_ranatra,
          family = poisson, offset = log(Effort),
          data = min.trap) 
##null model
m0 <- glm(Num_anura ~ 1,
          family = poisson, offset = log(Effort),
          data = min.trap) 

##check c-hat for global model
c_hat(m1) #uses Pearson's chi-square/df

##likelihood ratio test corrected for overdispersion
anovaOD(mod.simple = m0, mod.complex = m1, c.hat = c_hat(m1))
##compare without overdispersion correction
anovaOD(mod.simple = m0, mod.complex = m1)


##example with occupancy model
\dontrun{
##load unmarked package
if(require(unmarked)){
   
   data(bullfrog)
     
   ##detection data
   detections <- bullfrog[, 3:9]

   ##assemble in unmarkedFrameOccu
   bfrog <- unmarkedFrameOccu(y = detections)
     
   ##run model
   fm <- occu(~ 1 ~ Reed.presence, data = bfrog)
   ##null model
   fm0 <- occu(~ 1 ~ 1, data = bfrog)

   ##check GOF
   ##GOF <- mb.gof.test(fm, nsim = 1000)
   ##estimate of c-hat:  1.89

   ##display results after overdispersion adjustment
   anovaOD(fm0, fm, c.hat = 1.89)

   detach(package:unmarked)
}
}
}

\keyword{models}

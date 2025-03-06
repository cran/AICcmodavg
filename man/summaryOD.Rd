\name{summaryOD}
\Rdversion{1.1}
\alias{summaryOD}
\alias{summaryOD.glm}
\alias{summaryOD.unmarkedFitOccu}
\alias{summaryOD.unmarkedFitColExt}
\alias{summaryOD.unmarkedFitOccuRN}
\alias{summaryOD.unmarkedFitPCount}
\alias{summaryOD.unmarkedFitPCO}
\alias{summaryOD.unmarkedFitDS}
\alias{summaryOD.unmarkedFitGDS}
\alias{summaryOD.unmarkedFitOccuFP}
\alias{summaryOD.unmarkedFitMPois}
\alias{summaryOD.unmarkedFitGMM}
\alias{summaryOD.unmarkedFitGPC}
\alias{summaryOD.unmarkedFitOccuMulti}
\alias{summaryOD.unmarkedFitOccuMS}
\alias{summaryOD.unmarkedFitOccuTTD}
\alias{summaryOD.unmarkedFitMMO}
\alias{summaryOD.unmarkedFitDSO}
\alias{summaryOD.unmarkedFitGOccu}
\alias{summaryOD.unmarkedFitOccuComm}
\alias{summaryOD.glmerMod}
\alias{summaryOD.maxlikeFit}
\alias{summaryOD.multinom}
\alias{summaryOD.vglm}
\alias{print.summaryOD}

\title{
  Display Model Summary Corrected for Overdispersion
}
\description{
This function displays the estimates of a model with standard errors
corrected for overdispersion for a variety of model classes.  The output
includes either confidence intervals based on the normal approximation
or Wald hypothesis tests corrected for overdispersion.
}
\usage{
summaryOD(mod, c.hat = 1, conf.level = 0.95, 
          out.type = "confint", \dots)

\method{summaryOD}{glm}(mod, c.hat = 1, conf.level = 0.95, 
     out.type = "confint", \dots)

\method{summaryOD}{unmarkedFitOccu}(mod, c.hat = 1, conf.level = 0.95,
     out.type = "confint", \dots)

\method{summaryOD}{unmarkedFitColExt}(mod, c.hat = 1, conf.level = 0.95, 
     out.type = "confint", \dots)

\method{summaryOD}{unmarkedFitOccuRN}(mod, c.hat = 1, conf.level = 0.95, 
     out.type = "confint", \dots)

\method{summaryOD}{unmarkedFitPCount}(mod, c.hat = 1, conf.level = 0.95, 
     out.type = "confint", \dots)

\method{summaryOD}{unmarkedFitPCO}(mod, c.hat = 1, conf.level = 0.95, 
     out.type = "confint", \dots)

\method{summaryOD}{unmarkedFitDS}(mod, c.hat = 1, conf.level = 0.95, 
     out.type = "confint", \dots)

\method{summaryOD}{unmarkedFitGDS}(mod, c.hat = 1, conf.level = 0.95, 
     out.type = "confint", \dots)

\method{summaryOD}{unmarkedFitOccuFP}(mod, c.hat = 1, conf.level = 0.95, 
     out.type = "confint", \dots)

\method{summaryOD}{unmarkedFitMPois}(mod, c.hat = 1, conf.level = 0.95, 
     out.type = "confint", \dots)

\method{summaryOD}{unmarkedFitGMM}(mod, c.hat = 1, conf.level = 0.95, 
     out.type = "confint", \dots)

\method{summaryOD}{unmarkedFitGPC}(mod, c.hat = 1, conf.level = 0.95, 
     out.type = "confint", \dots)

\method{summaryOD}{unmarkedFitOccuMulti}(mod, c.hat = 1, conf.level = 0.95,
     out.type = "confint", \dots)

\method{summaryOD}{unmarkedFitOccuMS}(mod, c.hat = 1, conf.level = 0.95,
     out.type = "confint", \dots)

\method{summaryOD}{unmarkedFitOccuTTD}(mod, c.hat = 1, conf.level = 0.95,
     out.type = "confint", \dots)

\method{summaryOD}{unmarkedFitMMO}(mod, c.hat = 1, conf.level = 0.95,
     out.type = "confint", \dots)

\method{summaryOD}{unmarkedFitDSO}(mod, c.hat = 1, conf.level = 0.95,
     out.type = "confint", \dots)

\method{summaryOD}{unmarkedFitGOccu}(mod, c.hat = 1, conf.level = 0.95,
     out.type = "confint", \dots)

\method{summaryOD}{unmarkedFitOccuComm}(mod, c.hat = 1, conf.level = 0.95,
     out.type = "confint", \dots)

\method{summaryOD}{glmerMod}(mod, c.hat = 1, conf.level = 0.95, 
     out.type = "confint", \dots)

\method{summaryOD}{maxlikeFit}(mod, c.hat = 1, conf.level = 0.95, 
     out.type = "confint", \dots)

\method{summaryOD}{multinom}(mod, c.hat = 1, conf.level = 0.95, 
     out.type = "confint", \dots)

\method{summaryOD}{vglm}(mod, c.hat = 1, conf.level = 0.95, 
     out.type = "confint", \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mod}{
    an object of class \code{glm}, \code{glmmTMB},
    \code{maxlikeFit}, \code{mer}, \code{merMod}, \code{multinom},  
    \code{vglm}, and various \code{unmarkedFit} classes containing 
    the output of a model.
  }
  
  \item{c.hat}{
    value of overdispersion parameter (i.e., variance inflation factor)
    such as that obtained from \code{c_hat}, \code{mb.gof.test}, or
    \code{Nmix.gof.test}.
  }

  \item{conf.level}{
    the confidence level (\eqn{1 - \alpha}) requested for the computation of
    confidence intervals.
  }

  \item{out.type}{
    the type of summary requested for each parameter estimate. If
    \code{out.type = "confint"}, computes confidence intervals corrected
    for overdispersion, whereas \code{out.type = "nhst"} conducts
    null-hypothesis statistical testing corrected for overdispersion.
  }

  \item{\dots}{
    additional arguments passed to the function.
  }

}  
\details{
  Overdispersion occurs when the variance in the data exceeds that
  expected from a theoretical distribution such as the Poisson or
  binomial (McCullagh and Nelder 1989, Burnham and Anderson 2002).
  When the model is correct, small values of c-hat (1 < c-hat < 4) can
  reflect minor deviations from model assumptions (Burnham and Anderson
  2002).  In such cases, it is possible to adjust standard errors of
  parameter estimates by multiplying with \code{sqrt(c.hat)} (McCullagh
  and Nelder 1989).  This is the correction applied by
  \code{summaryOD}.

  Depending on the type of summary requested, i.e.,
  \code{out.type = "confint"} or \code{out.type = "nhst"},
  \code{summaryOD} will return either confidence intervals based on the
  normal approximation or Wald tests for each parameter estimate
  (Agresti 1990).

  For binomial distributions, note that values of \code{c.hat} > 1 are
  only appropriate with trials > 1 (i.e., \code{success/trial} or
  \code{cbind(success, failure)} syntax).  The function supports
  different model types such as Poisson GLM's and GLMM's, single-season
  occupancy models (MacKenzie et al. 2002), dynamic occupancy models
  (MacKenzie et al. 2003), or \emph{N}-mixture models (Royle 2004, Dail
  and Madsen 2011).
}
\value{
\code{summaryOD} returns an object of class \code{summaryOD} as a list with
the following components:

\item{out.type}{the type of output requested by the user.}
\item{c.hat}{the \code{c.hat} estimate used to adjust standard errors.} 
\item{conf.level}{the confidence level used to compute confidence
  intervals around the estimates.}
\item{outMat}{the output of the model corrected for overdispersion
  organized in a matrix.}
}

\references{
  Agresti, A. (2002) \emph{Categorical Data Analysis}. Second
  edition. John Wiley and Sons: New Jersey.

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

  Mazerolle, M. J. (2006) Improving data analysis in herpetology: using
  Akaike's Information Criterion (AIC) to assess the strength of
  biological hypotheses. \emph{Amphibia-Reptilia} \bold{27}, 169--180. 
  
  McCullagh, P., Nelder, J. A. (1989) \emph{Generalized Linear
  Models}. Second edition. Chapman and Hall: New York. 

  Royle, J. A. (2004) \emph{N}-mixture models for estimating population
  size from spatially replicated counts. \emph{Biometrics} \bold{60}, 
  108--115.

}
\author{
Marc J. Mazerolle
}

\seealso{
  \code{\link{c_hat}}, \code{\link{mb.gof.test}},
  \code{\link{Nmix.gof.test}}, \code{\link{anovaOD}} 
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

##check c-hat for global model
c_hat(m1) #uses Pearson's chi-square/df

##display results corrected for overdispersion
summaryOD(m1, c_hat(m1))
summaryOD(m1, c_hat(m1), out.type = "nhst")

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
   fm <- occu(~ 1 ~ 1, data = bfrog)

   ##check GOF
   ##GOF <- mb.gof.test(fm, nsim = 1000)
   ##estimate of c-hat:  1.89

   ##display results after overdispersion adjustment
   summaryOD(fm, c.hat = 1.89)
   summaryOD(fm, c.hat = 1.89, out.type = "nhst")

   detach(package:unmarked)
}
}
}

\keyword{models}

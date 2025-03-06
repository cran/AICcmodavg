\name{AICcCustom}
\Rdversion{1.1}
\alias{AICcCustom}
\title{
  Compute AIC, AICc, QAIC, and QAICc from User-supplied Input
}
\description{
  This function computes Akaike's information criterion (AIC), the
  second-order AIC (AICc), as well as their quasi-likelihood
  counterparts (QAIC, QAICc) from user-supplied input instead of
  extracting the values automatically from a model object.  This
  function is particularly useful for output imported from other
  software or for model classes that are not currently supported by
  \code{AICc}.
}
\usage{
AICcCustom(logL, K, return.K = FALSE, second.ord = TRUE, nobs = NULL,
           c.hat = 1)
}

\arguments{
  \item{logL}{
    the value of the model log-likelihood.
}
\item{K}{
  the number of estimated parameters in the model.
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
  the sample size required to compute the AICc or QAICc.
}
 \item{c.hat}{
   value of overdispersion parameter (i.e., variance inflation factor)
   such as that obtained from \code{c_hat}.  Note that values of
   \code{c.hat} different from 1 are only appropriate for binomial GLM's
   with trials > 1 (i.e., success/trial or cbind(success, failure)
   syntax), with Poisson GLM's, single-season or dynamic occupancy
   models (MacKenzie et al. 2002, 2003), \emph{N}-mixture models (Royle
   2004, Dail and Madsen 2011), or capture-mark-recapture models (e.g.,
   Lebreton et al. 1992).  If c.hat > 1, \code{AICcCustom} will return
   the quasi-likelihood analogue of the information criterion requested.
}
}
\details{
  \code{AICcCustom} computes one of the following four information criteria:
  
  Akaike's information criterion (AIC, Akaike 1973), the  second-order
  or small sample AIC (AICc, Sugiura   1978, Hurvich and Tsai 1989, 1991), the
  quasi-likelihood AIC (QAIC, Burnham and Anderson 2002), and the
  quasi-likelihood AICc (QAICc, Burnham and Anderson 2002).
}
\value{
  \code{AICcCustom} returns the AIC, AICc, QAIC, or QAICc, or the number
  of estimated parameters, depending on the values of the arguments.
}

\note{
  The actual (Q)AIC(c) values are not really interesting in themselves,
  as they depend directly on the data, parameters estimated, and
  likelihood function.  Furthermore, a single value does not tell much
  about model fit.  Information criteria become relevant when compared
  to one another for a given data set and set of candidate models.
}

\references{
  Akaike, H. (1973) Information theory as an extension of the maximum
  likelihood principle. In:  \emph{Second International Symposium on
  Information Theory}, pp. 267--281. Petrov, B.N., Csaki, F., Eds,
  Akademiai Kiado, Budapest.

  Burnham, K. P., Anderson, D. R. (2002) \emph{Model Selection and
  Multimodel Inference: a practical information-theoretic
  approach}. Second edition. Springer: New York.

  Dail, D., Madsen, L. (2011) Models for estimating abundance from
  repeated counts of an open population. \emph{Biometrics} \bold{67},
  577--587. 
  
  Hurvich, C. M., Tsai, C.-L. (1989) Regression and time series model
  selection in small samples. \emph{Biometrika} \bold{76}, 297--307.  

  Hurvich, C. M., Tsai, C.-L. (1991) Bias of the corrected AIC
  criterion for underfitted regression and time series
  models. \emph{Biometrika} \bold{78}, 499--509. 

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

  Sugiura, N. (1978) Further analysis of the data by Akaike's
  information criterion and the finite corrections. \emph{Communications
  in Statistics:  Theory and Methods} \bold{A7}, 13--26.
}
\author{
  Marc J. Mazerolle
}

\seealso{
\code{\link{AICc}}, \code{\link{aictabCustom}}, \code{\link{confset}}, 
\code{\link{evidence}}, \code{\link{c_hat}}, \code{\link{modavgCustom}}
}


\examples{
##cement data from Burnham and Anderson (2002, p. 101)
data(cement)
##run multiple regression - the global model in Table 3.2
glob.mod <- lm(y ~ x1 + x2 + x3 + x4, data = cement)

##extract log-likelihood
LL <- logLik(glob.mod)[1]

##extract number of parameters
K.mod <- coef(glob.mod) + 1

##compute AICc with full likelihood
AICcCustom(LL, K.mod, nobs = nrow(cement))
}
\keyword{models}

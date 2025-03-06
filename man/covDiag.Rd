\name{covDiag}
\Rdversion{1.1}
\alias{covDiag}
\alias{covDiag.default}
\alias{covDiag.unmarkedFitPCount}
\alias{covDiag.unmarkedFramePCount}
\alias{print.covDiag}
\title{
Compute Covariance Diagnostic for Lambda in \emph{N}-mixture Models
}
\description{
This function extracts the covariance diagnostic of Dennis et al. (2015) 
for lambda in \emph{N}-mixture models (Royle 2004) of the
\code{unmarkedFitPCount} class as well as in data frames of the
\code{unmarkedFramePcount} class.
}
\usage{
covDiag(object, \dots)

\method{covDiag}{unmarkedFitPCount}(object, \dots)

\method{covDiag}{unmarkedFramePCount}(object, \dots)

}

\arguments{
  \item{object}{
    an object of class \code{unmarkedFitPCount} or \code{unmarkedFramePCount}.
  }

  \item{\dots}{
    additional arguments passed to the function.
  }
}

\details{
  This function extracts the covariance diagnostic developed by Dennis 
  et al. (2015) for lambda in \emph{N}-mixture models.  Values <= 0
  suggest sparse data and potential problems during model fitting.
  \code{covDiag} can take data frames of the \code{unmarkedFramePcount}
  class as input.  For convenience, the function also takes the repeated
  count model object as input, extracts the raw data, and computes the
  covariance diagnostic.  Thus, different models on the same data set
  will have identical values for this covariance diagnostic.   
}
\value{
  \code{covDiag} returns a list with the following components:
  
  \item{cov.diag}{the value of the covariance diagnostic.}
  \item{message}{a string indicating whether a warning was issued (i.e.,
  \code{"Warning: lambda is infinite, data too sparse"}) or not (i.e.,
  \code{NULL}).}
  }

  \references{
    Dennis, E. B., Morgan, B. J. T., Ridout, M. S. (2015)
    Computational aspects of \emph{N}-mixture models. \emph{Biometrics}
    \bold{71}, 237--246.
    
    Royle, J. A. (2004) \emph{N}-mixture models for estimating population
    size from spatially replicated counts. \emph{Biometrics} \bold{60},
    108--115.
  }
  
\author{
Marc J. Mazerolle
}

\seealso{
\code{\link{modavg}}, \code{\link{modavgPred}},
\code{\link{Nmix.chisq}}, \code{\link{Nmix.gof.test}},
\code{\link{predictSE}}, \code{\link[unmarked]{pcount}} 
}
\examples{
##modified example from ?pcount
\dontrun{
if(require(unmarked)){
##Simulate data
set.seed(3)
nSites <- 100
nVisits <- 3
##covariate
x <- rnorm(nSites)               
beta0 <- 0
beta1 <- 1
##expected counts
lambda <- exp(beta0 + beta1*x)   
N <- rpois(nSites, lambda)      
y <- matrix(NA, nSites, nVisits)
p <- c(0.3, 0.6, 0.8)           
for(j in 1:nVisits) {
  y[,j] <- rbinom(nSites, N, p[j])
}
## Organize data
visitMat <- matrix(as.character(1:nVisits),
                   nSites, nVisits, byrow=TRUE)
     
umf <- unmarkedFramePCount(y=y, siteCovs=data.frame(x=x),
                           obsCovs=list(visit=visitMat))
## Fit model
fm1 <- pcount(~ visit ~ 1, umf, K=50)
covDiag(fm1)

##sparser data
p <- c(0.01, 0.001, 0.01)           
for(j in 1:nVisits) {
  y[,j] <- rbinom(nSites, N, p[j])
}
## Organize data
visitMat <- matrix(as.character(1:nVisits),
                   nSites, nVisits, byrow=TRUE)
     
umf <- unmarkedFramePCount(y=y, siteCovs=data.frame(x=x),
                           obsCovs=list(visit=visitMat))
## Fit model
fm.sparse <- pcount(~ visit ~ 1, umf, K=50)
covDiag(fm.sparse)
}
}
}
\keyword{models}

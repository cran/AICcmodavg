\name{extractLL}
\Rdversion{1.1}
\alias{extractLL}
\alias{extractLL.default}
\alias{extractLL.coxph}
\alias{extractLL.coxme}
\alias{extractLL.lmekin}
\alias{extractLL.maxlikeFit}
\alias{extractLL.unmarkedFit}
\alias{extractLL.vglm}

\title{
Extract Log-Likelihood of Model
}
\description{
This function extracts the log-likelihood from an object of \code{coxme},
\code{coxph}, \code{lmekin}, \code{maxlikeFit}, \code{vglm}, or various
\code{unmarkedFit} classes.
}
\usage{
extractLL(mod, \dots)

\method{extractLL}{coxme}(mod, type = "Integrated", \dots)

\method{extractLL}{coxph}(mod, \dots)

\method{extractLL}{lmekin}(mod, \dots)

\method{extractLL}{maxlikeFit}(mod, \dots)

\method{extractLL}{unmarkedFit}(mod, \dots)

\method{extractLL}{vglm}(mod, \dots)
}

\arguments{
  \item{mod}{
    an object of \code{coxme}, \code{coxph}, \code{lmekin},
  \code{maxlikeFit}, \code{vglm}, or \code{unmarkedFit} class resulting
  from the fit of \code{distsamp}, \code{gdistsamp}, \code{gmultmix},
  \code{multinomPois}, \code{gpcount}, \code{occu}, \code{occuRN},
  \code{colext}, \code{pcount}, or \code{pcountOpen}.} 
\item{\dots}{additional arguments passed to the function.}
 \item{type}{
   a character string indicating whether the integrated partial
   likelihood ("Integrated") or penalized likelihood ("Penalized") is to
   be used for a \code{coxme} object.}
}
\details{
This utility function extracts the information from a \code{coxme},
  \code{coxph}, \code{lmekin}, \code{maxlikeFit}, \code{vglm}, or
  \code{unmarkedFit} object resulting from \code{distsamp},
  \code{gdistsamp}, \code{gmultmix}, \code{multinomPois},
  \code{gpcount}, \code{occu}, \code{occuRN}, \code{colext},
  \code{pcount}, or \code{pcountOpen}. 
}
\value{
  These functions return the value of the log-likelihood of the model
  and associated degrees of freedom.
}

\author{
Marc J. Mazerolle
}

\seealso{
\code{\link{AICc}}, \code{\link{aictab}}, \code{\link[coxme]{coxme}},
\code{\link[survival]{coxph}}, \code{\link[coxme]{lmekin}},
\code{\link[maxlike]{maxlike}}, \code{\link[unmarked]{distsamp}},
\code{\link[unmarked]{gdistsamp}}, \code{\link[unmarked]{occu}},
\code{\link[unmarked]{occuRN}}, \code{\link[unmarked]{colext}},
\code{\link[unmarked]{pcount}}, \code{\link[unmarked]{pcountOpen}}
}

\examples{
##single-season occupancy model example modified from ?occu
\dontrun{
require(unmarked)
##single season
data(frogs)
pferUMF <- unmarkedFrameOccu(pfer.bin)
## add some fake covariates for illustration
siteCovs(pferUMF) <- data.frame(sitevar1 = rnorm(numSites(pferUMF)),
                                sitevar2 = rnorm(numSites(pferUMF))) 
     
## observation covariates are in site-major, observation-minor order
obsCovs(pferUMF) <- data.frame(obsvar1 = rnorm(numSites(pferUMF) *
                                 obsNum(pferUMF))) 

##run model set
fm1 <- occu(~ obsvar1 ~ sitevar1, pferUMF)

##extract log-likelihood
extractLL(fm1)
detach(package:unmarked)
}
}
\keyword{models}


\name{extract.LL}
\Rdversion{1.1}
\alias{extract.LL}
\alias{extract.LL.coxph}
\alias{extract.LL.unmarked}
\title{
Extract Log-Likelihood of Model
}
\description{
This function extracts the log-likelihood from an object of 'coxph' or various
'unmarkedFit' classes.
}
\usage{
extract.LL(mod)
extract.LL.coxph(mod)
extract.LL.unmarked(mod)
}

\arguments{
  \item{mod}{
    an object of 'coxph' or 'unmarkedFit' class resulting from the fit
    of 'distsamp', 'gdistsamp', 'occu', 'occuRN', 'colext', 'pcount', or
    'pcountOpen'

}
}
\details{
This utility function extracts the information from a 'coxph' or
'unmarkedFit' object resulting from 'distsamp', 'gdistsamp', 'occu',
'occuRN', 'colext', 'pcount', or 'pcountOpen'.
}
\value{
  'extract.LL.coxph' and 'extract.LL.unmarked' return the value of the
  log-likelihood of the model and associated degrees of freedom. 
}

\author{
Marc J. Mazerolle
}

\seealso{
\code{\link{AICc}}, \code{\link{aictab}}, \code{\link[survival]{coxph}},
\code{\link[unmarked]{distsamp}}, \code{\link[unmarked]{gdistsamp}},
\code{\link[unmarked]{occu}}, \code{\link[unmarked]{occuRN}},
\code{\link[unmarked]{colext}}, \code{\link[unmarked]{pcount}},
\code{\link[unmarked]{pcountOpen}} 
}

\examples{
##single-season occupancy model example modified from ?occu
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
extract.LL(fm1)
detach(package:unmarked)



##modified example of Cox regression from ?coxph
require(survival)
##Create a simple data set for a time-dependent model
test2 <- list(start=c(1, 2, 5, 2, 1, 7, 3, 4, 8, 8),
              stop =c(2, 3, 6, 7, 8, 9, 9, 9,14,17),
              event=c(1, 1, 1, 1, 1, 1, 1, 0, 0, 0),
              x    =c(1, 0, 0, 1, 0, 1, 1, 1, 0, 0) )
m.cox <- coxph(Surv(start, stop, event) ~ x, test2)
extract.LL(m.cox)
detach(package:survival)

}
\keyword{models}


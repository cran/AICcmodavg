\name{checkConv}
\Rdversion{1.1}
\alias{checkConv}
\alias{checkConv.default}
\alias{checkConv.betareg}
\alias{checkConv.clm}
\alias{checkConv.clmm}
\alias{checkConv.glm}
\alias{checkConv.glmmTMB}
\alias{checkConv.hurdle}
\alias{checkConv.lavaan}
\alias{checkConv.maxlikeFit}
\alias{checkConv.merMod}
\alias{checkConv.lmerModLmerTest}
\alias{checkConv.multinom}
\alias{checkConv.nls}
\alias{checkConv.polr}
\alias{checkConv.unmarkedFit}
\alias{checkConv.zeroinfl}
\alias{print.checkConv}
\title{
Check Convergence of Fitted Model
}
\description{
This function checks the convergence information contained in models of
various classes.
}
\usage{
checkConv(mod, \dots)

\method{checkConv}{betareg}(mod, \dots)

\method{checkConv}{clm}(mod, \dots)

\method{checkConv}{clmm}(mod, \dots)

\method{checkConv}{glm}(mod, \dots)

\method{checkConv}{glmmTMB}(mod, \dots)

\method{checkConv}{hurdle}(mod, \dots)

\method{checkConv}{lavaan}(mod, \dots)

\method{checkConv}{maxlikeFit}(mod, \dots)

\method{checkConv}{merMod}(mod, \dots)

\method{checkConv}{lmerModLmerTest}(mod, \dots)

\method{checkConv}{multinom}(mod, \dots)

\method{checkConv}{nls}(mod, \dots)

\method{checkConv}{polr}(mod, \dots)

\method{checkConv}{unmarkedFit}(mod, \dots)

\method{checkConv}{zeroinfl}(mod, \dots)
}

\arguments{
  \item{mod}{
    an object containing the output of a model of the classes mentioned
  above. 
  }

  \item{\dots}{
    additional arguments passed to the function.
  }
}

\details{
  This function checks the element of a model object that contains the
  convergence information from the optimization function.  The function
  is currently implemented for models of classes \code{betareg},
  \code{clm}, \code{clmm}, \code{glm}, \code{glmmTMB}, \code{hurdle},
  \code{lavaan}, \code{maxlikeFit}, \code{merMod},
  \code{lmerModLmerTest}, \code{multinom}, \code{nls}, \code{polr},
  \code{unmarkedFit}, and \code{zeroinfl}.  The function is particularly
  useful for functions with several groups of parameters, such as those
  of the \code{unmarked} package (Fiske and Chandler, 2011).
}
\value{
  \code{checkConv} returns a list with the following components:
  
  \item{converged}{a logical value indicating whether the algorithm
  converged or not.}
  \item{message}{a string containing the message from the optimization
  function.}
  }

  \references{
    Fiske, I., Chandler, R. (2011) unmarked: An R Package for
    fitting hierarchical models of wildlife occurrence and
    abundance. \emph{Journal of Statistical Software} \bold{43}, 1--23.
  }
  
\author{
Marc J. Mazerolle
}

\seealso{
\code{\link{checkParms}}, \code{\link{covDiag}},
\code{\link{mb.gof.test}}, \code{\link{Nmix.gof.test}}
}
\examples{
##example modified from ?pcount
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
checkConv(fm1)
}
}
}
\keyword{models}

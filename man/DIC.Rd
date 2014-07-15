\name{DIC}
\Rdversion{1.1}
\alias{DIC}
\alias{DIC.default}
\alias{DIC.bugs}
\alias{DIC.rjags}
\title{
Computing DIC
}
\description{
Functions to extract deviance information criterion (DIC).
}
\usage{
DIC(mod, return.pD = FALSE, \dots) 

\method{DIC}{bugs}(mod, return.pD = FALSE, \dots)

\method{DIC}{rjags}(mod, return.pD = FALSE, \dots)
}

\arguments{
  \item{mod}{
    an object of class \code{bugs} or \code{rjags} containing the output
    of a model.
  }
  
  \item{return.pD}{
    logical.  If \code{FALSE}, the function returns the DIC.  If
    \code{TRUE}, the function returns the effective number of estimated
    parameters (pD) for a given model.
  }
  \item{\dots}{
  additional arguments passed to the function.
  }
}
\details{
\code{DIC} is implemented for \code{bugs} and \code{rjags} classes.  The 
    function extracts the deviance information criterion (DIC,
    Spiegelhalter et al. 2002) or the effective number of parameters
    (pD).
}
\value{
\code{DIC} the DIC or pD depending on the values of the arguments.
}

\note{
The actual DIC values are not really interesting in themselves, as
they depend directly on the data, parameters estimated, and
likelihood function.  Furthermore, a single value does not tell much
about model fit.  Information criteria become relevant when compared to
Yone another for a given data set and set of candidate models.  Model
selection with hierarchical models is problematic as the classic DIC is
not appropriate for such types of models (Millar 2009).
}

\references{
Millar, R. B. (2009) Comparison of hierarchical Bayesian models for
overdispersed count data using DIC and Bayes'
factors. \emph{Biometrics}, \bold{65}, 962--969.

Spiegelhalter, D. J., Best, N. G., Carlin, B. P., van der Linde,
A. (2002). Bayesian measures of complexity and fit. \emph{Journal of the
  Royal Statistical Society, Series B} \bold{64}, 583--639.
}
\author{
Marc J. Mazerolle
}

\seealso{
  \code{\link{aictab}}, \code{\link{dictab}}, \code{\link{confset}}, \code{\link{evidence}}
}



\examples{
##from ?jags example in R2jags package
if(require(R2jags)) {
##example model file
model.file <- system.file(package="R2jags", "model", "schools.txt")
file.show(model.file)     

##data
J <- 8.0
y <- c(28.4,7.9,-2.8,6.8,-0.6,0.6,18.0,12.2)
sd <- c(14.9,10.2,16.3,11.0,9.4,11.4,10.4,17.6)
     
##arrange data in list     
jags.data <- list (J = J, y = y, sd = sd)

##initial values
jags.inits <- function(){
  list(theta=rnorm(J, 0, 100), mu=rnorm(1, 0, 100),
       sigma=runif(1, 0, 100))
}

##parameters to be monitored
jags.parameters <- c("theta", "mu", "sigma")
     
##run model
schools.sim <- jags(data = jags.data, inits = jags.inits,
                    parameters = jags.parameters,
                    model.file = model.file,
                    n.chains = 3, n.iter = 10)
##note that n.iter should be higher

##extract DIC
DIC(schools.sim)
##extract pD
DIC(schools.sim, return.pD = TRUE)
detach(package:R2jags)
}
}   
\keyword{models}

\name{DIC}
\Rdversion{1.1}
\alias{DIC}
\alias{DIC.bugs}
\alias{DIC.rjags}
\title{
Computing DIC
}
\description{
Functions to extract deviance information criterion (DIC).
}
\usage{
DIC(mod, return.pD = FALSE) 

DIC.bugs(mod, return.pD = FALSE)

DIC.rjags(mod, return.pD = FALSE)
}

\arguments{
  \item{mod}{
    an object of class 'bugs' or 'rjags' containing the output of a
  model.
  }
  
  \item{return.pD}{
    logical.  If FALSE, the function returns the DIC.  If TRUE, the
    function returns K (effective number of estimated parameters, pD)
    for a given model.
  }  
}
\details{
'DIC' is a function that calls 'DIC.bugs' or 'DIC.rjags'. The current
    function is implemented for 'bugs' and 'rjags' classes.  The
    function extracts the deviance information criterion (DIC,
    Spiegelhalter et al. 2002). 
}
\value{
%'AICc' selects one of the functions below based on the class of the 
%object:

%'AICc.clm' returns the AIC or AICc depending on the values of the
%arguments.

'DIC', 'DIC.bugs', and 'DIC.rjags' return the DIC or pD depending on the
values of the arguments.
}

\note{
The actual DIC values are not really interesting in themselves, as
they depend directly on the data, parameters estimated, and
likelihood function.  Furthermore, a single value does not tell much
about model fit.  Information criteria become relevant when compared to
one another for a given data set and set of candidate models.
}

\references{
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
require(R2jags)
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
}   
\keyword{models}
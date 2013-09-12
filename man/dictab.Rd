\name{dictab}
\Rdversion{1.1}
\alias{dictab}
\alias{dictab.bugs}
\alias{dictab.rjags}
\alias{print.dictab}
\title{
Create Model Selection Tables from Bayesian Analyses
}
\description{
This function creates a model selection table based on the deviance
information criterion (DIC).  The table ranks the models based on the
DIC and also provides delta DIC and DIC weights.  'dictab' selects the
appropriate function to create the model selection table based on the 
object class.  The current version works with objects of 'bugs' and
'rjags' classes.
}
\usage{
dictab(cand.set, modnames, sort = TRUE)

dictab.bugs(cand.set, modnames, sort = TRUE)

dictab.rjags(cand.set, modnames, sort = TRUE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cand.set}{
a list storing each of the models in the candidate model set.
}
  \item{modnames}{
a character vector of model names to facilitate the identification of
each model in the model selection table.
}
  \item{sort}{
logical.  If true, the model selection table is ranked according to the
DIC values.
}

}
\details{
'dictab' is a function that calls 'dictab.bugs' and 'dictab.rjags'.  The
current function is implemented for 'bugs' classes. The function
constructs a model selection table based on the DIC (Spiegelhalter et
al. 2002).  Note that DIC might not be appropriate to select among a set
of hierarchical models and that modifications to the information
criterion have been proposed (Millar 2009).
}
\value{
'dictab', 'dictab.bugs', and 'dictab.rjags' create an object of class
'dictab' with the following components: 

\item{Modname}{the names of each model of the candidate model set.}
\item{pD}{the effective number of estimated parameters for each
  model.}
\item{DIC}{the deviance information criterion for each model.}
\item{Delta_DIC}{the delta DIC of each model, measuring the difference
  in DIC between each model and the top-ranked model.}
\item{ModelLik}{the relative likelihood of the model given the
  data (exp(-0.5*delta[i])).  This is not to be confused with the
  likelihood of the parameters given the data.  The relative likelihood
  can then be normalized across all models to get the model probabilities.} 
\item{DICWt}{the DIC weights, sensu Burnham and Anderson (2002) and
  Anderson (2008). These measures indicate the level of support (i.e.,
  weight of evidence) in favor of any given model being the most
  parsimonious among the candidate model set.}
\item{Cum.Wt}{the cumulative DIC weights.  These are only meaningful
if results in table are sorted in decreasing order of DIC weights
(i.e., sort = TRUE).}
\item{Deviance}{the deviance of each model.}
}
\references{
Anderson, D. R. (2008) \emph{Model-based Inference in the Life Sciences:
a primer on evidence}. Springer: New York.

Burnham, K. P., Anderson, D. R. (2002) \emph{Model Selection and
Multimodel Inference: a practical information-theoretic
approach}. Second edition. Springer: New York.

Spiegelhalter, D. J., Best, N. G., Carlin, B. P., van der Linde,
A. (2002). Bayesian measures of complexity and fit. \emph{Journal of the
  Royal Statistical Society, Series B} \bold{64}, 583--639.
}
\author{
Marc J. Mazerolle
}
\seealso{
\code{\link{aictab}}, \code{\link{confset}}, \code{\link{DIC}},
\code{\link{evidence}}
} 
\examples{
##from ?jags example in R2jags package
if(require(R2jags)) {
model.file <- system.file(package="R2jags", "model", "schools.txt")
file.show(model.file)     

##data
J <- 8.0
y <- c(28.4,7.9,-2.8,6.8,-0.6,0.6,18.0,12.2)
sd <- c(14.9,10.2,16.3,11.0,9.4,11.4,10.4,17.6)
     
     
jags.data <- list (J = J, y = y, sd = sd)
jags.inits <- function(){
  list(theta=rnorm(J, 0, 100), mu=rnorm(1, 0, 100),
       sigma=runif(1, 0, 100))
}
jags.parameters <- c("theta", "mu", "sigma")
     
##run model
schools.sim <- jags(data = jags.data, inits = jags.inits,
                    parameters = jags.parameters,
                    model.file = model.file,
                    n.chains = 3, n.iter = 10)
#note that n.iter should be higher

##set up in list
Cand.mods <- list(schools.sim)
Model.names <- "hierarchical model"
##other models can be added to Cand.mods
##to compare them to the top model

##model selection table
dictab(cand.set = Cand.mods, modnames = Model.names)
detach(package:R2jags)
}
}
\keyword{models}

\name{dictab}
\Rdversion{1.1}
\alias{dictab}
\alias{dictab.default}
\alias{dictab.AICbugs}
\alias{dictab.AICrjags}
\alias{print.dictab}
\title{
Create Model Selection Tables from Bayesian Analyses
}
\description{
This function creates a model selection table based on the deviance
information criterion (DIC).  The table ranks the models based on the
DIC and also provides delta DIC and DIC weights.  \code{dictab} selects
the appropriate function to create the model selection table based on
the object class.  The current version works with objects of \code{bugs}
and \code{rjags} classes.
}
\usage{
dictab(cand.set, modnames = NULL, sort = TRUE, \dots)

\method{dictab}{AICbugs}(cand.set, modnames = NULL, sort = TRUE, \dots)

\method{dictab}{AICrjags}(cand.set, modnames = NULL, sort = TRUE, \dots)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cand.set}{
a list storing each of the models in the candidate model set.
}
  \item{modnames}{
a character vector of model names to facilitate the identification of
each model in the model selection table. If \code{NULL}, the function
uses the names in the cand.set list of candidate models. If no names
appear in the list, generic names (e.g., \code{Mod1}, \code{Mod2}) are
supplied in the table in the same order as in the list of candidate
models.
}
  \item{sort}{
logical.  If \code{TRUE}, the model selection table is ranked according
to the DIC values.
}
\item{\dots}{
  additional arguments passed to the function.
}
}
\details{
\code{dictab} internally creates a new class for the \code{cand.set}
list of candidate models, according to the contents of the list.  The
current function is implemented for \code{bugs} and \code{jags}
classes. The function constructs a model selection table based on the
DIC (Spiegelhalter et al. 2002).  Note that DIC might not be appropriate
to select among a set of hierarchical models and that modifications to
the information criterion have been proposed (Millar 2009). 
}
\value{
\code{dictab} creates an object of class \code{dictab} with the
following components:

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
(i.e., \code{sort = TRUE}).}
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

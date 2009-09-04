\name{importance}
\Rdversion{1.1}
\alias{importance}
\alias{print.importance}
\title{
Compute Importance Values of Variable
}
\description{
This function calculates the relative importance of variables (w+) based
on the sum of Akaike weights (model probabilities) of the models that
include the variable.  Note that this measure of evidence is only
appropriate when the variable appears in the same number of models as
those that do not include the variable.
}
\usage{
importance(cand.set, parm, modnames, c.hat = 1, second.ord = TRUE,
nobs = NULL)  
}

\arguments{
  \item{cand.set}{
a list storing each of the models in the candidate model set.
}
  \item{parm}{
the parameter of interest for which a measure of relative importance is
  required. 
}
  \item{modnames}{
a character vector of model names to facilitate the identification of
each model in the model selection table.
}
  \item{c.hat}{
value of overdispersion parameter (i.e., variance inflation factor) such
  as that obtained from 'c_hat'.  Note that values of c.hat different 
  from 1 are only appropriate for binomial GLM's with trials > 1 (i.e.,
  success/trial or cbind(success, failure) syntax) or with Poisson
  GLM's. If c.hat > 1, 'aictab' will return the quasi-likelihood
  analogue of the information criteria requested.
}
  \item{second.ord}{
logical. If TRUE, the function returns the second-order
Akaike information criterion (i.e., AICc).
}
  \item{nobs}{
this argument allows to specify a numeric value other than total sample
size to compute the AICc.  This is relevant only for linear mixed models
where sample size is not straightforward.  In such cases, one might use
total number of observations or number of independent clusters as the
value of 'nobs'. }
}
\value{
'importance' returns an object of class 'importance' consisting of the
following components:
\item{parm}{the parameter for which an importance value is required.} 
\item{w.plus}{the parameter for which an importance value is required.} 
\item{w.minus}{the sum of Akaike weights for the models that exclude the
  parameter of interest} 


}
\references{
Burnham, K. P., and Anderson, D. R. (2002) \emph{Model Selection and
  Multimodel Inference: a practical information-theoretic
  approach}. Second edition. Springer: New York.
}
\author{
Marc J. Mazerolle
}

\seealso{
\code{\link{AICc}}, \code{\link{aictab}}, \code{\link{c_hat}},
\code{\link{modavg}}, \code{\link{evidence}}, \code{\link{confset}},
\code{\link{modavgpred}} 
}
\examples{
#example on Orthodont data set in nlme
require(nlme)

#set up candidate model list
Cand.models <- list()
Cand.models[[1]] <- lme(distance ~ age, data = Orthodont, method = "ML")
#random is ~ age | Subject
Cand.models[[2]] <- lme(distance ~ age + Sex, data = Orthodont, random =
~ 1, method = "ML")
Cand.models[[3]] <- lme(distance ~ 1, data = Orthodont, random = ~ 1,
method = "ML") 
Cand.models[[4]] <- lme(distance ~ Sex, data = Orthodont, random = ~ 1,
method = "ML") 

#create a vector of model names
Modnames <- NULL
for (i in 1:length(Cand.models)) {
Modnames[i] <- paste("mod", i, sep = "")
}

importance(cand.set = Cand.models, parm = "age", modnames = Modnames,
second.ord = TRUE, nobs = NULL)

}
\keyword{models}
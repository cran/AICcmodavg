\name{evidence}
\Rdversion{1.1}
\alias{evidence}
\alias{print.evidence}
\title{
Compute Evidence Ratio Between Two Models
}
\description{
This function compares two models of a candidate model set based on
their evidence ratio (i.e., ratio of Akaike weights).  The default
computes the evidence ratio of the Akaike weights between the top-ranked
model and a lower-ranked model.  You must supply a model selection table
of class 'aictab' as the first argument.
}
\usage{
evidence(aic.table, model.high = "top", model.low)
}
\arguments{
  \item{aic.table}{
a model selection table of class 'aictab' such as that produced by 'aictab'.
}
  \item{model.high}{
the top-ranked model (default), or alternatively, the name of another
model as it appears in the model selection table.
}
  \item{model.low}{
the name of a lower-ranked model such as it appears in the model
selection table.
}
}
\details{
The default compares the Akaike weights of the top-ranked model to
another model of the candidate model set.  The evidence ratio can be
interpreted as the number of times a given model is more parsimonious
than a lower-ranked model.  If one desires an evidence ratio that does not
involve a comparison with the top-ranking model, the name of the
required model must be specified in the model.high argument. 
}
\value{
'evidence' produces an object of class 'evidence' with the following
components:
\item{Model.high}{  
the top-ranked model among the two compared.
}
\item{Model.low}{
  the lower-ranked model among the two compared.
  }
\item{Ev.ratio}{
  the evidence ratio between the two models compared.    
}
}
\references{
Burnham, K. P., Anderson, D. R. (2002) \emph{Model Selection and
  Multimodel Inference: a practical information-theoretic
  approach}. Second edition. Springer: New York. 
}
\author{
Marc J. Mazerolle
}

\seealso{
\code{\link{AICc}}, \code{\link{aictab}}, \code{\link{c_hat}},
\code{\link{modavg}}, \code{\link{importance}}, \code{\link{confset}}, \code{\link{modavgpred}}
}
\examples{
#run example from Burnham and Anderson (2002, p. 183) with two
#non-nested models
data(pine)
Cand.set <- list( )
Cand.set[[1]] <- lm(y ~ x, data = pine)
Cand.set[[2]] <- lm(y ~ z, data = pine)

#assign model names
Modnames <- c("raw density", "density corrected for resin content")

#compute model selection table
aicctable.out <- aictab(cand.set = Cand.set, modnames = Modnames)

#compute evidence ratio
evidence(aic.table = aicctable.out, model.low = "density corrected for resin content")          




#run models for the Orthodont data set in nlme
require(nlme)

#set up candidate model list
Cand.models<-list()
Cand.models[[1]]<- lme(distance ~ age, data = Orthodont, method="ML") # random is ~ age
Cand.models[[2]]<- lme(distance ~ age + Sex, data = Orthodont, random = ~ 1, method="ML")
Cand.models[[3]]<-lme(distance ~ 1, data = Orthodont, random = ~ 1, method="ML")

#create a vector of model names
Modnames<-NULL
for (i in 1:length(Cand.models)) {
Modnames[i]<-paste("mod", i, sep="")
}

#compute AICc table
aic.table.1 <- aictab(cand.set = Cand.models, modnames = Modnames,
second.ord = TRUE)

#compute evidence ratio between best model and second-ranked model
evidence(aic.table = aic.table.1, model.high = "top", model.low = "mod1") 

#compute evidence ratio between second-best model and third-ranked model
evidence(aic.table = aic.table.1, model.high = "mod1", model.low =
"mod3")
}
\keyword{models}


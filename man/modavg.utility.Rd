\name{modavg.utility}
\Rdversion{1.1}
\alias{reverse.parm}
\alias{reverse.exclude}
\title{
Accomodate Different Specifications of Interaction Terms
}
\description{
These utility functions enable the user to specify differently
interaction terms (e.g., A:B, B:A) across models for model averaging.
Both functions are called internally by 'modavg'. 
}
\usage{
reverse.parm(parm)
reverse.exclude(exclude)
}

\arguments{
  \item{parm}{
a parameter to be model-averaged, enclosed between quotes, as it appears
  in the output of some models.
}
  \item{exclude}{
a list of interaction or polynomial terms appearing in some models, as
  they would appear in the call to the model function (i.e., A*B, A:B).
  Models containing elements from the list will be excluded to obtain a
  model-averaged estimate.
}
}
\details{
Both functions have been added to avoid problems when users are not
  consistent in the specification of interaction terms across
  models.
}
\value{
'reverse.parm' returns all possible combinations of an interaction term to
  identify models that include the 'parm' of interest and find the
  corresponding estimate and standard error in the model object.

'reverse.exclude' returns a list of all possible combinations of 'exclude' 
  to identify models that should be excluded when computing a
  model-averaged estimate.  

}

\author{
Marc J. Mazerolle
}

\seealso{
\code{\link{modavg}}, \code{\link{modavg.shrink}}, \code{\link{modavgpred}}
}
\examples{
##a main effect
reverse.parm(parm = "Ageyoung") #does not return anything


##an interaction term as it might appear in the output
reverse.parm(parm = "Ageyoung:time") #returns the reverse


##exclude two interaction terms
reverse.exclude(exclude = list("Age*time", "A:B")) #returns all combinations
reverse.exclude(exclude = list("Age:time", "A*B")) #returns all combinations
}
\keyword{models}

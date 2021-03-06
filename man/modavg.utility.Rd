\name{modavg.utility}
\Rdversion{1.1}
\alias{reverse.parm}
\alias{reverse.exclude}
\alias{formatCands}
\alias{formulaShort}
\title{
Various Utility Functions
}
\description{
\code{reverse.parm} and \code{reverse.exclude} reverse the order of
variables in an interaction term.

\code{formatCands} creates new classes for lists containing candidate
models.

\code{formulaShort} prints a succinct formula from an \code{unmarkedFit} object.

}
\usage{
reverse.parm(parm)
reverse.exclude(exclude)
formatCands(cand.set)
formulaShort(mod, unmarked.type = NULL)
}

\arguments{
  \item{parm}{
a parameter to be model-averaged, enclosed between quotes, as it appears
  in the output of some models.
}
  \item{exclude}{
a list of interaction or polynomial terms appearing in some models, as
  they would appear in the call to the model function (i.e., \code{A*B},
  \code{A:B}).  Models containing elements from the list will be
  excluded to obtain a model-averaged estimate.
}
\item{cand.set}{a list storing each of the models in the candidate model set.
}
\item{mod}{an object storing the result of an \code{unmarkedFit}.
}
\item{unmarked.type}{a character string specifying the type of parameter in
  an \code{unmarkedFit} for which a formula is requested.  This argument
  uses the character string for each parameter group as defined by
  \code{unmarked}.
}
}
\details{
  These utility functions are used internally by \code{aictab},
  \code{modavg}, and other related functions.

  \code{reverse.parm} and \code{reverse.exclude} enable the user to
  specify differently interaction terms (e.g., \code{A:B}, \code{B:A})
  across models for model averaging.  These functions have been added to
  avoid problems when users are not consistent in the specification of
  interaction terms across models.

  \code{formatCands} creates new classes for the list of candidate
  models based on the contents of the list.  These new classes are used
  for method dispatch.

  \code{formulaShort} is used by \code{anovaOD}.
}
\value{
  \code{reverse.parm} returns all possible combinations of an interaction
  term to identify models that include the \code{parm} of interest and
  find the corresponding estimate and standard error in the model object.

  \code{reverse.exclude} returns a list of all possible combinations of
  \code{exclude} to identify models that should be excluded when
  computing a model-averaged estimate.

  \code{formatCands} adds a new class to the list of candidate
  models based on the classes of the models.

  \code{formulaShort} creates a character string for the formula related
  to a given parameter type from an \code{unmarkedFit} object.
}

\author{
Marc J. Mazerolle
}

\seealso{
\code{\link{aictab}}, \code{\link{anovaOD}}, \code{\link{modavg}},
\code{\link{modavgShrink}}, \code{\link{modavgPred}} 
}
\examples{
##a main effect
reverse.parm(parm = "Ageyoung") #does not return anything


##an interaction term as it might appear in the output
reverse.parm(parm = "Ageyoung:time") #returns the reverse


##exclude two interaction terms
reverse.exclude(exclude = list("Age*time", "A:B"))
##returns all combinations
reverse.exclude(exclude = list("Age:time", "A*B"))
##returns all combinations


##Mazerolle (2006) frog water loss example
data(dry.frog)

##setup a subset of models of Table 1
Cand.models <- list( )
Cand.models[[1]] <- lm(log_Mass_lost ~ Shade + Substrate +
                       cent_Initial_mass + Initial_mass2,
                       data = dry.frog)
Cand.models[[2]] <- lm(log_Mass_lost ~ Shade + Substrate +
                       cent_Initial_mass + Initial_mass2 +
                       Shade:Substrate, data = dry.frog)
Cand.models[[3]] <- lm(log_Mass_lost ~ cent_Initial_mass +
                       Initial_mass2, data = dry.frog)

formatCands(Cand.models)

\dontrun{
require(unmarked)
data(bullfrog)
bfrog <- unmarkedFrameOccu(y = bullfrog[, c("V1", "V2", "V3", "V4")],
                           siteCovs = bullfrog[, 1:2])
fm1 <- occu(~ 1 ~ Reed.presence, data = bfrog)
formulaShort(fm1, unmarked.type = "state")
formulaShort(fm1, unmarked.type = "det")
}
}
\keyword{models}

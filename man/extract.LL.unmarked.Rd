\name{extract.LL.unmarked}
\Rdversion{1.1}
\alias{extract.LL.unmarked}
\title{
Extract Log-Likelihood of Model
}
\description{
This function extracts the log-likelihood from an object of various
'unmarkedFit' classes.
}
\usage{
extract.LL.unmarked(mod)
}

\arguments{
  \item{mod}{
an object of 'unmarkedFit' class resulting from the fit of 'occu',
  'occuRN', 'colext', 'pcount', or 'pcountOpen'.
}
}
\details{
This utility function extracts the information from an 'unmarkedFit'
  object resulting from 'occu', 'occuRN', 'colext', 'pcount', or
  'pcountOpen' The function is called by 'AICc.unmarked' and
  'aictab.unmarked'.
}
\value{
  'extract.LL.unmarked' returns the value of the log-likelihood of the model.
}

\author{
Marc J. Mazerolle
}

\seealso{
\code{\link{AICc}}, \code{\link{aictab}},
\code{\link[unmarked]{occu}}, \code{\link[unmarked]{occuRN}},
\code{\link[unmarked]{colext}}, \code{\link[unmarked]{pcount}},
\code{\link[unmarked]{pcountOpen}} 
}
\examples{
##modified example of single-season model with heterogeneity from ?occuRN
require(unmarked)
data(birds)
woodthrushUMF <- unmarkedFrameOccu(woodthrush.bin)

##survey occasion-specific detection probabilities
fm.wood.rn <- occuRN(~ obsNum ~ 1, woodthrushUMF)

##extract log-likelihood
extract.LL.unmarked(fm.wood.rn)

detach(package:unmarked)
}
\keyword{models}


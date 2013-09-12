\name{beetle}
\Rdversion{1.1}
\alias{beetle}
\docType{data}
\title{
Flour Beetle Data
}
\description{
This data set illustrates the acute mortality of flour beetles
(\emph{Tribolium confusum}) following 5 hour exposure to carbon disulfide
gas. 
}
\usage{data(beetle)}
\format{
  A data frame with 8 rows and 4 variables.
  \describe{
    \item{\code{Dose}}{dose of carbon disulfide in mg/L.}
    \item{\code{Number_tested}}{number of beetles exposed to given dose
      of carbon disulfide.} 
    \item{\code{Number_killed}}{number of beetles dead after 5 hour
      exposure to given dose of carbon disulfide.}
    \item{\code{Mortality_rate}}{proportion of total beetles found dead
      after 5 hour exposure.} 
  }
}
\details{
Burnham and Anderson (2002, p. 195) use this data set originally from
Young and Young (1998) to show model selection for binomial models with
different link functions (logit, probit, cloglog). 
}
\source{
Burnham, K. P., Anderson, D. R. (2002) \emph{Model Selection and
  Multimodel Inference: a practical information-theoretic
  approach}. Second edition. Springer: New York.

Young, L. J., Young, J. H. (1998) \emph{Statistical Ecology}. Kluwer
Academic Publishers: London.    
}
\examples{
data(beetle)
## maybe str(beetle) ; plot(beetle) ...
}
\keyword{datasets}

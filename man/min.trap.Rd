\name{min.trap}
\Rdversion{1.1}
\alias{min.trap}
\docType{data}
\title{
Anuran Larvae Counts in Minnow Traps Across Pond Type
}
\description{
This data set consists of counts of anuran larvae as a function of pond
type, pond perimeter, and presence of water scorpions (\emph{Ranatra} sp.).
}
\usage{data(min.trap)}
\format{
  A data frame with 24 observations on the following 6 variables.
  \describe{
    \item{\code{Type}}{pond type, denotes the location of ponds in
      either bog or upland environment}
    \item{\code{Num_anura}}{number of anuran larvae in minnow traps}
    \item{\code{Effort}}{number of trap nights (i.e., number of traps x
      days of trapping) in each pond}
    \item{\code{Perimeter}}{pond perimeter in meters}
    \item{\code{Num_ranatra}}{number of water scorpions trapped in
      minnow traps}
    \item{\code{log.Perimeter}}{natural log of perimeter}
  }
}
\details{
Mazerolle (2006) uses this data set to illustrate model selection for
Poisson regression with low overdispersion. 
}
\source{
Mazerolle, M. J. (2006) Improving data analysis in herpetology: using
Akaike's Information Criterion (AIC) to assess the strength of
biological hypotheses. \emph{Amphibia-Reptilia} \bold{27}, 169--180. 
}
\examples{
data(min.trap)
## maybe str(min.trap) ; plot(min.trap) ...
}
\keyword{datasets}

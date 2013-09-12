\name{pine}
\Rdversion{1.1}
\alias{pine}
\docType{data}
\title{
Strength of Pine Wood Based on the Density Adjusted for Resin Content
}
\description{
This data set consists of the strength of pine wood as a function of
density or density adjusted for resin content.
}
\usage{data(pine)}
\format{
  A data frame with 42 observations on the following 3 variables.
  \describe{
    \item{\code{y}}{pine wood strength.}
    \item{\code{x}}{pine wood density.}
    \item{\code{z}}{pine wood density adjusted for resin content.}
  }
}
\details{
Burnham and Anderson (2002, p. 183) use this data set originally from
Carlin and Chib (1995) to illustrate model selection for two competing
and non-nested models. 
}
\source{
Burnham, K. P., Anderson, D. R. (2002) \emph{Model Selection and
  Multimodel Inference: a practical information-theoretic
  approach}. Second edition. Springer: New York. 

Carlin, B. P., Chib, S. (1995) Bayesian model choice via Markov chain
  Monte Carlo methods. \emph{Journal of the Royal Statistical Society, 
  Series B} \bold{57}, 473--484. 
}
\examples{
data(pine)
## maybe str(pine) ; plot(pine) ...
}
\keyword{datasets}

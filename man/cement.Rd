\name{cement}
\Rdversion{1.1}
\alias{cement}
\docType{data}
\title{
Heat Expended Following Hardening of Portland Cement
}
\description{
This data set illustrates the heat expended (calories) from mixtures of four
different ingredients of Portland cement expressed as a percentage by weight.}
\usage{data(cement)}
\format{
  A data frame with 13 observations on the following 5 variables.
  \describe{
    \item{\code{x1}}{calcium aluminate.}
    \item{\code{x2}}{tricalcium silicate.}
    \item{\code{x3}}{tetracalcium alumino ferrite.}
    \item{\code{x4}}{dicalcium silicate.}
    \item{\code{y}}{calories of heat per gram of cement following 180
      days of hardening.}
  }
}
\details{
Burnham and Anderson (2002, p. 101) use this data set originally from Woods et
al. (1932) to select among a set of multiple regression models.
}
\source{
Burnham, K. P., Anderson, D. R. (2002) \emph{Model Selection and
  Multimodel Inference: a practical information-theoretic
  approach}. Second edition. Springer: New York. 

Woods, H., Steinour, H. H., Starke, H. R. (1932) Effect of composition
  of Portland cement on heat evolved during hardening. \emph{Industrial
  and Engineering Chemistry} \bold{24}, 1207--1214.   
}

\examples{
data(cement)
## maybe str(cement) ; plot(cement) ...
}
\keyword{datasets}

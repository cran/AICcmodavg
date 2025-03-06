\name{turkey}
\Rdversion{1.1}
\alias{turkey}
\docType{data}
\title{
Turkey Weight Gain
}
\description{
This one-way ANOVA data set presents turkey weight gain in pounds across
five diets. 
}
\usage{data(turkey)}
\format{
  A data frame with 30 rows and 2 variables.
  \describe{
    \item{\code{Diet}}{diet factor with 5 levels.}
    \item{\code{Weight.gain}}{weight gain in pounds.} 
  }
}
\details{
Heiberger and Holland (2004) and Ott (1993) analyze
this data set to illustrate one-way ANOVA.
}
\source{
Heiberger, R. M., Holland, B. (2004) \emph{Statistical Analysis and Data
  Display:  an intermediate course with examples in S-Plus, R, and
  SAS}. Springer: New York.

Ott, R. L. (1993) \emph{An Introduction to Statistical Methods and Data
  Analysis}. Fourth edition. Duxbury: Pacific Grove, CA.
}
\examples{
data(turkey)
str(turkey)
}
\keyword{datasets}

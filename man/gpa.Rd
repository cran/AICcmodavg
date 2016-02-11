\name{gpa}
\Rdversion{1.1}
\alias{gpa}
\docType{data}
\title{
GPA Data and Standardized Test Scores
}
\description{
This data set features the first-year college GPA and four standardized
tests conducted before matriculation. 
}
\usage{data(gpa)}
\format{
  A data frame with 20 rows and 5 variables.
  \describe{
    \item{\code{gpa.y}}{first-year GPA.}
    \item{\code{sat.math.x1}}{SAT math score.} 
    \item{\code{sat.verb.x2}}{SAT verbal score.}
    \item{\code{hs.math.x3}}{high school math score.}
    \item{\code{hs.engl.x4}}{high school English score.} 
  }
}
\details{
Burnham and Anderson (2002, p. 225) use this data set originally from
Graybill and Iyer (1994) to show model selection for all subsets
regression.
}
\source{
Burnham, K. P., Anderson, D. R. (2002) \emph{Model Selection and
  Multimodel Inference: a practical information-theoretic
  approach}. Second edition. Springer: New York.

Graybill, F. A., Iyer, H. K. (1994) \emph{Regression analysis: concepts
  and applications.} Duxbury Press: Belmont.
}
\examples{
data(gpa)
str(gpa)
}
\keyword{datasets}
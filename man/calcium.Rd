\name{calcium}
\Rdversion{1.1}
\alias{calcium}
\docType{data}
\title{
Blood Calcium Concentration in Birds
}
\description{
This data set features calcium concentration in the plasma of birds
of both sexes following a hormonal treatment.
}
\usage{data(calcium)}
\format{
  A data frame with 20 rows and 3 variables.
  \describe{
    \item{\code{Calcium}}{calcium concentration in mg/100 ml in the
      blood of birds.}
    \item{\code{Hormone}}{a factor with two levels indicating whether
      the bird received a hormonal treatment or not.}
    \item{\code{Sex}}{a factor with two levels coding for the sex of
      birds.}
  }
}
\details{
Zar (1984, p. 206) illustrates a two-way ANOVA with interaction with
this data set.
}
\source{
Zar, J. H. (1984) \emph{Biostatistical analysis}. Second
edition. Prentice Hall: Englewood Cliffs, New Jersey.
}
\examples{
data(calcium)
str(calcium)
}
\keyword{datasets}

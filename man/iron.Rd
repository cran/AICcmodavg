\name{iron}
\Rdversion{1.1}
\alias{iron}
\docType{data}
\title{
Iron Content in Food
}
\description{
This data set, originally from Adish et al. (1999), describes the iron
content of food cooked in different pot types.
}
\usage{data(iron)}
\format{
  A data frame with 36 rows and 3 variables.
  \describe{
    \item{\code{Pot}}{pot type, one of "aluminium", "clay", or "iron".} 
    \item{\code{Food}}{food type, one of "legumes", "meat", or
      "vegetables".}
    \item{\code{Iron}}{iron content measured in mg/100 g of food.}
  }
}
\details{
Heiberger and Holland (2004, p. 378) use these data as an exercise on 
two-way ANOVA with interaction.
}
\source{
Heiberger, R. M., Holland, B. (2004) \emph{Statistical Analysis and Data
  Display:  an intermediate course with examples in S-Plus, R, and
  SAS}. Springer: New York.

Adish, A. A., Esrey, S. A., Gyorkos, T. W., Jean-Baptiste, J., Rojhani,
A. (1999) Effect of consumption of food cooked in iron pots on iron
status and growth of young children: a randomised trial. \emph{The
  Lancet} \bold{353}, 712--716.
}
\examples{
data(iron)
str(iron)
}
\keyword{datasets}
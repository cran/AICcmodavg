\name{lizards}
\Rdversion{1.1}
\alias{lizards}
\docType{data}
\title{
Habitat Preference of Lizards
}
\description{
This data set describes the habitat preference of two species of
lizards, \emph{Anolis grahami} and \emph{A}. \emph{opalinus}, on the
island of Jamaica and is originally from Schoener (1970). McCullagh and
Nelder (1989) and Burnham and Anderson (2002) reanalyzed the data. Note
that a typo occurs in table 3.11 of Burnham and Anderson (2002).
}
\usage{data(lizards)}
\format{
  A data frame with 48 rows and 6 variables.
  \describe{
    \item{\code{Insolation}}{position of perch, either \code{shaded} or
      \code{sunny}.}
    \item{\code{Diameter}}{diameter of the perch, either \code{< 2 in}
      or \code{>= 2 in}.}
    \item{\code{Height}}{perch height, either \code{< 5} or
      \code{>= 5}.}
    \item{\code{Time}}{time of day, either \code{morning},
      \code{midday}, or \code{afternoon}.} 
    \item{\code{Species}}{species observed, either \code{grahami} or
      \code{opalinus}.}
    \item{\code{Counts}}{number of individuals observed.}
  }
}
\details{
Burnham and Anderson (2002, p. 137) use this data set originally from
Schoener (1970) to illustrate model selection for log-linear models. 
}
\source{
Burnham, K. P., Anderson, D. R. (2002) \emph{Model Selection and
  Multimodel Inference: a practical information-theoretic
  approach}. Second edition. Springer: New York.

McCullagh, P., Nelder, J. A. (1989) \emph{Generalized Linear
  Models}. Second edition. Chapman and Hall: New York.

Schoener, T. W. (1970) Nonsynchronous spatial overlap of lizards in
patchy habitats. \emph{Ecology} \bold{51}, 408--418. 
}
\examples{
data(lizards)
\dontrun{
##log-linear model as in Burnham and Anderson 2002, p. 137
##main effects
m1 <- glm(Counts ~ Insolation + Diameter + Height + Time + Species, 
          family = poisson, data = lizards)


##main effects and all second order interactions = base
m2 <- glm(Counts ~ Insolation + Diameter + Height + Time + Species +
          Insolation:Diameter + Insolation:Height + Insolation:Time +
          Insolation:Species + Diameter:Height + Diameter:Time +
          Diameter:Species + Height:Time + Height:Species +
          Time:Species, family = poisson, data = lizards)

##base - DT
m3 <- glm(Counts ~ Insolation + Diameter + Height + Time + Species +
          Insolation:Diameter + Insolation:Height + Insolation:Time +
          Insolation:Species + Diameter:Height + Diameter:Species +
          Height:Time + Height:Species + Time:Species, 
          family = poisson, data = lizards)

##base + HDI + HDT + HDS
m4 <- glm(Counts ~ Insolation + Diameter + Height + Time + Species +
          Insolation:Diameter + Insolation:Height + Insolation:Time +
          Insolation:Species + Diameter:Height + Diameter:Time +
          Diameter:Species + Height:Time + Height:Species +
          Time:Species + Height:Diameter:Insolation +
          Height:Diameter:Time + Height:Diameter:Species, 
          family = poisson, data = lizards)

##base + HDI + HDS + HIT + HIS + HTS + ITS
m5 <- glm(Counts ~ Insolation + Diameter + Height + Time + Species + 
          Insolation:Diameter + Insolation:Height + Insolation:Time +
          Insolation:Species + Diameter:Height + Diameter:Time +
          Diameter:Species + Height:Time + Height:Species +
          Time:Species + Height:Diameter:Insolation +
          Height:Diameter:Species + Height:Insolation:Time +
          Height:Insolation:Species + Height:Time:Species +
          Insolation:Time:Species, family = poisson, data = lizards)

##base + HIT + HIS + HTS + ITS
m6 <- glm(Counts ~ Insolation + Diameter + Height + Time + Species + 
          Insolation:Diameter + Insolation:Height + Insolation:Time +
          Insolation:Species + Diameter:Height + Diameter:Time +
          Diameter:Species + Height:Time + Height:Species +
          Time:Species + Height:Insolation:Time +
          Height:Insolation:Species +  Height:Time:Species +
          Insolation:Time:Species,  family = poisson, data = lizards)

##base + HIS + HTS + ITS
m7 <- glm(Counts ~ Insolation + Diameter + Height + Time + Species +
          Insolation:Diameter + Insolation:Height + Insolation:Time +
          Insolation:Species + Diameter:Height + Diameter:Time +
          Diameter:Species + Height:Time + Height:Species +
          Time:Species + Height:Insolation:Species +
          Height:Time:Species + Insolation:Time:Species, 
          family = poisson, data = lizards)

##base + HIT + HIS + HTS + ITS - DT
m8 <- glm(Counts ~ Insolation + Diameter + Height + Time + Species +
          Insolation:Diameter + Insolation:Height + Insolation:Time +
          Insolation:Species + Diameter:Height + Diameter:Species +
          Height:Time + Height:Species + Time:Species +
          Height:Insolation:Time + Height:Insolation:Species +
          Height:Time:Species + Insolation:Time:Species, 
          family = poisson, data = lizards)

##base + HIT + HIS + ITS - DT
m9 <- glm(Counts ~ Insolation + Diameter + Height + Time + Species +
          Insolation:Diameter + Insolation:Height + Insolation:Time +
          Insolation:Species + Diameter:Height + Diameter:Species +
          Height:Time + Height:Species + Time:Species +
          Height:Insolation:Time + Height:Insolation:Species +
          Insolation:Time:Species, 
          family = poisson, data = lizards)

##base + HIT + HIS - DT
m10 <- glm(Counts ~ Insolation + Diameter + Height + Time + Species +
          Insolation:Diameter + Insolation:Height + Insolation:Time +
          Insolation:Species + Diameter:Height + Diameter:Species +
          Height:Time + Height:Species + Time:Species +
          Height:Insolation:Time + Height:Insolation:Species,
          family = poisson, data = lizards)

##set up in list
Cands <- list(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10)
Modnames <- paste("m", 1:length(Cands), sep = "")

##model selection
library(AICcmodavg)
aictab(Cands, Modnames)
}
}
\keyword{datasets}

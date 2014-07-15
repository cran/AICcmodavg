\name{salamander}
\Rdversion{1.1}
\alias{salamander}
\docType{data}
\title{
Salamander capture-mark-recapture data
}
\description{
This is a capture-mark-recapture data set on male and female Spotted
Salamanders (\emph{Ambystoma maculatum}) recorded by Husting (1965). A
total of 1244 unique individuals were captured in pitfall traps at a
breeding site between 1959 and 1963.
}
\usage{data(salamander)}
\format{
  A data frame with 36 observations on the following 7 variables.
  \describe{
    \item{\code{T1959}}{a binary variable, either 1 (captured) or 0 (not
      captured) during the 1959 breeding season.}
    \item{\code{T1960}}{a binary variable, either 1 (captured) or 0 (not
      captured) during the 1960 breeding season.}
    \item{\code{T1961}}{a binary variable, either 1 (captured) or 0 (not
      captured) during the 1961 breeding season.}
    \item{\code{T1962}}{a binary variable, either 1 (captured) or 0 (not
      captured) during the 1962 breeding season.}
    \item{\code{T1963}}{a binary variable, either 1 (captured) or 0 (not
      captured) during the 1963 breeding season.}
    \item{\code{Males}}{a numeric variable indicating the total number
      of males with a given capture history. Negative values
      indicate losses on capture (animals not released on last captured). }
    \item{\code{Females}}{a numeric variable indicating the total number
      of females with a given capture history. Negative values
      indicate losses on capture (animals not released on last captured). }
  }
}
\details{
This data set is used to illustrate classic Cormack-Jolly-Seber and
related models (Cormack 1964, Jolly 1965, Seber 1965, Lebreton et
al. 1992). 
}
\source{
  Cormack, R. M. (1964) Estimates of survival from the sighting of
  marked animals. \emph{Biometrika} \bold{51}, 429--438.

  Husting, E. L. (1965) Survival and breeding structure in a population
  of \emph{Ambystoma maculatum}. \emph{Copeia} \bold{1965}, 352--362. 
  
  Jolly, G. M. (1965) Explicit estimates from capture-recapture data
  with both death and immigration: stochastic model. \emph{Biometrika} 
  \bold{52}, 225--247.

  Laake, J. L. (2013) \emph{RMark: an R interface for analysis of
  capture-recapture data with MARK.} Alaska Fisheries Science Center
  (AFSC), National Oceanic and Atmospheric Administration, National
  Marine Fisheries Service, AFSC Report 2013-01.
  
  Lebreton, J.-D., Burnham, K. P., Clobert, J., Anderson, D. R. (1992)
  Modeling survival and testing biological hypotheses using marked
  animals: a unified approach with case-studies. \emph{Ecological
    Monographs} \bold{62}, 67-118.

  Seber, G. A. F. (1965) A note on the multiple-recapture
  census. \emph{Biometrika} \bold{52}, 249--259.
}
\examples{
data(salamander)
str(salamander)

##convert raw capture data to capture histories
captures <- salamander[, c("T1959", "T1960", "T1961", "T1962", "T1963")]
salam.ch <- apply(captures, MARGIN = 1, FUN = function(i)
                  paste(i, collapse = ""))

##organize as a data frame readable by RMark package (Laake 2013)
##RMark requires at least one column called "ch" 
##and another "freq" if summarized captures are provided
salam.full <- data.frame(ch = rep(salam.ch, 2),
                         freq = c(salamander$Males, salamander$Females),
                         Sex = c(rep("male", length(salam.ch)),
                           rep("female", length(salam.ch))))
str(salam.full)
salam.full$ch <- as.character(salam.full$ch)

##delete rows with 0 freqs
salam.full.orig <- salam.full[which(salam.full$freq != 0), ]

\dontrun{
if(require(unmarked)) {
##process.data: format data, define time intervals and groups
salam.process <- process.data(salam.full.orig, model = "CJS",
                              begin.time = 1959, groups = "Sex",
                              nocc = 5)
##default dot model phi(.)p(.)
mark(salam.process, model = "CJS")
}
}
}
\keyword{datasets}

\name{newt}
\Rdversion{1.1}
\alias{newt}
\docType{data}
\title{
Newt Capture-mark-recapture Data
}
\description{
This is a capture-mark-recapture data set on adult male and female
Red-spotted Newts (\emph{Notophthalmus viridescens}) recorded by Gill
(1985).  A total of 1079 unique individuals were captured in pitfall
traps at a breeding site (White Oak Flat pond, Virginia) between 1975
and 1983.

}
\usage{data(newt)}
\format{
  A data frame with 78 observations on the following 11 variables.
  \describe{
    \item{\code{T1975}}{a binary variable, either 1 (captured) or 0 (not
      captured) during the 1975 breeding season.}
    \item{\code{T1976}}{a binary variable, either 1 (captured) or 0 (not
      captured) during the 1976 breeding season.}
    \item{\code{T1977}}{a binary variable, either 1 (captured) or 0 (not
      captured) during the 1977 breeding season.}
    \item{\code{T1978}}{a binary variable, either 1 (captured) or 0 (not
      captured) during the 1978 breeding season.}
    \item{\code{T1979}}{a binary variable, either 1 (captured) or 0 (not
      captured) during the 1979 breeding season.}
    \item{\code{T1980}}{a binary variable, either 1 (captured) or 0 (not
      captured) during the 1980 breeding season.}
    \item{\code{T1981}}{a binary variable, either 1 (captured) or 0 (not
      captured) during the 1981 breeding season.}
    \item{\code{T1982}}{a binary variable, either 1 (captured) or 0 (not
      captured) during the 1982 breeding season.}
    \item{\code{T1983}}{a binary variable, either 1 (captured) or 0 (not
      captured) during the 1983 breeding season.}
    \item{\code{Males}}{a numeric variable indicating the total number
      of males with a given capture history.} 
    \item{\code{Females}}{a numeric variable indicating the total number
      of females with a given capture history.}
  }
}
\details{
A single cohort of individuals was followed throughout the study, as all
individuals were marked in 1975 and no new individuals were added during
the subsequent years.  This data set is used to illustrate classic
Cormack-Jolly-Seber and related models (Cormack 1964, Jolly 1965, Seber
1965, Lebreton et al. 1992, Mazerolle 2015).
}
\source{
  Cormack, R. M. (1964) Estimates of survival from the sighting of
  marked animals. \emph{Biometrika} \bold{51}, 429--438.

  Gill, D. E. (1985) Interpreting breeding patterns from census data: a
  solution to the Husting dilemma. \emph{Ecology} \bold{66}, 344--354. 
  
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

  Mazerolle, M. J. (2015) Estimating detectability and biological
  parameters of interest with the use of the R
  environment. \emph{Journal of Herpetology} \bold{49}, 541--559. 
  
  Seber, G. A. F. (1965) A note on the multiple-recapture
  census. \emph{Biometrika} \bold{52}, 249--259.
}

\examples{
data(newt)
str(newt)

##convert raw capture data to capture histories
captures <- newt[, c("T1975", "T1976", "T1977", "T1978", "T1979",
                      "T1980", "T1981", "T1982", "T1983")]
newt.ch <- apply(captures, MARGIN = 1, FUN = function(i)
                  paste(i, collapse = ""))

##organize as a data frame readable by RMark package (Laake 2013)
##RMark requires at least one column called "ch" 
##and another "freq" if summarized captures are provided
newt.full <- data.frame(ch = rep(newt.ch, 2),
                        freq = c(newt$Males, newt$Females),
                        Sex = c(rep("male", length(newt.ch)),
                        rep("female", length(newt.ch))))
str(newt.full)
newt.full$ch <- as.character(newt.full$ch)

##delete rows with 0 freqs
newt.full.orig <- newt.full[which(newt.full$freq != 0), ]
}
\keyword{datasets}

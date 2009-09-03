\name{dry.frog}
\Rdversion{1.1}
\alias{dry.frog}
\docType{data}
\title{
Frog dehydration experiment on three different substrate types.
}
\description{
This is a data set modified from Mazerolle and Desrochers (2005) on the mass lost by frogs after spending two hours on one of three substrates that are encountered in some landscape types.
}
\usage{data(dry.frog)}
\format{
  A data frame with 121 observations on the following 16 variables.
  \describe{
    \item{\code{Individual}}{a numeric vector}
    \item{\code{Species}}{a factor with levels \code{Racla}}
    \item{\code{Shade}}{a numeric vector}
    \item{\code{SVL}}{a numeric vector}
    \item{\code{Substrate}}{a factor with levels \code{PEAT} \code{SOIL} \code{SPHAGNUM}}
    \item{\code{Initial_mass}}{a numeric vector}
    \item{\code{Mass_lost}}{a numeric vector}
    \item{\code{Airtemp}}{a numeric vector}
    \item{\code{Wind_cat}}{a numeric vector}
    \item{\code{Cloud}}{a numeric vector}
    \item{\code{cent_Initial_mass}}{a numeric vector}
    \item{\code{Initial_mass2}}{a numeric vector}
    \item{\code{cent_Air}}{a numeric vector}
    \item{\code{Perc.cloud}}{a numeric vector}
    \item{\code{Wind}}{a numeric vector}
    \item{\code{log_Mass_lost}}{a numeric vector}
  }
}
\details{
Note that the original analysis in Mazerolle and Desrochers (2005)
consisted of generalized estimating equations for three mass
measurements:  mass at time 0, 1 hour, and 2 hours following exposure on
the substrate. 
}
\source{
Mazerolle, M. J., Desrochers, A. (2005) Landscape resistance to frog
movements. \emph{Canadian Journal of Zoology} \bold{83}, 455--464.
}
\examples{
data(dry.frog)
## maybe str(dry.frog) ; plot(dry.frog) ...
}
\keyword{datasets}

\name{dry.frog}
\Rdversion{1.1}
\alias{dry.frog}
\docType{data}
\title{
Frog Dehydration Experiment on Three Substrate Types
}
\description{
This is a data set modified from Mazerolle and Desrochers (2005) on the
mass lost by frogs after spending two hours on one of three substrates
that are encountered in some landscape types.
}
\usage{data(dry.frog)}
\format{
  A data frame with 121 observations on the following 16 variables.
  \describe{
    \item{\code{Individual}}{a numeric identifier unique to each
      individual.} 
    \item{\code{Species}}{a factor with levels \code{Racla}.}
    \item{\code{Shade}}{a numeric vector, either 1 (shade) or 0 (no
      shade).} 
    \item{\code{SVL}}{the snout-vent length of the individual.}
    \item{\code{Substrate}}{the substrate type, a factor with levels
      \code{PEAT}, \code{SOIL}, and \code{SPHAGNUM}.}
    \item{\code{Initial_mass}}{the initial mass of individuals.}
    \item{\code{Mass_lost}}{the mass lost in g.}
    \item{\code{Airtemp}}{the air temperature in degrees C.}
    \item{\code{Wind_cat}}{the wind intensity, either 0 (no wind), 1
      (low wind), 2 (moderate wind), or 3 (strong wind).}
    \item{\code{Cloud}}{cloud cover expressed as a percentage.}
    \item{\code{cent_Initial_mass}}{centered inital mass.}
    \item{\code{Initial_mass2}}{initial mass squared.}
    \item{\code{cent_Air}}{centered air temperature.}
    \item{\code{Perc.cloud}}{proportion of cloud cover}
    \item{\code{Wind}}{wind intensity, either 1 (no or low wind) or 1
      (moderate to strong wind).}
    \item{\code{log_Mass_lost}}{log of mass lost.}
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

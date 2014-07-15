\name{tortoise}
\Rdversion{1.1}
\alias{tortoise}
\docType{data}
\title{
Gopher Tortoise distance sampling data
}
\description{
This simulated data set is based on the biological parameters for the
Gopher Tortoise (\emph{Gopherus polyphemus}) reported by Smith et
al. (2009). A half-normal distribution with a scale of 10 and without an 
adjustment factor was used to simulate the distance data for a study
area of 120 km^2. An effort of 500 m in 300 line transects was
deployed. A density of 72 individuals per km^2 was used in the
simulation using the approach outlined in Buckland et al. (2001).
}
\usage{data(tortoise)}
\format{
  A data frame with 410 observations on the following 5 variables.
  \describe{
    \item{\code{Region.Label}}{a numeric identifier for the study area.} 
    \item{\code{Area}}{a numeric variable for the surface area of the
      study area in square meters.}
    \item{\code{Sample.Label}}{a numeric identifier for each line
      transect relating each observation to its corresponding transect.}
    \item{\code{Effort}}{Effort in meters expended in each line transect.}
    \item{\code{distance}}{a numeric variable for the perpendicular
      distances in meters relative to the transect line for each of the
      individuals detected during the survey. Note that transects
      without detections have a value of NA for this variable.}
  }
}
\details{
This data set is used to illustrate classic distance sampling (Buckland
et al. 2001). 
}
\source{
  Buckland, S. T., Anderson, D. R., Burnham, K. P., Laake, J. L.,
  Borchers, D. L., Thomas, L. (2001) \emph{Introduction to distance
    sampling: estimating abundance of biological populations.} Oxford
  University Press: Oxford.

  Smith, L. L., Linehan, J. M., Stober, J. M., Elliott, M. J., Jensen,
  J. B. (2009) An evaluation of distance sampling for large-scale gopher
  tortoise surveys in Georgia, USA. \emph{Applied Herpetology} \bold{6},
  355--368.
}
\examples{
data(tortoise)
str(tortoise)

##plot distance data to determine if truncation is required
##(Buckland et al. 2001, pp. 15--17)
hist(tortoise$distance)
}
\keyword{datasets}

\name{countDist}
\Rdversion{1.1}
\alias{countDist}
\alias{countDist.default}
\alias{countDist.unmarkedFrameDS}
\alias{countDist.unmarkedFitDS}
\alias{countDist.unmarkedFrameGDS}
\alias{countDist.unmarkedFitGDS}
\alias{countDist.unmarkedFrameDSO}
\alias{countDist.unmarkedFitDSO}
\alias{print.countDist}
\title{
Compute Summary Statistics from Distance Sampling Data
}
\description{
This function extracts various summary statistics from distance sampling
data of various \code{unmarkedFrame} and \code{unmarkedFit} classes.
}
\usage{
countDist(object, plot.freq = TRUE, plot.distance = TRUE, cex.axis = 1,
          cex.lab = 1, cex.main = 1, \dots)

\method{countDist}{unmarkedFrameDS}(object, plot.freq = TRUE,
                   plot.distance = TRUE, cex.axis = 1, cex.lab = 1,
                   cex.main = 1, \dots)

\method{countDist}{unmarkedFitDS}(object, plot.freq = TRUE,
                   plot.distance = TRUE, cex.axis = 1, cex.lab = 1,
                   cex.main = 1, \dots)

\method{countDist}{unmarkedFrameGDS}(object, plot.freq = TRUE,
                   plot.distance = TRUE, cex.axis = 1,
                   cex.lab = 1, cex.main = 1, \dots)

\method{countDist}{unmarkedFitGDS}(object, plot.freq = TRUE,
                   plot.distance = TRUE, cex.axis = 1,
                   cex.lab = 1, cex.main = 1, \dots)

\method{countDist}{unmarkedFrameDSO}(object, plot.freq = TRUE,
                   plot.distance = TRUE, cex.axis = 1,
                   cex.lab = 1, cex.main = 1, plot.seasons = FALSE,
                   \dots)

\method{countDist}{unmarkedFitDSO}(object, plot.freq = TRUE,
                   plot.distance = TRUE, cex.axis = 1,
                   cex.lab = 1, cex.main = 1, plot.seasons = FALSE,
                   \dots)
}

\arguments{
  \item{object}{
    an object of various \code{unmarkedFrame} or \code{unmarkedFit}
    classes containing distance sampling data. 
  }
 \item{plot.freq}{
   logical.  Specifies if the count data (pooled across seasons and
   distance classes) should  be plotted.
 }
 \item{plot.distance}{
   logical.  Specifies if the counts in each distance class (pooled
   across seasons) should be plotted.  
 }
 \item{cex.axis}{expansion factor influencing the size of axis
   annotations on plots produced by the function.
 }
 \item{cex.lab}{expansion factor influencing the size of axis labels on
   plots produced by the function.
 }
 \item{cex.main}{expansion factor influencing the size of the main title
   above plots produced by the function.
 } 
 \item{plot.seasons}{
   logical.  Specifies if the count data should be plotted for each
   distance class and season separately.  This argument is only relevant
   for data collected across more than a single season.
 }
 \item{\dots}{
   additional arguments passed to the function.
 }
}

\details{
  This function computes a number of summary statistics in data sets
  used for the distance sampling models of Royle et al. (2004),
  Chandler et al. (2011), and distance-sampling versions of models of Dail
  and Madsen (2011) and Hostetler and Chandler (2015) based on Sollmann et
  al. (2015).
  
  \code{countDist} can take data frames of the
  \code{unmarkedFrameDS}, \code{unmarkedFrameGDS},
  \code{unmarkedFrameDSO} classes as input.  For convenience, the
  function can also extract the raw data from model objects of classes
  \code{unmarkedFitDS}, \code{unmarkedFitGDS}, and
  \code{unmarkedFitDSO}.  Note that different model objects using the
  same data set will have identical values.
}
\value{
  \code{countDist} returns a list with the following components:
  
  \item{count.table.full}{a table with the frequency of each observed
    count pooled across distances classes.}
  
  \item{count.table.seasons}{a list of tables with the frequency of
    each season-specific count pooled across distance classes.}

  \item{dist.sums.full}{a table with the frequency of counts in each
    distance class across the entire sampling seasons.}
  
  \item{dist.table.seasons}{a list of tables with the frequency of
    counts in each distance class for each primary period.}

  \item{dist.names}{a character string of labels for the distance
    classes.}

  \item{n.dist.classes}{the number of distance classes.}
  
  \item{out.freqs}{a matrix where the rows correspond to each sampling
    season and where columns consist of the number of sites sampled in
    season \eqn{t} (\code{sampled}) and the number of sites with at
    least one detection in season \eqn{t} (\code{detected}).  For
    multiseason data, the matrix includes the number of sites sampled in
    season \eqn{t - 1} with colonizations observed in season \eqn{t}
    (\code{colonized}), the number of sites sampled in season \eqn{t -
    1} with extinctions observed in season \eqn{t} (\code{extinct}), the
    number of sites sampled in season \eqn{t - 1} without changes
    observed in season \eqn{t} (\code{static}), and the number of sites
    sampled in season \eqn{t} that were also sampled in season \eqn{t -
    1} (\code{common}).
  }

  \item{out.props}{a matrix where the rows correspond to each sampling
    season and where columns consist of the proportion of sites in
    season \emph{t} with at least one detection (\code{naive.occ}).  For
    multiseason data, the matrix includes the proportion of sites
    sampled in season \eqn{t - 1} with colonizations observed in season
    \eqn{t} (\code{naive.colonization}), the proportion of sites
    sampled in season \eqn{t - 1} with extinctions observed in season
    \eqn{t} (\code{naive.extinction}), and the proportion of sites
    sampled in season \eqn{t - 1} with no changes observed in season
    \eqn{t}.
  }

  \item{n.seasons}{the number of seasons (primary periods) in the data
  set.
  }

  \item{n.visits.season}{the maximum number of visits per season in the
  data set.
}

  \item{missing.seasons}{logical vector indicating whether data were
  collected or not during a given season (primary period), where
  \code{TRUE} indicates that no data were collected during the season. 
}

}


\references{
  Chandler, R. B., Royle, J. A., King, D. I. (2011) Inference about
  density and temporary emigration in unmarked
  populations. \emph{Ecology} \bold{92}, 1429--1435.

  Dail, D., Madsen, L. (2011) Models for estimating abundance from
  repeated counts of an open population. \emph{Biometrics} \bold{67},
  577--587.

  Hostetler, J. A., Chandler, R. B. (2015) Improved state-space models
  for inference about spatial and temporal variation in abundance from
  count data. \emph{Ecology} \bold{96}, 1713--1723.

  Royle, J. A., Dawson, D. K., Bates, S. (2004) Modeling abundance
  effects in distance sampling. \emph{Ecology} \bold{85}, 1591--1597.

  Sollmann, R., Gardner, B., Chandler, R. B., Royle, J. A., Sillett,
  T. S. (2015) An open-population hierarchical distance sampling
  model. \bold{Ecology} \bold{96}, 325--331.
    
}
  
\author{
Marc J. Mazerolle
}

\seealso{
\code{\link{covDiag}}, \code{\link{detHist}}, \code{\link{detTime}},
\code{\link{countHist}}, \code{\link{Nmix.chisq}},
\code{\link{Nmix.gof.test}}
}
\examples{
##modified example from ?distsamp
\dontrun{
if(require(unmarked)){
data(linetran)
##format data
ltUMF <- with(linetran, {
              unmarkedFrameDS(y = cbind(dc1, dc2, dc3, dc4),
              siteCovs = data.frame(Length, area, habitat),
              dist.breaks = c(0, 5, 10, 15, 20),
              tlength = linetran$Length * 1000, survey = "line",
              unitsIn = "m") 
             })

##compute descriptive stats from data object
countDist(ltUMF)
     
##Half-normal detection function
fm1 <- distsamp(~ 1 ~ 1, ltUMF)
##compute descriptive stats from model object
countDist(fm1)
}
}
}
\keyword{models}

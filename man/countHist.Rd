\name{countHist}
\Rdversion{1.1}
\alias{countHist}
\alias{countHist.default}
\alias{countHist.unmarkedFramePCount}
\alias{countHist.unmarkedFitPCount}
\alias{countHist.unmarkedFrameGPC}
\alias{countHist.unmarkedFitGPC}
\alias{countHist.unmarkedFrameMPois}
\alias{countHist.unmarkedFitMPois}
\alias{countHist.unmarkedFramePCO}
\alias{countHist.unmarkedFitPCO}
\alias{countHist.unmarkedFrameGMM}
\alias{countHist.unmarkedFitGMM}
\alias{countHist.unmarkedFrameMMO}
\alias{countHist.unmarkedFitMMO}
\alias{print.countHist}
\title{
Compute Summary Statistics from Count Histories
}
\description{
This function extracts various summary statistics from count data of
various \code{unmarkedFrame} and \code{unmarkedFit} classes.
}
\usage{
countHist(object, plot.freq = TRUE, cex.axis = 1, cex.lab = 1,
          cex.main = 1, \dots) 

\method{countHist}{unmarkedFramePCount}(object, plot.freq = TRUE,
                   cex.axis = 1, cex.lab = 1, cex.main = 1, \dots)

\method{countHist}{unmarkedFitPCount}(object, plot.freq = TRUE,
                   cex.axis = 1, cex.lab = 1, cex.main = 1, \dots)

\method{countHist}{unmarkedFrameGPC}(object, plot.freq = TRUE,
                   cex.axis = 1, cex.lab = 1, cex.main = 1, \dots)

\method{countHist}{unmarkedFitGPC}(object, plot.freq = TRUE,
                   cex.axis = 1, cex.lab = 1, cex.main = 1, \dots)

\method{countHist}{unmarkedFrameMPois}(object, plot.freq = TRUE,
                   cex.axis = 1, cex.lab = 1, cex.main = 1, \dots)

\method{countHist}{unmarkedFitMPois}(object, plot.freq = TRUE,
                   cex.axis = 1, cex.lab = 1, cex.main = 1, \dots)

\method{countHist}{unmarkedFramePCO}(object, plot.freq = TRUE,
                   cex.axis = 1, cex.lab = 1, cex.main = 1,
                   plot.seasons = FALSE, \dots)

\method{countHist}{unmarkedFitPCO}(object, plot.freq = TRUE,
                   cex.axis = 1, cex.lab = 1, cex.main = 1,
                   plot.seasons = FALSE, \dots)

\method{countHist}{unmarkedFrameGMM}(object, plot.freq = TRUE,
                   cex.axis = 1, cex.lab = 1, cex.main = 1,
                   plot.seasons = FALSE, \dots)

\method{countHist}{unmarkedFitGMM}(object, plot.freq = TRUE,
                   cex.axis = 1, cex.lab = 1, cex.main = 1,
                   plot.seasons = FALSE, \dots)

\method{countHist}{unmarkedFrameMMO}(object, plot.freq = TRUE,
                   cex.axis = 1, cex.lab = 1, cex.main = 1,
                   plot.seasons = FALSE, \dots)

\method{countHist}{unmarkedFitMMO}(object, plot.freq = TRUE,
                   cex.axis = 1, cex.lab = 1, cex.main = 1,
                   plot.seasons = FALSE, \dots)

}

\arguments{
  \item{object}{
    an object of various \code{unmarkedFrame} or \code{unmarkedFit}
    classes containing count history data. 
  }
 \item{plot.freq}{
   logical.  Specifies if the count data (pooled across seasons) should
   be plotted.
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
   season separately.  This argument is only relevant for data collected
   across more than a single season.
 }
 \item{\dots}{
   additional arguments passed to the function.
 }
}

\details{
  This function computes a number of summary statistics in data sets
  used for various \emph{N}-mixture models including those of Royle
  (2004a, b), Dail and Madsen (2011), and Chandler et al. (2011).
  
  \code{countHist} can take data frames of the
  \code{unmarkedFramePCount}, \code{unmarkedFrameGPC},
  \code{unmarkedFrameMPois}, \code{unmarkedFramePCO},
  \code{unmarkedFrameGMM}, \code{unmarkedFrameMMO} classes as input.
  For convenience, the function can also extract the raw data from model
  objects of classes \code{unmarkedFitPCount}, \code{unmarkedFitGPC},
  \code{unmarkedFitMPois}, \code{unmarkedFitPCO}, \code{unmarkedFitGMM},
  and \code{unmarkedFitMMO}. Note that different model objects using the
  same data set will have identical values.
}
\value{
  \code{countHist} returns a list with the following components:
  
  \item{count.table.full}{a table with the frequency of each observed
    count.}
  
  \item{count.table.seasons}{a list of tables with the frequency of
    each season-specific count.}

  \item{hist.table.full}{a table with the frequency of each count
    history across the entire sampling period.}
  
  \item{hist.table.seasons}{a list of tables with the frequency of
    each count history for each primary period (season).}
  
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

  \item{n.seasons}{the number of seasons (primary periods) in the data set.
  }

  \item{n.visits.season}{the maximum number of visits per season in the data set.
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
  
  Royle, J. A. (2004a) \emph{N}-mixture models for estimating population
  size from spatially replicated counts. \emph{Biometrics} \bold{60},
  108--115.
  
  Royle, J. A. (2004b) Generalized estimators of avian abundance from
  count survey data. \emph{Animal Biodiversity and Conservation}
  \bold{27}, 375--386. 
}
  
\author{
Marc J. Mazerolle
}

\seealso{
\code{\link{covDiag}}, \code{\link{detHist}}, \code{\link{detTime}},
\code{\link{countDist}}, \code{\link{Nmix.chisq}},
\code{\link{Nmix.gof.test}}
}
\examples{
##modified example from ?pcount
\dontrun{
if(require(unmarked)){
data(mallard)
mallardUMF <- unmarkedFramePCount(mallard.y, siteCovs = mallard.site,
                                  obsCovs = mallard.obs)
##compute descriptive stats from data object
countHist(mallardUMF)

##run single season model
fm.mallard <- pcount(~ ivel+ date + I(date^2) ~ length + elev +
                     forest, mallardUMF, K=30)
##compute descriptive stats from model object
countHist(fm.mallard)
}
}
}
\keyword{models}

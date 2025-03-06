\name{detHist}
\Rdversion{1.1}
\alias{detHist}
\alias{detHist.default}
\alias{detHist.unmarkedFitColExt}
\alias{detHist.unmarkedFitOccu}
\alias{detHist.unmarkedFitOccuFP}
\alias{detHist.unmarkedFitOccuRN}
\alias{detHist.unmarkedFitOccuMulti}
\alias{detHist.unmarkedFitOccuMS}
\alias{detHist.unmarkedFrameOccu}
\alias{detHist.unmarkedFrameOccuFP}
\alias{detHist.unmarkedFrameOccuMulti}
\alias{detHist.unmarkedMultFrame}
\alias{detHist.unmarkedFrameOccuMS}
\alias{detHist.unmarkedFitGOccu}
\alias{detHist.unmarkedFrameOccuComm}
\alias{detHist.unmarkedFitOccuComm}
\alias{print.detHist}
\title{
Compute Summary Statistics from Detection Histories
}
\description{
This function extracts various summary statistics from detection history
data of various \code{unmarkedFrame} and \code{unmarkedFit} classes.
}
\usage{
detHist(object, \dots)

\method{detHist}{unmarkedFitColExt}(object, \dots)

\method{detHist}{unmarkedFitOccu}(object, \dots)

\method{detHist}{unmarkedFitOccuFP}(object, \dots)

\method{detHist}{unmarkedFitOccuRN}(object, \dots)

\method{detHist}{unmarkedFitOccuMulti}(object, \dots)

\method{detHist}{unmarkedFitOccuMS}(object, \dots)

\method{detHist}{unmarkedFrameOccu}(object, \dots)

\method{detHist}{unmarkedFrameOccuFP}(object, \dots)

\method{detHist}{unmarkedMultFrame}(object, \dots)

\method{detHist}{unmarkedFrameOccuMulti}(object, \dots)

\method{detHist}{unmarkedFrameOccuMS}(object, \dots)

\method{detHist}{unmarkedFitGOccu}(object, \dots)

\method{detHist}{unmarkedFrameOccuComm}(object, \dots)

\method{detHist}{unmarkedFitOccuComm}(object, \dots)

}

\arguments{
  \item{object}{
    an object of various \code{unmarkedFrame} or \code{unmarkedFit}
    classes containing detection history data. 
  }

  \item{\dots}{
    additional arguments passed to the function.
  }
}

\details{
  This function computes a number of summary statistics in data sets
  used for single-season occupancy models (MacKenzie et al. 2002),
  dynamic occupancy models (MacKenzie et al. 2003), Royle-Nichols models
  (Royle and Nichols 2003), false-positive occupancy models (Royle
  and Link 2006, Miller et al. 2011), multispecies occupancy models
  (Rota et al. 2016), and multistate occupancy models (Nichols et
  al. 2007, MacKenzie et al. 2009).
  
  \code{detHist} can take data frames of the \code{unmarkedFrameOccu},
  \code{unmarkedFrameOccuFP}, \code{unmarkedMultFrame},
  \code{unmarkedFrameOccuMulti}, \code{unmarkedFrameOccuMS} classes as
  input.  For convenience, the  function can also extract the raw data
  from model objects of classes \code{unmarkedFitColExt},
  \code{unmarkedFitOccu},  \code{unmarkedFitOccuFP},
  \code{unmarkedFitOccuRN}, \code{unmarkedFitOccuMulti},
  \code{unmarkedFitOccuMS}, \code{unmarkedFitGOccu}, and
  \code{unmarkedFitOccuComm}.  Note that different model objects using
  the same data set will have identical values. 
}
\value{
  For objects of classes \code{unmarkedFitOccu}, \code{unmarkedFitOccuRN},
  \code{unmarkedFitOccuFP}, \code{unmarkedFitColExt},
  \code{unmarkedFitOccuMS}, \code{unmarkedFitGOccu},
  \code{unmarkedFitOccuComm}, \code{unmarkedFrameOccu},
  \code{unmarkedFrameOccuFP}, \code{unmarkedMultFrame},
  \code{unmarkedFrameOccuMS},  \code{unmarkedFrameOccuComm},
  \code{detHist} returns a list with the following components: 
  
  \item{hist.table.full}{a table with the frequency of each observed
    detection history. For \code{unmarkedFitGOccu} objects, this element
    consists in a summary of the detections pooled across devices
    (secondary periods) within each primary period. This element is
    \code{NULL} for \code{unmarkedFrameOccuComm} and
    \code{unmarkedFitOccuComm} objects.}
  
  \item{hist.table.seasons}{a list of tables with the frequency of
    each season-specific detection history.}

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
    1} (\code{common}).  For multispecies data, \code{out.freqs}
    presents, for each species, the number of sites sampled and the number
    of sites with at least one detection.
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
    \eqn{t}. For multispecies data, \code{out.props} presents the
    proportion of sites with a least one detection for each species.
  }

  \item{n.seasons}{the number of seasons (primary periods) in the data set.
  }

  \item{n.visits.season}{the maximum number of visits per season in the data set.
  }

  \item{n.species}{the number of species in the data set.
  }

  \item{missing.seasons}{logical vector indicating whether data were
  collected or not during a given season (primary period), where
  \code{TRUE} indicates that no data were collected during the season. 
}


  For objects of classes \code{unmarkedFitOccuMulti} and
  \code{unmarkedFrameOccuMulti}, \code{detHist} returns a list with the
  following components: 
  
  \item{hist.table.full}{a table with the frequency of each observed
    detection history.  The species are coded with letters and follow
    the same order of presentation as in the other parts of the output.}
  
  \item{hist.table.species}{a list of tables with the frequency of
    each species-specific detection history.  The last element of
    \code{hist.table.species} features the number of sites with
    co-occurrence of the different species (\code{coOcc}).} 

  \item{out.freqs}{a matrix where the rows correspond to each species
    and where columns consist of the number of sites sampled during the
    season (\code{sampled}) and the number of sites with at least one
    detection (\code{detected}).}

  \item{out.props}{a matrix where the rows correspond to each species
    and where columns consist of the proportion of sites with at least
    one detection during the season (\code{naive.occ}).
  }

  \item{n.seasons}{the number of seasons (primary periods) in the data set.
  }

  \item{n.visits.season}{the maximum number of visits per season in the data set.
  }

  \item{n.species}{the number of species in the data set.
  }

  \item{missing.seasons}{logical vector indicating whether data were
  collected or not during a given season (primary period), where
  \code{TRUE} indicates that no data were collected during the season. 
}
}


\references{
  MacKenzie, D. I., Nichols, J. D., Lachman, G. B., Droege, S., Royle, 
  J. A., Langtimm, C. A. (2002) Estimating site occupancy rates when 
  detection probabilities are less than one. \emph{Ecology} \bold{83},
  2248--2255.

  MacKenzie, D. I., Nichols, J. D., Hines, J. E., Knutson, M. G.,
  Franklin, A. B. (2003) Estimating site occupancy, colonization, and
  local extinction when a species is detected imperfectly. \emph{Ecology}
  \bold{84}, 2200--2207.  

  MacKenzie, D. I., Nichols, J. D., Seamans, M. E., Gutierrez,
  R. J. (2009) Modeling species occurrence dynamics with multiple states
  and imperfect detection. \emph{Ecology} \bold{90}, 823--835.

  Mazerolle, M. J. (2015) Estimating detectability and biological
  parameters of interest with the use of the R
  environment. \emph{Journal of Herpetology} \bold{49}, 541--559.
  
  Miller, D. A. W., Nichols, J. D., McClintock, B. T., Campbell
  Grant, E. H., Bailey, L. L. (2011) Improving occupancy estimation when
  two types of observational error occur: non-detection and species
  misidentification. \emph{Ecology} \bold{92}, 1422--1428.

  Nichols, J. D., Hines, J. E., Mackenzie, D. I., Seamans, M. E.,
  Gutierrez, R. J. (2007) Occupancy estimation and modeling with
  multiple states and state uncertainty. \emph{Ecology} \bold{88},
  1395--1400.

  Rota, C. T., Ferreira, M. A. R., Kays, R. W., Forrester, T. D.,
  Kalies, E. L., McShea, W. J., Parsons, A. W., Millspaugh, J. J. (2016)
  A multispecies occupancy model for two or more interacting
  species. \emph{Methods in Ecology and Evolution} \bold{7},
  1164--1173.

  Royle, J. A., Link, W. A. (2006) Generalized site occupancy models
  allowing for false positive and false negative errors. \emph{Ecology}
  \bold{87}, 835--841.
  
  Royle, J. A., Nichols, J. D. (2003) Estimating abundance from
  repeated presence-absence data or point counts. \emph{Ecology}
  \bold{84}, 777--790.

}
  
\author{
Marc J. Mazerolle
}

\seealso{
\code{\link{covDiag}}, \code{\link{countHist}}, \code{\link{countDist}},
\code{\link{detTime}}, \code{\link{mb.chisq}}, \code{\link{mb.gof.test}}
}
\examples{
##data from Mazerolle (2015)
\dontrun{
data(bullfrog)

##detection data
detections <- bullfrog[, 3:9]

##load unmarked package
if(require(unmarked)){

##assemble in unmarkedFrameOccu
bfrog <- unmarkedFrameOccu(y = detections)

##compute descriptive stats from data object
detHist(bfrog)

##run model
fm <- occu(~ 1 ~ 1, data = bfrog)
##compute descriptive stats from model object
detHist(fm)
}
}
}
\keyword{models}

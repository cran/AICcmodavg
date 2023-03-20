\name{detTime}
\Rdversion{1.1}
\alias{detTime}
\alias{detTime.default}
\alias{detTime.unmarkedFrameOccuTTD}
\alias{detTime.unmarkedFitOccuTTD}
\alias{print.detTime}
\title{
Compute Summary Statistics from Time to Detection Data
}
\description{
This function extracts various summary statistics from time to detection
data of various \code{unmarkedFrame} and \code{unmarkedFit} classes.
}
\usage{
detTime(object, plot.time = TRUE, plot.seasons = FALSE,
        cex.axis = 1, cex.lab = 1, cex.main = 1, \dots)

\method{detTime}{unmarkedFrameOccuTTD}(object, plot.time = TRUE,
                 plot.seasons = FALSE, cex.axis = 1, cex.lab = 1,
                 cex.main = 1, \dots)

\method{detTime}{unmarkedFitOccuTTD}(object, plot.time = TRUE,
                   plot.seasons = FALSE, cex.axis = 1, cex.lab = 1,
                   cex.main = 1, \dots)
}

\arguments{
  \item{object}{
    an object of various \code{unmarkedFrame} or \code{unmarkedFit}
    classes containing time to detection data. 
  }
 \item{plot.time}{
   logical.  Specifies if the time to detection data (pooled across
   seasons) should  be plotted.
 }
 \item{plot.seasons}{
   logical.  Specifies if the time to detection data should be plotted
   for each season separately.  This argument is only relevant
   for data collected across more than a single season.
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
\item{\dots}{
  additional arguments passed to the function.
}
}

\details{
  This function computes a number of summary statistics in data sets
  used for the time to detection models of Garrard et al. (2008,
  2013).
  
  \code{detTime} can take data frames of the \code{unmarkedFrameOccuTTD}
  class as input, or can also extract the raw data from model objects of
  the \code{unmarkedFitOccuTTD} class.  Note that different model
  objects using the same data set will have identical values.
}
\value{
  \code{detTime} returns a list with the following components:
  
  \item{time.table.full}{a table with the quantiles of time to
    detection data pooled across seasons, but excluding censored
    observations.}
  
  \item{time.table.seasons}{a list of tables with the quantiles of
    season-specific time to detection data, but excluding censored
    observations.}

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
  Garrard, G. E., Bekessy, S. A., McCarthy, M. A., Wintle, B. A. (2008)
  When have we looked hard enough? A novel method for setting minimum
  survey effort protocols for flora surveys. \emph{Austral Ecology}
  \bold{33}, 986--998.
  
  Garrard, G. E., McCarthy, M. A., Williams, N. S., Bekessy, S. A.,
  Wintle, B. A. (2013) A general model of detectability using species
  traits. \emph{Methods in Ecology and Evolution} \bold{4}, 45--52.
}
  
\author{
Marc J. Mazerolle
}

\seealso{
\code{\link{countDist}}, \code{\link{countHist}}, \code{\link{detHist}},
\code{\link{Nmix.chisq}}, \code{\link{Nmix.gof.test}}
}
\examples{
##example from ?occuTTD
\dontrun{
if(require(unmarked)){
N <- 500; J <- 1
     
##Simulate occupancy
scovs <- data.frame(elev=c(scale(runif(N, 0,100))),
                    forest=runif(N,0,1),
                    wind=runif(N,0,1))
     
beta_psi <- c(-0.69, 0.71, -0.5)
psi <- plogis(cbind(1, scovs$elev, scovs$forest) %*% beta_psi)
z <- rbinom(N, 1, psi)
     
##Simulate detection
Tmax <- 10 #Same survey length for all observations
beta_lam <- c(-2, -0.2, 0.7)
rate <- exp(cbind(1, scovs$elev, scovs$wind) %*% beta_lam)
ttd <- rexp(N, rate)
ttd[z==0] <- Tmax #Censor unoccupied sites
ttd[ttd>Tmax] <- Tmax #Censor when ttd was greater than survey length
     
##Build unmarkedFrame
umf <- unmarkedFrameOccuTTD(y=ttd, surveyLength=Tmax, siteCovs=scovs)

##compute descriptive stats from data object
detTime(umf)
     
##Fit model
fit.occuTTD <- occuTTD(psiformula=~elev+forest, detformula=~elev+wind, data=umf)
##extract info from model object
detTime(fit.occuTTD)     
     
}
}
}
\keyword{models}

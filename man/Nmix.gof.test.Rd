\name{Nmix.gof.test}
\Rdversion{1.1}
\alias{Nmix.chisq}
\alias{Nmix.chisq.default}
\alias{Nmix.chisq.unmarkedFitPCount}
\alias{Nmix.chisq.unmarkedFitPCO}
\alias{Nmix.chisq.unmarkedFitMPois}
\alias{Nmix.chisq.unmarkedFitDS}
\alias{Nmix.chisq.unmarkedFitGDS}
\alias{Nmix.chisq.unmarkedFitGPC}
\alias{Nmix.chisq.unmarkedFitGMM}
\alias{Nmix.chisq.unmarkedFitMMO}
\alias{Nmix.chisq.unmarkedFitDSO}
\alias{Nmix.gof.test}
\alias{Nmix.gof.test.default}
\alias{Nmix.gof.test.unmarkedFitPCount}
\alias{Nmix.gof.test.unmarkedFitPCO}
\alias{Nmix.gof.test.unmarkedFitDS}
\alias{Nmix.gof.test.unmarkedFitGDS}
\alias{Nmix.gof.test.unmarkedFitGMM}
\alias{Nmix.gof.test.unmarkedFitGPC}
\alias{Nmix.gof.test.unmarkedFitMPois}
\alias{Nmix.gof.test.unmarkedFitMMO}
\alias{Nmix.gof.test.unmarkedFitDSO}
\alias{print.Nmix.chisq}
\title{
Compute Chi-square Goodness-of-fit Test for N-mixture Models
}
\description{
These functions compute a goodness-of-fit test for N-mixture models
based on Pearson's chi-square.
}
\usage{
##methods for 'unmarkedFitPCount', 'unmarkedFitPCO', 
##'unmarkedFitDS', 'unmarkedFitGDS', 'unmarkedFitGMM',
##'unmarkedFitGPC', and 'unmarkedFitMPois' classes
Nmix.chisq(mod, \dots)

Nmix.gof.test(mod, nsim = 5, plot.hist = TRUE, report = NULL,
parallel = TRUE, ncores, cex.axis = 1, cex.lab = 1, cex.main = 1,
lwd = 1, \dots)

}
\arguments{
  \item{mod}{
    the \emph{N}-mixture model of \code{unmarkedFitPCount},
    \code{unmarkedFitPCO}, \code{unmarkedFitDS}, \code{unmarkedFitGDS},
    \code{unmarkedFitGMM}, \code{unmarkedFitGPC},
    \code{unmarkedFitMPois}, \code{unmarkedFitMMO}, or
    \code{unmarkedFitDSO} classes for which a goodness-of-fit test is
    required. 
  }
\item{nsim}{the number of bootstrapped samples.
}
\item{plot.hist}{logical.  Specifies that a histogram of the
  bootstrapped test statistic is to be included in the output.
}
\item{report}{If \code{NULL}, the test statistic for each iteration is
  not printed in the terminal.  Otherwise, an integer indicating the
  number of values of the test statistic that should be printed on
  the same line.  For example, if \code{report = 3}, the values of the
  test statistic for three iterations are reported on each line.
}
\item{parallel}{logical.  If \code{TRUE}, requests that \code{parboot}
  use multiple cores to accelerate computations of the bootstrap.
}
\item{ncores}{integer indicating the number of cores to use when
  bootstrapping in parallel during the analysis of simulated data sets.
  If \code{ncores} is not specified, one less than the number of
  available cores on the computer is used.
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
\item{lwd}{expansion factor of line width on plots produced by the
  function.
}
\item{\dots}{
  additional arguments passed to the function.
}
}
\details{
The Pearson chi-square can be used to assess the fit of N-mixture
models.  Instead of relying on the theoretical distribution of the
chi-square, a parametric bootstrap approach is implemented to obtain
\emph{P}-values with the \code{parboot} function of the \code{unmarked}
package.  \code{Nmix.chisq} computes the observed chi-square statistic
based on the observed and expected counts from the model.
\code{Nmix.gof.test} calls internally \code{Nmix.chisq} and
\code{parboot} to generate simulated data sets based on the model and
compute the chi-square test statistic.

It is also possible to obtain an estimate of the overdispersion
parameter (c-hat) for the model at hand by dividing the observed
chi-square statistic by the mean of the statistics obtained from
simulation (MacKenzie and Bailey 2004, McKenny et al. 2006).  This
method of estimating c-hat is similar to the one implemented for
capture-mark-recapture models in program MARK (White and Burnham 1999).

Note that values of c-hat > 1 indicate overdispersion (variance >
mean).  Values much higher than 1 (i.e., > 4) probably indicate
lack-of-fit.  In cases of moderate overdispersion, one can multiply the
variance-covariance matrix of the estimates by c-hat.  As a result, the
SE's of the estimates are inflated (c-hat is also known as a variance
inflation factor).

In model selection, c-hat should be estimated from the global model and
the same value of c-hat applied to the entire model set.  Specifically,
a global model is the most complex model which can be simplified to
yield all the other (nested) models of the set.  When no single global
model exists in the set of models considered, such as when sample size
does not allow a complex model, one can estimate c-hat from 'subglobal'
models.  Here, 'subglobal' models denote models from which only a subset
of the models of the candidate set can be derived.  In such cases, one
can use the smallest value of c-hat for model selection (Burnham and
Anderson 2002).

Note that c-hat counts as an additional parameter estimated and should 
be added to \emph{K}.  All functions in package \code{AICcmodavg}
automatically add 1 when the \code{c.hat} argument > 1 and apply the
same value of c-hat for the entire model set.  When c-hat > 1, functions
compute quasi-likelihood information criteria (either QAICc or QAIC,
depending on the value of the \code{second.ord} argument) by scaling the
log-likelihood of the model by c-hat.  The value of c-hat can influence
the ranking of the models:  as c-hat increases, QAIC or QAICc will favor
models with fewer parameters.  As an additional check against this
potential problem, one can generate several model selection tables by
incrementing values of c-hat to assess the model selection uncertainty.
If ranking changes only slightly up to the c-hat value observed, one can
be confident in making inference.

In cases of underdispersion (c-hat < 1), it is recommended to keep the
value of c-hat to 1.  However, note that values of c-hat << 1 can also
indicate lack-of-fit and that an alternative model should be investigated. 

}
\value{
  \code{Nmix.chisq} returns two value:
\item{chi.square}{the Pearson chi-square statistic.
}
\item{model.type}{the class of the fitted model.
}


\code{Nmix.gof.test} returns the following components:

\item{model.type}{the class of the fitted model.
}
\item{chi.square}{the Pearson chi-square statistic.
}
\item{t.star}{the bootstrapped chi-square test statistics (i.e.,
  obtained for each of the simulated data sets).
}
\item{p.value}{the \emph{P}-value assessed from the parametric bootstrap,
  computed as the proportion of the simulated test statistics greater
  than or equal to the observed test statistic.
}
\item{c.hat.est}{the estimate of the overdispersion parameter, c-hat,
  computed as the observed test statistic divided by the mean of the
  simulated test statistics.
}
\item{nsim}{the number of bootstrap samples. The recommended number of
  samples varies with the data set, but should be on the order of 1000
  or 5000, and in cases with a large number of visits, even 10 000
  samples, namely to reduce the effect of unusually small values of the
  test statistics.
}
}
\references{
  Burnham, K. P., Anderson, D. R. (2002) \emph{Model Selection and
    Multimodel Inference: a practical information-theoretic
    approach}. Second edition. Springer: New York.

  MacKenzie, D. I., Bailey, L. L. (2004) Assessing the fit of
  site-occupancy models. \emph{Journal of Agricultural, Biological, and
    Environmental Statistics} \bold{9}, 300--318.

  McKenny, H. C., Keeton, W. S., Donovan, T. M. (2006). Effects of
  structural complexity enhancement on eastern red-backed salamander
  (\emph{Plethodon cinereus}) populations in northern hardwood
  forests. \emph{Forest Ecology and Management} \bold{230}, 186--196.
  
  White, G. C., Burnham, K. P. (1999). Program MARK: Survival estimation
  from populations of marked animals. \emph{Bird Study} \bold{46
    (Supplement)}, 120--138. 
}
\author{
Marc J. Mazerolle
}

\seealso{
\code{\link{AICc}}, \code{\link{c_hat}}, \code{\link{evidence}},
\code{\link{modavg}}, \code{\link{importance}},
\code{\link{mb.gof.test}}, \code{\link{modavgPred}},
\code{\link[unmarked]{pcount}}, \code{\link[unmarked]{pcountOpen}},
\code{\link[unmarked]{parboot}}
}

\examples{
##N-mixture model example modified from ?pcount
\dontrun{
require(unmarked)
##single season
data(mallard)
mallardUMF <- unmarkedFramePCount(mallard.y, siteCovs = mallard.site,
                                  obsCovs = mallard.obs)
##run model
fm.mallard <- pcount(~ ivel+ date + I(date^2) ~ length + elev + forest,
                     mallardUMF, K=30)

##compute observed chi-square
obs <- Nmix.chisq(fm.mallard)
obs

##round to 4 digits after decimal point
print(obs, digits.vals = 4)

##compute observed chi-square, assess significance, and estimate c-hat
obs.boot <- Nmix.gof.test(fm.mallard, nsim = 10)
##note that more bootstrap samples are recommended
##(e.g., 1000, 5000, or 10 000)
obs.boot
print(obs.boot, digits.vals = 4, digits.chisq = 4)
detach(package:unmarked)
}
}

\keyword{models}

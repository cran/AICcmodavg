\name{mb.gof.test}
\Rdversion{1.1}
\alias{mb.gof.test}
\alias{mb.chisq}
\alias{print.mb.chisq}
\title{
Compute MacKenzie and Bailey Goodness-of-fit Test for Single Season Occupancy Models
}
\description{
These functions compute the MacKenzie and Bailey (2004) goodness-of-fit
test for single season occupancy models based on Pearson's chi-square.
}
\usage{
mb.chisq(mod, print.table = TRUE)
mb.gof.test(mod, nsim = 5, plot.hist = TRUE)
}
\arguments{
  \item{mod}{
the model for which a goodness-of-fit test is required.
}
\item{print.table}{logical.  Specifies if the detailed table of observed
  and expected values is to be included in the output.
}
\item{nsim}{the number of bootstrapped samples.
}
\item{plot.hist}{logical.  Specifies that a histogram of the
  bootstrapped test statistic is to be included in the output.
}
}
\details{
MacKenzie and Bailey (2004) and MacKenzie et al. (2006) suggest using
the Pearson chi-square to assess the fit of single season occupancy
models.  Given low expected frequencies, the chi-square statistic will
deviate from the theoretical distribution and it is recommended to use a
parametric bootstrap approach to obtain P-values with the \code{parboot}
function of the \code{unmarked} package.  \code{mb.chisq} computes the
table of observed and expected values based on the detection histories
and single season occupancy model used. \code{mb.gof.test} calls
internally \code{mb.chisq} and \code{parboot} to generate simulated data
sets based on the model and compute the MacKenzie and Bailey test
statistic.  Missing values are accomodated by creating cohorts for each
pattern of missing values.

It is also possible to obtain an estimate of the overdispersion
parameter (c-hat) for the model at hand by dividing the observed
chi-square statistic by the mean of the statistics obtained from
simulation. 

Note that values of c-hat > 1 indicate overdispersion (variance > mean),
but that values much higher than 1 (i.e., > 4) probably indicate
lack-of-fit.  In cases of moderate overdispersion, one usually
multiplies the variance-covariance matrix of the estimates by c-hat.  As
a result, the SE's of the estimates are inflated (c-hat is also known as
a variance inflation factor). 

In model selection, c-hat should be estimated from the global model and
the same value of c-hat applied to the entire model set.  Specifically,
a global model is the most complex model from which all the other models
of the set are simpler versions (nested).  When no single global model
exists in the set of models considered, such as when sample size does
not allow a complex model, one can estimate c-hat from 'subglobal'
models.  Here, 'subglobal' models denote models from which only a subset
of the models of the candidate set can be derived.  In such cases, one
can use the smallest value of c-hat for model selection (Burnham and
Anderson 2002).

Note that c-hat counts as an additional parameter estimated and should
be added to K.  All functions in package \code{AICcmodavg} automatically
add 1 when the \code{c.hat} argument > 1 and apply the same value of
c-hat for the entire model set.  When c-hat > 1, functions compute
quasi-likelihood information criteria (either QAICc or QAIC, depending
on the value of the \code{second.ord} argument) by scaling the
log-likelihood of the model by c-hat.  The value of c-hat can influence
the ranking of the models:  as c-hat increases, QAIC or QAICc will favor
models with fewer parameters.  As an additional check against this
potential problem, one can generate several model selection tables by
incrementing values of c-hat to assess the model selection uncertainty.
If ranking changes little up to the c-hat value observed, one can be
confident in making inference.

In cases of underdispersion (c-hat < 1), it is recommended to keep the
value of c-hat to 1.  However, note that values of c-hat << 1 can also
indicate lack-of-fit and that an alternative model should be investigated. 

}
\value{
  \code{mb.chisq} returns the following components:
  
  \item{chisq.table}{the table of observed and expected values for each
  detection history and its chi-square component (if \code{print.table = TRUE}).
}
\item{chi.square}{the Pearson chi-square statistic.
}


  \code{mb.gof.test} returns the following components:

  \item{chisq.table}{the table of observed and expected values for each
  detection history and its chi-square component (if \code{print.table = TRUE}).
}
\item{chi.square}{the Pearson chi-square statistic.
}
\item{t.star}{the bootstrapped chi-square test statistics (i.e.,
  obtained for each of the simulated data sets).
}
\item{p.value}{the P-value assessed from the parametric bootstrap,
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

MacKenzie, D. I., Nichols, J. D., Royle, J. A., Pollock, K. H., Bailey,
  L. L., Hines, J. E. (2006) \emph{Occupancy estimation and
  modeling: inferring patterns and dynamics of species
  occurrence}. Academic Press: New York.
}
\author{
Marc J. Mazerolle
}

\seealso{
\code{\link{AICc}}, \code{\link{c_hat}}, \code{\link{evidence}},
\code{\link{modavg}}, \code{\link{importance}},
\code{\link{modavgPred}}, \code{\link{Nmix.gof.test}},
\code{\link[unmarked]{occu}}, \code{\link[unmarked]{parboot}} 
}

\examples{
\dontrun{
##single-season occupancy model example modified from ?occu
if(require(unmarked)) {
##single season
data(frogs)
pferUMF <- unmarkedFrameOccu(pfer.bin)
## add some fake covariates for illustration
siteCovs(pferUMF) <- data.frame(sitevar1 = rnorm(numSites(pferUMF)),
                                sitevar2 = rnorm(numSites(pferUMF))) 
     
## observation covariates are in site-major, observation-minor order
obsCovs(pferUMF) <- data.frame(obsvar1 = rnorm(numSites(pferUMF) *
                                 obsNum(pferUMF))) 

##run model
fm1 <- occu(~ obsvar1 ~ sitevar1, pferUMF)

##compute observed chi-square
obs <- mb.chisq(fm1)
obs
##round to 4 digits after decimal point
print(obs, digits.vals = 4)
}

\dontrun{
##compute observed chi-square, assess significance, and estimate c-hat
obs.boot <- mb.gof.test(fm1, nsim = 3)
##note that more bootstrap samples are recommended
##(e.g., 1000, 5000, or 10 000)
obs.boot
print(obs.boot, digits.vals = 4, digits.chisq = 4)
}


\dontrun{
##data with missing values
mat1 <- matrix(c(0, 0, 0), nrow = 120, ncol = 3, byrow = TRUE)
mat2 <- matrix(c(0, 0, 1), nrow = 23, ncol = 3, byrow = TRUE)
mat3 <- matrix(c(1, NA, NA), nrow = 42, ncol = 3, byrow = TRUE)
mat4 <- matrix(c(0, 1, NA), nrow = 33, ncol = 3, byrow = TRUE)
y.mat <- rbind(mat1, mat2, mat3, mat4)
y.sim.data <- unmarkedFrameOccu(y = y.mat)
m1 <- occu(~ 1 ~ 1, data = y.sim.data)

mb.gof.test(m1, nsim = 3)
##note that more bootstrap samples are recommended
##(e.g., 1000, 5000, or 10 000) 
}

detach(package:unmarked)
}
}

\keyword{models}

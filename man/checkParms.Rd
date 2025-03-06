\name{checkParms}
\Rdversion{1.1}
\encoding{utf-8}
\alias{checkParms}
\alias{checkParms.default}
\alias{checkParms.betareg}
\alias{checkParms.clm}
\alias{checkParms.clmm}
\alias{checkParms.coxme}
\alias{checkParms.coxph}
\alias{checkParms.glm}
\alias{checkParms.glmmTMB}
\alias{checkParms.gls}
\alias{checkParms.gnls}
\alias{checkParms.hurdle}
\alias{checkParms.lm}
\alias{checkParms.lme}
\alias{checkParms.lmekin}
\alias{checkParms.maxlikeFit}
\alias{checkParms.merMod}
\alias{checkParms.lmerModLmerTest}
\alias{checkParms.multinom}
\alias{checkParms.nlme}
\alias{checkParms.nls}
\alias{checkParms.polr}
\alias{checkParms.rlm}
\alias{checkParms.survreg}
\alias{checkParms.unmarkedFit}
\alias{checkParms.vglm}
\alias{checkParms.zeroinfl}
\alias{print.checkParms}
\title{
Identify Parameters with Large Standard Errors
}
\description{
This function identifies parameter estimates with large standard errors
in a model.  It is particularly useful for complex models with different
parameter types such as those of \code{unmarkedFit} classes implemented
in package \code{unmarked} (Fiske and Chandler, 2011), as well as other
types of regression models.
}
\usage{
checkParms(mod, se.max = 25, simplify = TRUE, \dots)

\method{checkParms}{betareg}(mod, se.max = 25, \dots)

\method{checkParms}{clm}(mod, se.max = 25, \dots)

\method{checkParms}{clmm}(mod, se.max = 25, \dots)

\method{checkParms}{coxme}(mod, se.max = 25, \dots)

\method{checkParms}{coxph}(mod, se.max = 25, \dots)

\method{checkParms}{glm}(mod, se.max = 25, \dots)

\method{checkParms}{glmmTMB}(mod, se.max = 25, \dots)

\method{checkParms}{gls}(mod, se.max = 25, \dots)

\method{checkParms}{gnls}(mod, se.max = 25, \dots)

\method{checkParms}{hurdle}(mod, se.max = 25, \dots)

\method{checkParms}{lm}(mod, se.max = 25, \dots)

\method{checkParms}{lme}(mod, se.max = 25, \dots)

\method{checkParms}{lmekin}(mod, se.max = 25, \dots)

\method{checkParms}{maxlikeFit}(mod, se.max = 25, \dots)

\method{checkParms}{merMod}(mod, se.max = 25, \dots)

\method{checkParms}{lmerModLmerTest}(mod, se.max = 25, \dots)

\method{checkParms}{multinom}(mod, se.max = 25, \dots)

\method{checkParms}{nlme}(mod, se.max = 25, \dots)

\method{checkParms}{nls}(mod, se.max = 25, \dots)

\method{checkParms}{polr}(mod, se.max = 25, \dots)

\method{checkParms}{rlm}(mod, se.max = 25, \dots)

\method{checkParms}{survreg}(mod, se.max = 25, \dots)

\method{checkParms}{unmarkedFit}(mod, se.max = 25, simplify = TRUE,
\dots)

\method{checkParms}{vglm}(mod, se.max = 25, \dots)

\method{checkParms}{zeroinfl}(mod, se.max = 25, \dots)

}
\arguments{
  \item{mod}{
    a model of \code{unmarkedFit} classes or other regression model.
    This model is checked to determine the occurrence of large standard
    errors for parameter estimates.
}
\item{se.max}{
  specifies the value beyond which standard errors are deemed high for
  the model at hand.  The function will determine the number of
  estimates with standard errors that exceed \code{se.max}.
}
\item{simplify}{
  this argument is only valid for models of \code{unmarkedFit}
  classes which consist of several parameter types for detection
  probability and demographic parameters (e.g., abundance, occupancy,
  extinction).  If \code{TRUE}, the function returns a matrix with a
  single row identifying the parameter type and estimate with the
  highest standard error.  If \code{FALSE}, the function returns a
  matrix with as many rows as there are parameter types in the model.
  In the latter case, the estimate with the highest standard error for
  each parameter type is presented.
}
\item{\dots}{
  additional arguments passed to the function.
  }
}
\details{
In some complex models such as certain hierarchical models (Royle and
Dorazio 2008, Kéry and Royle 2015), issues in estimating parameters and
their standard errors can occur.  Large standard errors can be
indicative of problems in estimating certain parameters due to sparse
data, parameters on the boundary, or model misspecification.  The
\code{checkParms} function computes the number of parameter estimates
with standard errors larger than \code{se.max} and identifies the
parameter estimate with the largest standard error across all parameter
types (\code{simplify = TRUE}) or for each parameter type
(\code{simplify = FALSE}).

To help identify large standard errors, users can standardize numeric
explanatory variables to zero mean and unit variance.  The
\code{checkParms} function can also be useful to identify boundary
estimates in classic generalized models or their extensions (Venables
and Ripley 2002).
}
\value{
  \code{checkParms} returns a list of class \code{checkParms} with the
  following components:

  \item{model.class}{the class of the model for which diagnostics are
    requested.}

  \item{se.max}{the value of SE used as a threshold in diagnostics.  The
  function reports the number of parameter estimates with SE >
  \emph{se.max}.}
  
  \item{result}{a matrix consisting of three columns, namely, the
    identity of the parameter estimate with the highest SE
    (\code{variable}), its standard error (\code{max.se}), and the
    number of parameter estimates with SE larger than \code{se.max}
    (\code{n.high.se}).  For classical regression models with a single
    response variable, the row name is labeled \code{beta}.  For
    \code{unmarkedFit} models, the matrix either consists of a single
    row (\code{simplify = TRUE)} labeled with the name of the parameter
    type (e.g., psi, gam, eps, p) where the highest SE occurs, or
    consists of as many rows as there are parameter types
    (\code{simplify = FALSE}).
  }
}
\references{
  Agresti, A. (2002) \emph{Categorical data analysis}. John Wiley and
    Sons, Inc.: Hoboken.

  Fiske, I., Chandler, R. (2011) unmarked: An R Package for
  fitting hierarchical models of wildlife occurrence and
  abundance. \emph{Journal of Statistical Software} \bold{43}, 1--23.
    
  Kéry, M., Royle, J. A. (2015) \emph{Applied hierarchical modeling in
  ecology: analysis of distribution, abundance and species richness in R
  and BUGS}. Academic Press, New York, USA. 
  
  Royle, J. A., Dorazio, R. M. (2008) \emph{Hierarchical modeling and
  inference in ecology: the analysis of data from populations,
  metapopulations and communities}. Academic Press: New York.

Venables, W. N., Ripley, B. D. (2002) \emph{Modern applied statistics
  with S}, 2nd edition. Springer-Verlag: New York. 
}

\author{
Marc J. Mazerolle
}

\seealso{
  \code{\link{c_hat}}, \code{\link{detHist}}, \code{\link{checkConv}},
  \code{\link{countDist}}, \code{\link{countHist}},
  \code{\link{extractCN}}, \code{\link{mb.gof.test}},
  \code{\link{Nmix.gof.test}}, \code{\link[unmarked]{parboot}} 
}

\examples{
##example with multiple-season occupancy model modified from ?colext
\dontrun{
require(unmarked)
data(frogs)
umf <- formatMult(masspcru)
obsCovs(umf) <- scale(obsCovs(umf))
siteCovs(umf) <- rnorm(numSites(umf))
yearlySiteCovs(umf) <- data.frame(year = factor(rep(1:7,
                                    numSites(umf))))

##model with with year-dependent transition rates
fm.yearly <- colext(psiformula = ~ 1, gammaformula = ~ year,
                    epsilonformula = ~ year,
                    pformula = ~ JulianDate + I(JulianDate^2),
                    data = umf)

##check for high SE's and report highest
##across all parameter types
checkParms(fm.yearly, simplify = TRUE)

##check for high SE's and report highest
##for each parameter type
checkParms(fm.yearly, simplify = FALSE)
detach(package:unmarked)
}

##example from Agresti 2002 of logistic regression
##with parameters estimated at the boundary (complete separation)
\dontrun{
x <- c(10, 20, 30, 40, 60, 70, 80, 90)
y <- c(0, 0, 0, 0, 1, 1, 1, 1)

m1 <- glm(y ~ x, family = binomial)
checkParms(m1)
}
}


\keyword{models}

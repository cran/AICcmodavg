\name{aictabCustom}
\Rdversion{1.1}
\alias{aictabCustom}
\title{
Create Model Selection Tables from User-supplied Input Based on (Q)AIC(c)
}
\description{
  This function creates a model selection table from model input
  (log-likelihood, number of estimated parameters) supplied by the user
  instead of extracting the values automatically from a list of
  candidate models.  The models are ranked based on one of the following
  information criteria:  AIC, AICc, QAIC, QAICc.  The table ranks the
  models based on the selected information criteria and also provides
  delta AIC and Akaike weights.  
}
\usage{
aictabCustom(logL, K, modnames = NULL, second.ord = TRUE, nobs = NULL,
             sort = TRUE, c.hat = 1) 
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{logL}{
    a vector of log-likelihood values for the models in the candidate
    model set.
  }
  \item{K}{
    a vector containing the number of estimated parameters for each
    model in the candidate model set.
  }
  \item{modnames}{
    a character vector of model names to facilitate the identification of
    each model in the model selection table. If \code{NULL}, the function
    uses the names in the cand.set list of candidate models (i.e., a named
    list).  If no names appear in the list and no character vector is
    provided, generic names (e.g., \code{Mod1}, \code{Mod2}) are supplied in
    the table in the same order as in the list of candidate models.
  }
  \item{second.ord}{
    logical.  If \code{TRUE}, the function returns the second-order
    Akaike information criterion (i.e., AICc).
  }
  \item{nobs}{
    the sample size required to compute the AICc or QAICc.
  }
  \item{sort}{
    logical.  If \code{TRUE}, the model selection table is ranked according
    to the (Q)AIC(c) values.
  }
  \item{c.hat}{
   value of overdispersion parameter (i.e., variance inflation factor)
   such as that obtained from \code{c_hat}.  Note that values of
   \code{c.hat} different from 1 are only appropriate for binomial GLM's
   with trials > 1 (i.e., success/trial or cbind(success, failure)
   syntax), with Poisson GLM's, single-season or dynamic occupancy
   models (MacKenzie et al. 2002, 2003), \emph{N}-mixture models (Royle
   2004, Dail and Madsen 2011), or capture-mark-recapture models (e.g.,
   Lebreton et al. 1992).  If c.hat > 1, \code{aictabCustom} will return the
   quasi-likelihood analogue of the information criterion requested.
}
}
\details{
  \code{aictabCustom} constructs a model selection table based on one of
  the four information criteria: AIC, AICc, QAIC, and QAICc.  This
  function is most useful when model input is imported into R from other
  software (e.g., Program MARK, PRESENCE) or for model classes that are
  not yet supported by \code{aictab}.
}

  \value{
    \code{aictabCustom} creates an object of class \code{aictab} with the
    following components:

    \item{Modname}{the name of each model of the candidate model set.}
    \item{K }{the number of estimated parameters for each model.}
    \item{(Q)AIC(c) }{the information criteria requested for each model
      (AICc, AICc, QAIC, QAICc).}
    \item{Delta_(Q)AIC(c) }{the appropriate delta AIC component depending on
      the information criteria selected.}
    \item{ModelLik }{the relative likelihood of the model given the
      data (exp(-0.5*delta[i])).  This is not to be confused with the
      likelihood of the parameters given the data.  The relative likelihood
      can then be normalized across all models to get the model probabilities.} 
    \item{(Q)AIC(c)Wt}{the Akaike weights, also termed "model
      probabilities" sensu Burnham and Anderson (2002) and Anderson (2008).
      These measures indicate the level of support (i.e., weight of
      evidence) in favor of any given model being the most parsimonious
      among the candidate model set.}
    \item{Cum.Wt }{the cumulative Akaike weights.  These are only meaningful
      if results in table are sorted in decreasing order of Akaike weights
      (i.e., sort = TRUE).}
    \item{c.hat}{if c.hat was specified as an argument, it is included in
      the table.}
    \item{LL}{if c.hat = 1 and parameters estimated by maximum likelihood,
      the log-likelihood of each model.}
    \item{Quasi.LL}{if c.hat > 1, the quasi log-likelihood of each model.}
  }
\references{
Anderson, D. R. (2008) \emph{Model-based Inference in the Life Sciences:
a primer on evidence}. Springer: New York.

Burnham, K. P., Anderson, D. R. (2002) \emph{Model Selection and
Multimodel Inference: a practical information-theoretic
approach}. Second edition. Springer: New York.

Dail, D., Madsen, L. (2011) Models for estimating abundance from 
repeated counts of an open population. \emph{Biometrics} \bold{67},
577--587. 

Lebreton, J.-D., Burnham, K. P., Clobert, J., Anderson, D. R. (1992)
Modeling survival and testing biological hypotheses using marked
animals: a unified approach with case-studies. \emph{Ecological
  Monographs} \bold{62}, 67--118. 

MacKenzie, D. I., Nichols, J. D., Lachman, G. B., Droege, S., Royle,
J. A., Langtimm, C. A. (2002) Estimating site occupancy rates when
detection probabilities are less than one. \emph{Ecology} \bold{83},
2248--2255.

MacKenzie, D. I., Nichols, J. D., Hines, J. E., Knutson, M. G.,
Franklin, A. B. (2003) Estimating site occupancy, colonization, and
local extinction when a species is detected imperfectly. \emph{Ecology}
\bold{84}, 2200--2207.

Mazerolle, M. J. (2006) Improving data analysis in herpetology: using
Akaike's Information Criterion (AIC) to assess the strength of
biological hypotheses. \emph{Amphibia-Reptilia} \bold{27}, 169--180.

Royle, J. A. (2004) \emph{N}-mixture models for estimating population
size from spatially replicated counts. \emph{Biometrics} \bold{60},
108--115.
}
\author{
Marc J. Mazerolle
}
\seealso{
\code{\link{AICcCustom}}, \code{\link{bictabCustom}},
\code{\link{confset}}, \code{\link{c_hat}}, \code{\link{evidence}},
\code{\link{ictab}}, \code{\link{modavgCustom}}
}
\examples{
##vector with model LL's
LL <- c(-38.8876, -35.1783, -64.8970)

##vector with number of parameters
Ks <- c(7, 9, 4)

##create a vector of names to trace back models in set
Modnames <- c("Cm1", "Cm2", "Cm3")

##generate AICc table
aictabCustom(logL = LL, K = Ks, modnames = Modnames, nobs = 121,
             sort = TRUE) 
}

\keyword{models}

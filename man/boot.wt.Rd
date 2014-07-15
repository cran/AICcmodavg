\name{boot.wt}
\Rdversion{1.1}
\alias{boot.wt}
\alias{boot.wt.default}
\alias{boot.wt.AICaov.lm}
\alias{boot.wt.AICsclm.clm}
\alias{boot.wt.AICglm.lm}
\alias{boot.wt.AIClm}
\alias{boot.wt.AICmultinom.nnet}
\alias{boot.wt.AICpolr}
\alias{boot.wt.AICrlm.lm}
\alias{boot.wt.AICvglm}
\alias{boot.wt.AICzeroinfl}
\title{
Compute Model Selection Relative Frequencies
}
\description{
This function computes the model selection relative frequencies based on
the nonparametric bootstrap (Burnham and Anderson 2002).  Models are
ranked based on the AIC, AICc, QAIC, or QAICc.  The function currently
supports objects of \code{clm}, \code{glm}, \code{lm}, \code{multinom},
\code{polr}, \code{rlm}, \code{vglm}, and \code{zeroinfl} classes.
}
\usage{
boot.wt(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL,
        sort = TRUE, nsim = 100, \dots)

\method{boot.wt}{AICaov.lm}(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL,
       sort = TRUE, nsim = 100, \dots)

\method{boot.wt}{AICsclm.clm}(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL,
       sort = TRUE, nsim = 100, \dots)

\method{boot.wt}{AICglm.lm}(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL,
       sort = TRUE, nsim = 100, c.hat = 1, \dots)
       
\method{boot.wt}{AIClm}(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL,
       sort = TRUE, nsim = 100, \dots)

\method{boot.wt}{AICmultinom.nnet}(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL,
       sort = TRUE, nsim = 100, c.hat = 1, \dots)

\method{boot.wt}{AICpolr}(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL,
       sort = TRUE, nsim = 100, \dots)

\method{boot.wt}{AICrlm.lm}(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL,
       sort = TRUE, nsim = 100, \dots)

\method{boot.wt}{AICvglm}(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL,
       sort = TRUE, nsim = 100, c.hat = 1, \dots)

\method{boot.wt}{AICzeroinfl}(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL,
       sort = TRUE, nsim = 100, \dots)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cand.set}{
a list storing each of the models in the candidate model set.
}
  \item{modnames}{
a character vector of model names to facilitate the identification of
each model in the model selection table. If \code{NULL}, the function
uses the names in the cand.set list of candidate models. If no names
appear in the list, generic names (e.g., \code{Mod1}, \code{Mod2}) are
supplied in the table in the same order as in the list of candidate
models.
}
  \item{second.ord}{
logical.  If \code{TRUE}, the function returns the second-order
Akaike information criterion (i.e., AICc).
}
  \item{nobs}{
this argument allows to specify a numeric value other than total sample
size to compute the AICc (i.e., \code{nobs} defaults to total number of
observations).  This is relevant only for certain types of models such 
as mixed models where sample size is not straightforward.  In 
such cases, one might use total number of observations or number of 
independent clusters (e.g., sites) as the value of \code{nobs}.
}
\item{sort}{
logical.  If \code{TRUE}, the model selection table is ranked according
to the (Q)AIC(c) values.
}
  \item{c.hat}{
value of overdispersion parameter (i.e., variance inflation factor) such
as that obtained from \code{c_hat}.  Note that values of c.hat different
from 1 are only appropriate for binomial GLM's with trials > 1 (i.e.,
success/trial or cbind(success, failure) syntax) or with Poisson
GLM's.  If c.hat > 1, \code{AICc} will return the quasi-likelihood
analogue of the information criterion requested.
}

\item{nsim}{the number of bootstrap iterations. Burnham and Anderson
  (2002) recommend at least 1000 and up to 10 000 iterations for certain
  problems.
}
\item{\dots}{
  additional arguments passed to the function.
  }
}
\details{
\code{boot.wt} is implemented for \code{glm}, \code{lm},
  \code{multinom}, \code{polr}, \code{rlm}, \code{vglm}, and
  \code{zeroinfl} classes.  During each bootstrap iteration, the data
  are resampled with replacement, all the models specified in
  \code{cand.set} are updated with the new data set, and the top-ranked
  model is saved. When all iterations are completed, the relative
  frequency of selection is computed for each model appearing in the
  candidate model set.

Relative frequencies of the models are often similar to Akaike
weights, and the latter are often preferred due to their link with
a Bayesian perspective (Burnham and Anderson 2002).  \code{boot.wt} is
most useful for teaching purposes of sampling-theory based relative
frequencies of model selection. The current implementation is only
appropriate with completely randomized designs. For more complex data
structures (e.g., blocks or random effects), the bootstrap should be
modified accordingly.
}

\value{
\code{boot.wt} creates an object of class \code{aictab} with the
following components: 

\item{Modname}{the names of each model of the candidate model set.}
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
\item{PiWt }{the relative frequencies of model selection from the bootstrap.}
\item{c.hat}{if c.hat was specified as an argument, it is included in
  the table.}
}
\references{
Anderson, D. R. (2008) \emph{Model-based Inference in the Life Sciences:
a primer on evidence}. Springer: New York.

Burnham, K. P., Anderson, D. R. (2002) \emph{Model Selection and
Multimodel Inference: a practical information-theoretic
approach}. Second edition. Springer: New York.

Burnham, K. P., Anderson, D. R. (2004) Multimodel inference:
understanding AIC and BIC in model selection. \emph{Sociological Methods
and Research} \bold{33}, 261--304.

Mazerolle, M. J. (2006) Improving data analysis in herpetology: using
Akaike's Information Criterion (AIC) to assess the strength of
biological hypotheses. \emph{Amphibia-Reptilia} \bold{27}, 169--180.
}
\author{
Marc J. Mazerolle
}
\seealso{
\code{\link{AICc}}, \code{\link{confset}}, \code{\link{c_hat}},
\code{\link{evidence}}, \code{\link{importance}}, \code{\link{modavg}}, 
\code{\link{modavgShrink}}, \code{\link{modavgPred}} 
}
\examples{
##Mazerolle (2006) frog water loss example
data(dry.frog)

##setup a subset of models of Table 1
Cand.models <- list( )
Cand.models[[1]] <- lm(log_Mass_lost ~ Shade + Substrate +
                       cent_Initial_mass + Initial_mass2,
                       data = dry.frog)
Cand.models[[2]] <- lm(log_Mass_lost ~ Shade + Substrate +
                       cent_Initial_mass + Initial_mass2 +
                       Shade:Substrate, data = dry.frog)
Cand.models[[3]] <- lm(log_Mass_lost ~ cent_Initial_mass +
                       Initial_mass2, data = dry.frog)
Cand.models[[4]] <- lm(log_Mass_lost ~ Shade + cent_Initial_mass +
                       Initial_mass2, data = dry.frog)
Cand.models[[5]] <- lm(log_Mass_lost ~ Substrate + cent_Initial_mass +
                       Initial_mass2, data = dry.frog)

##create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(Cand.models), sep = " ")

##generate AICc table with bootstrapped relative
##frequencies of model selection
boot.wt(cand.set = Cand.models, modnames = Modnames, sort = TRUE)


\dontrun{
##Burnham and Anderson (2002) flour beetle data
data(beetle)
##models as suggested by Burnham and Anderson p. 198          
Cand.set <- list( )
Cand.set[[1]] <- glm(Mortality_rate ~ Dose, family =
                     binomial(link = "logit"), weights = Number_tested,
                     data = beetle)
Cand.set[[2]] <- glm(Mortality_rate ~ Dose, family =
                     binomial(link = "probit"), weights = Number_tested,
                     data = beetle)
Cand.set[[3]] <- glm(Mortality_rate ~ Dose, family =
                     binomial(link ="cloglog"), weights = Number_tested,
                     data = beetle)

##create a vector of names to trace back models in set
Modnames <- paste("Mod", 1:length(Cand.set), sep = " ")

##model selection table with bootstrapped
##relative frequencies
boot.wt(cand.set = Cand.set, modnames = Modnames)
}
}

\keyword{models}

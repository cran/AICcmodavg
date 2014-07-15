\name{confset}
\Rdversion{1.1}
\alias{confset}
\alias{print.confset}

\title{
Computing Confidence Set for the Kullback-Leibler Best Model
}
\description{
This function computes the confidence set on the best model given
the data and model set.  \code{confset} implements three different
methods proposed by Burnham and Anderson (2002).

}
\usage{
confset(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL,
        method = "raw", level = 0.95, delta = 6, c.hat = 1)  
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
logical.  If TRUE, the function returns the second-order
Akaike information criterion (i.e., AICc).
}
  \item{nobs}{
this argument allows to specify a numeric value other than total sample
size to compute the AICc (i.e., 'nobs' defaults to total number of
observations).  This is relevant only for mixed models or various models
of 'unmarkedFit' classes where sample size is not straightforward.  In
such cases, one might use total number of observations or number of
independent clusters (e.g., sites) as the value of 'nobs'.
}
  \item{method}{
a character value, either as 'raw', 'ordinal', or 'ratio', indicating
the method for determining the confidence set for the best model (see
'description' above for details).
}
  \item{level}{
the level of confidence (i.e., sum of model probabilities) used to
determine the confidence set on the best model when using the 'raw'
method.  Note that the argument is not used for the other methods of
determining the confidence set on the best model. 
}
  \item{delta}{
the delta (Q)AIC(c) value associated with the cutoff point to determine
the confidence set for the best model.  Note that the argument is only
used when method = 'ratio'.
}
  \item{c.hat}{
value of overdispersion parameter (i.e., variance inflation factor) such
as that obtained from \code{c_hat}.  Note that values of c.hat different
from 1 are only appropriate for binomial GLM's with trials > 1 (i.e.,
success/trial or cbind(success, failure) syntax), with Poisson
GLM's, or single-season occupancy models (MacKenzie et al. 2002).  If
c.hat > 1, \code{AICc} will return the quasi-likelihood analogue of the
information criterion requested.  This option is not supported for
generalized linear mixed models of the \code{mer} or \code{merMod}
classes.
}
}
\details{
The first and simplest (\code{method = 'raw'}), relies on summing the
Akaike weights (i.e., model probabilities) of the ranked models until we
reach a given cutpoint (e.g., 0.95 for a 95 percent set).

The second method (\code{method = 'ordinal'}) suggested is based on the
classification of the models on an ordinal scale based on the delta
(Q)AIC(c).  The models are grouped in different classes based on their
weight of support as determined by the delta (Q)AIC(c) values:
substantial support (delta (Q)AIC(c) <= 2), some support (2 < delta
(Q)AIC(c) <= 7), little support (7 < delta (Q)AIC(c) <= 10), no support
(delta (Q)AIC(c) > 10).

The third method (\code{method = 'ratio'}) is based on identifying the
ratios of model likelihoods (i.e., exp(-delta_(Q)AIC(c)/2) ) that exceed
a cutpoint, similar to the building of profile likelihood intervals.  An
evidence ratio of each model relative to the top-ranked model is
computed and the ratios exceeding the cutpoint determine which models
are included in the confidence set.  Note here that small cutoff points 
are suggested (e.g., 0.125, 0.050).  The cutoff point is linked to delta
(Q)AIC(c) by the following relationship:  cutoff =
exp(-1*delta_(Q)AIC(c)/2).  

  
}
\value{
\code{confset} returns an object of class \code{confset} as a list with
the following components, depending on which method is used:

when \code{method = 'raw'}:
\item{method}{identifies the method of determining the confidence set on
the best model.}
\item{level}{the confidence level used to determine the confidence set
  on the best model.} 
\item{table}{a reduced table with the models included in the confidence
  set.}

when \code{method = 'ordinal'}:
\item{method}{identifies the method of determining the confidence set on
the best model.}
\item{substantial}{a reduced table with the models included in the
  confidence set for which delta (Q)AIC(c) <= 2.}
\item{some}{a reduced table with the models included in the confidence
  set for which 2 < delta (Q)AIC(c) <= 7.}
\item{little}{a reduced table with the models included in the confidence
  set for which 7 < delta (Q)AIC(c) <= 10.}
\item{none}{a reduced table with the models included in the confidence
  set for which delta (Q)AIC(c) > 10.}

when \code{method = 'ratio'}:
\item{method}{identifies the method of determining the confidence set on
  the best model.}
\item{cutoff}{the cutoff value for the ratios used to determine the
  confidence set on the best model.}
\item{delta}{the delta (Q)AIC(c) used to compute the cutoff value for
  ratios to determine the confidence set on the best model.}
\item{table}{a reduced table with the models included in the confidence
  set.}

}
\references{
Burnham, K. P., Anderson, D. R. (2002) \emph{Model Selection and
  Multimodel Inference: a practical information-theoretic
  approach}. Second edition. Springer: New York.

MacKenzie, D. I., Nichols, J. D., Lachman, G. B., Droege, S., Royle,
J. A., Langtimm, C. A. (2002) Estimating site occupancy rates when
detection probabilities are less than one. \emph{Ecology} \bold{83},
2248--2255.
}
\author{
Marc J. Mazerolle
}

\seealso{
\code{\link{AICc}}, \code{\link{aictab}}, \code{\link{c_hat}},
\code{\link{evidence}}, \code{\link{importance}}, \code{\link{modavg}},
\code{\link{modavgShrink}}, \code{\link{modavgPred}}
}
\examples{
##anuran larvae example from Mazerolle (2006)
data(min.trap)
##assign "UPLAND" as the reference level as in Mazerolle (2006)          
min.trap$Type <- relevel(min.trap$Type, ref = "UPLAND") 

##set up candidate models          
Cand.mod <- list()
##global model          
Cand.mod[[1]] <- glm(Num_anura ~ Type + log.Perimeter + Num_ranatra,
                     family = poisson, offset = log(Effort),
                     data = min.trap) 
Cand.mod[[2]] <- glm(Num_anura ~ Type + log.Perimeter, family = poisson,
                     offset = log(Effort), data = min.trap) 
Cand.mod[[3]] <- glm(Num_anura ~ Type + Num_ranatra, family = poisson,
                     offset = log(Effort), data = min.trap) 
Cand.mod[[4]] <- glm(Num_anura ~ Type, family = poisson,
                     offset = log(Effort), data = min.trap) 
Cand.mod[[5]] <- glm(Num_anura ~ log.Perimeter + Num_ranatra,
                     family = poisson, offset = log(Effort),
                     data = min.trap) 
Cand.mod[[6]] <- glm(Num_anura ~ log.Perimeter, family = poisson,
                     offset = log(Effort), data = min.trap) 
Cand.mod[[7]] <- glm(Num_anura ~ Num_ranatra, family = poisson,
                     offset = log(Effort), data = min.trap) 
Cand.mod[[8]] <- glm(Num_anura ~ 1, family = poisson,
                     offset = log(Effort), data = min.trap) 
          
##check c-hat for global model
c_hat(Cand.mod[[1]]) #uses Pearson's chi-square/df
##note the very low overdispersion: in this case, the analysis could be
##conducted without correcting for c-hat as its value is reasonably close
##to 1  

##assign names to each model
Modnames <- c("type + logperim + invertpred", "type + logperim",
              "type + invertpred", "type", "logperim + invertpred",
              "logperim", "invertpred", "intercept only") 

##compute confidence set based on 'raw' method
confset(cand.set = Cand.mod, modnames = Modnames, second.ord = TRUE,
        method = "raw")  



##example with linear mixed model
require(nlme)

##set up candidate model list for Orthodont data set shown in Pinheiro
##and Bates (2000:  Mixed-effect models in S and S-PLUS. Springer Verlag:
##New York.)
Cand.models <- list()
Cand.models[[1]] <- lme(distance ~ age, random = ~age | Subject,
                        data = Orthodont, method = "ML")
Cand.models[[2]] <- lme(distance ~ age + Sex, data = Orthodont,
                        random = ~ 1 | Subject, method = "ML")
Cand.models[[3]] <- lme(distance ~ 1, data = Orthodont,
                        random = ~ 1 | Subject, method = "ML")

##create a vector of model names
Modnames <- paste("mod", 1:length(Cand.models), sep = "")

##compute confidence set based on 'raw' method
confset(cand.set = Cand.models, modnames = Modnames, second.ord = TRUE,
        method = "raw")
##round to 4 digits after decimal point
print(confset(cand.set = Cand.models, modnames = Modnames,
              second.ord = TRUE, method = "raw"), digits = 4)

confset(cand.set = Cand.models, modnames = Modnames, second.ord = TRUE,
        level = 0.9, method = "raw")

##compute confidence set based on 'ordinal' method
confset(cand.set = Cand.models, modnames = Modnames, second.ord = TRUE,
        method = "ordinal")

##compute confidence set based on 'ratio' method
confset(cand.set = Cand.models, modnames = Modnames, second.ord = TRUE,
        method = "ratio", delta = 4)

confset(cand.set = Cand.models, modnames = Modnames, second.ord = TRUE,
        method = "ratio", delta = 8)
}

\keyword{models}

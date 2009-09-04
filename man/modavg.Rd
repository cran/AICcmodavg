\name{modavg}
\Rdversion{1.1}
\alias{modavg}
\alias{modavg.glm}
\alias{modavg.lme}
\alias{print.modavg}

\title{
Model-average Parameter Estimate (Multimodel Inference)
}
\description{
This function model-averages the estimate of a parameter of interest
among a set of candidate models, computes the unconditional standard
error and unconditional confidence intervals as described in Buckland et
al. (1997) and Burnham and Anderson (2002).
}
\usage{
modavg(cand.set, parm, modnames, c.hat = 1, gamdisp = NULL,
conf.level = 0.95, second.ord =
TRUE, nobs = NULL)

modavg.glm(cand.set, parm, modnames, c.hat = 1, gamdisp = NULL,
conf.level = 0.95, second.ord = TRUE, nobs = NULL) 

modavg.lme(cand.set, parm, modnames, conf.level = 0.95, second.ord =
TRUE, nobs = NULL) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cand.set}{
a list storing each of the models in the candidate model set.
}
  \item{parm}{
the parameter of interest for which a model-averaged estimate is
required. 
}
  \item{modnames}{
a character vector of model names to facilitate the identification of
each model in the model selection table.
}
  \item{c.hat}{
value of overdispersion parameter (i.e., variance inflation factor) such
as that obtained from 'c_hat'.  Note that values of c.hat different from
1 are only appropriate either for binomial GLM's with trials > 1 (i.e.,
success/trial or cbind(success, failure) syntax) or for Poisson GLM's.
If c.hat > 1, 'aictab' will return the quasi-likelihood analogue of the
information criteria requested and multiply the variance-covariance
matrix of the estimates by this value (i.e., SE's are multiplied by
sqrt(c.hat)).
}
  \item{gamdisp}{
if gamma GLM is used, the dispersion parameter should be specified here
to apply the same value to each model.
}
  \item{conf.level}{
the confidence level requested for the computation of unconditional
confidence intervals.
}
  \item{second.ord}{
logical. If TRUE, the function returns the second-order
Akaike information criterion (i.e., AICc).
}
  \item{nobs}{
this argument allows to specify a numeric value other than total sample
size to compute the AICc.  This is relevant only for linear mixed models
where sample size is not straightforward.  In such cases, one might use
total number of observations or number of independent clusters as the
value of 'nobs'. 
}
}
\details{
The parameter must be specified as it appears in the model output. For
factors, one must specify the name of the variable and the level of
interest.

In cases where one desires a model-averaged estimate for a main
effect term that is also involved in an interaction in certain models,
modavg (correctly) excludes models with interactions involving the main
effect term when it finds that the variable also appears in an
interaction.  Care must be taken when some models include interaction or
polynomial terms, because main effect terms do not have the same
interpretation when an interaction term appears in the same model.  To
avoid problems, one should use the usual operators to build interaction
terms or polynomial terms.  If polynomial or interaction terms are built
manually, you should include the entire name of the variables involved
because modavg( ) uses regexpr( ) to search for duplicates in model
formula. 

'modavg' is a function that calls either 'modavg.glm' or 'modavg.lme',
depending on the class of the object.  The current function is
implemented for a list containing objects of 'lm','glm', and 'lme'
classes. 
}
\value{
'modavg', 'modavg.glm', and 'modavg.lme' create an object of
class 'modavg' with the following components:

\item{Parameter}{the parameter for which a model-averaged estimate was
  obtained} 
\item{Mod.avg.table}{the reduced model selection table based on models
  including the parameter of interest}
\item{Mod.avg.beta}{the model-averaged estimate based on all models
  including the parameter of interest (see details above regarding the
  exclusion of models where parameter of interest is involved in an
  interaction)}
\item{Uncond.SE}{the unconditional standard error for the model-averaged
  estimate (as opposed to the conditional SE based on a single model)}
\item{Conf.level}{the confidence level used to compute the confidence
  interval} 
\item{Lower.CL}{the lower confidence limit}
\item{Upper.CL}{the upper confidence limit}

}
\references{
Anderson, D. R. (2008) \emph{Model-based Inference in the Life Sciences:
  a primer on evidence}. Springer: New York.

Buckland, S. T., Burnham, K. P., Augustin, N. H. (1997) Model selection:
an integral part of inference. \emph{Biometrics} \bold{53}, 603--618.  

Burnham, K. P., Anderson, D. R. (2002) \emph{Model Selection and
  Multimodel Inference: a practical information-theoretic
  approach}. Second edition. Springer: New York.

Burnham, K. P., Anderson, D. R. (2004) Multimodel inference:
  understanding AIC and BIC in model selection. \emph{Sociological
  Methods and Research} \bold{33}, 261--304.

Mazerolle, M. J. (2006) Improving data analysis in herpetology: using
Akaike's Information Criterion (AIC) to assess the strength of
biological hypotheses. \emph{Amphibia-Reptilia} \bold{27}, 169--180. 

}
\author{
Marc J. Mazerolle
}

\seealso{
\code{\link{AICc}}, \code{\link{aictab}}, \code{\link{c_hat}},
\code{\link{importance}}, \code{\link{confset}}, \code{\link{evidence}},
\code{\link{modavgpred}}
}
\examples{
#anuran larvae example from Mazerolle (2006)
data(min.trap)
#assign "UPLAND" as the reference level as in Mazerolle (2006)          
min.trap$Type <- relevel(min.trap$Type, ref = "UPLAND") 

#set up candidate models          
Cand.mod <- list()
#global model          
Cand.mod[[1]] <- glm(Num_anura ~ Type + log.Perimeter + Num_ranatra,
family = poisson, offset = log(Effort), data = min.trap) 
Cand.mod[[2]] <- glm(Num_anura ~ Type + log.Perimeter, family = poisson,
offset = log(Effort), data = min.trap) 
Cand.mod[[3]] <- glm(Num_anura ~ Type + Num_ranatra, family = poisson,
offset = log(Effort), data = min.trap) 
Cand.mod[[4]] <- glm(Num_anura ~ Type, family = poisson,
offset = log(Effort), data = min.trap) 
Cand.mod[[5]] <- glm(Num_anura ~ log.Perimeter + Num_ranatra,
family = poisson, offset = log(Effort), data = min.trap) 
Cand.mod[[6]] <- glm(Num_anura ~ log.Perimeter, family = poisson,
offset = log(Effort), data = min.trap) 
Cand.mod[[7]] <- glm(Num_anura ~ Num_ranatra, family = poisson,
offset = log(Effort), data = min.trap) 
Cand.mod[[8]] <- glm(Num_anura ~ 1, family = poisson,
offset = log(Effort), data = min.trap) 
          
#check c-hat for global model
c_hat(Cand.mod[[1]]) #uses Pearson's chi-square/df
#note the very low overdispersion: in this case, the analysis could be
#conducted without correcting for c-hat as its value is reasonably close
#to 1  

#assign names to each model
Modnames <- c("type + logperim + invertpred", "type + logperim", "type +
invertpred", "type", "logperim + invertpred", "logperim", "invertpred",
"intercept only") 

#compute model-averaged estimate of TypeBOG
modavg(parm = "TypeBOG", cand.set = Cand.mod, modnames = Modnames)


#compute residual deviance as in Mazerolle (2006)          
Cand.mod[[1]]$deviance/Cand.mod[[1]]$df.residual

#compute model-averaged estimate of TypeBOG as in Table 4 of
#Mazerolle (2006) 
modavg(parm = "TypeBOG", cand.set = Cand.mod, modnames = Modnames,
c.hat = 1.11) 


}
\keyword{models}
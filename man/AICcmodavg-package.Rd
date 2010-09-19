\name{AICcmodavg-package}
\Rdversion{1.1}
\alias{AICcmodavg-package}
\alias{AICcmodavg}
\docType{package}
\title{
Model Selection and Multimodel Inference Based on (Q)AIC(c)
}
\description{
Description:  This package includes functions to create model selection
tables based on Akaike's information criterion (AIC) and the
second-order AIC (AICc), as well as their quasi-likelihood counterparts
(QAIC, QAICc).  Tables are printed with delta AIC and Akaike
weights.  The package also includes functions to conduct model averaging
(multimodel inference) for a given parameter of interest or predicted
values.  Other handy functions enable the computation of relative
variable importance, evidence ratios, and confidence sets for the best
model.  The present version works with linear models ('lm' class),
generalized linear models ('glm' class), linear models fit by
generalized least squares ('gls' class), linear mixed models ('lme'
class), generalized linear mixed models ('mer' class), multinomial and
ordinal logistic regressions ('multinom' and 'polr' classes).  
}
\details{
\tabular{ll}{
Package: \tab AICcmodavg\cr
Type: \tab Package\cr
Version: \tab 1.09\cr
Date: \tab 2010-09-16\cr
License: \tab GPL (>=2 )\cr
LazyLoad: \tab yes\cr
}
This package contains several useful functions for model selection and
multimodel inference:
\itemize{
\item{\code{\link{AICc}}}  {Computes AIC, AICc, and their quasi-likelihood
  counterparts (QAIC, QAICc).}
\item{\code{\link{aictab}}}  {Constructs model selection tables with
  number of parameters, AIC, delta AIC, Akaike weights or variants based
  on other AICc, QAIC, and QAICc for a set of candidate models.}
\item{\code{\link{confset}}}  {Determines the confidence set for the best
  model based on one of three criteria.}
\item{\code{\link{evidence}}}  {Computes the evidence ratio between the
  highest-ranked model based on the information criteria selected and a
  lower-ranked model.}
\item{\code{\link{extractSE.mer}}}  {Extracts standard errors of the
  fixed effects of a generalized linear mixed model of class 'mer' and
  adds the labels.}
\item{\code{\link{fam.link.mer}}}  {Extracts the distribution family and 
  link function from a generalized linear mixed model of class 'mer'.} 
\item{\code{\link{importance}}}  {Computes importance values (w+) for
  the support of a given parameter among set of candidate models.}
\item{\code{\link{modavg}}}  {Computes model-averaged estimate,
  unconditional standard error, and unconditional confidence interval of a
  parameter of interest among a set of candidate models.}
\item{\code{\link{modavgpred}}}  {Computes model-average predictions and
  unconditional SE's among entire set of candidate models.}
\item{\code{\link{c_hat}}}  {Computes an estimate of variance
  inflation factor for binomial or Poisson GLM's based on Pearson's
  chi-square.}
\item{\code{\link{predictSE.lme}}}  {Computes predictions and associated 
  standard errors for 'lme' object based on fixed effects (i.e., 
  population predictions).}
\item{\code{\link{predictSE.mer}}}  {Computes predictions and associated 
  standard errors for 'mer' object based on fixed effects (i.e., 
  population predictions).} 
}
}
\author{
  Marc J. Mazerolle <marc.mazerolle@uqat.ca>.  Special thanks to
  T. Ergon for the original idea of storing candidate models in a list.  
}
\references{
Anderson, D. R. (2008) \emph{Model-based inference in the life sciences:
  a primer on evidence}. Springer: New York.

Burnham, K. P., and Anderson, D. R. (2002) \emph{Model selection and
  multimodel inference: a practical information-theoretic approach}. Second
  edition. Springer: New York. 

Burnham, K. P., Anderson, D. R. (2004) Multimodel inference:
  understanding AIC and BIC in model selection. \emph{Sociological
  Methods and Research} \bold{33}, 261--304.

Mazerolle, M. J. (2006) Improving data analysis in herpetology: using
Akaike's Information Criterion (AIC) to assess the strength of
biological hypotheses. \emph{Amphibia-Reptilia} \bold{27}, 169--180. 
}
\keyword{ models }

\name{AICc}
\Rdversion{1.1}
\alias{AICc}
\alias{AICc.glm}
\alias{AICc.gls}
\alias{AICc.lme}
\alias{AICc.mer}
\alias{AICc.mult}
\alias{AICc.polr}
\title{
Computing AIC, AICc, QAIC, and QAICc
}
\description{
Functions to computes Akaike's information criterion (AIC), the
second-order AIC (AICc), as well as their quasi-likelihood counterparts
(QAIC, QAICc).
}
\usage{
AICc(mod, return.K = FALSE, c.hat = 1, second.ord = TRUE, nobs = NULL) 

AICc.glm(mod, return.K = FALSE, c.hat = 1, second.ord = TRUE, nobs = 
NULL)

AICc.gls(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL)

AICc.lme(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL)

AICc.mer(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL)

AICc.mult(mod, return.K = FALSE, c.hat = 1, second.ord = TRUE, nobs =
NULL) 

AICc.polr(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL)

}
\arguments{
  \item{mod}{
an object of class 'lm', 'glm', 'gls', 'lme', 'mer', 'multinom', or 'polr' containing
the output of a model. 
}
  \item{return.K}{
logical. If FALSE, the function returns the information criteria
specified. If TRUE, the function returns K (number of estimated
parameters) for a given model.  Using this argument to facilitate
computation of tables was an original idea from T. Ergon.
}
  \item{c.hat}{
value of overdispersion parameter (i.e., variance inflation factor) such
as that obtained from 'c_hat'.  Note that values of c.hat different from
1 are only appropriate for binomial GLM's with trials > 1 (i.e.,
success/trial or cbind(success, failure) syntax) or with Poisson
GLM's.  If c.hat > 1, 'AICc' will return the quasi-likelihood analogue
of the information criterion requested.  This option is not yet
supported for generalized linear mixed models of the 'mer' class. 
}

\item{second.ord}{ 
logical.  If TRUE, the function returns the second-order Akaike 
information criterion (i.e., AICc). 
}
\item{nobs}{
this argument allows to specify a numeric value other than total sample
size to compute the AICc (i.e., 'nobs' defaults to total number of
observations).  This is relevant only for linear mixed models where
sample size is not straightforward. In such cases, one might use total
number of observations or number of independent clusters as the value of
'nobs'.
}
}
\details{
'AICc' is a function that calls either 'AICc.glm', 'AICc.gls',
'AICc.lme', 'AICc.mer', 'AICc.mult', or 'AICc.polr', depending on the
class of the object.  The current function is implemented for
'lm','glm', 'gls', 'lme', 'mer', 'multinom', and 'polr' classes and
computes one of the following four information criteria:  Akaike's
information criterion (AIC, Akaike 1973), the second-order or small
sample AIC (AICc, Sugiura 1978, Hurvich and Tsai 1991), the
quasi-likelihood AIC (QAIC, Burnham and Anderson 2002), and the
quasi-likelihood AICc (QAICc, Burnham and Anderson 2002). Note that AIC
and AICc values are meaningful to select among 'gls' or 'lme' models fit
by maximum likelihood; AIC and AICc based on REML are valid to select
among different models that only differ in their random effects
(Pinheiro and Bates 2000). 

}
\value{
'AICc' selects one of the functions below based on the class of the 
object:

'AICc.glm' returns the AIC, AICc, QAIC, or QAICc depending on the values
of the arguments.

'AICc.gls' returns the AIC or AICc depending on the values of the
arguments. 

'AICc.lme' returns the AIC or AICc depending on the values of the
arguments. 

'AICc.mer' returns the AIC or AICc depending on the values of the
arguments. 

'AICc.mult' returns the AIC, AICc, QAIC, or QAICc depending on the
values of the arguments.

'AICc.polr' returns the AIC or AICc depending on the values of
the arguments.


}

\note{
The actual (Q)AIC(c) values are not really interesting in themselves, as
they depend directly on the data, parameters estimated, and
likelihood function.  Furthermore, a single value does not tell much
about model fit.  Information criteria become relevant when compared to
one another for a given data set and set of candidate models.
}

\references{
Akaike, H. (1973) Information theory as an extension of the maximum 
likelihood principle. In:  \emph{Second International Symposium on Information 
Theory}, pp. 267--281. Petrov, B.N., Csaki, F., Eds, Akademiai Kiado,
Budapest. 

Anderson, D. R. (2008) \emph{Model-based Inference in the Life Sciences: a
primer on evidence}. Springer: New York.

Burnham, K. P., Anderson, D. R. (2002) \emph{Model Selection and
Multimodel Inference: a practical information-theoretic
approach}. Second edition. Springer: New York.

Burnham, K. P., Anderson, D. R. (2004) Multimodel inference:
understanding AIC and BIC in model selection. \emph{Sociological Methods and Research} \bold{33},
261--304.

Hurvich, C. M., Tsai, C.-L. (1991) Bias of the corrected AIC
criterion for underfitted regression and time series
models. \emph{Biometrika} \bold{78}, 499--509. 

Pinheiro, J. C., Bates, D. M. (2000) \emph{Mixed-effect models in S and
S-PLUS}. Springer Verlag: New York.  

Sugiura, N. (1978) Further analysis of the data by Akaike's information 
criterion and the finite corrections. \emph{Communications in
Statistics:  Theory and Methods} \bold{A7}, 13--26. 
}
\author{
Marc J. Mazerolle
}

\seealso{
\code{\link{aictab}}, \code{\link{confset}}, \code{\link{importance}},
\code{\link{evidence}}, \code{\link{c_hat}}, \code{\link{modavg}},
\code{\link{modavgpred}}
}
\examples{
##cement data from Burnham and Anderson (2002, p. 101)
data(cement)
##run multiple regression - the global model in Table 3.2
glob.mod <- lm(y ~ x1 + x2 + x3 + x4, data = cement)

##compute AICc with full likelihood
AICc(glob.mod, return.K = FALSE)

##compute AIC with full likelihood 
AICc(glob.mod, return.K = FALSE, second.ord = FALSE)
##note that Burnham and Anderson (2002) did not use full likelihood
##in Table 3.2 and that the MLE estimate of the variance was
##rounded to 2 digits after decimal point  



##compute AICc for mixed model on Orthodont data set in Pinheiro and
##Bates (2000)
require(nlme)
m1 <- lme(distance ~ age, random = ~1 | Subject, data = Orthodont,
method= "ML")
AICc(m1, return.K = FALSE)

}
\keyword{models}
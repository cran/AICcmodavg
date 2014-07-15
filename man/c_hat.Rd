\name{c_hat}
\Rdversion{1.1}
\alias{c_hat}
\title{
Compute Estimate of Dispersion for Poisson and Binomial GLM's
}
\description{
This function computes an estimate of c-hat for binomial or Poisson
GLM's based on Pearson's chi-square divided by the residual degrees of freedom.
}
\usage{
c_hat(mod)
}
\arguments{
  \item{mod}{
the model for which a c-hat estimate is required.
}
}
\details{
Poisson and binomial GLM's do not have a parameter for the variance and
it is usually held fixed to 1 (i.e., mean = variance).  However, one must
check whether this assumption is appropriate by estimating the
overdispersion parameter (c-hat).  Though one can obtain an estimate of
c-hat by dividing the residual deviance by the residual degrees of
freedom, McCullagh and Nelder (1989) recommend using Pearson's
chi-square divided by the residual degrees of freedom, which performs
better.  The latter is the method implemented by \code{c_hat}.

Note that values of c-hat > 1 indicate overdispersion (variance > mean),
but that values much higher than 1 (i.e., > 4) probably indicate
lack-of-fit.  In cases of moderate overdispersion, one usually
multiplies the variance-covariance matrix of the estimates by c-hat.  As
a result, the SE's of the estimates are inflated (c-hat is also known as
a variance inflation factor). 

In model selection, c-hat should be estimated from the global model of
the candidate model set and the same value of c-hat applied to the
entire model set.  Specifically, a global model is the most complex model
from which all the other models of the set are simpler versions
(nested).  When no single global model exists in the set of models 
considered, such as when sample size does not allow a complex model, one
can estimate c-hat from 'subglobal' models.  Here, 'subglobal' models
denote models from which only a subset of the models of the candidate
set can be derived.  In such cases, one can use the smallest value of
c-hat for model selection (Burnham and Anderson 2002).

Note that c-hat counts as an additional parameter estimated and should
be added to K.  All functions in package \code{AICcmodavg} automatically
add 1 when the \code{c.hat} argument > 1 and apply the same value of
c-hat for the entire model set.  When \code{c.hat > 1}, functions
compute quasi-likelihood information criteria (either QAICc or QAIC,
depending on the value of the \code{second.ord} argument) by scaling the
log-likelihood of the model by \code{c.hat}.  The value of \code{c.hat}
can influence the ranking of the models:  as c-hat increases, QAIC or
QAICc will favor models with fewer parameters.  As an additional check
against this potential problem, one can create several model selection
tables by incrementing values of c-hat to assess the model selection
uncertainty.  If ranking changes little up to the c-hat value observed,
one can be confident in making inference. 

In cases of underdispersion (c-hat < 1), it is recommended to keep the
value of \code{c.hat} to 1.  However, note that values of c-hat << 1 can
also indicate lack-of-fit and that an alternative model (and distribution)
should be investigated. 

Note that it is only possible to estimate c-hat for binomial models with
trials > 1 (i.e., success/trial or cbind(success, failure) syntax) or
with Poisson GLM's.
}
\value{
\code{c_hat} returns the estimated c-hat value
}
\references{
Anderson, D. R. (2008) \emph{Model-based Inference in the Life Sciences:
  a primer on evidence}. Springer: New York. 

Burnham, K. P., Anderson, D. R. (2002) \emph{Model Selection and
  Multimodel Inference: a practical information-theoretic
  approach}. Second edition. Springer: New York.

Burnham, K. P., Anderson, D. R. (2004) Multimodel inference:
  understanding AIC and BIC in model selection. \emph{Sociological
  Methods and Research} \bold{33}, 261--304. 

Mazerolle, M. J. (2006) Improving data analysis in herpetology: using
Akaike's Information Criterion (AIC) to assess the strength of
biological hypotheses. \emph{Amphibia-Reptilia} \bold{27}, 169--180.

McCullagh, P., Nelder, J. A. (1989) \emph{Generalized Linear
  Models}. Second edition. Chapman and Hall: New York. 
}
\author{
Marc J. Mazerolle
}

\seealso{
\code{\link{AICc}}, \code{\link{confset}}, \code{\link{evidence}}, 
\code{\link{modavg}}, \code{\link{importance}},
\code{\link{modavgPred}}, \code{\link{mb.gof.test}}, \code{\link{Nmix.gof.test}}
}

\examples{
#binomial glm example
set.seed(seed = 10)
resp <- rbinom(n = 60, size = 1, prob = 0.5)
set.seed(seed = 10)
treat <- as.factor(sample(c(rep(x = "m", times = 30), rep(x = "f",
                                           times = 30))))
age <- as.factor(c(rep("young", 20), rep("med", 20), rep("old", 20)))
#each invidual has its own response (n = 1)
mod1 <- glm(resp ~ treat + age, family = binomial)
\dontrun{
c_hat(mod1) #gives an error because model not appropriate for
##computation of c-hat
}

##computing table to summarize successes
table(resp, treat, age)
dat2 <- as.data.frame(table(resp, treat, age)) #not quite what we need
data2 <- data.frame(success = c(9, 4, 2, 3, 5, 2),
                    sex = c("f", "m", "f", "m", "f", "m"),
                    age = c("med", "med", "old", "old", "young",
                      "young"), total = c(13, 7, 10, 10, 7, 13))
data2$prop <- data2$success/data2$total
data2$fail <- data2$total - data2$success

##run model with success/total syntax using weights argument
mod2 <- glm(prop ~ sex + age, family = binomial, weights = total,
            data = data2)
c_hat(mod2)

##run model with other syntax cbind(success, fail)
mod3 <- glm(cbind(success, fail) ~ sex + age, family = binomial,
            data = data2) 
c_hat(mod3)
}

\keyword{models}

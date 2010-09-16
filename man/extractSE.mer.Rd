\name{extractSE.mer}
\Rdversion{1.1}
\alias{extractSE.mer}
\title{
Extract SE of Fixed Effects of 'glmer' Fit
}
\description{
This function extracts the standard errors (SE) of the fixed effects of a
generalized linear mixed model fit with 'glmer' and adds the appropriate
labels.
}
\usage{
extractSE.mer(mod)
}

\arguments{
  \item{mod}{
an object of 'mer' class resulting from the fit of 'glmer'.
}
}
\details{
This is an extractor function that uses 'vcov.mer' from the 'lme4'
  package and it is called by 'modavg.mer'.
}
\value{
Returns the SE's of the fixed effects with the appropriate labels for each.
}

\author{
Marc J. Mazerolle
}

\seealso{
\code{\link{modavg}}, \code{\link[lme4]{glmer}}, \code{\link[lme4]{lmer}}
}
\examples{
##modified example from ?glmer
require(lme4)
##create proportion of incidence
cbpp$prop <- cbpp$incidence/cbpp$size
gm1 <- glmer(prop ~ period + (1 | herd), family = binomial, weights =
size, data = cbpp)
##print summary
summary(gm1)
##extract variance-covariance matrix of fixed effects
vcov(gm1)
##extract SE's of fixed effects - no labels
sqrt(diag(vcov(gm1))) #no labels
extractSE.mer(gm1)  #with labels
detach(package:lme4)
}
\keyword{models}


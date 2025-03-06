\name{extractSE}
\Rdversion{1.1}
\alias{extractSE}
\alias{extractSE.default}
\alias{extractSE.coxme}
\alias{extractSE.lmekin}
\alias{extractSE.mer}
\alias{extractSE.merMod}
\alias{extractSE.lmerModLmerTest}
\title{
Extract SE of Fixed Effects
}
\description{
This function extracts the standard errors (SE) of the fixed effects of
a mixed model fit with \code{coxme}, \code{glmer}, \code{lmer},
\code{lmerModLmerTest}, and \code{lmekin} and adds the appropriate labels.
}
\usage{

extractSE(mod, \dots)

\method{extractSE}{coxme}(mod, \dots)

\method{extractSE}{lmekin}(mod, \dots)

\method{extractSE}{mer}(mod, \dots)

\method{extractSE}{merMod}(mod, \dots)

\method{extractSE}{lmerModLmerTest}(mod, \dots)
}

\arguments{
  \item{mod}{
    an object of \code{coxme}, \code{lmekin}, \code{mer}, \code{merMod},
    or \code{lmerModTest} class.
}
\item{\dots}{additional arguments passed to the function.}
}
\details{
These extractor functions use \code{vcov.coxme}, \code{vcov.lmekin},
  \code{vcov.mer}, and \code{vcov.merMod}. Some of these functions are
  called by \code{modavg} and \code{modavgShrink}, depending on the
  class of the objects. 
}
\value{
Returns the SE's of the fixed effects with the appropriate labels for each.
}

\author{
Marc J. Mazerolle
}

\seealso{
\code{\link{modavg}}, \code{\link[lme4]{glmer}},
\code{\link[lme4]{lmer}}, \code{\link[coxme]{coxme}},
\code{\link[coxme]{lmekin}} 
}
\examples{
##modified example from ?glmer
\dontrun{
if(require(lme4)) {
##create proportion of incidence
cbpp$prop <- cbpp$incidence/cbpp$size
gm1 <- glmer(prop ~ period + (1 | herd), family = binomial,
             weights = size, data = cbpp)
##print summary
summary(gm1)
##extract variance-covariance matrix of fixed effects
vcov(gm1)
##extract SE's of fixed effects - no labels
sqrt(diag(vcov(gm1))) #no labels
extractSE(gm1)  #with labels
detach(package:lme4)
}
}
}

\keyword{models}

\name{fam.link.mer}
\Rdversion{1.1}
\alias{fam.link.mer}
\title{
Extract Distribution Family and Link Function
}
\description{
This function extracts the distribution family and link function a
generalized linear mixed model fit with 'glmer' or 'lmer'.
}
\usage{
fam.link.mer(mod)
}

\arguments{
  \item{mod}{
an object of 'mer' class resulting from the fit of 'glmer' or 'lmer'.
}
}
\details{
This is an extractor function that extracts the information from the 
  call of 'glmer' or 'lmer' and it is called by 'modavg.mer',
  'modavgpred.mer', and  'predictSE.mer'.
}
\value{
  'fam.link.mer' returns a list with the following components:
  
  \item{family}{the family of the distribution of the model.}
  \item{link}{the link function of the model.}
  \item{supp.link}{a character value indicating whether the link
  function used is supported by 'predictSE.mer' and 'modavgpred.mer'.}
}

\author{
Marc J. Mazerolle
}

\seealso{
\code{\link{modavg}}, \code{\link{modavgpred}},
\code{\link{predictSE.mer}}, \code{\link[lme4]{glmer}}, \code{\link[lme4]{lmer}}
}
\examples{
##modified example from glmer
require(lme4)
##create proportion of incidence
cbpp$prop <- cbpp$incidence/cbpp$size
gm1 <- glmer(prop ~ period + (1 | herd), family = binomial, weights =
size, data = cbpp)
fam.link.mer(gm1)
gm2 <- glmer(prop ~ period + (1 | herd), family = binomial(link =
"cloglog"), weights = size, data = cbpp)
fam.link.mer(gm2)

##example with linear mixed model with Orthodont data from
##Pinheiro and Bates (2000)
data(Orthodont, package = "nlme")
m1 <- lmer(distance ~ Sex + (1 | Subject), data = Orthodont, REML =
FALSE)
fam.link.mer(m1)
m2 <- glmer(distance ~ Sex + (1 | Subject), family = gaussian(link =
"log"), data = Orthodont, REML = FALSE)
fam.link.mer(m2)

detach(package:lme4)
}
\keyword{models}

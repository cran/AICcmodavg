\name{mult.comp}
\Rdversion{1.1}
\alias{mult.comp}
\alias{print.mult.comp}
\title{
Create Model Selection Tables based on Multiple Comparisons
}
\description{
This function is an alternative to traditional multiple comparison
tests in designed experiments. It creates a model selection table based
on different grouping patterns of a factor and computes model-averaged
predictions for each of the factor levels.  The current version works
with objects of 'glm', 'gls', 'lm','lme', 'mer', 'merMod', and 'rlm'
classes.
}

\usage{
mult.comp(mod, factor.id, letter.labels = TRUE, sort = TRUE, c.hat = 1,
          second.ord = TRUE, nobs = NULL, newdata = NULL,
          type = "response", uncond.se = "revised", gamdisp = NULL) 
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mod}{
a model of one of the above-mentioned classes that includes at least one
factor as an explanatory variable.
}
  \item{factor.id}{
the factor of interest, on which the groupings (multiple comparisons)
are based. The user must supply the name of the categorical variable
between quotes as it appears in the model formula.
}
  \item{letter.labels}{
logical. If TRUE, letters are used as labels to denote the grouping
structure.  If FALSE, numbers are used as group labels.
}
  \item{sort}{
logical.  If true, the model selection table is ranked according to the
(Q)AIC(c) values.
}
  \item{c.hat}{
value of overdispersion parameter (i.e., variance inflation factor) such
as that obtained from 'c_hat'.  Note that values of c.hat different from
1 are only appropriate for binomial GLM's with trials > 1 (i.e.,
success/trial or cbind(success, failure) syntax) or with Poisson
GLM's.  If c.hat > 1, 'AICc' will return the quasi-likelihood analogue
of the information criterion requested.  This option is not supported
for generalized linear mixed models of the 'mer' class.
}
  \item{second.ord}{
logical.  If TRUE, the function returns the second-order
Akaike information criterion (i.e., AICc).
}
  \item{nobs}{
this argument allows to specify a numeric value other than total sample
size to compute the AICc (i.e., 'nobs' defaults to total number of
observations).  This is relevant only for certain types of models such
as mixed models where sample size is not straightforward.  In
such cases, one might use total number of observations or number of
independent clusters (e.g., sites) as the value of 'nobs'.
}
\item{newdata}{
  a data frame with the same structure as that of the original data
  frame for which we want to make predictions.  This data frame should 
  hold all variables constant other than the 'factor.id' variable.  All
  levels of the 'factor.id' variables should be included in the
  'newdata' data frame to get model-averaged predictions for each
  level.  If NULL, model-averaged predictions are computed for each level
  of the 'factor.id' variable while the values of the other explanatory
  variables are taken from the first row of the original data set.
}
\item{type}{the scale of prediction requested, one of "response" or
  "link". The latter is only relevant for 'glm' and 'mer' classes.  Note
  that the value "terms" is not defined for 'mult.comp'.
}
\item{uncond.se}{either, "old", or "revised", specifying the equation
used to compute the unconditional standard error of a model-averaged
estimate.  With uncond.se = "old", computations are based on equation
4.9 of Burnham and Anderson (2002), which was the former way to compute
unconditional standard errors.  With uncond.se = "revised", equation
6.12 of Burnham and Anderson (2002) is used.  Anderson (2008, p. 111)
recommends use of the revised version for the computation of
unconditional standard errors and it is now the default.  Note that
versions of package AICcmodavg < 1.04 used the old method to compute
unconditional standard errors.
}
\item{gamdisp}{
  the value of the gamma dispersion parameter.
}
}
\details{
A number of pairwise comparison tests are available for traditional
experimental designs, some controlling for the experiment-wise error and
others for comparison-wise errors (Day and Quinn 1991).  With the advent
of information-theoretic approaches, there has been a need for methods
analogous to multiple comparison tests in a model selection
framework.  Dayton (1998) and Burnham et al. (2011) suggested using
different parameterizations or grouping patterns of a factor to perform
multiple comparisons with model selection.  As such, it is possible to
assess the support in favor of certain grouping patterns based on a
factor.

For example, a factor with three levels has four possible grouping
patterns: {abc} (all groups are different), {abb} (the first group
differs from the other two), {aab} (the first two groups differ from the
third), and {aaa} (all groups are equal).  'mult.comp' implements such
an approach by pooling groups of the factor variable in a model and
updating the model, for each grouping pattern possible.  The models are
ranked according to one of four information criteria (AIC, AICc, QAIC,
and QAICc), and the labels in the table correspond to the grouping
pattern.  Note that the factor levels are sorted according to their means
for the response variable before being assigned to a group.  The
function also returns model-averaged predictions and unconditional
standard errors for each level of the 'factor.id' variable based on the
support in favor of each model (i.e., grouping pattern).

The number of grouping patterns increases substantially with the number
of factor levels, as 2^(k - 1), where k is the number of factor levels.
'mult.comp' supports factors with a maximum of 6 levels.  Also note that
'mult.comp' does not handle models where the 'factor.id' variable is
involved in an interaction.  In such cases, one should create the
interaction variable manually before fitting the model (see examples).
}
\value{
'mult.comp' creates a list of class 'mult.comp' with the following
components: 

\item{factor.id}{the factor for which grouping patterns are investigated.}
\item{models}{a list with the output of each model representing a
  different grouping pattern for the factor of interest.}
\item{model.names}{a vector of model names denoting the grouping pattern
  for each level of the factor.}
\item{model.table}{the model selection table for the models
  corresponding to each grouping pattern for the factor of interest.}
\item{ordered.levels}{the levels of the factor ordered according to the
  mean of the response variable. The grouping patterns (and model names)
  in the model selection table are based on the same order.}
\item{model.avg.est}{a matrix with the model-averaged prediction and
  unconditional standard error for each level of the factor.}
}
\references{
Burnham, K. P., Anderson, D. R., Huyvaert, K. P. (2011) AIC model
selection and multimodel inference in behaviorial ecology: some
background, observations and comparisons. \emph{Behavioral Ecology and
  Sociobiology} \bold{65} 23--25. 

Day, R. W., Quinn, G. P. (1989) Comparisons of treatments after an
analysis of variance in ecology. \emph{Ecological Monographs}
\bold{59} 433--463. 

Dayton, C. M. (1998) Information criteria for the paired-comparisons
problem. \emph{American Statistician} \bold{52} 144--151.
}
\author{
Marc J. Mazerolle
}
\seealso{
\code{\link{aictab}}, \code{\link{confset}}, \code{\link{c_hat}},
\code{\link{evidence}}, \code{\link[multcomp]{glht}}, \code{\link[gmodels]{fit.contrast}} 
}
\examples{
##one-way ANOVA example
data(turkey)

##convert diet to factor
turkey$Diet <- as.factor(turkey$Diet)
##run one-way ANOVA
m.aov <- lm(Weight.gain ~ Diet, data = turkey)

##compute models with different grouping patterns
##and also compute model-averaged group means
out <- mult.comp(m.aov, factor.id = "Diet")
##look at results
out

##look at grouping structure of a given model
##and compare with original variable
cbind(model.frame(out$models[[2]]), turkey$Diet)

##evidence ratio
evidence(out$model.table)


\dontrun{
##two-way ANOVA with interaction
data(calcium)

m.aov2 <- lm(Calcium ~ Hormone + Sex + Hormone:Sex, data = calcium)

##multiple comparisons
mult.comp(m.aov2, factor.id = "Hormone")
##returns an error because 'Hormone' factor is
##involved in an interaction

##create interaction variable
calcium$inter <- interaction(calcium$Hormone, calcium$Sex)

##run model with interaction
m.aov.inter <- lm(Calcium ~ inter, data = calcium)

##compare both
logLik(m.aov2)
logLik(m.aov.inter)
##both are identical

##multiple comparisons
mult.comp(m.aov.inter, factor.id = "inter")
}


\dontrun{
##Poisson regression
##example from ?glm
##Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
print(d.AD <- data.frame(treatment, outcome, counts))
glm.D93 <- glm(counts ~ outcome + treatment, data = d.AD, family = poisson)

mult.comp(mod = glm.D93, factor.id = "outcome")
}

\dontrun{
##Poisson regression
##example from ?glm
##Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
print(d.AD <- data.frame(treatment, outcome, counts))
glm.D93 <- glm(counts ~ outcome + treatment, data = d.AD, family = poisson)

mult.comp(mod = glm.D93, factor.id = "outcome")
}



\dontrun{
##example specifying 'newdata'
data(dry.frog)
m1 <- lm(log_Mass_lost ~ Shade + Substrate +
      cent_Initial_mass + Initial_mass2,
      data = dry.frog)

mult.comp(m1, factor.id = "Substrate",
          newdata = data.frame(
            Substrate = c("PEAT", "SOIL", "SPHAGNUM"),
            Shade = 0, cent_Initial_mass = 0,
            Initial_mass2 = 0))
}




}
\keyword{models}

\name{mult.comp}
\Rdversion{1.1}
\alias{mult.comp}
\alias{print.mult.comp}
\title{
Create Model Selection Tables based on Multiple Comparisons
}
\description{
This function is an alternative to traditional multiple comparison
tests in designed experiments. It creates a list of models with
parameterizations based on different groupings of the levels of a
factor. The function then creates a model selection table and ranks the
models according to one of the following information criteria: AIC,
AICc, QAIC, QAICc.  The table also provides delta AIC and Akaike
weights. The current version works with objects of 'glm', 'gls', 'lm',
'lme', 'mer', 'merMod', and 'rlm' classes.
}
\usage{
mult.comp(mod, factor.id, letter.labels = TRUE, sort = TRUE, c.hat = 1,
          second.ord = TRUE, nobs = NULL)
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
}
\details{
A number of pairwise comparison tests are available for traditional
experimental designs, some controlling for the experiment-wise error and
others for comparison-wise errors (Day and Quinn 1991).  With the advent
of information-theoretic approaches, there has been a need for methods
analogous to multiple comparison tests in a model selection
framework. Dayton (1998) and Burnham et al. (2011) suggested using
different parameterizations or grouping patterns of a factor to perform
multiple comparisons with model selection. As such, it is possible to
assess the support in favor of certain grouping patterns based on a
factor.

For example, a factor with three levels has four possible grouping
patterns: {abc} (all groups are different), {abb} (the first group
differs from the other two), {aab} (the first two groups differ from the
third), and {aaa} (all groups are equal).  'mult.comp' implements such
an approach by pooling groups of the factor variable in a model and
updating the model, for each grouping pattern possible. The models are
ranked according to one of four information criteria (AIC, AICc, QAIC,
and QAICc), and the labels in the table correspond to the grouping
pattern.  Note that the factor levels are sorted according to their means
for the response variable before being assigned to a group.
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
}
\references{
Burnham, K. P., Anderson, D. R., Huyvaert, K. P. (2011) AIC model
selection and multimodel inference in behaviorial ecology: some
background, observations and comparisons. \emph{Behavioral Ecology and
  Sociobiology}, \bold{65}, 23--25. 

Day, R. W., Quinn, G. P. (1989) Comparisons of treatments after an
analysis of variance in ecology. \emph{Ecological Monographs},
\bold{59}, 433--463. 

Dayton, C. M. (1998) Information criteria for the paired-comparisons
problem. \emph{American Statistician}, \bold{52}, 144--151.
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
set.seed(seed = 10)
ctl <- rnorm(n = 15, 2.3, sd = 0.3)
trt1 <- rnorm(n = 15, 2.5, sd = 0.3)
trt2 <- rnorm(n = 15, 2.1, sd = 0.3)

##combine data
dats <- data.frame(Y = c(ctl, trt1, trt2),
                   Treat = c(rep("ctl", 15),
                     rep("trt1", 15),
                     rep("trt2", 15)))

##basic model
m1 <- lm(Y ~ Treat, data = dats)

##do multiple comparisons based on grouping structures
out <- mult.comp(mod = m1, factor.id = "Treat")

##look at results
out

##look at grouping structure of a given model
##and compare with original variable
cbind(model.frame(out$models[[2]]), dats$Treat)

##evidence ratio
evidence(out$model.table)


##example from ?glm
##Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
print(d.AD <- data.frame(treatment, outcome, counts))
glm.D93 <- glm(counts ~ outcome + treatment, data = d.AD, family = poisson)

mult.comp(mod = glm.D93, factor.id = "outcome")
}
\keyword{models}

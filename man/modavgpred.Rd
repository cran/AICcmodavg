\name{modavgpred}
\Rdversion{1.1}
\alias{modavgpred}
\alias{print.modavgpred}
\title{
Compute Model-averaged Predictions
}
\description{
This function computes the model-averaged predictions and unconditional
standard errors based on the entire candidate model set.  The function
is currently implemented for 'lm' and 'glm' object classes that are
stored in a list. 
}
\usage{
modavgpred(cand.set, modnames, newdata, type = "response", c.hat = 1,
gamdisp = NULL, second.ord = TRUE, nobs = NULL, uncond.se = "revised")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cand.set}{
a list storing each of the models in the candidate model set
}
  \item{modnames}{
a character vector of model names to facilitate the identification of
each model in the model selection table
}
  \item{newdata}{
a data.frame with the same structure as that of the original data.frame
for which we want to make predictions
}
  \item{type}{
the scale of prediction requested, one of "response" or "link". Note
that the value "terms" is not defined for 'modavgpred').
}
  \item{c.hat}{
value of overdispersion parameter (i.e., variance inflation factor) such
as that obtained from 'c_hat'.  Note that values of c.hat different from
1 are only appropriate for binomial GLM's with trials > 1 (i.e.,
success/trial or cbind(success, failure) syntax) or with Poisson GLM's.
If c.hat > 1, 'aictab' will return the quasi-likelihood analogue of the
information criteria requested and multiply the variance-covariance
matrix of the estimates by this value (i.e., SE's are multiplied by
sqrt(c.hat)). 
}
\item{gamdisp}{
  the value of the gamma dispersion parameter
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
}
\details{
Unconditional confidence intervals computed by 'modavgpred' assume
asymptotic normality of the estimator, and are appropriate to estimate
confidence intervals of beta estimates (or estimates on the linear 
predictor scale).  For predictions of some types of response variables
(e.g., discrete values such as counts), the normal approximation may be
inappropriate.  In such cases, it is probably better to compute the
confidence intervals on the linear predictor scale and then
back-transform the limits to the scale of the response variable.
Burnham and Anderson (2002, p. 164) suggest alternative methods of
computing unconditional confidence intervals for small degrees of
freedom with profile likelihood intervals or bootstrapping.  
}
\value{
'modavgpred' returns an object of class 'modavgpred' with the
following components:
\item{type}{
  the scale of predicted values (response or link).
}
\item{mod.avg.pred}{
  the model-averaged prediction over the entire candidate model set.
}
\item{uncond.se}{
  the unconditional standard error of each model-averaged prediction.
}

}
\references{
Anderson, D. R. (2008) \emph{Model-based Inference in the Life Sciences:
  a primer on evidence}. Springer: New York.

Burnham, K. P., Anderson, D. R. (2002) \emph{Model Selection and
  Multimodel Inference: a practical information-theoretic
  approach}. Second edition. Springer: New York.  

}
\author{
Marc J. Mazerolle
}

\seealso{
\code{\link{AICc}}, \code{\link{modavg}},  \code{\link{predict}},
\code{\link{aictab}}, \code{\link{importance}},  \code{\link{c_hat}},
\code{\link{evidence}}, \code{\link{confset}} 
}

\examples{
#example from subset of models in Table 1 in Mazerolle (2006)
data(dry.frog)

Cand.models <- list( )
Cand.models[[1]] <- lm(log_Mass_lost ~ Shade + Substrate +
cent_Initial_mass + Initial_mass2, data = dry.frog)
Cand.models[[2]] <- lm(log_Mass_lost ~ Shade + Substrate +
cent_Initial_mass + Initial_mass2 + Shade:Substrate, data = dry.frog)
Cand.models[[3]] <- lm(log_Mass_lost ~ cent_Initial_mass +
Initial_mass2, data = dry.frog)
Cand.models[[4]] <- lm(log_Mass_lost ~ Shade + cent_Initial_mass +
Initial_mass2, data = dry.frog)
Cand.models[[4]] <- lm(log_Mass_lost ~ Shade + cent_Initial_mass +
Initial_mass2, data = dry.frog)
Cand.models[[5]] <- lm(log_Mass_lost ~ Substrate + cent_Initial_mass +
Initial_mass2, data = dry.frog)

#setup model names
Modnames <- NULL
for (i in 1:length(Cand.models)) {
Modnames[i] <- paste("mod", i, sep = "")
}

#compute model-averaged value and unconditional SE of predicted log of
#mass lost for frogs of average mass in shade for each substrate type

#first create data set to use for predictions
new.dat <- data.frame(Shade = c(1, 1, 1), cent_Initial_mass = c(0, 0,
0), Initial_mass2 = c(0, 0, 0), Substrate = c("SOIL", "SPHAGNUM",
"PEAT")) 

#compare unconditional SE's using both methods
modavgpred(cand.set = Cand.models, modnames = Modnames, newdata =
new.dat, type = "response", uncond.se = "old")
modavgpred(cand.set = Cand.models, modnames = Modnames, newdata =
new.dat, type = "response", uncond.se = "revised")
#round to 4 digits after decimal point
print(modavgpred(cand.set = Cand.models, modnames = Modnames, newdata =
new.dat, type = "response", uncond.se = "revised"), digits = 4)

#Gamma glm
#clotting data example from gamma.shape( ) in MASS package of
#Venables and Ripley (2002, Modern applied statistics with
#S. Springer-Verlag: New York.)
clotting <- data.frame(u = c(5, 10, 15, 20, 30, 40, 60, 80, 100),
         lot1 = c(118, 58, 42, 35, 27, 25, 21, 19, 18),
         lot2 = c(69, 35, 26, 21, 18, 16, 13, 12, 12))
clot1 <- glm(lot1 ~ log(u), data = clotting, family = Gamma)

library(MASS)
gamma.dispersion(clot1) #dispersion parameter
gamma.shape(clot1) #reciprocal of dispersion parameter ==
#shape parameter 
summary(clot1, dispersion = gamma.dispersion(clot1))  #better

#create list with models
Cand <- list()
Cand[[1]] <- glm(lot1 ~ log(u), data = clotting, family = Gamma)
Cand[[2]] <- glm(lot1 ~ 1, data = clotting, family = Gamma)

#create vector of model names
Modnames <- NULL
for (i in 1:length(Cand)) {
Modnames[i] <- paste("mod", i, sep = "")
}

#compute model-averaged predictions on scale of response variable for
#all observations
modavgpred(cand.set = Cand, modnames = Modnames, newdata = clotting,
gamdisp = gamma.dispersion(clot1), type = "response") 

#compute model-averaged predictions on scale of linear predictor
modavgpred(cand.set = Cand, modnames = Modnames, newdata = clotting,
gamdisp = gamma.dispersion(clot1), type = "link")

#compute model-averaged predictions on scale of linear predictor
\dontrun{
modavgpred(cand.set = Cand, modnames = Modnames, newdata = clotting,
gamdisp = gamma.dispersion(clot1), type = "terms") #returns an error
#because it is not defined for modavgpred
}
}
\keyword{models}
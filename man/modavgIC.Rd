\name{modavgIC}
\Rdversion{1.1}
\alias{modavgIC}
\alias{print.modavgIC}

\title{
Compute Model-averaged Parameter Estimate from User-supplied Information
Criterion
}
\description{
This function model-averages the estimate of a parameter of interest
among a set of candidate models, and computes the unconditional standard
error and unconditional confidence intervals as described in Buckland et
al. (1997) and Burnham and Anderson (2002).  Computations are based on
the values of the information criterion supplied manually by the user.
}
\usage{
modavgIC(ic, K, modnames = NULL, estimate, se, uncond.se = "revised",
         conf.level = 0.95, ic.name = NULL) 
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ic}{
    a vector of information criterion values for each model in the candidate
    model set.
  }
  \item{K}{
    a vector containing the number of estimated parameters for each
    model in the candidate model set.
  }
  \item{modnames}{
    a character vector of model names to identify each model in the
    model selection table.  If \code{NULL}, generic names (e.g.,
    \code{Mod1}, \code{Mod2}) are supplied in the table in the same
    order as the information criterion values.
  }
  \item{estimate}{
    a vector of estimates for each of the models in the candidate model
    set.  Estimates can be either beta estimates for a parameter of
    interest or a single prediction from each model.
  }
  \item{se}{
    a vector of standard errors for each of the estimates appearing in the
    \code{estimate} vector.
  }  
  \item{uncond.se}{
    either, \code{"old"}, or \code{"revised"}, specifying the equation
    used to compute the unconditional standard error of a model-averaged
    estimate.  With \code{uncond.se = "old"}, computations are based on
    equation 4.9 of Burnham and Anderson (2002), which was the former way
    to compute unconditional standard errors.  With \code{uncond.se =
      "revised"}, equation 6.12 of Burnham and Anderson (2002) is used.
    Anderson (2008, p. 111) recommends use of the revised version for the
    computation of unconditional standard errors and it is now the
    default.
  }
  \item{conf.level}{
    the confidence level (\eqn{1 - \alpha}) requested for the computation of
    unconditional confidence intervals.
  }
  \item{ic.name}{
   a character string denoting the name of the information criterion
   input by the user.  This character string will appear in certain
   column labels of the model selection table.
 }
}

\details{
  \code{modavgIC} computes a model-averaged estimate from the vector
  of parameter estimates specified in \code{estimate}.  Estimates and
  their associated standard errors must be specified in the same order
  as the values of the information criterion, the number of estimated
  parameters, and the model names.  Estimates provided may be for a
  parameter of interest (i.e., beta estimates) or predictions from each
  model.  This function is most useful for information criterion other
  than AIC, AICc, QAIC, and QAICc (e.g., WAIC: Watanabe 2010) or for
  classes not supported by \code{modavg}, \code{modavgCustom}, or
  \code{modavgPred}.
}

\value{
\code{modavgIC} creates an object of class \code{modavgIC} with
the following components: 

\item{Mod.avg.table}{the model selection table.}
\item{Mod.avg.est}{the model-averaged estimate.}
\item{Uncond.SE}{the unconditional standard error for the model-averaged
  estimate.}
\item{Conf.level}{the confidence level used to compute the confidence
  interval.} 
\item{Lower.CL}{the lower confidence limit.}
\item{Upper.CL}{the upper confidence limit.}
}
\references{
Anderson, D. R. (2008) \emph{Model-based Inference in the Life Sciences:
a primer on evidence}. Springer: New York.

Buckland, S. T., Burnham, K. P., Augustin, N. H. (1997) Model selection:
an integral part of inference. \emph{Biometrics} \bold{53}, 603--618.  

Burnham, K. P., Anderson, D. R. (2002) \emph{Model Selection and
Multimodel Inference: a practical information-theoretic
approach}. Second edition. Springer: New York.

Watanabe, S. (2010) Asymptotic equivalence of Bayes cross validation and
widely applicable information criterion in singular learning
theory. \emph{Journal of Machine Learning Research} \bold{11},
3571--3594. 
}
\author{
Marc J. Mazerolle
}

\seealso{
\code{\link{aictabCustom}}, \code{\link{ictab}},
\code{\link{modavg}}, \code{\link{modavgCustom}},
\code{\link{modavgShrink}}, \code{\link{modavgPred}}
}
\examples{
\dontrun{
##model averaging parameter estimate based on WAIC
##create a vector of names to trace back models in set
Modnames <- c("global model", "interactive model",
              "additive model", "invertpred model")

##WAIC values
waic <- c(105.74, 107.36, 108.24, 100.57)
##number of effective parameters
effK <- c(7.45, 5.61, 6.14, 6.05)

##vector of predictions
Preds <- c(0.106, 0.137, 0.067, 0.050)
##vector of SE's for prediction
Ses <- c(0.128, 0.159, 0.054, 0.039)

##compute model-averaged estimate and unconditional SE based on WAIC
modavgIC(ic = waic, K = effK, modnames = Modnames, 
         estimate = Preds, se = Ses,
         ic.name = "WAIC")
}
}
\keyword{models}

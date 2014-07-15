\name{predictSE}
\Rdversion{1.1}
\alias{predictSE}
\alias{predictSE.default}
\alias{predictSE.gls}
\alias{predictSE.lme}
\alias{predictSE.mer}
\alias{predictSE.merMod}
\alias{predictSE.unmarkedFitPCount}
\alias{predictSE.unmarkedFitPCO}

\title{
Computing Predicted Values and Standard Errors
}
\description{
Function to compute predicted values based on linear predictor and
associated standard errors from various fitted models.
}
\usage{
predictSE(mod, newdata, se.fit = TRUE, print.matrix = FALSE, \dots)

\method{predictSE}{gls}(mod, newdata, se.fit = TRUE, print.matrix = 
        FALSE, \dots) 

\method{predictSE}{lme}(mod, newdata, se.fit = TRUE, print.matrix =
        FALSE, level = 0, \dots)

\method{predictSE}{mer}(mod, newdata, se.fit = TRUE, print.matrix =
        FALSE, level = 0, type = "response", \dots)

\method{predictSE}{merMod}(mod, newdata, se.fit = TRUE, print.matrix =
        FALSE, level = 0, type = "response", \dots)

\method{predictSE}{unmarkedFitPCount}(mod, newdata, se.fit = TRUE,
        print.matrix = FALSE, type = "response", c.hat = 1, parm.type =
        "lambda", \dots)

\method{predictSE}{unmarkedFitPCO}(mod, newdata, se.fit = TRUE, 
        print.matrix = FALSE, type = "response", c.hat = 1,
        parm.type = "lambda", \dots)
}
\arguments{
  \item{mod}{
    an object of class \code{gls}, \code{lme}, \code{mer},
        \code{merMod}, \code{unmarkedFitPCount}, or
        \code{unmarkedFitPCO} containing the output of a model.
  }
  \item{newdata}{
    a data frame with the same structure as that of the original data
    frame for which we want to make predictions.
  }
  \item{se.fit}{
    logical.  If \code{TRUE}, compute standard errors on predictions.
  }

  \item{print.matrix}{
    logical.  If \code{TRUE}, the output is returned as a matrix, with
    predicted values and standard errors in columns.  If \code{FALSE},
    the output is returned as a list.
  }
  \item{level}{ 
    the level for which predicted values and standard errors are to be
    computed.  The current version of the function only supports
    predictions for the populations excluding random effects (i.e.,
    \code{level = 0}).
  }
  \item{type}{
    specifies the type of prediction requested.  This argument can take
    the value \code{response} or \code{link}, for predictions on the
    scale of the response variable or on the scale of the linear
    predictor, respectively.
  }
  \item{c.hat}{
    value of overdispersion parameter (i.e., variance inflation factor)
    such as that obtained from \code{Nmix.gof.test}.  If \code{c.hat >
    1}, \code{predictSE} will multiply the variance-covariance matrix of
    the predictions by this value (i.e., SE's are multiplied by
    \code{sqrt(c.hat)}).  High values of \code{c.hat} (e.g., \code{c.hat
    > 4}) may indicate that model structure is inappropriate.
  }
  \item{parm.type}{ 
    the parameter for which predictions are made based on the
    \emph{N}-mixture model of class \code{unmarkedFitPCount} or
    \code{unmarkedFitPCO} classes.
  }
  \item{\dots}{
    additional arguments passed to the function.
  }
}
\details{
  \code{predictSE} computes predicted values and associated standard
    errors.  Standard errors are approximated using the delta method
    (Oehlert 1992).  Predictions and standard errors for objects of
    \code{gls} class and mixed models of \code{lme}, \code{mer},
    \code{merMod} classes exclude the correlation or variance structure
    of the model.

  \code{predictSE} computes predicted values on abundance and standard
  errors based on the estimates from an \code{unmarkedFitPCount} or
  \code{unmarkedFitPCO} object.  Currently, only predictions on
  abundance (i.e., \code{parm.type = "lambda"}) with the zero-inflated
  Poisson distribution is supported.  For other parameters or
  distributions for models of \code{unmarkedFit} classes, use
  \code{predict} from the \code{unmarked} package. 
}
\value{
  \code{predictSE} returns requested values either as a matrix
  (\code{print.matrix = TRUE}) or list (\code{print.matrix = FALSE})
  with components:
  
    \item{fit}{
      the predicted values.
    }
    \item{se.fit}{
      the standard errors of the predicted values (if \code{se.fit = TRUE}). 
    }
  }

\note{
  For standard errors with better properties, especially for small
  samples, one can opt for simulations (see Gelman and Hill 2007), or
  nonparametric bootstrap (Efron and Tibshirani 1998).
}

\references{
  Efron, B., Tibshirani, R. J. (1998) \emph{An Introduction to the
    Bootstrap}. Chapman & Hall/CRC: New York. 
  
  Gelman, A., Hill, J. (2007) \emph{Data Analysis Using Regression and
    Multilevel/Hierarchical Models}. Cambridge University Press: New York. 
  
  Oehlert, G. W. (1992) A note on the delta method. \emph{American
    Statistician} \bold{46}, 27--29.

}
\author{
Marc J. Mazerolle
}

\seealso{
\code{\link[nlme]{gls}}, \code{\link[nlme]{lme}}, \code{\link[lme4]{glmer}},
\code{\link[lme4]{simulate.merMod}}, \code{\link[boot]{boot}},
\code{\link[unmarked]{parboot}}, \code{\link[unmarked]{nonparboot}},
\code{\link[unmarked]{pcount}}, \code{\link[unmarked]{pcountOpen}},
\code{\link[unmarked]{unmarkedFit-class}} 
}
\examples{
##Orthodont data from Pinheiro and Bates (2000) revisited
require(nlme)
m1 <- gls(distance ~ age, correlation = corCompSymm(value = 0.5, form = ~ 1 | Subject),
          data = Orthodont, method= "ML")
\dontrun{
##compare against lme fit
logLik(m1)
logLik(lme(distance ~ age, random = ~1 | Subject, data = Orthodont,
          method= "ML"))
##both are identical
}

##compute predictions and SE's for different ages
predictSE(m1, newdata = data.frame(age = c(8, 10, 12, 14)))

\dontrun{
##contagious bovine pleuropneumonia example modified from lme4
require(lme4)
data(cbpp)
##create proportion of incidence
cbpp$prop <- cbpp$incidence/cbpp$size
gm1 <- glmer(prop ~ period + (1 | herd), family = binomial,
             weights = size, data = cbpp)

##create a data set to make predictions
newherd<- data.frame(period = as.factor(c("1", "2", "3", "4")))

##predictions on logit link scale
predictSE(mod = gm1, newdata = newherd, se.fit = TRUE,
              type = "link", level = 0, print.matrix = FALSE)

##predictions on scale of original response variable
predictSE(mod = gm1, newdata = newherd, se.fit = TRUE,
          print.matrix = TRUE, level = 0, type = "response")
}

\dontrun{
if(require(unmarked)) {
##example with mallard data set from unmarked package
data(mallard)
mallardUMF <- unmarkedFramePCount(mallard.y, siteCovs = mallard.site,
                                  obsCovs = mallard.obs)
##run model with zero-inflated Poisson abundance
fm.mall.one <- pcount(~ ivel + date  ~ length + forest, mallardUMF, K=30,
                      mixture = "ZIP")
##make prediction
predictSE(fm.mall.one, type = "response", parm.type = "lambda",
          newdata = data.frame(length = 0, forest = 0, elev = 0))
##compare against predict
predict(fm.mall.one, type = "state", backTransform = TRUE,
        newdata = data.frame(length = 0, forest = 0, elev = 0))

##add offset in model to scale abundance per transect length
fm.mall.off <- pcount(~ ivel + date  ~ forest + offset(length), mallardUMF, K=30,
                      mixture = "ZIP")
##make prediction
predictSE(fm.mall.off, type = "response", parm.type = "lambda",
          newdata = data.frame(length = 10, forest = 0, elev = 0))
##compare against predict
predict(fm.mall.off, type = "state", backTransform = TRUE,
        newdata = data.frame(length = 10, forest = 0, elev = 0))
detach(package:unmarked)
}
}

}
\keyword{models}
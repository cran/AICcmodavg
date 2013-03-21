\name{predictSE.zip}
\Rdversion{1.1}
\alias{predictSE.zip}
\title{
Computing Predicted Values and Standard Errors
}
\description{
Function to compute predicted values and standard errors from output of
N-mixture models for count data based on the zero-inflated Poisson of
'unmarkedFitPCount' or 'unmarkedFitPCO' classes.
}
\usage{
predictSE.zip(mod, newdata, se.fit = TRUE, parm.type = "lambda",
              type = "response", print.matrix = FALSE)
}
\arguments{
  \item{mod}{
    an object of class 'unmarkedFitPCount' or 'unmarkedFitPCO'
    containing the output of a model. 
  }
  \item{newdata}{
    a data frame with the same structure as that of the original data
    frame for which we want to make predictions.
  }
  \item{se.fit}{
    logical.  If TRUE, compute standard errors on predictions.
  }

  \item{parm.type}{ 
    the parameter for which predictions are made based on the
    model. Currently, only predictions on abundance
    (i.e., parm.type = "lambda") is supported.
  }
  \item{type}{
    the scale of predicted values (response or link). Currently, only
    type = "response" is supported.
}
  
  \item{print.matrix}{
    logical.  If TRUE, the output is returned as a matrix, with
    predicted values and standard errors in columns.  If FALSE, the
    output is returned as a list.
  }
}
\details{
  'predictSE.zip' computes predicted values on abundance based on the
    estimates from an 'unmarkedFitPCount' or 'unmarkedFitPCO' object as
    well as the associated standard errors.  Standard errors are
    approximated using the delta method (Oehlert 1992).
}
\value{
  'predictSE.zip' returns requested values either as a matrix
    ('print.matrix = TRUE') or list ('print.matrix = FALSE') with
    components:
    \item{fit}{
      the predicted values.
    }
    \item{se.fit}{
      the standard errors of the predicted values (if 'se.fit = TRUE'). 
    }
  }

\note{
  Alternatively, one can get standard errors on predictions by using a
  bootstrap approach. 
}

\references{

Oehlert, G. W. (1992) A note on the delta method. \emph{American
  Statistician} \bold{46}, 27--29.

}
\author{
Marc J. Mazerolle
}

\seealso{
\code{\link{modavgpred}}, \code{\link[unmarked]{pcount}},
\code{\link[unmarked]{pcountOpen}}, \code{\link[unmarked]{unmarkedFit-class}}
}
\examples{
library(unmarked)
##example with mallard data set from unmarked package
data(mallard)
mallardUMF <- unmarkedFramePCount(mallard.y, siteCovs = mallard.site,
                                  obsCovs = mallard.obs)
##run model with zero-inflated Poisson abundance
fm.mall.one <- pcount(~ ivel + date  ~ length + forest, mallardUMF, K=30,
                      mixture = "ZIP")
##make prediction
predictSE.zip(fm.mall.one, type = "response", parm.type = "lambda",
              newdata = data.frame(length = 0, forest = 0, elev = 0))
##compare against predict
predict(fm.mall.one, type = "state", backTransform = TRUE,
        newdata = data.frame(length = 0, forest = 0, elev = 0))
}
\keyword{models}
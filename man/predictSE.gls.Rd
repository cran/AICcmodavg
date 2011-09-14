\name{predictSE.gls}
\Rdversion{1.1}
\alias{predictSE.gls}
\title{
Computing Predicted Values and Standard Errors
}
\description{
Function to compute predicted values based on linear predictor and
associated standard errors from a generalized least squares model.
}
\usage{
predictSE.gls(mod, newdata, se.fit = TRUE, print.matrix = FALSE)
}
\arguments{
  \item{mod}{
    an object of class 'gls' containing the output of a model.
  }
  \item{newdata}{
    a data frame with the same structure as that of the original data
    frame for which we want to make predictions.
  }
  \item{se.fit}{
    logical.  If TRUE, compute standard errors on predictions.
  }

  \item{print.matrix}{
    logical.  If TRUE, the output is returned as a matrix, with
    predicted values and standard errors in columns.  If FALSE, the
    output is returned as a list.
  }
}
\details{
  'predictSE.gls' computes predicted values based on the linear
    predictor and associated standard errors, excluding the correlation
    or variance structure of the model.  Standard errors are
    approximated using the delta method (Oehlert 1992). 
}
\value{
  'predictSE.gls' returns requested values either as a matrix
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
  For standard errors with better properties, especially for small
  samples, one can opt for simulations (see Gelman and Hill 2007).
}

\references{
Gelman, A., Hill, J. (2007) \emph{Data Analysis Using Regression and
  Multilevel/Hierarchical Models}. Cambridge University Press: New York. 
  
Oehlert, G. W. (1992) A note on the delta method. \emph{American
  Statistician} \bold{46}, 27--29.

}
\author{
Marc J. Mazerolle
}

\seealso{
\code{\link{gls}}, \code{\link{predict.gls}}
}
\examples{
##Orthodont data from Pinheiro and Bates (2000) revisited
require(nlme)
m1 <- gls(distance ~ age, correlation = corCompSymm(value = 0.5, form = ~ 1 | Subject), data = Orthodont,
          method= "ML")
\dontrun{
##compare against lme fit
  identical(logLik(m1),
  logLik(lme(distance ~ age, random = ~1 | Subject, data = Orthodont,
          method= "ML")))
  ##both are identical
}

##compute predictions and SE's for different ages
predictSE.gls(m1, newdata = data.frame(age = c(8, 10, 12, 14)))

}
\keyword{models}
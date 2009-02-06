\name{predictSE.lme}
\Rdversion{1.1}
\alias{predictSE.lme}
\title{
Computing Predicted Values and Standard Errors
}
\description{
Function to compute predicted values based on fixed effects (population
predictions) and associated standard errors from a linear mixed effect
model.
}
\usage{
predictSE.lme(mod, newdata, se.fit = TRUE, level = 0,
print.matrix = FALSE) 
}
\arguments{
  \item{mod}{
    an object of class 'lme' containing the output of a model.
  }
  \item{newdata}{
    a data frame with the same structure as that of the original data
    frame for which we want to make predictions.
  }
  \item{se.fit}{
    logical.  If TRUE, compute standard errors on predictions.
  }

  \item{level}{ 
    the level for which predicted values and standard errors are to be 
    computed.  The current version of the function only supports
    predictions for the populations excluding random effects (i.e.,
    level = 0).
  }
  \item{print.matrix}{
    logical.  If TRUE, the output is returned as a matrix, with
    predicted values and standard errors in columns.  If FALSE, the
    output is returned as a list.
  }
}
\details{
  'predictSE.lme' computes predicted values based on fixed effects and 
  associated standard errors.  Standard errors are approximated using
  the delta method (Oehlert 1992).
}
\value{
  'predictSE.lme' returns requested values either as a matrix
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
  To get the standard errors on predictions that include random effects,
  one can opt for simulations (see Gelman and Hill 2007).
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
\code{\link{lme}}, \code{\link{predict.lme}}
}
\examples{
##Orthodont data set in Pinheiro and Bates (2000)
require(nlme)
m1 <- lme(distance ~ age, random = ~1 | Subject, data = Orthodont,
method= "ML")
##create a new data frame to make predictions
newOrth <- data.frame(Sex = c("Male","Male","Female","Female"),
                           age = c(15, 20, 15, 20),
                           Subject = c("M01","M01","F30","F30"))
predictSE.lme(m1, newdata = newOrth, level = 0) 
predict(m1, newdata = newOrth, level = 0) #compare against predict.lme( )

##because only 'level = 0' is supported, the grouping levels are not
##necessary for prediction  
newd <- data.frame(Sex = c("Male","Male","Female","Female"),
                           age = c(15, 20, 15, 20))
predictSE.lme(m1, newdata = newd, level = 0)
predictSE.lme(m1, newdata = newd, level = 0, print.matrix = TRUE)

\dontrun{
predictSE.lme(m1, newdata = newd, level = 1, print.matrix = TRUE)
##generates an error
}
}
\keyword{models}
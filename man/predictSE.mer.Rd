\name{predictSE.mer}
\Rdversion{1.1}
\alias{predictSE.mer}
\title{
Computing Predicted Values and Standard Errors
}
\description{
Function to compute predicted values based on fixed effects (population
predictions) and associated standard errors from a generalized linear mixed effect
model.
}
\usage{
predictSE.mer(mod, newdata, se.fit = TRUE, type = "response", level = 0, 
print.matrix = FALSE) 
}
\arguments{
  \item{mod}{
    an object of class 'mer' containing the output of a model.
  }
  \item{newdata}{
    a data frame with the same structure as that of the original data
    frame for which we want to make predictions.
  }
  \item{se.fit}{
    logical.  If TRUE, compute standard errors on predictions.
  }
  \item{type}{
    specifies the type of prediction requested.  This argument can take
    the value 'response' or 'link', for predictions on the scale of the
    response variable or on the scale of the linear predictor,
    respectively.
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
  'predictSE.mer' computes predicted values based on fixed effects and
  associated standard errors.  Standard errors are approximated using
  the delta method (Oehlert 1992). 
}
\value{
  'predictSE.mer' returns requested values either as a matrix
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
\code{\link{predictSE.lme}}, \code{\link[lme4]{mer-class}}, \code{\link[lme4]{lmer}}, \code{\link[lme4]{glmer}}
}
\examples{
##contagious bovine pleuropneumonia example modified from lme4
require(lme4)
data(cbpp)
##create proportion of incidence
cbpp$prop <- cbpp$incidence/cbpp$size
gm1 <- glmer(prop ~ period + (1 | herd), family = binomial, weights =
size, data = cbpp)

##create a data set to make predictions
newherd<- data.frame(period = as.factor(c("1", "2", "3", "4")))

##predictions on logit link scale
predictSE.mer(mod = gm1, newdata = newherd, se.fit = TRUE, type =
"link", level = 0, print.matrix = FALSE)

##predictions on scale of original response variable
predictSE.mer(mod = gm1, newdata = newherd, se.fit = TRUE, type =
"response", level = 0, print.matrix = TRUE)



##example with linear mixed model with Orthodont data from
##Pinheiro and Bates (2000)
data(Orthodont, package = "nlme")
m2 <- lmer(distance ~ Sex + (1 | Subject), data = Orthodont, REML =
FALSE)

##create data set for predictions for all combinations of Sex and 2 ages
neworth <- expand.grid(Sex = c("Male", "Female"), age = c(8, 10))

##compute predicted values
predictSE.mer(m2, newdata = neworth)
##the following yields the same answer because the model
##uses an identity link
predictSE.mer(m2, newdata = neworth, type = "link") 

detach(package:lme4)

\dontrun{
predictSE.mer(m2, newdata = neworth, level = 1, print.matrix = TRUE)
##generates an error
}
}
\keyword{models}
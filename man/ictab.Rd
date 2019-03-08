\name{ictab}
\Rdversion{1.1}
\alias{ictab}
\alias{print.ictab}

\title{
Create Model Selection Tables from User-supplied Information
Criterion
}
\description{
  This function creates a model selection table from information
  criterion values supplied by the user.  The table ranks the models
  based on the values of the information criterion and also displays
  delta values and information criterion weights.
}
\usage{
ictab(ic, K, modnames = NULL, sort = TRUE, ic.name = NULL)
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
  \item{sort}{
    logical.  If \code{TRUE}, the model selection table is ranked according
    to the values of the information criterion.
  }
  \item{ic.name}{
   a character string denoting the name of the information criterion
   input by the user.  This character string will appear in certain
   column labels of the model selection table.
 }
}
\details{
  \code{ictab} constructs a model selection table based on the
  information criterion values supplied by the user.  This function is
  most useful for information criterion other than AIC, AICc, QAIC, and
  QAICc (e.g., WAIC: Watanabe 2010) or for classes not supported by
  \code{aictab} or \code{bictab}. 
}

  \value{
    \code{ictab} creates an object of class \code{ictab} with the
    following components:

    \item{Modname}{the name of each model of the candidate model set.}
    \item{K }{the number of estimated parameters for each model.}
    \item{IC}{the values of the information criterion input by the
    user.  If a value for \code{ic.name} is provided, the table modifies
    the labels of the table.}
    \item{Delta_IC}{the delta information criterion component comparing
    each model to the top-ranked model.}
    \item{ModelLik }{the relative likelihood of the model given the
      data (exp(-0.5*delta[i])).  This is not to be confused with the
      likelihood of the parameters given the data.  The relative likelihood
      can then be normalized across all models to get the model probabilities.} 
    \item{ICWt}{the information criterion weights, also termed "model
      probabilities" sensu Burnham and Anderson (2002) and Anderson (2008).
      These measures indicate the level of support (i.e., weight of
      evidence) in favor of any given model being the most parsimonious
      among the candidate model set.}
    \item{Cum.Wt }{the cumulative information criterion weights.  These
      are only meaningful if results in table are sorted in decreasing
      order of Akaike weights (i.e., \code{sort = TRUE}).}
  }
\references{
Anderson, D. R. (2008) \emph{Model-based Inference in the Life Sciences:
a primer on evidence}. Springer: New York.

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
\code{\link{aictabCustom}}, \code{\link{confset}}, \code{\link{evidence}},
\code{\link{modavgCustom}}, \code{\link{modavgIC}}
}
\examples{
##create a vector of names to trace back models in set
Modnames <- c("global model", "interactive model",
              "additive model", "invertpred model")

##WAIC values
waic <- c(105.74, 107.36, 108.24, 100.57)
##number of effective parameters
effK <- c(7.45, 5.61, 6.14, 6.05)

##generate WAIC table
ictab(ic = waic, K = effK, modnames = Modnames,
      sort = TRUE, ic.name = "WAIC")  
}

\keyword{models}

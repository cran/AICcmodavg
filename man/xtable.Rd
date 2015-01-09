\name{xtable}
\Rdversion{1.1}
\alias{xtable.aictab}
\alias{xtable.boot.wt}
\alias{xtable.dictab}
\alias{xtable.mb.chisq}
\alias{xtable.modavg}
\alias{xtable.modavgShrink}
\alias{xtable.modavgPred}
\alias{xtable.modavgEffect}
\alias{xtable.multComp}

\title{
  Format Objects to LaTeX or HTML
}
\description{
  Functions to format various objects following model selection and
  multimodel inference to LaTeX or HTML tables.  These functions extend the
  methods from the \code{xtable} package (Dahl 2014).
}
\usage{
\method{xtable}{aictab}(x, caption = NULL, label = NULL, align = NULL, 
        digits = NULL, display = NULL, nice.names = TRUE,
        include.AICc = TRUE, include.LL = TRUE, include.Cum.Wt = FALSE,
        \dots)

\method{xtable}{boot.wt}(x, caption = NULL, label = NULL, align = NULL,
        digits = NULL, display = NULL, nice.names = TRUE,
        include.AICc = TRUE, include.AICcWt = FALSE, \dots)

\method{xtable}{dictab}(x, caption = NULL, label = NULL, align = NULL, 
        digits = NULL, display = NULL, nice.names = TRUE,
        include.DIC = TRUE, include.Cum.Wt = FALSE, \dots)

\method{xtable}{mb.chisq}(x, caption = NULL, label = NULL, align = NULL,
        digits = NULL, display = NULL, nice.names = TRUE,
        include.detection.histories = TRUE, \dots)

\method{xtable}{modavg}(x, caption = NULL, label = NULL, align = NULL,
        digits = NULL, display = NULL, nice.names = TRUE,
        print.table = FALSE, \dots)

\method{xtable}{modavgEffect}(x, caption = NULL, label = NULL,
        align = NULL, digits = NULL, display = NULL, nice.names = TRUE, 
        print.table = FALSE, \dots)

\method{xtable}{modavgPred}(x, caption = NULL, label = NULL,
        align = NULL, digits = NULL, display = NULL, nice.names = TRUE,  
        \dots)

\method{xtable}{modavgShrink}(x, caption = NULL, label = NULL,
        align = NULL, digits = NULL, display = NULL, nice.names = TRUE, 
        print.table = FALSE, \dots)

\method{xtable}{multComp}(x, caption = NULL, label = NULL,
        align = NULL, digits = NULL, display = NULL, nice.names = TRUE, 
        print.table = FALSE, \dots)
}

\arguments{
  \item{x}{
    an object of class \code{aictab}, \code{boot.wt}, \code{dictab}, 
    \code{mb.chisq}, \code{modavg}, \code{modavgEffect},
    \code{modavgPred}, \code{modavgShrink}, or \code{multComp} resulting
    from model selection or multimodel inference.
  }
  
  \item{caption}{
    a character vector of length 1 or 2 storing the caption or title of
    the table.  If the vector is of length 2, the second item is the
    short caption used when LaTeX generates a list of tables.  The
    default value is \code{NULL} and suppresses the caption.
  }
  
  \item{label}{ 
    a character vector storing the LaTeX label or HTML anchor.  The default
    value is \code{NULL} and suppresses the label. 
  }

  \item{align}{
    a character vector of length equal to the number of columns of the
    table specifying the alignment of the elements.  Note that the rownames
    are considered as an additional column and require an alignment value.
  }

  \item{digits}{
    a numeric vector of length one or equal to the number of columns in
    the table (including the rownames) specifying the number of digits to
    display in each column. 
  }
 
  \item{display}{
    a character vector of length equal to the number of columns
    (including the rownames) specifying the format of each column.  For
    example, use \code{s} for strings, \code{f} for numbers in the
    regular format, or \code{d} for integers. See \code{formatC} for
    additional possible values.
  }

  \item{nice.names}{
    logical.  If \code{TRUE}, column labels are modified to improve their
    appearance in the table.  If \code{FALSE}, simpler labels are used,
    or the ones supplied directly by the user in the object storing the
    output.
  }

  \item{include.AICc}{logical.  If \code{TRUE}, the column containing
    the information criterion of each model is printed in the table.
    If \code{FALSE}, the column is suppressed.
  }

  \item{include.DIC}{logical.  If \code{TRUE}, the column containing
    the deviance information criterion (DIC) of each model is printed in
    the table.  If \code{FALSE}, the column is suppressed.
  }
  
  \item{include.LL}{logical.  If \code{TRUE}, the column containing the
    log-likelihood of each model is printed in the table.  If
    \code{FALSE}, the column is suppressed.
  }

  \item{include.Cum.Wt}{logical.  If \code{TRUE}, the column containing
    the cumulative Akaike weights is printed in the table.  If
    \code{FALSE}, the column is suppressed. 
  }

  \item{include.AICcWt}{logical.  If \code{TRUE}, the column containing
    the Akaike weight of each model is printed in the table.  If
    \code{FALSE}, the column is suppressed.
  }

  \item{include.detection.histories}{logical.  If \code{TRUE}, the
    column containing detection histories is printed in the table.  If
    \code{FALSE}, the column is suppressed. 
  }
    
  \item{print.table}{logical.  If \code{TRUE}, the model selection table
    is printed and other sections of the output are suppressed (e.g.,
    model-averaged estimates).  If \code{FALSE}, the model selection
    table is suppressed and only the other portion of the output is
    printed in the table.
  }
  
\item{\dots}{
  additional arguments passed to the function.
  }
}

\details{
  \code{xtable} creates an object of the \code{xtable} class inheriting
  from the \code{data.frame} class.  This object can then be used with
  \code{print.xtable} for added flexibility such as suppressing row names,
  modifying caption placement, and format tables in LaTeX or HTML
  format.
}

\references{
  Dahl, D. B. (2014) xtable: Export tables to LaTeX or HTML. R package
  version 1.7-3. \url{http://CRAN.R-project.org/package=xtable}.
}

\author{
  Marc J. Mazerolle
}

\seealso{
\code{\link{aictab}}, \code{\link{boot.wt}}, \code{\link{dictab}},
\code{\link{formatC}}, \code{\link{mb.chisq}}, \code{\link{modavg}},
\code{\link{modavgEffect}}, \code{\link{modavgPred}},
\code{\link{modavgShrink}}, \code{\link{multComp}},
\code{\link[xtable]{xtable}}, \code{\link[xtable]{print.xtable}}
}


\examples{
if(require(xtable)) {
##model selection
data(dry.frog)
##setup candidate models
Cand.models <- list( )
Cand.models[[1]] <- lm(log_Mass_lost ~ Shade + Substrate +
                       cent_Initial_mass + Initial_mass2,
                       data = dry.frog)
Cand.models[[2]] <- lm(log_Mass_lost ~ Shade + Substrate +
                       cent_Initial_mass + Initial_mass2 +
                       Shade:Substrate, data = dry.frog)
Cand.models[[3]] <- lm(log_Mass_lost ~ cent_Initial_mass +
                       Initial_mass2, data = dry.frog)
Model.names <- c("additive", "interaction", "no shade")

##model selection table
out <- aictab(cand.set = Cand.models, modnames = Model.names)

xtable(out)
##exclude AICc and LL
xtable(out, include.AICc = FALSE, include.LL = FALSE)
##remove row names and add caption
print(xtable(out, caption = "Model selection based on AICc"),
      include.rownames = FALSE, caption.placement = "top")


##model-averaged estimate of Initial_mass2
mavg.mass <- modavg(cand.set = Cand.models, parm = "Initial_mass2",
                    modnames = Model.names)
#model-averaged estimate
xtable(mavg.mass, print.table = FALSE)
#table with contribution of each model
xtable(mavg.mass, print.table = TRUE)  


##model-averaged predictions for first 10 observations
preds <- modavgPred(cand.set = Cand.models, modnames = Model.names,
                    newdata = dry.frog[1:10, ])
xtable(preds)
}
}
\keyword{models}

\name{xtable}
\Rdversion{1.1}
\alias{xtable.aictab}
\alias{xtable.anovaOD}
\alias{xtable.bictab}
\alias{xtable.boot.wt}
\alias{xtable.checkParms}
\alias{xtable.countDist}
\alias{xtable.countHist}
\alias{xtable.detHist}
\alias{xtable.detTime}
\alias{xtable.dictab}
\alias{xtable.ictab}
\alias{xtable.mb.chisq}
\alias{xtable.modavg}
\alias{xtable.modavgCustom}
\alias{xtable.modavgEffect}
\alias{xtable.modavgIC}
\alias{xtable.modavgPred}
\alias{xtable.modavgShrink}
\alias{xtable.multComp}
\alias{xtable.summaryOD}

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
        digits = NULL, display = NULL, auto = FALSE,
        nice.names = TRUE, include.AICc = TRUE,
        include.LL = TRUE, include.Cum.Wt = FALSE, \dots)

\method{xtable}{anovaOD}(x, caption = NULL, label = NULL, align = NULL,
       digits = NULL, display = NULL, auto = FALSE, 
       nice.names = TRUE, \dots)

\method{xtable}{bictab}(x, caption = NULL, label = NULL, align = NULL, 
        digits = NULL, display = NULL, auto = FALSE,
        nice.names = TRUE, include.BIC = TRUE, include.LL = TRUE,
        include.Cum.Wt = FALSE, \dots)

\method{xtable}{boot.wt}(x, caption = NULL, label = NULL, align = NULL,
        digits = NULL, display = NULL, auto = FALSE,
        nice.names = TRUE, include.AICc = TRUE, include.AICcWt = FALSE, \dots)

\method{xtable}{countDist}(x, caption = NULL, label = NULL,
        align = NULL, digits = NULL, display = NULL, auto = FALSE,
        nice.names = TRUE, table.countDist = "distance", \dots)

\method{xtable}{checkParms}(x, caption = NULL, label = NULL,
        align = NULL, digits = NULL, display = NULL, auto = FALSE,
        nice.names = TRUE, include.variable = TRUE, include.max.se =
        TRUE, include.n.high.se = TRUE, \dots)

\method{xtable}{countHist}(x, caption = NULL, label = NULL,
        align = NULL, digits = NULL, display = NULL, auto = FALSE,
        nice.names = TRUE, table.countHist = "count", \dots)

\method{xtable}{detHist}(x, caption = NULL, label = NULL,
        align = NULL, digits = NULL, display = NULL, auto = FALSE,
        nice.names = TRUE, table.detHist = "freq", \dots)

\method{xtable}{detTime}(x, caption = NULL, label = NULL,
        align = NULL, digits = NULL, display = NULL, auto = FALSE,
        nice.names = TRUE, table.detTime = "freq", \dots)

\method{xtable}{dictab}(x, caption = NULL, label = NULL, align = NULL, 
        digits = NULL, display = NULL, auto = FALSE, 
        nice.names = TRUE, include.DIC = TRUE, include.Cum.Wt = FALSE, \dots)

\method{xtable}{ictab}(x, caption = NULL, label = NULL, align = NULL,
        digits = NULL, display = NULL, auto = FALSE, nice.names = TRUE,
        include.IC = TRUE, include.Cum.Wt = FALSE, \dots)

\method{xtable}{mb.chisq}(x, caption = NULL, label = NULL, align = NULL,
        digits = NULL, display = NULL, auto = FALSE, nice.names = TRUE,
        include.detection.histories = TRUE, \dots)

\method{xtable}{modavg}(x, caption = NULL, label = NULL, align = NULL,
        digits = NULL, display = NULL, auto = FALSE, nice.names = TRUE,
        print.table = FALSE, \dots)

\method{xtable}{modavgCustom}(x, caption = NULL, label = NULL,
        align = NULL, digits = NULL, display = NULL, auto = FALSE,
        nice.names = TRUE, print.table = FALSE, \dots)

\method{xtable}{modavgEffect}(x, caption = NULL, label = NULL,
        align = NULL, digits = NULL, display = NULL, auto = FALSE,
        nice.names = TRUE, print.table = FALSE, \dots)

\method{xtable}{modavgIC}(x, caption = NULL, label = NULL, align = NULL,
        digits = NULL, display = NULL, auto = FALSE, nice.names = TRUE,
        print.table = FALSE, \dots)

\method{xtable}{modavgPred}(x, caption = NULL, label = NULL,
        align = NULL, digits = NULL, display = NULL, auto = FALSE,
        nice.names = TRUE, \dots)

\method{xtable}{modavgShrink}(x, caption = NULL, label = NULL,
        align = NULL, digits = NULL, display = NULL, auto = FALSE,
        nice.names = TRUE, print.table = FALSE, \dots)

\method{xtable}{multComp}(x, caption = NULL, label = NULL,
        align = NULL, digits = NULL, display = NULL, auto = FALSE,
        nice.names = TRUE, print.table = FALSE, \dots)

\method{xtable}{summaryOD}(x, caption = NULL, label = NULL, align = NULL,
                           digits = NULL, display = NULL, auto = FALSE, 
                           nice.names = TRUE, \dots)

}

\arguments{
  \item{x}{
    an object of class \code{aictab}, \code{anovaOD}, \code{bictab},
    \code{boot.wt}, \code{checkParms}, \code{countDist},
    \code{countHist}, \code{detHist}, \code{detTime}, \code{dictab},
    \code{ictab}, \code{mb.chisq}, \code{modavg}, \code{modavgEffect},
    \code{modavgCustom}, \code{modavgIC}, \code{modavgPred},
    \code{modavgShrink}, \code{multComp}, or \code{summaryOD}.
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
    regular format, or \code{d} for integers.  See \code{formatC} for
    additional possible values.
  }


  \item{auto}{
    Logical, indicating whether to apply automatic format when no value
    is passed to \code{align}, \code{digits}, or \code{display}. This
    \code{autoformat} (based on \code{xalign}, \code{xdigits}, and \code{xdisplay})
    can be useful to quickly format a typical \code{matrix} or
    \code{data.frame}. Default value is \code{FALSE}.
  }
  
  \item{nice.names}{
    logical.  If \code{TRUE}, column labels are modified to improve their
    appearance in the table.  If \code{FALSE}, simpler labels are used,
    or the ones supplied directly by the user in the object storing the
    output.
  }

  \item{include.AICc}{logical.  If \code{TRUE}, the column containing
    the information criterion (AIC, AICc, QAIC, or QAICc) of each model
    is printed in the table.  If \code{FALSE}, the column is
    suppressed.
  }

  \item{include.BIC}{logical.  If \code{TRUE}, the column containing
    the information criterion (BIC or QBIC) of each model is printed in
    the table.  If \code{FALSE}, the column is suppressed.
  }

  \item{include.DIC}{logical.  If \code{TRUE}, the column containing
    the deviance information criterion (DIC) of each model is printed in
    the table.  If \code{FALSE}, the column is suppressed.
  }
  
  \item{include.IC}{logical.  If \code{TRUE}, the column containing
    the information criterion of each model is printed in the table.  If
    \code{FALSE}, the column is suppressed.
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

  \item{include.variable}{logical.  If \code{TRUE}, the
    column containing the variable name is printed in the table.  If
    \code{FALSE}, the column is suppressed.
  }

    \item{include.max.se}{logical.  If \code{TRUE}, the
      column containing the maximum SE in the model is printed in the
      table.  If \code{FALSE}, the column is suppressed.
    }

    \item{include.n.high.se}{logical.  If \code{TRUE}, the
      column containing the number of SE's greater than the threshold 
      specified by the user is printed in the table.  If \code{FALSE},
      the column is suppressed.
    }
  
  \item{print.table}{logical.  If \code{TRUE}, the model selection table
    is printed and other sections of the output are suppressed (e.g.,
    model-averaged estimates).  If \code{FALSE}, the model selection
    table is suppressed and only the other portion of the output is
    printed in the table.
  }
  
  \item{table.detHist}{character string specifying, either \code{"freq"},
  \code{"prop"}, or \code{"hist"}.  If \code{table.type = "freq"}, the
  function returns a table of frequencies of sites sampled, of sites with
  at least one detection, and for data with multiple primary periods,
  the frequencies of sites with observed extinctions and colonizations.
  If \code{table.type = "prop"}, the table returns the proportion of
  sites with at least one detection, and for data with multiple periods,
  the proportion of sites with observed extinctions and colonizations.
  If \code{table.type = "hist"}, the function returns the frequencies of
  each observed detection history.
}

\item{table.detTime}{character string specifying, either \code{"freq"},
  \code{"prop"}, or \code{"dist"}.  If \code{table.type = "freq"}, the
  function returns a table of frequencies of sites sampled, of sites with
  at least one detection, and for data with multiple primary periods,
  the frequencies of sites with observed extinctions and colonizations.
  If \code{table.type = "prop"}, the table returns the proportion of
  sites with at least one detection, and for data with multiple periods,
  the proportion of sites with observed extinctions and colonizations.
  If \code{table.type = "dist"}, the function returns the quantiles of
  detection times.
}

\item{table.countDist}{character string specifying, either
  \code{"distance"}, \code{"count"}, \code{"freq"}, or \code{"prop"}.  If
  \code{table.type = "distance"}, the function returns a table of counts
  summarized for each distance class.  If \code{table.type = "count"},
  the function returns the table of frequencies of counts observed
  across sites.  If \code{table.type = "freq"}, the function returns a
  table of frequencies of sites sampled, of sites with at least one
  detection, and for data with multiple primary periods, the frequencies
  of sites with observed extinctions and colonizations.  If
  \code{table.type = "prop"}, the table returns the proportion of sites
  with at least one detection, and for data with multiple periods, the
  proportion of sites with observed extinctions and colonizations.
}

\item{table.countHist}{character string specifying, either
  \code{"count"}, \code{"freq"}, \code{"prop"}, or \code{"hist"}.  If
  \code{table.type = "count"}, the function returns the table of
  frequencies of counts observed across sites.  If \code{table.type =
  "freq"}, the function returns a table of frequencies of sites sampled,
  of sites with at least one detection, and for data with multiple
  primary periods, the frequencies of sites with observed extinctions
  and colonizations.  If \code{table.type = "prop"}, the table returns
  the proportion of sites with at least one detection, and for data with
  multiple periods, the proportion of sites with observed extinctions
  and colonizations.  If \code{table.type = "hist"}, the function
  returns the frequencies of each observed count history.
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
  version 1.7-3. \url{https://cran.r-project.org/package=xtable}.
}

\author{
  Marc J. Mazerolle
}

\seealso{
\code{\link{aictab}}, \code{\link{boot.wt}}, \code{\link{dictab}},
\code{\link{formatC}}, \code{\link{ictab}}, \code{\link{mb.chisq}},
\code{\link{modavg}}, \code{\link{modavgCustom}},
\code{\link{modavgIC}}, \code{\link{modavgEffect}},
\code{\link{modavgPred}}, \code{\link{modavgShrink}},
\code{\link{multComp}}, \code{\link{summaryOD}}, \code{\link{anovaOD}},
\code{\link[xtable]{xtable}}, \code{\link[xtable]{print.xtable}}
}


\examples{
if(require(xtable)) {
##model selection example
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

##model selection table - AICc
out <- aictab(cand.set = Cand.models, modnames = Model.names)

xtable(out)
##exclude AICc and LL
xtable(out, include.AICc = FALSE, include.LL = FALSE)
##remove row names and add caption
print(xtable(out, caption = "Model selection based on AICc"),
      include.rownames = FALSE, caption.placement = "top")


##model selection table - BIC
out2 <- bictab(cand.set = Cand.models, modnames = Model.names)

xtable(out2)
##exclude AICc and LL
xtable(out2, include.BIC = FALSE, include.LL = FALSE)
##remove row names and add caption
print(xtable(out2, caption = "Model selection based on BIC"),
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


##example of diagnostics
\dontrun{
if(require(unmarked)){
##distance sampling example from ?distsamp
data(linetran)
ltUMF <- with(linetran, {
              unmarkedFrameDS(y = cbind(dc1, dc2, dc3, dc4),
              siteCovs = data.frame(Length, area, habitat),
              dist.breaks = c(0, 5, 10, 15, 20),
              tlength = linetran$Length * 1000, survey = "line",
              unitsIn = "m")
              })

##summarize counts across distance classes
xtable(countDist(ltUMF), table.countDist = "distance")
##summarize counts across all sites
xtable(countDist(ltUMF), table.countDist = "count")

##Half-normal detection function
fm1 <- distsamp(~ 1 ~ 1, ltUMF)
##determine parameters with highest SE's
xtable(checkParms(fm1))
}
}
}
\keyword{models}

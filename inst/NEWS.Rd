\name{NEWS}
\title{AICcmodavg News}
\encoding{UTF-8}

\section{Changes in version 2.3-1 (2020-08-21)}{
  \itemize{
    \item corrected bug in margin width when plotting results of
    'mb.gof.test' with dynamic occupancy model models (thanks to
    M. Meredith)
    \item added 'ncores' argument to 'mb.gof.test' and 'Nmix.gof.test'
    to control number of cores used parametric bootstrap (suggested by
    M. Meredith)
    \item added arguments 'cex.axis', 'cex.lab', 'cex.main', and 'lwd'
    to add flexibility to various functions producing plots (i.e.,
    'mb.gof.test', 'Nmix.gof.test', 'detTime', 'countDist', 'countHist')
    \item created an NEWS.Rd file
    \item created a vignette for a tour of main functions of the package
  }
}



\section{Changes in version 2.3-0 (2020-07-03)}{
  \itemize{
    \item added support for models of classes 'unmarkedFitOccuMS',
    'unmarkedFitOccuTTD', 'unmarkedFitMMO, 'unmarkedFitDSO'
    \item created new function 'detTime' to summarize time to detection data
    for models of class 'unmarkedFitOccuTTD'
    \item addedd 'modavgPred' method for 'unmarkedFitOccuMulti' and
    'unmarkedFitOccuMS' objects
    \item modified default plotting parameters in functions 'countHist',
    'countDist', 'detTime', 'mb.gof.test', and 'Nmix.gof.test'
    \item modified 'mb.chisq' methods to allow for new default behavior
    of R to use stringsAsFactors = FALSE for characters in data frames
    \item correction bug for GOF Royle-Nichols model (thanks to Ken
    Kellner)
  }
}


\section{Changes in version 2.2-2 (2019-05-05)}{
  \itemize{
    \item added methods for 'lmerTest' objects
    \item corrected bug when model names are NULL in 'modavgEffect' for
    certain model types
  }
}



\section{Changes in version 2.2-1 (2019-03-07)}{
  \itemize{
    \item added 'DIC' and 'dictab' methods for objects of classes
    'jagsUI'
    \item added new functions 'ictab' and 'modavgIC' for generic
    information criteria, and added 'print' and 'xtable' methods for
    these classes
    \item modified 'evidence' to include 'ictab' objects
    \item modified sorting of models in 'aictab', 'bictab', 'dictab',
    'aictabCustom', and 'bictabCustom' based on values of delta AIC
    instead of Akaike weights to avoid rounding error with low values
    \item modified 'xtable' for nicer printing of K for 'aictab', 'bictab',
    and 'anovaOD'
    \item modified value of 'out.type' argument to 'confint' in
    'summaryOD'
    \item added 'xtable' method for 'modavgCustom' class and modified
    'xtable' method for 'modavg' when working with models of 'multinom'
    class
    \item added examples in 'AICcmodavg-package' Rd file
  }
}



\section{Changes in version 2.2-0 (2019-02-25)}{
  \itemize{
    \item added 'summaryOD' to compute confidence intervals around
    estimates or hypothesis tests adjusted for overdispersion (request
    by M. Kéry)
    \item added methods for objects of 'glm.nb' class
    \item added methods for glmmTMB models (request by Y. Paillet)
    \item added 'parallel' argument in 'Nmix.gof.test' and 'mb.gof.test'
    for new 'parallel' argument in 'parboot'
    \item included support of goodness-of-fit statistics of 'mb.chisq'
    and 'mb.gof.test' for Royle-Nichols models fit with 'occuRN' (thanks
    to Dan Linden)
    \item included support for multispecies occupancy data and models of
    classes 'unmarkedFrameOccuMulti' and 'unmarkedFitOccuMulti'
    \item added 'anovaOD' function to compute likelihood ratio tests
    corrected for overdispersion for models of various classes that use
    the Poisson or binomial distributions
    \item modified 'covDiag' to handle missing values in data set
    \item updated man page of 'modavgShrink' to include method for
    'coxme' class
    \item updated man page of 'importance' to correct typo (thanks to
    J. Carlisle)
    \item corrected bug in modavg for parameter "b" in false-positive
    occupancy models
    \item modified sort in 'aictab' and 'bictab' based on delta AIC
    instead of Akaike weight
    \item modified printing of NA to "." to save space and improve
    alignment in 'detHist', 'countDist', 'xtable', and 'mb.chisq'
    methods
    \item modified 'modavgPred' and 'modavgEffect' to use 'predict' for
    ZIP models fit with 'pcount'
  }
}



\section{Changes in version 2.1-1 (2017-06-19)}{
  \itemize{
    \item created 'extractX' generic function to extract all predictors
    that appear in a given list of candidate models and provide info on
    these variables (added imports of getData from nlme and unmarked,
    and of formula for lme4)
    \item corrected minor bug for unmarkedFitPCO model class with
    'modavgShrink' to modify label of gamma depending on H0 (gamConst,
    gamAR, gamTrend)
    \item cleaned up modavgEffect
    \item added type = "certain" and "falsepos" as values for the
    'parm.type' argument of 'modavgPred', 'modavgEffect',
    'modavgShrink', and 'modavg' methods for objects of
    'unmarkedFitOccuFP' class and added check for consistency of parameter
    occurrence in different models
    \item added "iota" as a value for the 'parm.type' argument for
    'modavgPred', 'modavgEffect', 'modavgShrink', and 'modavg' methods
    for objects of 'unmarkedFitPCO' class and added check for
    consistency of parameter occurrence in different models
    \item added "detect" as a value for the 'parm.type' argument of
    'modavgPred', 'modavgEffect', 'modavgShrink', and 'modavg' methods
    for objects of 'unmarkedFitDS' class
    \item added "detect" as a value for the 'parm.type' argument of
    'modavgPred', 'modavgEffect', 'modavgShrink', and 'modavg' methods
    for objects of 'unmarkedFitGDS' class
  }
}



\section{Changes in version 2.1-0 (2016-11-17)}{
  \itemize{
    \item corrected minor bug in storing lower.CL in 'modavgPred' result
    \item modified 'fam.link.mer' to account for negative binomial with
    'glmer.nb'
    \item corrected minor bug in rounding output of 'c_hat'
    \item corrected bug in 'countDist' with GDS objects including NA's
    \item modified ifelse( ) test when single site includes NA's in
    'mb.chisq.unmarkedFitOccu' and 'mb.chisq.unmarkedFitColExt' (thanks to
    D. Linden)
    \item added drop = FALSE to preserve matrix in
    'mb.chisq.unmarkedFitOccu' and 'mb.chisq.unmarkedFitColExt'
    \item added 'report' argument to print test statistic in terminal
    for functions 'mb.gof.test' and 'Nmix.gof.test'
    \item added 'checkConv' function to check convergence status of
    algorithm at termination for different model classes
    \item added 'checkParms' function to identify parameters with high
    SE's in different model classes
    \item modified gls method for 'predictSE' to extract model formula
    following suggestion by C. R. Andersen
    \item corrected bug in creating final plot of theoretical
    distribution of sums of chi-squares in mb.gof.test for dynamic
    occupancy models (thanks to M. Kéry)
    \item corrected bug in 'covDiag'
    \item added 'useBIC' and 'bictab' methods for different model
    objects based on Bayesian Information Criterion
    \item added 'useBICCustom' function to compute BIC or QBIC from
    user-supplied inputs (log-likelihood, K, n, c.hat)
    \item added 'bictabCustom' function to create model selection table
    using BIC or QBIC from user-supplied inputs
    \item added methods for 'xtable' for objects of 'checkParms' and
    'bictab' classes
    \item added functionality of 'evidence' with dictab and bictab
    objects
    \item added new 'useBIC' argument to 'modavgCustom' to conduct
    multimodel inference based on BIC from user-supplied inputs
    \item removed functions that were moved to defunct during previous
    cycle of package
  }
}



\section{Changes in version 2.0-4 (2016-02-10)}{
  \itemize{
    \item updated email address
    \item modified 'mb.gof.test' to deal with dynamic occupancy data for
    entire visits without data, and seasons with only a single visit
    (Thanks to G. Sadoti for comments)
    \item corrected bug for negative binomial models and dispersion
    argument in 'modavg', 'modavgShrink', 'modavgPred', 'modavgEffect',
    'multComp'
    \item corrected bug on parameter count when negative binomial
    formated to glm object with glm.convert( )
    \item updated functionality for overdispersion correction of ZIP
    N-mixture models in unmarked with 'modavgPred'
    \item added new estimators of overdispersion in 'c_hat' for GLM's
    (deviance, Farrington, and Fletcher) and print method
    \item added comment on help page of 'mb.gof.test' to clarify that
    the table only shows the observed detection histories, but that
    unobserved detection histories are included in the computations, and
    that the chi-square test statistic should be compared against
    bootstrap distribution because low expected frequencies invalidate
    chi-square assumption
    \item modified 'modavgEffect' to allow multiple variables to vary
    and print generic group names
    \item modified 'modavgShrink' to allow cases where variables do not
    appear with the same frequency in models, but issues a warning
    \item corrected bug in 'modavgShrink' to deal with 'unmarkedFitGDS'
    objects on lambda
    \item added support for model selection of 'lavaan' class objects
    \item added support for model selection of fitdistr class objects
    (package MASS)
    \item added support for model selection of fitdist class objects
    (package fitdistrplus)
    \item added functionality to request CI's (back-transformed when
    using type = "response") with 'modavgPred' and modified
    'xtable.modavgPred'
    \item added covariance diagnostic for lambda in N-mixture models
    based on Dennis et al. 2015 Biometrics 71: 237-246
    \item added 'detHist' function to summarize information on detection
    histories
    \item added 'countHist' function to summarize information on count
    data
    \item added 'countDist' function to summarize information on count
    data collected by distance sampling
    \item added 'xtable' methods for objects of 'detHist' and
    'countHist' classes
    \item added fat and body measurement data set
    \item added GPA data set and test scores
    \item moved deprecated functions to defunct
    \item cleaned up print method of functions
  }
}



\section{Changes in version 2.0-3 (2015-01-15)}{
  \itemize{
    \item modified CITATION file to use meta$Date and meta$Version as
    requested by CRAN
    \item added the Red-spotted newt CMR data set (Gill 1985)
  }
}



\section{Changes in version 2.0-2 (2015-01-05)}{
  \itemize{
    \item added support for c-hat argument for dynamic occupancy models
    fit with 'colext' and additional mixture models fit from 'distsamp',
    'gdistsamp', 'gmultmix', 'gpcount', 'multinomPois', and 'occuRN'
    \item made 'Nmix.chisq' and 'Nmix.gof.test' generic and added
    support for different N-mixture models and updated the Rd file
    \item made 'mb.chisq' and 'mb.gof.test' generic for unmarkedFitOccu
    and unmarkedFitColExt models, modified structure of output, added
    new arguments 'print.seasons' for 'mb.gof.test.unmarkedFitColExt',
    and updated the Rd file
    \item added new function to estimate c-hat in lme4 models and made
    'c_hat' function generic and updated the Rd file
    \item clarified examples in Rd file of 'modavg' to clearly
    illustrate how to exclude models when model-averaging the main
    effect of a variable also involved in an interaction
    \item corrected bug in 'boot.wt' when sorting models in ascending
    order of Akaike weight
    \item created a separate class for 'boot.wt' objects and print
    method to facilitate LaTeX formatting
    \item created methods for 'xtable' to format objects from 'aictab',
    'boot.wt', 'dictab', 'modavg', 'modavgShrink', 'modavgPred',
    'modavgEffect', 'multComp', 'mb.chisq' to LaTeX 
  }
}



\section{Changes in version 2.01 (2014-10-09)}{
  \itemize{
    \item corrected bugs in 'modavg' and 'modavgShrink' for models of
    class 'unmarkedFitGMM' and 'unmarkedFitGPC'
    \item corrected bugs in 'dictab' to return deviance of all models in
    model selection table
    \item added logical test to identify cases with only a single cohort
    without missing values for 'mb.chisq'
    \item added custom versions of 'AICc', 'aictabCustom', and
    'modavgCustom' when user manually supplies log-Likelihood, K,
    estimates, and SE for model selection and multimodel inference when
    model type is not supported or SE's are not available from predict
    methods (convenient for importing output from other software such as
    PRESENCE)
    \item added support for hurdle models for count data for model
    selection and model averaging parameter estimates
    \item added support for model selection of models of the 'gnls'
    class
    \item added support for model selection and multimodel inference of
    models of the 'betareg' class
    \item added support for model selection and multimodel inference of
    models of the 'survreg' class
    \item modified 'multComp' for computation of critical values with z
    value instead of t for certain model types
    \item corrected minor bug of 'multComp' for some models with a
    factor consisting of two levels
    \item added 'extractCN' function for models of 'unmarkedFit' classes
    \item modified namespace
  }
}



\section{Changes in version 2.00 (2014-07-14)}{
  \itemize{
    \item modified 'predictSE.lme' to handle the use of formulas
    specified externally from model call (suggested by Clark Andersen)
    \item added support for named lists of candidate models and add
    generic model names when modnames vector is NULL in model selection
    and model averaging functions
    \item replaced functions 'extract.LL', 'modavg.effect',
    'modavgpred', 'modavg.shrink', 'mult.comp',  with 'extractLL',
    'modavgEffect', 'modavgPred', 'modavgShrink', and 'multComp'
    respectively, to avoid problems with generic functions 
    \item 'extractLL', 'extractSE', 'AICc', 'aictab', 'modavg',
    'modavgEffect', 'modavgShrink', 'modavgPred', 'multComp', and
    predictSE made generic
    \item changed order of some arguments to present the most often used
    first
    \item added a check in 'multComp' for missing values in data set,
    and added options to correct confidence intervals for multiple
    comparisons
    \item added support for presence-only models of maxlikeFit class
    \item added support for models of classes 'coxme' and 'lmekin'
    \item added support for models of classes 'vglm' (package VGAM) and
    'zeroinfl' (package pscl)
    \item added support for models of 'unmarkedFitGMM',
    'unmarkedFitMPois', and unmarkedFitGPC classes
    \item added support for model averaging availability parameter (phi)
    for models of 'unmarkedFitGMM', 'unmarkedFitMPois', and
    unmarkedFitGPC classes
    \item corrected minor bug in some cases of missing values in
    'Nmix.chisq'
    \item corrected minor bug in 'mb.gof.test' that produced warning( )
    \item added 'aov' method for all generic functions of the package
    \item simplified Rd files AICc.rd, aictab.Rd, boot.wt.Rd, DIC.Rd,
    dictab.Rd, importance.Rd, modavg.Rd, modavg.effect.Rd
    \item added Bullfrog occupancy data set
    \item added Spotted salamander capture-mark-recapture data
    \item added Gopher Tortoise distance sampling data set
  }
}



\section{Changes in version 1.35 (2013-11-17)}{
  \itemize{
    \item corrected minor bug in chi-square for N-mixture models in
    'Nmix.chisq'
    \item corrected minor bug in title display, now uses
    as.expression(substitute())
  }
}



\section{Changes in version 1.34 (2013-11-13)}{
  \itemize{
    \item added functionality to compute model-averaged estimates
    (means) with 'mult.comp'
    \item added functionality for clogit models
    \item added one-way ANOVA turkey data set
    \item added two-way ANOVA iron data set
    \item modified 'predictSE.zip' to handle offsets and added c.hat
    argument
    \item added functionality of overdispersion correction for N-mixture
    models in 'aictab', 'modavg', 'modavgpred', 'modavg.effect',
    'modavg.shrink', 'predictSE.zip'
    \item added functions 'Nmix.gof.test' and 'Nmix.chisq' to assess GOF
    of N-mixture models
  }
}



\section{Changes in version 1.33 (2013-09-12)}{
  \itemize{
    \item added 'mult.comp' function to compute multiple comparisons on
    factors
    \item added 'boot.wt' function to compute bootstrapped relative
    frequencies
    \item corrected minor bug in 'mb.chisq' for certain patterns of
    missing values
    \item added conditional test to avoid problems on CRAN with some
    packages unavailable on certain platforms
    \item added the lizards habitat preference data set
  }
}



\section{Changes in version 1.32 (2013-07-05)}{
  \itemize{
    \item corrected minor bug in 'modavg' when model-averaging with
    models of classes 'clm' and 'clmm' when the 'exclude' argument was
    used
    \item corrected bug with certain types of missing values in
    'mb.chisq'
    \item added function to extract DIC and create model selection
    tables based on DIC (functions 'DIC' and 'dictab')
  }
}



\section{Changes in version 1.31 (2013-06-26)}{
  \itemize{
    \item added functions to work with new 'merMod' classes of package
    lme4
  }
}



\section{Changes in version 1.30 (2013-05-16)}{
  \itemize{
    \item modified example code of predictSE.zip( ) to avoid errors on
    some platforms during package building
  }
}



\section{Changes in version 1.29 (2013-04-22)}{
  \itemize{
    \item added 'mb.chisq' to implement MacKenzie and Bailey (2004)
    goodness-of-fit test on single-season occupancy models and
    'mb.gof.test' to compute P-values, c-hat, and histogram of
    bootstrapped test statistics
    \item corrected bug for offsets in 'predictSE.mer'
  }
}



\section{Changes in version 1.28 (2013-03-21)}{
  \itemize{
    \item added functionality for zip models to compute SE's of
    predictions on abundance of N-mixture models from unmarked
    \item added functionality for 'occuFP' models (detections with
    false-positives)
  }
}



\section{Changes in version 1.27 (2013-02-07)}{
  \itemize{
    \item implemented model-averaging estimates and predictions for
    'gdistsamp' models
    \item added support for nonlinear mixed models of class 'nlme'
  }
}



\section{Changes in version 1.26 (2012-08-30)}{
  \itemize{
    \item added support for robust regression models of class 'rlm' -
    'AICc', 'aictab', 'importance', 'modavg', 'modavg.effect',
    'modavg.shrink', 'modavgpred'
    \item added support for ordinal logistic regression models of
    classes 'clm' and 'clmm' - 'AICc', 'aictab', 'importance', 'modavg',
    'modavg.effect', 'modavg.shrink'
    \item modified modavg.unmarked( ) to enable treatment of gamma
    parameter in pcountOpen models with different formulations (Const,
    Trend, AR)
    \item added reference of Luckacs et al. 2010 about model averaging
    with shrinkage in modavg.shrink.Rd file
    \item added references Burnham et al. 2011 for using model selection
    as alternative to multiple comparisons in modavg.effect.Rd file
    \item removed white space from formulas to reduce chance of reaching
    character limit in 'parse' in 'predictSE.mer', 'predictSE.lme', and
    'predictSE.gls'
    \item changed treatment of parm with modavg.unmarked to add object
    parm.strip (line 35) to enable functionality with 'regexpr'
  }
}



\section{Changes in version 1.25 (2012-06-07)}{
  \itemize{
    \item added Cox proportional hazards model support
    \item added modavg.effect functions to compute model-averaged effect
    sizes for true experiments as an information-theoretic analogue to
    multiple comparisons
  }
}



\section{Changes in version 1.24 (2012-02-06)}{
  \itemize{
    \item corrected minor bug in 'AICc.nls' to determine number of
    observations automatically when "second.ord = TRUE"
  }
}



\section{Changes in version 1.23 (2012-01-27)}{
  \itemize{
    \item corrected minor bug following revision 1.22 on 'predictSE.mer'
    to avoid problems with factors
  }
}



\section{Changes in version 1.22 (January 2012-01-19)}{
  \itemize{
    \item modified example in 'modavg' help file
    \item added error message for models without intercept with 'modavg'
    \item added support for offset variable in 'predictSE.mer' for
    Poisson GLMM
    \item added support for 'distsamp' and 'gdistsamp' models from
    unmarked package (Poisson or negative binomial abundance with log
    link)
    \item added example in 'modavg' and 'modavgpred' for 'distsamp'
    models
    \item added example for Poisson GLMM with log link and offset
  }
}



\section{Changes in version 1.21 (2011-11-04)}{
  \itemize{
    \item corrected error in computing sample size for binomial model
    with events/trials syntax in 'AICc.glm' and 'AICc.mer'
  }
}



\section{Changes in version 1.20 (2011-10-03)}{
  \itemize{
    \item added functionality of 'nobs' argument for objects created in
    unmarked package
  }
}



\section{Changes in version 1.19 (2011-09-14)}{
  \itemize{
    \item created function 'predictSE.gls' to compute predicted values
    and SE's from 'gls' model fits and added functionality of 'modavgpred'
    with 'gls' model fit
    \item added documentation with examples for 'gls'
  }
}



\section{Changes in version 1.18 (2011-08-26)}{
  \itemize{
    \item added NAMESPACE
    \item modified functions 'extractSE.mer' and streamlined extraction
    of estimate labels of 'mer' objects in 'modavg.mer',
    'modavg.shrink.mer', 'importance'
  }
}



\section{Changes in version 1.17 (2011-06-17)}{
  \itemize{
    \item fixed minor bug for external calls to 'AICc' to assign
    'unmarkedFit' classes
    \item added functionality of 'AICc', 'aictab', 'evidence',
    'importance', 'confset' with 'nls' fits
  }
}



\section{Changes in version 1.16 (2011-05-27)}{
  \itemize{
    \item added model selection and multimodel inference functionality
    with unmarked package for objects of classes unmarkedFitOccu,
    unmarkedFitColExt, unmarkedFitPCount, unmarkedFitPCO,
    unmarkedFitOccuRN for models incorporating detectability
    \item modifications of functions 'aictab', 'importance', 'confset',
    'evidence', 'modavg', 'modavgpred', 'modavg.shrink', and
    extract.LL.unmarked utility function added
    \item added new argument parm.type, to specify class of parameters
    for objects of unmarkedFit classes
    \item added several examples for dealing with objects of
    'unmarkedFit' classes
    \item modified 'modavg', 'modavg.shrink', 'importance', and
    'reverse.exclude' functions to strip white spaces within elements
    entered inadvertently by user in the parm and exclude arguments
  }
}



\section{Changes in version 1.15 (2011-04-21)}{
  \itemize{
    \item in 'modavg' and related functions, added functionality of
    dealing with specification of interaction terms (e.g., A:B, B:A, A*B,
    B*A) by adding two new utility functions to recognize the terms when
    user is not consistent in his model syntax
    \item created function 'modavg.shrink' for the shrinkage version of
    model-averaging estimates to reduce model selection bias when the same
    number of models include/exclude a variable
  }
}



\section{Changes in version 1.14 (2011-03-11)}{
  \itemize{
    \item fixed minor bug in 'modavg.mer' regarding exclusion of
    interaction terms as 'formula' includes spaces within interaction
    terms
  }
}



\section{Changes in version 1.13 (2011-01-10)}{
  \itemize{
    \item fixed minor bug in 'AICc.lme' on the computation of n in
    certain lme's with several levels
  }
}



\section{Changes in version 1.12 (2010-11-22)}{
  \itemize{
    \item fixed minor bug in 'modavgpred.lme' that generated error
    message when second.ord = FALSE
    \item fixed minor bug by dispatching 'mer' class objects to
    'AICc.mer'
    \item fixed minor bug with incorrect computation of n for 'mer'
    models with binomial distribution using cbind(incidence, size -
    incidence) or binary syntax
  }
}
  


\section{Changes in version 1.11 (2010-10-08)}{
  \itemize{
    \item fixed minor bug in 'modavgpred.glm' that generated an error
    message for Poisson, binomial, and negative binomial distributions
    when 'c.hat' = 1
    \item added additional example in 'predictSE.mer' with gaussian
    distribution but with log link
  }
}



\section{Changes in version 1.10 (2010-10-07)}{
  \itemize{
    \item fixed minor bug in 'modavgpred.glm' for SE's of predictions
    from gaussian model
    \item fixed minor bug in 'modavg.glm' for SE's of beta estimates
    from negative binomial glm
  }
}



\section{Changes in version 1.09 (2010-09-16)}{
  \itemize{
    \item added function 'predictSE.mer' to compute predictions and SE's
    on predictions accounting only for fixed effects (mean population
    prediction) from objects of 'mer' class
    \item added functionality to 'modavgpred' to include computation of
    model-averaged predictions from objects of class 'mer'
    \item added a check for link function to ensure same link function
    was used for all models in 'modavg' when model averaging beta
    estimates
    \item created a utility function to determine family and link
    function for objects of 'mer' class called by 'modavg',
    'modavgpred', and 'predictSE.mer'
    \item added a warning when some models have an identical (Q)AIC(c)
    and are possibly redundant in 'aictab'
    \item corrected minor bug in number of observations in 'modavgpred'
    \item added a check for binomial family with events/trials syntax to
    use sum of weights as number of observations in 'AICc.mer' when nobs =
    NULL
  }
}



\section{Changes in version 1.08 (2010-08-16)}{
  \itemize{
    \item added Suggests "lme4" in DESCRIPTION file
  }
}



\section{Changes in version 1.07 (2010-05-24)}{
  \itemize{
    \item modified modavg to include 'nobs' argument in call to
    'modavg.lme' and 'modavg.gls'
    \item added support for generalized linear mixed models of the mer
    class for 'AICc', 'aictab', 'modavg', 'evidence', 'importance',
    'confset' functions
    \item updated information on support for 'mer' class in help files
    \item added examples in some of the help files with 'mer'
  }
}



\section{Changes in version 1.06 (2010-02-06)}{
  \itemize{
    \item added a 'polr' class check to 'AICc' function
    \item added function 'predictSE.lme' to make predictions based on
    fixed effects for objects of class 'lme' (mean population predictions)
    and compute associated SE's
    \item added function 'modvagpred.lme' to compute model-averaged
    predictions and unconditional standard errors for 'lme' models
    \item added statement in Rd files for 'aictab', 'confset',
    'evidence', 'importance', 'modavgpred', 'modavg',  and generics
    'nobs = NULL' defaults to total number of observations in
    computations of information criteria
    \item added NEWS file to document changes across versions of
    AICcmodavg package
    \item added check that same response variable is used throughout
    model set for 'aictab' and generics with lapply(X = cand.set, FUN =
    function(b) formula(b)[2])
    \item added functionality of 'AICc', 'aictab', 'confset',
    'evidence', 'importance', 'modavg' with objects of 'gls' class
    \item added additional flexibility to 'evidence' function when model
    selection table provided is not sorted ('evidence' now sorts table
    internally) and set default comparison to top-ranked model vs
    second-ranked model
    \item added check/error message when no dispersion parameter (i.e.,
    'gamdisp') is specified for 'modavg' and 'modavgpred' and generics
    with gamma regression models
    \item added additional check/error message when c.hat < 1 in
    'AICc.glm'
  }
}



\section{Changes in version 1.05 (2009-12-03)}{
  \itemize{
    \item modified number of digits to print with 'importance' to be
    consistent with revision 1.04
    \item added parenthesis where one was missing in aictab.Rd under
    value at line for (Q)AIC(c)
    \item corrected example in evidence.Rd for faulty name of model.low
    - should have been 'model.low = raw density'
    \item deleted CumWt when 'sort = FALSE' in 'aictab' - by modifying
    argument in 'print.aictab'
    \item for 'AICc.glm', 'AICc.mult', c.hat > 1 && c.hat <= 4 now gives
    warning
    \item added a column for logLik, quasi-logLik, or restricted logLik
    in 'aictab'
    \item added options in 'print.aictab' to select items to display
    \item added a check for 'exclude =' in 'modavg' and generics to warn
    if exclude is not a list
    \item documented 'aictab' with new LL, Quasi.LL, Res.LL elements
    \item used 'vcov' and argument 'dispersion =' to scale SE's in a
    single step instead of adding more code to scale it on separate lines
    - works with Poisson, Binomial, Gamma
    \item in 'modavg.polr' and 'modavg.mult', added a logical test to
    distinguish between 'coef' or zeta which are the intercepts
    \item modified 'print.aictab' so that c_hat is printed only when c_hat
    appears in AICc table (i.e, when any(names(x)=="c_hat"))
    \item modified 'print.modavg' so that c_hat is printed only when c_hat
    appears in AICc table (i.e, when any(names(oldtab)=="c_hat")) AND
    deleted Cum.Wt from x$Mod.avg.table
    \item added a check for the order for aictab and modavg and generics
    so that LL and other measures are displayed properly in table and
    sorted at the very end (placed line of sort and cumsum at very end of
    code) with Results[rev(order(Results[,5])), ]
    \item added Cum.Wt in last column of Results
    \item changed order of elements with colnames so that prints nicely -
    this affects 'importance', 'confset', and 'evidence' - number of
    columns used might differ depending on what is requested
    \item in 'modavg' and generics, changed 'regexpr' to determine whether
    parameter appears alone in the model and added 'fixed = TRUE', to find
    matches "as is" to avoid problems with 'mult' (e.g., "well|mild")
  }
}



\section{Changes in version 1.04 (2009-11-11)}{
  \itemize{
    \item added the new equation for unconditional SE based on equation
    6.12 of Burnham and Anderson (2002) in 'modavg' and 'modavgpred' and
    new argument 'uncond.se' with default value "revised" -  'uncond.se =
    old' uses equation 4.19
    \item added a check to determine whether all link functions are
    identical for 'modavg' as averaging over beta estimates under
    different link functions is not appropriate and returns error)
    \item added a check to determine whether all link functions are
    identical for 'modavgpred' when linear predictor is used
    \item changed default to display only 2 digits after decimal point
    for all output 'print' generic functions
  }
}



\section{Changes in version 1.03 (2009-10-21)}{
  \itemize{
    \item replaced 'regexpr' in 'modavg.glm', 'modavg.lme',
    'modavg.mult' with 'match'
    \item -added functionality for 'polr'
    \item cleaned up code with indentations
    \item added 'exclude' argument to manually exclude models with
    interaction terms/quadratic terms in 'modavg'
    \item added 'warn' argument (to suppress warning that variables with
    similar names appear in a given model) to 'modavg' function.  This is
    a fundamental difference with the previous versions of the function as
    now the user must specify the variables (as a list) to exclude models
    for the computation of model-averaged estimates (case when
    interaction/polynomial terms are present).  This modication also
    solved the problem of variables with similar names
  }
}



\section{Changes in version 1.02 (2009-09-29)}{
  \itemize{
    \item corrected small bug for computation of AICc (in 'AICc.glm')
    for binomial models with events/trials syntax to use the correct
    samplesize
    \item added function to compute (Q)AIC(c) for multinomial models and
    updated 'aictab' function accordingly
    \item added citation file
    \item changed 'Depends' to 'Suggests' for packages MASS, nnet, nlme
    \item modified 'importance' to include 'multinom' objects
    \item removed dependencies with nlme to avoid potential problems with
    lme4
  }
}



\section{Changes in version 1.01 (2009-09-03)}{
  \itemize{
    \item cleaned up documentation file, added a few more sentences and
    added space within formulas for the example to display properly in pdf
  }
}

##approximate F-test in presence of overdispersion
##generic
anovaOD <- function(mod.simple, mod.complex, c.hat = 1, nobs = NULL, ...){
    UseMethod("anovaOD", mod.simple)
}



anovaOD.default <- function(mod.simple, mod.complex, c.hat = 1, nobs = NULL, ...){
    stop("\nFunction not yet defined for this object class\n")
}



##glm
anovaOD.glm <- function(mod.simple, mod.complex, c.hat = 1, nobs = NULL, ...){

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##check for distributions
    modFamily1 <- family(mod.simple)$family
    modFamily2 <- family(mod.complex)$family
    if(!identical(modFamily1, modFamily2)) stop("\nComparisons only appropriate for models using the same distribution\n")

    if(!identical(modFamily1, "poisson") && !identical(modFamily1, "binomial")) {
        if(c.hat > 1) stop("\ndistribution not appropriate for overdispersion correction\n\n")
    }

    ##for binomial, check that number of trials > 1
    if(identical(modFamily1, "binomial")) {
        if(!any(mod.simple$prior.weights > 1)) stop("\nOverdispersion correction only appropriate for success/trials syntax\n\n")
    }

    ##response variable
    y1 <- mod.simple$y
    y2 <- mod.complex$y

    ##number of observations
    if(is.null(nobs)) {
        nobs <- length(y1)
    }
    
    ##check that data are the same
    if(!identical(y1, y2)) stop("\nData set should be identical to compare models\n")
  
    ##extract log-likelihood
    LL1 <- logLik(mod.simple)
    LL2 <- logLik(mod.complex)
    LL.simple <- LL1[1]
    LL.complex <- LL2[1]

    ##extract number of estimated parameters
    K.simple <- attr(LL1, "df")
    K.complex <- attr(LL2, "df")
    
    ##residual df of complex model
    df.complex <- nobs - K.complex

    ##extract model formula
    simpleForm <- formula(mod.simple)
    form.simple <- paste(simpleForm[2], "~", simpleForm[3])
    complexForm <- formula(mod.complex)
    form.complex <- paste(complexForm[2], "~", complexForm[3])
    
    ##- 2 * (logLik(simple) - logLik(complex))
    LR <- -2 * (LL.simple - LL.complex)
    ##difference in number of estimated parameters
    K.diff <- K.complex - K.simple

    ##use chi-square if no overdispersion
    if(c.hat == 1) {
        Chistat <- LR
        df <- K.diff
        pval <- 1 - pchisq(Chistat, df = df)

        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Chistat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "Chistat", "pval")
    
    } else {
        ##compute F statistic
        Fstat <- LR/((K.diff)*c.hat)

        ##compute df
        df.num <- K.diff
        df.denom <- df.complex

        ##P value
        pval <- 1 - pf(Fstat, df1 = df.num, df2 = df.denom)
        
        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Fstat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "F", "pval")
    }
    
    ##assemble in list
    outList <- list(form.simple = form.simple,
                    form.complex = form.complex,
                    c.hat = c.hat,
                    devMat = devMat)
    class(outList) <- c("anovaOD", "list")
    return(outList)
}



##vglm
anovaOD.vglm <- function(mod.simple, mod.complex, c.hat = 1, nobs = NULL, ...){

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##check family of vglm to avoid problems
    fam.type1 <- mod.simple@family@vfamily[1]
    fam.type2 <- mod.complex@family@vfamily[1]
    if(!identical(fam.type1, fam.type2)) stop("\nComparisons only appropriate for models using the same distribution\n")
    
    if(!(fam.type1 == "poissonff" || fam.type1 == "binomialff" || fam.type1 == "multinomial")) stop("\nDistribution not supported by function\n")
    if(fam.type1 == "binomialff") {
        if(!any(mod.simple@prior.weights > 1)) stop("\nOverdispersion correction only appropriate for success/trials syntax\n\n")
    }

    ##number of observations
    y1 <- mod.simple@y
    y2 <- mod.complex@y

    ##number of observations
    if(is.null(nobs)) {
        nobs <- length(y1)
    }
    
    ##check that data are the same
    if(!identical(y1, y2)) stop("\nData set should be identical to compare models\n")

    ##extract log-likelihood
    LL1 <- extractLL(mod.simple)
    LL2 <- extractLL(mod.complex)
    LL.simple <- LL1[1]
    LL.complex <- LL2[1]

    ##extract number of estimated parameters
    K.simple <- attr(LL1, "df")
    K.complex <- attr(LL2, "df")
    
    ##extract model formula
    simpleForm <- formula(mod.simple)
    form.simple <- paste(simpleForm[2], "~", simpleForm[3])
    complexForm <- formula(mod.complex)
    form.complex <- paste(complexForm[2], "~", complexForm[3])
    
    ##residual df of complex model
    df.complex <- nobs - K.complex

    ##- 2 * (logLik(simple) - logLik(complex))
    LR <- -2 * (LL.simple - LL.complex)
    ##difference in number of estimated parameters
    K.diff <- K.complex - K.simple
    
    ##use chi-square if no overdispersion
    if(c.hat == 1) {
        Chistat <- LR
        df <- K.diff
        pval <- 1 - pchisq(Chistat, df = df)

        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Chistat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "Chistat", "pval")
    
    } else {

        ##compute F statistic
        Fstat <- LR/((K.diff)*c.hat)

        ##compute df
        df.num <- K.diff
        df.denom <- df.complex

        ##P value
        pval <- 1 - pf(Fstat, df1 = df.num, df2 = df.denom)
        
        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Fstat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "F", "pval")
    }
    
    ##assemble in list
    outList <- list(form.simple = form.simple,
                    form.complex = form.complex,
                    c.hat = c.hat,
                    devMat = devMat)
    class(outList) <- c("anovaOD", "list")
    return(outList)
}



##multinom
anovaOD.multinom <- function(mod.simple, mod.complex, c.hat = 1, nobs = NULL, ...){

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##number of observations
    nobs1 <- nrow(mod.simple$fitted.values)
    nobs2 <- nrow(mod.complex$fitted.values)

    ##check that sample size is the same
    if(!identical(nobs1, nobs2)) stop("\nData set should be identical to compare models\n")

    ##number of observations
    if(is.null(nobs)) {
        nobs <- nobs1
    }
    
    ##extract log-likelihood
    LL1 <- logLik(mod.simple)
    LL2 <- logLik(mod.complex)
    LL.simple <- LL1[1]
    LL.complex <- LL2[1]

    ##extract number of estimated parameters
    K.simple <- attr(LL1, "df")
    K.complex <- attr(LL2, "df")

    ##extract model formula
    simpleForm <- formula(mod.simple)
    form.simple <- paste(simpleForm[2], "~", simpleForm[3])
    complexForm <- formula(mod.complex)
    form.complex <- paste(complexForm[2], "~", complexForm[3])
    
    ##residual df of complex model
    df.complex <- nobs - K.complex

    ##- 2 * (logLik(simple) - logLik(complex))
    LR <- -2 * (LL.simple - LL.complex)
    ##difference in number of estimated parameters
    K.diff <- K.complex - K.simple

    ##use chi-square if no overdispersion
    if(c.hat == 1) {
        Chistat <- LR
        df <- K.diff
        pval <- 1 - pchisq(Chistat, df = df)

        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Chistat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "Chistat", "pval")
    } else {

        ##compute F statistic
        Fstat <- LR/((K.diff)*c.hat)

        ##compute df
        df.num <- K.diff
        df.denom <- df.complex

        ##P value
        pval <- 1 - pf(Fstat, df1 = df.num, df2 = df.denom)
        
        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Fstat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "F", "pval")
    }
    
    ##assemble in list
    outList <- list(form.simple = form.simple,
                    form.complex = form.complex,
                    c.hat = c.hat,
                    devMat = devMat)
    class(outList) <- c("anovaOD", "list")
    return(outList)
}


###############################
##residual DF for mixed models is difficult
##glmerMod
anovaOD.glmerMod <- function(mod.simple, mod.complex, c.hat = 1, nobs = NULL, ...){

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }
   
    ##check for distributions
    modFamily1 <- family(mod.simple)$family
    modFamily2 <- family(mod.complex)$family
    if(!identical(modFamily1, modFamily2)) stop("\nComparisons only appropriate for models using the same distribution\n")

    if(!identical(modFamily1, "poisson") && !identical(modFamily1, "binomial")) {
        if(c.hat > 1) stop("\ndistribution not appropriate for overdispersion correction\n\n")
    }

    ##number of observations
    y1 <- mod.simple@resp$y
    y2 <- mod.complex@resp$y

    ##check that data are the same
    if(!identical(y1, y2)) stop("\nData set should be identical to compare models\n")

    ##number of observations
    if(is.null(nobs)) {
        nobs <- length(y1)
    }
    
    ##for binomial, check that number of trials > 1
    if(identical(modFamily1, "binomial")) {
        if(!any(mod.simple@resp$weights > 1)) stop("\nOverdispersion correction only appropriate for success/trials syntax\n\n")
    }

        
    ##extract LL
    LL.simple <- logLik(mod.simple)
    LL.complex <- logLik(mod.complex)

    ##extract number of estimated parameters
    K.simple <- attr(logLik(mod.simple), "df")
    K.complex <- attr(logLik(mod.complex), "df")

    ##residual df of complex model
    df.complex <- nobs - K.complex

    ##extract model formula
    simpleForm <- formula(mod.simple)
    form.simple <- paste(simpleForm[2], "~", simpleForm[3])
    complexForm <- formula(mod.complex)
    form.complex <- paste(complexForm[2], "~", complexForm[3])

    ##- 2 * (logLik(simple) - logLik(complex))
    LR <- -2 * (LL.simple - LL.complex)
    ##difference in number of estimated parameters
    K.diff <- K.complex - K.simple

    ##use chi-square if no overdispersion
    if(c.hat == 1) {
        Chistat <- LR
        df <- K.diff
        pval <- 1 - pchisq(Chistat, df = df)

        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Chistat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "Chistat", "pval")
    
    } else {
        ##compute F statistic
        Fstat <- LR/((K.diff)*c.hat)

        ##compute df
        df.num <- K.diff
        df.denom <- df.complex

        ##P value
        pval <- 1 - pf(Fstat, df1 = df.num, df2 = df.denom)
        
        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Fstat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "F", "pval")
    }
    
    ##assemble in list
    outList <- list(form.simple = form.simple,
                    form.complex = form.complex,
                    c.hat = c.hat,
                    devMat = devMat)
    class(outList) <- c("anovaOD", "list")
    return(outList)
}



##glmmTMB
anovaOD.glmmTMB <- function(mod.simple, mod.complex, c.hat = 1, nobs = NULL, ...) {

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##check for distributions
    ##check for distributions
    modFamily1 <- family(mod.simple)$family
    modFamily2 <- family(mod.complex)$family
    if(!identical(modFamily1, modFamily2)) stop("\nComparisons only appropriate for models using the same distribution\n")

    if(!identical(modFamily1, "poisson") && !identical(modFamily1, "binomial")) {
        if(c.hat > 1) stop("\ndistribution not appropriate for overdispersion correction\n\n")
    }
    
    ##extract response
    y1Name <- mod.simple$modelInfo$respCol
    y2Name <- mod.simple$modelInfo$respCol
    y1 <- mod.simple$frame[, y1Name, drop = FALSE]
    y2 <- mod.simple$frame[, y2Name, drop = FALSE]

    ##check that data are the same
    if(!identical(y1, y2)) stop("\nData set should be identical to compare models\n")
    
    ##number of observations
    if(is.null(nobs)) {
        nobs <- nrow(y1)
    }    

    ##if binomial, check if n > 1 for each case
    if(modFamily1 == "binomial") {
        resp <- mod.simple$frame[, mod.simple$modelInfo$respCol]
        if(!is.matrix(resp)) {
            if(!any(names(mod.simple$frame) == "(weights)")) {
                stop("\nOverdispersion correction only appropriate for success/trials syntax\n\n")
            }
        }
    }
    

    ##extract LL
    LL.simple <- logLik(mod.simple)
    LL.complex <- logLik(mod.complex)

    ##extract number of estimated parameters
    K.simple <- attr(logLik(mod.simple), "df")
    K.complex <- attr(logLik(mod.complex), "df")

    ##residual df of complex model
    df.complex <- nobs - K.complex

    ##extract model formula
    simpleForm <- formula(mod.simple)
    form.simple <- paste(simpleForm[2], "~", simpleForm[3])
    complexForm <- formula(mod.complex)
    form.complex <- paste(complexForm[2], "~", complexForm[3])

    ##- 2 * (logLik(simple) - logLik(complex))
    LR <- -2 * (LL.simple - LL.complex)
    ##difference in number of estimated parameters
    K.diff <- K.complex - K.simple

    ##use chi-square if no overdispersion
    if(c.hat == 1) {
        Chistat <- LR
        df <- K.diff
        pval <- 1 - pchisq(Chistat, df = df)

        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Chistat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "Chistat", "pval")
    
    } else {
        ##compute F statistic
        Fstat <- LR/((K.diff)*c.hat)

        ##compute df
        df.num <- K.diff
        df.denom <- df.complex

        ##P value
        pval <- 1 - pf(Fstat, df1 = df.num, df2 = df.denom)
        
        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Fstat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "F", "pval")
    }
    
    ##assemble in list
    outList <- list(form.simple = form.simple,
                    form.complex = form.complex,
                    c.hat = c.hat,
                    devMat = devMat)

    class(outList) <- c("anovaOD", "list")
    return(outList)
}



##occu
anovaOD.unmarkedFitOccu <- function(mod.simple, mod.complex, c.hat = 1, nobs = NULL, ...) {

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##extract response
    y1 <- mod.simple@data@y
    y2 <- mod.complex@data@y

    ##check that data are the same
    if(!identical(y1, y2)) stop("\nData set should be identical to compare models\n")
    
    ##number of observations
    if(is.null(nobs)) {
        nobs <- nrow(y1)
    }    

    ##extract LL
    LL.simple <- logLik(mod.simple)
    LL.complex <- logLik(mod.complex)

    ##extract number of estimated parameters
    K.simple <- length(coef(mod.simple))
    K.complex <- length(coef(mod.complex))

    ##residual df of complex model
    df.complex <- nobs - K.complex

    ##extract model formula
    ##simple model
    form.simplePsi <- formulaShort(mod.simple, unmarked.type = "state")
    form.simplePsi2 <- paste("psi(", form.simplePsi, ")", sep = "")

    form.simpleDet <- formulaShort(mod.simple, unmarked.type = "det")
    form.simpleDet2 <- paste("p(", form.simpleDet, ")", sep = "")
    form.simple <- paste(form.simplePsi2, form.simpleDet2, sep = "")
    
    ##complex model
    form.complexPsi <- formulaShort(mod.complex, unmarked.type = "state")
    form.complexPsi2 <- paste("psi(", form.complexPsi, ")", sep = "")
    form.complexDet <- formulaShort(mod.complex, unmarked.type = "det")
    form.complexDet2 <- paste("p(", form.complexDet, ")", sep = "")
    form.complex <- paste(form.complexPsi2, form.complexDet2, sep = "")
    
    
    ##- 2 * (logLik(simple) - logLik(complex))
    LR <- -2 * (LL.simple - LL.complex)
    ##difference in number of estimated parameters
    K.diff <- K.complex - K.simple

    ##use chi-square if no overdispersion
    if(c.hat == 1) {
        Chistat <- LR
        df <- K.diff
        pval <- 1 - pchisq(Chistat, df = df)

        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Chistat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "Chistat", "pval")
    
    } else {
        ##compute F statistic
        Fstat <- LR/((K.diff)*c.hat)

        ##compute df
        df.num <- K.diff
        df.denom <- df.complex

        ##P value
        pval <- 1 - pf(Fstat, df1 = df.num, df2 = df.denom)
        
        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Fstat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "F", "pval")
    }
    
    ##assemble in list
    outList <- list(form.simple = form.simple,
                    form.complex = form.complex,
                    c.hat = c.hat,
                    devMat = devMat)
    
    class(outList) <- c("anovaOD", "list")
    return(outList)
}



##colext
anovaOD.unmarkedFitColExt <- function(mod.simple, mod.complex, c.hat = 1, nobs = NULL, ...) {

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##extract response
    y1 <- mod.simple@data@y
    y2 <- mod.complex@data@y

    ##check that data are the same
    if(!identical(y1, y2)) stop("\nData set should be identical to compare models\n")
    
    ##number of observations
    if(is.null(nobs)) {
        nobs <- nrow(y1)
    }    

    ##extract LL
    LL.simple <- logLik(mod.simple)
    LL.complex <- logLik(mod.complex)

    ##extract number of estimated parameters
    K.simple <- length(coef(mod.simple))
    K.complex <- length(coef(mod.complex))

    ##residual df of complex model
    df.complex <- nobs - K.complex

    ##extract model formula
    ##simple model
    form.simplePsi <- formulaShort(mod.simple, unmarked.type = "psi")
    form.simplePsi2 <- paste("psi(", form.simplePsi, ")", sep = "")

    form.simpleGam <- formulaShort(mod.simple, unmarked.type = "col")
    form.simpleGam2 <- paste("gam(", form.simpleGam, ")", sep = "")

    form.simpleEps <- formulaShort(mod.simple, unmarked.type = "ext")
    form.simpleEps2 <- paste("eps(", form.simpleEps, ")", sep = "")

    form.simpleDet <- formulaShort(mod.simple, unmarked.type = "det")
    form.simpleDet2 <- paste("p(", form.simpleDet, ")", sep = "")
    form.simple <- paste(form.simplePsi2, form.simpleGam2,
                         form.simpleEps2, form.simpleDet2, sep = "")

    ##complex model
    form.complexPsi <- formulaShort(mod.complex, unmarked.type = "psi")
    form.complexPsi2 <- paste("psi(", form.complexPsi, ")", sep = "")

    form.complexGam <- formulaShort(mod.complex, unmarked.type = "col")
    form.complexGam2 <- paste("gam(", form.complexGam, ")", sep = "")

    form.complexEps <- formulaShort(mod.complex, unmarked.type = "ext")
    form.complexEps2 <- paste("eps(", form.complexEps, ")", sep = "")

    form.complexDet <- formulaShort(mod.complex, unmarked.type = "det")
    form.complexDet2 <- paste("p(", form.complexDet, ")", sep = "")
    form.complex <- paste(form.complexPsi2, form.complexGam2,
                          form.complexEps2, form.complexDet2, sep = "")

    
    ##- 2 * (logLik(simple) - logLik(complex))
    LR <- -2 * (LL.simple - LL.complex)
    ##difference in number of estimated parameters
    K.diff <- K.complex - K.simple

    ##use chi-square if no overdispersion
    if(c.hat == 1) {
        Chistat <- LR
        df <- K.diff
        pval <- 1 - pchisq(Chistat, df = df)

        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Chistat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "Chistat", "pval")
    
    } else {
        ##compute F statistic
        Fstat <- LR/((K.diff)*c.hat)

        ##compute df
        df.num <- K.diff
        df.denom <- df.complex

        ##P value
        pval <- 1 - pf(Fstat, df1 = df.num, df2 = df.denom)
        
        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Fstat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "F", "pval")
    }
    
    ##assemble in list
    outList <- list(form.simple = form.simple,
                    form.complex = form.complex,
                    c.hat = c.hat,
                    devMat = devMat)
    
    class(outList) <- c("anovaOD", "list")
    return(outList)
}



##occuRN
anovaOD.unmarkedFitOccuRN <- function(mod.simple, mod.complex, c.hat = 1, nobs = NULL, ...) {

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##extract response
    y1 <- mod.simple@data@y
    y2 <- mod.complex@data@y
    
    ##check that data are the same
    if(!identical(y1, y2)) stop("\nData set should be identical to compare models\n")
    
    ##number of observations
    if(is.null(nobs)) {
        nobs <- nrow(y1)
    }    

    ##extract LL
    LL.simple <- logLik(mod.simple)
    LL.complex <- logLik(mod.complex)

    ##extract number of estimated parameters
    K.simple <- length(coef(mod.simple))
    K.complex <- length(coef(mod.complex))

    ##residual df of complex model
    df.complex <- nobs - K.complex

    ##extract model formula
    ##simple model
    form.simpleLam <- formulaShort(mod.simple, unmarked.type = "state")
    form.simpleLam2 <- paste("lam(", form.simpleLam, ")", sep = "")

    form.simpleDet <- formulaShort(mod.simple, unmarked.type = "det")
    form.simpleDet2 <- paste("p(", form.simpleDet, ")", sep = "")
    form.simple <- paste(form.simpleLam2, form.simpleDet2, sep = "")
    
    ##complex model
    form.complexLam <- formulaShort(mod.complex, unmarked.type = "state")
    form.complexLam2 <- paste("lam(", form.complexLam, ")", sep = "")
    form.complexDet <- formulaShort(mod.complex, unmarked.type = "det")
    form.complexDet2 <- paste("p(", form.complexDet, ")", sep = "")
    form.complex <- paste(form.complexLam2, form.complexDet2, sep = "")
    
    
    ##- 2 * (logLik(simple) - logLik(complex))
    LR <- -2 * (LL.simple - LL.complex)
    ##difference in number of estimated parameters
    K.diff <- K.complex - K.simple

    ##use chi-square if no overdispersion
    if(c.hat == 1) {
        Chistat <- LR
        df <- K.diff
        pval <- 1 - pchisq(Chistat, df = df)

        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Chistat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "Chistat", "pval")
    
    } else {
        ##compute F statistic
        Fstat <- LR/((K.diff)*c.hat)

        ##compute df
        df.num <- K.diff
        df.denom <- df.complex

        ##P value
        pval <- 1 - pf(Fstat, df1 = df.num, df2 = df.denom)
        
        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Fstat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "F", "pval")
    }
    
    ##assemble in list
    outList <- list(form.simple = form.simple,
                    form.complex = form.complex,
                    c.hat = c.hat,
                    devMat = devMat)
    
    class(outList) <- c("anovaOD", "list")
    return(outList)
}



##occuFP
anovaOD.unmarkedFitOccuFP <- function(mod.simple, mod.complex, c.hat = 1, nobs = NULL, ...) {

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##extract response
    y1 <- mod.simple@data@y
    y2 <- mod.complex@data@y

    ##check that data are the same
    if(!identical(y1, y2)) stop("\nData set should be identical to compare models\n")
    
    ##number of observations
    if(is.null(nobs)) {
        nobs <- nrow(y1)
    }    

    ##extract LL
    LL.simple <- logLik(mod.simple)
    LL.complex <- logLik(mod.complex)

    ##extract number of estimated parameters
    K.simple <- length(coef(mod.simple))
    K.complex <- length(coef(mod.complex))

    ##residual df of complex model
    df.complex <- nobs - K.complex

    ##extract model formula
    ##simple model
    ##determine if certain detections (b) occur
    form.simplePsi <- formulaShort(mod.simple, unmarked.type = "state")
    form.simplePsi2 <- paste("psi(", form.simplePsi, ")", sep = "")

    form.simpleFp <- formulaShort(mod.simple, unmarked.type = "fp")
    form.simpleFp2 <- paste("fp(", form.simpleFp, ")", sep = "")

    form.simpleDet <- formulaShort(mod.simple, unmarked.type = "det")
    form.simpleDet2 <- paste("p(", form.simpleDet, ")", sep = "")

    if(exists("b", mod.simple@estimates@estimates)){
        form.simpleB <- formulaShort(mod.simple, unmarked.type = "b")
        form.simpleB2 <- paste("b(", form.simpleB, ")", sep = "")

        form.simple <- paste(form.simplePsi2, form.simpleFp2,
                             form.simpleB2, form.simpleDet2, sep = "")
    } else {
        form.simple <- paste(form.simplePsi2, form.simpleFp2,
                             form.simpleDet2, sep = "")
    }

    ##complex model
    ##determine if certain detections (b) occur
    form.complexPsi <- formulaShort(mod.complex, unmarked.type = "state")
    form.complexPsi2 <- paste("psi(", form.complexPsi, ")", sep = "")

    form.complexFp <- formulaShort(mod.complex, unmarked.type = "fp")
    form.complexFp2 <- paste("fp(", form.complexFp, ")", sep = "")

    form.complexDet <- formulaShort(mod.complex, unmarked.type = "det")
    form.complexDet2 <- paste("p(", form.complexDet, ")", sep = "")

    if(exists("b", mod.complex@estimates@estimates)){
        form.complexB <- formulaShort(mod.complex, unmarked.type = "b")
        form.complexB2 <- paste("b(", form.complexB, ")", sep = "")

        form.complex <- paste(form.complexPsi2, form.complexFp2,
                             form.complexB2, form.complexDet2, sep = "")
    } else {
        form.complex <- paste(form.complexPsi2, form.complexFp2,
                             form.complexDet2, sep = "")
    }
    
    ##- 2 * (logLik(simple) - logLik(complex))
    LR <- -2 * (LL.simple - LL.complex)
    ##difference in number of estimated parameters
    K.diff <- K.complex - K.simple

    ##use chi-square if no overdispersion
    if(c.hat == 1) {
        Chistat <- LR
        df <- K.diff
        pval <- 1 - pchisq(Chistat, df = df)

        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Chistat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "Chistat", "pval")
    
    } else {
        ##compute F statistic
        Fstat <- LR/((K.diff)*c.hat)

        ##compute df
        df.num <- K.diff
        df.denom <- df.complex

        ##P value
        pval <- 1 - pf(Fstat, df1 = df.num, df2 = df.denom)
        
        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Fstat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "F", "pval")
    }
    
    ##assemble in list
    outList <- list(form.simple = form.simple,
                    form.complex = form.complex,
                    c.hat = c.hat,
                    devMat = devMat)
    
    class(outList) <- c("anovaOD", "list")
    return(outList)
}



##pcount
anovaOD.unmarkedFitPCount <- function(mod.simple, mod.complex, c.hat = 1, nobs = NULL, ...) {

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##check for distributions
    modFamily1 <- mod.simple@mixture
    modFamily2 <- mod.complex@mixture
    if(!identical(modFamily1, modFamily2)) stop("\nComparisons only appropriate for models using the same mixture distribution\n")

    if(!identical(modFamily1, "P") && !identical(modFamily1, "ZIP")) {
        if(c.hat > 1) stop("\ndistribution not appropriate for overdispersion correction\n\n")
    }
    
    ##extract response
    y1 <- mod.simple@data@y
    y2 <- mod.complex@data@y

    ##check that data are the same
    if(!identical(y1, y2)) stop("\nData set should be identical to compare models\n")
    
    ##number of observations
    if(is.null(nobs)) {
        nobs <- nrow(y1)
    }    

    ##extract LL
    LL.simple <- logLik(mod.simple)
    LL.complex <- logLik(mod.complex)

    ##extract number of estimated parameters
    K.simple <- length(coef(mod.simple))
    K.complex <- length(coef(mod.complex))

    ##residual df of complex model
    df.complex <- nobs - K.complex

    ##extract model formula
    ##simple model
    form.simpleLam <- formulaShort(mod.simple, unmarked.type = "state")
    form.simpleLam2 <- paste("lam(", form.simpleLam, ")", sep = "")

    form.simpleDet <- formulaShort(mod.simple, unmarked.type = "det")
    form.simpleDet2 <- paste("p(", form.simpleDet, ")", sep = "")
    form.simple <- paste(form.simpleLam2, form.simpleDet2, sep = "")
    
    ##complex model
    form.complexLam <- formulaShort(mod.complex, unmarked.type = "state")
    form.complexLam2 <- paste("lam(", form.complexLam, ")", sep = "")
    form.complexDet <- formulaShort(mod.complex, unmarked.type = "det")
    form.complexDet2 <- paste("p(", form.complexDet, ")", sep = "")
    form.complex <- paste(form.complexLam2, form.complexDet2, sep = "")
    
    
    ##- 2 * (logLik(simple) - logLik(complex))
    LR <- -2 * (LL.simple - LL.complex)
    ##difference in number of estimated parameters
    K.diff <- K.complex - K.simple

    ##use chi-square if no overdispersion
    if(c.hat == 1) {
        Chistat <- LR
        df <- K.diff
        pval <- 1 - pchisq(Chistat, df = df)

        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Chistat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "Chistat", "pval")
    
    } else {
        ##compute F statistic
        Fstat <- LR/((K.diff)*c.hat)

        ##compute df
        df.num <- K.diff
        df.denom <- df.complex

        ##P value
        pval <- 1 - pf(Fstat, df1 = df.num, df2 = df.denom)
        
        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Fstat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "F", "pval")
    }
    
    ##assemble in list
    outList <- list(form.simple = form.simple,
                    form.complex = form.complex,
                    c.hat = c.hat,
                    devMat = devMat)
    
    class(outList) <- c("anovaOD", "list")
    return(outList)
}



##pcountOpen
anovaOD.unmarkedFitPCO <- function(mod.simple, mod.complex, c.hat = 1, nobs = NULL, ...) {

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##check for distributions
    modFamily1 <- mod.simple@mixture
    modFamily2 <- mod.complex@mixture
    if(!identical(modFamily1, modFamily2)) stop("\nComparisons only appropriate for models using the same mixture distribution\n")

    if(!identical(modFamily1, "P") && !identical(modFamily1, "ZIP")) {
        if(c.hat > 1) stop("\ndistribution not appropriate for overdispersion correction\n")
    }

    ##extract response
    y1 <- mod.simple@data@y
    y2 <- mod.complex@data@y

    ##check that data are the same
    if(!identical(y1, y2)) stop("\nData set should be identical to compare models\n")
    
    ##number of observations
    if(is.null(nobs)) {
        nobs <- nrow(y1)
    }    

    ##extract LL
    LL.simple <- logLik(mod.simple)
    LL.complex <- logLik(mod.complex)

    ##extract number of estimated parameters
    K.simple <- length(coef(mod.simple))
    K.complex <- length(coef(mod.complex))

    ##residual df of complex model
    df.complex <- nobs - K.complex

    ##extract model formula
    ##simple model
    form.simpleLam <- formulaShort(mod.simple, unmarked.type = "lambda")
    form.simpleLam2 <- paste("lam(", form.simpleLam, ")", sep = "")

    form.simpleGam <- formulaShort(mod.simple, unmarked.type = "gamma")
    form.simpleGam2 <- paste("gam(", form.simpleGam, ")", sep = "")

    form.simpleOmega <- formulaShort(mod.simple, unmarked.type = "omega")
    form.simpleOmega2 <- paste("omega(", form.simpleOmega, ")", sep = "")

    form.simpleDet <- formulaShort(mod.simple, unmarked.type = "det")
    form.simpleDet2 <- paste("p(", form.simpleDet, ")", sep = "")

    if(exists("iota", mod.simple@estimates@estimates)){
        form.simpleIota <- formulaShort(mod.simple, unmarked.type = "iota")
        form.simpleIota2 <- paste("iota(", form.simpleIota, ")", sep = "")

        form.simple <- paste(form.simpleLam2, form.simpleGam2,
                             form.simpleOmega2, form.simpleIota2,
                             form.simpleDet2, sep = "")
    } else {
        form.simple <- paste(form.simpleLam2, form.simpleGam2,
                             form.simpleOmega2, 
                             form.simpleDet2, sep = "")
    }

    
    ##complex model
    form.complexLam <- formulaShort(mod.complex, unmarked.type = "lambda")
    form.complexLam2 <- paste("lam(", form.complexLam, ")", sep = "")

    form.complexGam <- formulaShort(mod.complex, unmarked.type = "gamma")
    form.complexGam2 <- paste("gam(", form.complexGam, ")", sep = "")

    form.complexOmega <- formulaShort(mod.complex, unmarked.type = "omega")
    form.complexOmega2 <- paste("omega(", form.complexOmega, ")", sep = "")

    form.complexDet <- formulaShort(mod.complex, unmarked.type = "det")
    form.complexDet2 <- paste("p(", form.complexDet, ")", sep = "")

    if(exists("iota", mod.complex@estimates@estimates)){
        form.complexIota <- formulaShort(mod.complex, unmarked.type = "iota")
        form.complexIota2 <- paste("iota(", form.complexIota, ")", sep = "")

        form.complex <- paste(form.complexLam2, form.complexGam2,
                             form.complexOmega2, form.complexIota2,
                             form.complexDet2, sep = "")
    } else {
        form.complex <- paste(form.complexLam2, form.complexGam2,
                             form.complexOmega2, 
                             form.complexDet2, sep = "")
    }
    
    
    ##- 2 * (logLik(simple) - logLik(complex))
    LR <- -2 * (LL.simple - LL.complex)
    ##difference in number of estimated parameters
    K.diff <- K.complex - K.simple

    ##use chi-square if no overdispersion
    if(c.hat == 1) {
        Chistat <- LR
        df <- K.diff
        pval <- 1 - pchisq(Chistat, df = df)

        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Chistat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "Chistat", "pval")
    
    } else {
        ##compute F statistic
        Fstat <- LR/((K.diff)*c.hat)

        ##compute df
        df.num <- K.diff
        df.denom <- df.complex

        ##P value
        pval <- 1 - pf(Fstat, df1 = df.num, df2 = df.denom)
        
        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Fstat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "F", "pval")
    }
    
    ##assemble in list
    outList <- list(form.simple = form.simple,
                    form.complex = form.complex,
                    c.hat = c.hat,
                    devMat = devMat)
    
    class(outList) <- c("anovaOD", "list")
    return(outList)
}



##distsamp
anovaOD.unmarkedFitDS <- function(mod.simple, mod.complex, c.hat = 1, nobs = NULL, ...) {

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##extract response
    y1 <- mod.simple@data@y
    y2 <- mod.complex@data@y

    ##check that data are the same
    if(!identical(y1, y2)) stop("\nData set should be identical to compare models\n")
    
    ##number of observations
    if(is.null(nobs)) {
        nobs <- nrow(y1)
    }    

    ##extract LL
    LL.simple <- logLik(mod.simple)
    LL.complex <- logLik(mod.complex)

    ##extract number of estimated parameters
    K.simple <- length(coef(mod.simple))
    K.complex <- length(coef(mod.complex))

    ##residual df of complex model
    df.complex <- nobs - K.complex

    ##extract model formula
    ##simple model
    form.simpleLam <- formulaShort(mod.simple, unmarked.type = "state")
    form.simpleLam2 <- paste("lam(", form.simpleLam, ")", sep = "")

    form.simpleDet <- formulaShort(mod.simple, unmarked.type = "det")
    form.simpleDet2 <- paste("p(", form.simpleDet, ")", sep = "")
    key.simple <- mod.simple@keyfun
    form.simple <- paste(form.simpleLam2, form.simpleDet2, "  (", key.simple, ")", sep = "")
    
    ##complex model
    form.complexLam <- formulaShort(mod.complex, unmarked.type = "state")
    form.complexLam2 <- paste("lam(", form.complexLam, ")", sep = "")
    form.complexDet <- formulaShort(mod.complex, unmarked.type = "det")
    form.complexDet2 <- paste("p(", form.complexDet, ")", sep = "")
    key.complex <- mod.complex@keyfun
    form.complex <- paste(form.complexLam2, form.complexDet2, "  (", key.complex, ")", sep = "")
    
    
    ##- 2 * (logLik(simple) - logLik(complex))
    LR <- -2 * (LL.simple - LL.complex)
    ##difference in number of estimated parameters
    K.diff <- K.complex - K.simple

    ##use chi-square if no overdispersion
    if(c.hat == 1) {
        Chistat <- LR
        df <- K.diff
        pval <- 1 - pchisq(Chistat, df = df)

        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Chistat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "Chistat", "pval")
    
    } else {
        ##compute F statistic
        Fstat <- LR/((K.diff)*c.hat)

        ##compute df
        df.num <- K.diff
        df.denom <- df.complex

        ##P value
        pval <- 1 - pf(Fstat, df1 = df.num, df2 = df.denom)
        
        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Fstat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "F", "pval")
    }
    
    ##assemble in list
    outList <- list(form.simple = form.simple,
                    form.complex = form.complex,
                    c.hat = c.hat,
                    devMat = devMat)
    
    class(outList) <- c("anovaOD", "list")
    return(outList)
}



##gdistsamp
anovaOD.unmarkedFitGDS <- function(mod.simple, mod.complex, c.hat = 1, nobs = NULL, ...) {

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##check for distributions
    modFamily1 <- mod.simple@mixture
    modFamily2 <- mod.complex@mixture
    if(!identical(modFamily1, modFamily2)) stop("\nComparisons only appropriate for models using the same mixture distribution\n")

    if(!identical(modFamily1, "P") && !identical(modFamily1, "ZIP")) {
        if(c.hat > 1) stop("\ndistribution not appropriate for overdispersion correction\n")
    }

 
    ##extract response
    y1 <- mod.simple@data@y
    y2 <- mod.complex@data@y

    ##check that data are the same
    if(!identical(y1, y2)) stop("\nData set should be identical to compare models\n")
    
    ##number of observations
    if(is.null(nobs)) {
        nobs <- nrow(y1)
    }    

    ##extract LL
    LL.simple <- logLik(mod.simple)
    LL.complex <- logLik(mod.complex)

    ##extract number of estimated parameters
    K.simple <- length(coef(mod.simple))
    K.complex <- length(coef(mod.complex))

    ##residual df of complex model
    df.complex <- nobs - K.complex

    ##extract model formula
    ##simple model
    form.simpleLam <- formulaShort(mod.simple, unmarked.type = "lambda")
    form.simpleLam2 <- paste("lam(", form.simpleLam, ")", sep = "")

    form.simplePhi <- formulaShort(mod.simple, unmarked.type = "phi")
    form.simplePhi2 <- paste("phi(", form.simplePhi, ")", sep = "")

    form.simpleDet <- formulaShort(mod.simple, unmarked.type = "det")
    form.simpleDet2 <- paste("p(", form.simpleDet, ")", sep = "")
    key.simple <- mod.simple@keyfun
    form.simple <- paste(form.simpleLam2, form.simplePhi2,
                         form.simpleDet2, "  (", key.simple, ")",
                         sep = "")
    
    ##complex model
    form.complexLam <- formulaShort(mod.complex, unmarked.type = "lambda")
    form.complexLam2 <- paste("lam(", form.complexLam, ")", sep = "")

    form.complexPhi <- formulaShort(mod.complex, unmarked.type = "phi")
    form.complexPhi2 <- paste("phi(", form.complexPhi, ")", sep = "")

    form.complexDet <- formulaShort(mod.complex, unmarked.type = "det")
    form.complexDet2 <- paste("p(", form.complexDet, ")", sep = "")
    key.complex <- mod.complex@keyfun
    form.complex <- paste(form.complexLam2, form.complexPhi2,
                          form.complexDet2, "  (", key.complex, ")",
                          sep = "")
    
    
    ##- 2 * (logLik(simple) - logLik(complex))
    LR <- -2 * (LL.simple - LL.complex)
    ##difference in number of estimated parameters
    K.diff <- K.complex - K.simple

    ##use chi-square if no overdispersion
    if(c.hat == 1) {
        Chistat <- LR
        df <- K.diff
        pval <- 1 - pchisq(Chistat, df = df)

        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Chistat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "Chistat", "pval")
    
    } else {
        ##compute F statistic
        Fstat <- LR/((K.diff)*c.hat)

        ##compute df
        df.num <- K.diff
        df.denom <- df.complex

        ##P value
        pval <- 1 - pf(Fstat, df1 = df.num, df2 = df.denom)
        
        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Fstat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "F", "pval")
    }
    
    ##assemble in list
    outList <- list(form.simple = form.simple,
                    form.complex = form.complex,
                    c.hat = c.hat,
                    devMat = devMat)
    
    class(outList) <- c("anovaOD", "list")
    return(outList)
}



##multinomPois
anovaOD.unmarkedFitMPois <- function(mod.simple, mod.complex, c.hat = 1, nobs = NULL, ...) {

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##extract response
    y1 <- mod.simple@data@y
    y2 <- mod.complex@data@y

    ##check that data are the same
    if(!identical(y1, y2)) stop("\nData set should be identical to compare models\n")
    
    ##number of observations
    if(is.null(nobs)) {
        nobs <- nrow(y1)
    }    

    ##extract LL
    LL.simple <- logLik(mod.simple)
    LL.complex <- logLik(mod.complex)

    ##extract number of estimated parameters
    K.simple <- length(coef(mod.simple))
    K.complex <- length(coef(mod.complex))

    ##residual df of complex model
    df.complex <- nobs - K.complex

    ##extract model formula
    ##simple model
    form.simpleLam <- formulaShort(mod.simple, unmarked.type = "state")
    form.simpleLam2 <- paste("lam(", form.simpleLam, ")", sep = "")

    form.simpleDet <- formulaShort(mod.simple, unmarked.type = "det")
    form.simpleDet2 <- paste("p(", form.simpleDet, ")", sep = "")
    form.simple <- paste(form.simpleLam2, form.simpleDet2, sep = "")
    
    ##complex model
    form.complexLam <- formulaShort(mod.complex, unmarked.type = "state")
    form.complexLam2 <- paste("lam(", form.complexLam, ")", sep = "")
    form.complexDet <- formulaShort(mod.complex, unmarked.type = "det")
    form.complexDet2 <- paste("p(", form.complexDet, ")", sep = "")
    form.complex <- paste(form.complexLam2, form.complexDet2, sep = "")
    
    
    ##- 2 * (logLik(simple) - logLik(complex))
    LR <- -2 * (LL.simple - LL.complex)
    ##difference in number of estimated parameters
    K.diff <- K.complex - K.simple

    ##use chi-square if no overdispersion
    if(c.hat == 1) {
        Chistat <- LR
        df <- K.diff
        pval <- 1 - pchisq(Chistat, df = df)

        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Chistat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "Chistat", "pval")
    
    } else {
        ##compute F statistic
        Fstat <- LR/((K.diff)*c.hat)

        ##compute df
        df.num <- K.diff
        df.denom <- df.complex

        ##P value
        pval <- 1 - pf(Fstat, df1 = df.num, df2 = df.denom)
        
        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Fstat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "F", "pval")
    }
    
    ##assemble in list
    outList <- list(form.simple = form.simple,
                    form.complex = form.complex,
                    c.hat = c.hat,
                    devMat = devMat)
    
    class(outList) <- c("anovaOD", "list")
    return(outList)
}



##gmultmix
anovaOD.unmarkedFitGMM <- function(mod.simple, mod.complex, c.hat = 1, nobs = NULL, ...) {

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##check for distributions
    modFamily1 <- mod.simple@mixture
    modFamily2 <- mod.complex@mixture
    if(!identical(modFamily1, modFamily2)) stop("\nComparisons only appropriate for models using the same mixture distribution\n")

    if(!identical(modFamily1, "P") && !identical(modFamily1, "ZIP")) {
        if(c.hat > 1) stop("\ndistribution not appropriate for overdispersion correction\n\n")
    }
    
    ##extract response
    y1 <- mod.simple@data@y
    y2 <- mod.complex@data@y

    ##check that data are the same
    if(!identical(y1, y2)) stop("\nData set should be identical to compare models\n")
    
    ##number of observations
    if(is.null(nobs)) {
        nobs <- nrow(y1)
    }    

    ##extract LL
    LL.simple <- logLik(mod.simple)
    LL.complex <- logLik(mod.complex)

    ##extract number of estimated parameters
    K.simple <- length(coef(mod.simple))
    K.complex <- length(coef(mod.complex))

    ##residual df of complex model
    df.complex <- nobs - K.complex

    ##extract model formula
    ##simple model
    form.simpleLam <- formulaShort(mod.simple, unmarked.type = "lambda")
    form.simpleLam2 <- paste("lam(", form.simpleLam, ")", sep = "")

    form.simplePhi <- formulaShort(mod.simple, unmarked.type = "phi")
    form.simplePhi2 <- paste("phi(", form.simplePhi, ")", sep = "")

    form.simpleDet <- formulaShort(mod.simple, unmarked.type = "det")
    form.simpleDet2 <- paste("p(", form.simpleDet, ")", sep = "")
    form.simple <- paste(form.simpleLam2, form.simplePhi2,
                         form.simpleDet2, sep = "")
    
    ##complex model
    form.complexLam <- formulaShort(mod.complex, unmarked.type = "lambda")
    form.complexLam2 <- paste("lam(", form.complexLam, ")", sep = "")

    form.complexPhi <- formulaShort(mod.complex, unmarked.type = "phi")
    form.complexPhi2 <- paste("phi(", form.complexPhi, ")", sep = "")

    form.complexDet <- formulaShort(mod.complex, unmarked.type = "det")
    form.complexDet2 <- paste("p(", form.complexDet, ")", sep = "")
    form.complex <- paste(form.complexLam2, form.complexPhi2,
                          form.complexDet2, sep = "")
    
    
    ##- 2 * (logLik(simple) - logLik(complex))
    LR <- -2 * (LL.simple - LL.complex)
    ##difference in number of estimated parameters
    K.diff <- K.complex - K.simple

    ##use chi-square if no overdispersion
    if(c.hat == 1) {
        Chistat <- LR
        df <- K.diff
        pval <- 1 - pchisq(Chistat, df = df)

        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Chistat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "Chistat", "pval")
    
    } else {
        ##compute F statistic
        Fstat <- LR/((K.diff)*c.hat)

        ##compute df
        df.num <- K.diff
        df.denom <- df.complex

        ##P value
        pval <- 1 - pf(Fstat, df1 = df.num, df2 = df.denom)
        
        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Fstat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "F", "pval")
    }
    
    ##assemble in list
    outList <- list(form.simple = form.simple,
                    form.complex = form.complex,
                    c.hat = c.hat,
                    devMat = devMat)
    
    class(outList) <- c("anovaOD", "list")
    return(outList)
}



##gpcount
anovaOD.unmarkedFitGPC <- function(mod.simple, mod.complex, c.hat = 1, nobs = NULL, ...) {

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##check for distributions
    modFamily1 <- mod.simple@mixture
    modFamily2 <- mod.complex@mixture
    if(!identical(modFamily1, modFamily2)) stop("\nComparisons only appropriate for models using the same mixture distribution\n")

    if(!identical(modFamily1, "P") && !identical(modFamily1, "ZIP")) {
        if(c.hat > 1) stop("\ndistribution not appropriate for overdispersion correction\n\n")
    }
    
    ##extract response
    y1 <- mod.simple@data@y
    y2 <- mod.complex@data@y

    ##check that data are the same
    if(!identical(y1, y2)) stop("\nData set should be identical to compare models\n")
    
    ##number of observations
    if(is.null(nobs)) {
        nobs <- nrow(y1)
    }    

    ##extract LL
    LL.simple <- logLik(mod.simple)
    LL.complex <- logLik(mod.complex)

    ##extract number of estimated parameters
    K.simple <- length(coef(mod.simple))
    K.complex <- length(coef(mod.complex))

    ##residual df of complex model
    df.complex <- nobs - K.complex

    ##extract model formula
    ##simple model
    form.simpleLam <- formulaShort(mod.simple, unmarked.type = "lambda")
    form.simpleLam2 <- paste("lam(", form.simpleLam, ")", sep = "")

    form.simplePhi <- formulaShort(mod.simple, unmarked.type = "phi")
    form.simplePhi2 <- paste("phi(", form.simplePhi, ")", sep = "")

    form.simpleDet <- formulaShort(mod.simple, unmarked.type = "det")
    form.simpleDet2 <- paste("p(", form.simpleDet, ")", sep = "")
    form.simple <- paste(form.simpleLam2, form.simplePhi2,
                         form.simpleDet2, sep = "")
    
    ##complex model
    form.complexLam <- formulaShort(mod.complex, unmarked.type = "lambda")
    form.complexLam2 <- paste("lam(", form.complexLam, ")", sep = "")

    form.complexPhi <- formulaShort(mod.complex, unmarked.type = "phi")
    form.complexPhi2 <- paste("phi(", form.complexPhi, ")", sep = "")

    form.complexDet <- formulaShort(mod.complex, unmarked.type = "det")
    form.complexDet2 <- paste("p(", form.complexDet, ")", sep = "")
    form.complex <- paste(form.complexLam2, form.complexPhi2,
                          form.complexDet2, sep = "")
    
    
    ##- 2 * (logLik(simple) - logLik(complex))
    LR <- -2 * (LL.simple - LL.complex)
    ##difference in number of estimated parameters
    K.diff <- K.complex - K.simple

    ##use chi-square if no overdispersion
    if(c.hat == 1) {
        Chistat <- LR
        df <- K.diff
        pval <- 1 - pchisq(Chistat, df = df)

        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Chistat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "Chistat", "pval")
    
    } else {
        ##compute F statistic
        Fstat <- LR/((K.diff)*c.hat)

        ##compute df
        df.num <- K.diff
        df.denom <- df.complex

        ##P value
        pval <- 1 - pf(Fstat, df1 = df.num, df2 = df.denom)
        
        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Fstat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "F", "pval")
    }
    
    ##assemble in list
    outList <- list(form.simple = form.simple,
                    form.complex = form.complex,
                    c.hat = c.hat,
                    devMat = devMat)
    
    class(outList) <- c("anovaOD", "list")
    return(outList)
}



##occuMulti
anovaOD.unmarkedFitOccuMulti <- function(mod.simple, mod.complex, c.hat = 1, nobs = NULL, ...) {

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##extract response
    y1 <- mod.simple@data@y
    y2 <- mod.complex@data@y

    ##check that data are the same
    if(!identical(y1, y2)) stop("\nData set should be identical to compare models\n")
    
    ##number of observations
    if(is.null(nobs)) {
        nobs <- nrow(y1)
    }    

    ##number of species
    nspecies <- length(mod.simple@data@ylist)
    genericNames <- paste("sp", 1:nspecies, sep = "")
    
    ##extract LL
    LL.simple <- logLik(mod.simple)
    LL.complex <- logLik(mod.complex)

    ##extract number of estimated parameters
    K.simple <- length(coef(mod.simple))
    K.complex <- length(coef(mod.complex))

    ##residual df of complex model
    df.complex <- nobs - K.complex

    ##extract model formula
    ##simple model
    ##extract labels of fDesign
    simple.fDesign <- mod.simple@data@fDesign
    simple.colNames <- substr(x = colnames(simple.fDesign),
                              start = 1, stop = 2)
    ##extract state formulas
    simpleState <- mod.simple@stateformulas
    ##replace ~ 1 by "."
    simplerState <- gsub(pattern = "~1", replacement = ".", x = simpleState)
    form.simplePsi <- paste(simple.colNames, "(", simplerState, ")", sep = "")
    form.simplePsi2 <- paste(form.simplePsi, collapse = "")
    form.simplePsi3 <- paste("psi[", form.simplePsi2, "]", collapse = "")

    ##extract detection formulas
    simpleDet <- mod.simple@detformulas
    ##replace ~ 1 by "."
    simplerDet <- gsub(pattern = "~1", replacement = ".", x = simpleDet)
    form.simpleDet <- paste(genericNames, "(", simplerDet, ")", sep = "")
    form.simpleDet2 <- paste(form.simpleDet, collapse = "")
    form.simpleDet3 <- paste("p[", form.simpleDet2, "]", collapse = "")

    form.simple <- paste(form.simplePsi3, form.simpleDet3, sep = "")

    ##complex model
    ##extract labels of fDesign
    complex.fDesign <- mod.complex@data@fDesign
    complex.colNames <- substr(x = colnames(complex.fDesign),
                              start = 1, stop = 2)
    ##extract state formulas
    complexState <- mod.complex@stateformulas
    ##replace ~ 1 by "."
    complexrState <- gsub(pattern = "~1", replacement = ".", x = complexState)
    form.complexPsi <- paste(complex.colNames, "(", complexrState, ")", sep = "")
    form.complexPsi2 <- paste(form.complexPsi, collapse = "")
    form.complexPsi3 <- paste("psi[", form.complexPsi2, "]", collapse = "")

    ##extract detection formulas
    complexDet <- mod.complex@detformulas
    ##replace ~ 1 by "."
    complexrDet <- gsub(pattern = "~1", replacement = ".", x = complexDet)
    form.complexDet <- paste(genericNames, "(", complexrDet, ")", sep = "")
    form.complexDet2 <- paste(form.complexDet, collapse = "")
    form.complexDet3 <- paste("p[", form.complexDet2, "]", collapse = "")

    form.complex <- paste(form.complexPsi3, form.complexDet3, sep = "")
    
    
    ##- 2 * (logLik(simple) - logLik(complex))
    LR <- -2 * (LL.simple - LL.complex)
    ##difference in number of estimated parameters
    K.diff <- K.complex - K.simple

    ##use chi-square if no overdispersion
    if(c.hat == 1) {
        Chistat <- LR
        df <- K.diff
        pval <- 1 - pchisq(Chistat, df = df)

        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Chistat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "Chistat", "pval")
    
    } else {
        ##compute F statistic
        Fstat <- LR/((K.diff)*c.hat)

        ##compute df
        df.num <- K.diff
        df.denom <- df.complex

        ##P value
        pval <- 1 - pf(Fstat, df1 = df.num, df2 = df.denom)
        
        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Fstat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "F", "pval")
    }
    
    ##assemble in list
    outList <- list(form.simple = form.simple,
                    form.complex = form.complex,
                    c.hat = c.hat,
                    devMat = devMat)
    
    class(outList) <- c("anovaOD", "list")
    return(outList)
}



##occuMS
anovaOD.unmarkedFitOccuMS <- function(mod.simple, mod.complex, c.hat = 1, nobs = NULL, ...) {

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##extract response
    y1 <- mod.simple@data@y
    y2 <- mod.complex@data@y

    ##check for each model single season vs dynamic
    nseason.simple <- mod.simple@data@numPrimary
    nseason.complex <- mod.complex@data@numPrimary
   
    ##check that data are the same
    if(!identical(y1, y2)) stop("\nData set should be identical to compare models\n")
    
    ##number of observations
    if(is.null(nobs)) {
        nobs <- nrow(y1)
    }    

    ##extract LL
    LL.simple <- logLik(mod.simple)
    LL.complex <- logLik(mod.complex)

    ##extract number of estimated parameters
    K.simple <- length(coef(mod.simple))
    K.complex <- length(coef(mod.complex))

    ##residual df of complex model
    df.complex <- nobs - K.complex

    ##extract model formula
    ##simple model
    if(nseason.simple == 1) {
    form.simplePsi <- formulaShort(mod.simple, unmarked.type = "state")
    form.simplePsi2 <- paste("psi(", form.simplePsi, ")", sep = "")

    form.simpleDet <- formulaShort(mod.simple, unmarked.type = "det")
    form.simpleDet2 <- paste("p(", form.simpleDet, ")", sep = "")
    form.simple <- paste(form.simplePsi2, form.simpleDet2, sep = "")
    } else {
        
        form.simplePsi <- formulaShort(mod.simple, unmarked.type = "state")
        form.simplePsi2 <- paste("psi(", form.simplePsi, ")", sep = "")

        form.simplePhi <- formulaShort(mod.simple, unmarked.type = "transition")
        form.simplePhi2 <- paste("phi(", form.simplePhi, ")", sep = "")
        
        form.simpleDet <- formulaShort(mod.simple, unmarked.type = "det")
        form.simpleDet2 <- paste("p(", form.simpleDet, ")", sep = "")
        form.simple <- paste(form.simplePsi2, form.simplePhi2, form.simpleDet2, sep = "")
        
    }

    
    ##complex model
    if(nseason.complex == 1) {
        form.complexPsi <- formulaShort(mod.complex, unmarked.type = "state")
        form.complexPsi2 <- paste("psi(", form.complexPsi, ")", sep = "")
        
        form.complexDet <- formulaShort(mod.complex, unmarked.type = "det")
        form.complexDet2 <- paste("p(", form.complexDet, ")", sep = "")
        form.complex <- paste(form.complexPsi2, form.complexDet2, sep = "")
    } else {

        form.complexPsi <- formulaShort(mod.complex, unmarked.type = "state")
        form.complexPsi2 <- paste("psi(", form.complexPsi, ")", sep = "")

        form.complexPhi <- formulaShort(mod.complex, unmarked.type = "transition")
        form.complexPhi2 <- paste("phi(", form.complexPhi, ")", sep = "")

        form.complexDet <- formulaShort(mod.complex, unmarked.type = "det")
        form.complexDet2 <- paste("p(", form.complexDet, ")", sep = "")
        form.complex <- paste(form.complexPsi2, form.complexPhi2, form.complexDet2, sep = "")
        
    }

        
    ##- 2 * (logLik(simple) - logLik(complex))
    LR <- -2 * (LL.simple - LL.complex)
    ##difference in number of estimated parameters
    K.diff <- K.complex - K.simple

    ##use chi-square if no overdispersion
    if(c.hat == 1) {
        Chistat <- LR
        df <- K.diff
        pval <- 1 - pchisq(Chistat, df = df)

        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Chistat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "Chistat", "pval")
    
    } else {
        ##compute F statistic
        Fstat <- LR/((K.diff)*c.hat)

        ##compute df
        df.num <- K.diff
        df.denom <- df.complex

        ##P value
        pval <- 1 - pf(Fstat, df1 = df.num, df2 = df.denom)
        
        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Fstat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "F", "pval")
    }
    
    ##assemble in list
    outList <- list(form.simple = form.simple,
                    form.complex = form.complex,
                    c.hat = c.hat,
                    devMat = devMat)
    
    class(outList) <- c("anovaOD", "list")
    return(outList)
}



##occuTTD
anovaOD.unmarkedFitOccuTTD <- function(mod.simple, mod.complex, c.hat = 1, nobs = NULL, ...) {

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##extract response
    y1 <- mod.simple@data@y
    y2 <- mod.complex@data@y

    ##check for each model single season vs dynamic
    nseason.simple <- mod.simple@data@numPrimary
    nseason.complex <- mod.complex@data@numPrimary
    
    ##check that data are the same
    if(!identical(y1, y2)) stop("\nData set should be identical to compare models\n")
    
    ##number of observations
    if(is.null(nobs)) {
        nobs <- nrow(y1)
    }    

    ##extract LL
    LL.simple <- logLik(mod.simple)
    LL.complex <- logLik(mod.complex)

    ##extract number of estimated parameters
    K.simple <- length(coef(mod.simple))
    K.complex <- length(coef(mod.complex))

    ##residual df of complex model
    df.complex <- nobs - K.complex

    ##extract model formula
    ##simple model
    if(nseason.simple == 1) { 
        form.simplePsi <- formulaShort(mod.simple, unmarked.type = "psi")
        form.simplePsi2 <- paste("psi(", form.simplePsi, ")", sep = "")

        form.simpleDet <- formulaShort(mod.simple, unmarked.type = "det")
        form.simpleDet2 <- paste("p(", form.simpleDet, ")", sep = "")
        form.simple <- paste(form.simplePsi2, form.simpleDet2, sep = "")
    } else {

        form.simplePsi <- formulaShort(mod.simple, unmarked.type = "psi")
        form.simplePsi2 <- paste("psi(", form.simplePsi, ")", sep = "")

        form.simpleGam <- formulaShort(mod.simple, unmarked.type = "col")
        form.simpleGam2 <- paste("gam(", form.simpleGam, ")", sep = "")

        form.simpleEps <- formulaShort(mod.simple, unmarked.type = "ext")
        form.simpleEps2 <- paste("eps(", form.simpleEps, ")", sep = "")

        form.simpleDet <- formulaShort(mod.simple, unmarked.type = "det")
        form.simpleDet2 <- paste("p(", form.simpleDet, ")", sep = "")
        form.simple <- paste(form.simplePsi2, form.simpleGam2,
                             form.simpleEps2, form.simpleDet2, sep = "")
    }

    ##complex model
    if(nseason.complex == 1) {
        form.complexPsi <- formulaShort(mod.complex, unmarked.type = "psi")
        form.complexPsi2 <- paste("psi(", form.complexPsi, ")", sep = "")
        form.complexDet <- formulaShort(mod.complex, unmarked.type = "det")
        form.complexDet2 <- paste("p(", form.complexDet, ")", sep = "")
        form.complex <- paste(form.complexPsi2, form.complexDet2, sep = "")
    } else {

        form.complexPsi <- formulaShort(mod.complex, unmarked.type = "psi")
        form.complexPsi2 <- paste("psi(", form.complexPsi, ")", sep = "")

        form.complexGam <- formulaShort(mod.complex, unmarked.type = "col")
        form.complexGam2 <- paste("gam(", form.complexGam, ")", sep = "")

        form.complexEps <- formulaShort(mod.complex, unmarked.type = "ext")
        form.complexEps2 <- paste("eps(", form.complexEps, ")", sep = "")

        form.complexDet <- formulaShort(mod.complex, unmarked.type = "det")
        form.complexDet2 <- paste("p(", form.complexDet, ")", sep = "")
        form.complex <- paste(form.complexPsi2, form.complexGam2,
                              form.complexEps2, form.complexDet2, sep = "")

    }
    
    
    ##- 2 * (logLik(simple) - logLik(complex))
    LR <- -2 * (LL.simple - LL.complex)
    ##difference in number of estimated parameters
    K.diff <- K.complex - K.simple

    ##use chi-square if no overdispersion
    if(c.hat == 1) {
        Chistat <- LR
        df <- K.diff
        pval <- 1 - pchisq(Chistat, df = df)

        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Chistat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "Chistat", "pval")
    
    } else {
        ##compute F statistic
        Fstat <- LR/((K.diff)*c.hat)

        ##compute df
        df.num <- K.diff
        df.denom <- df.complex

        ##P value
        pval <- 1 - pf(Fstat, df1 = df.num, df2 = df.denom)
        
        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Fstat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "F", "pval")
    }
    
    ##assemble in list
    outList <- list(form.simple = form.simple,
                    form.complex = form.complex,
                    c.hat = c.hat,
                    devMat = devMat)
    
    class(outList) <- c("anovaOD", "list")
    return(outList)
}


##multmixOpen
anovaOD.unmarkedFitMMO <- function(mod.simple, mod.complex, c.hat = 1, nobs = NULL, ...) {

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##check for distributions
    modFamily1 <- mod.simple@mixture
    modFamily2 <- mod.complex@mixture
    if(!identical(modFamily1, modFamily2)) stop("\nComparisons only appropriate for models using the same mixture distribution\n")

    if(!identical(modFamily1, "P") && !identical(modFamily1, "ZIP")) {
        if(c.hat > 1) stop("\ndistribution not appropriate for overdispersion correction\n")
    }

    ##extract response
    y1 <- mod.simple@data@y
    y2 <- mod.complex@data@y

    ##check that data are the same
    if(!identical(y1, y2)) stop("\nData set should be identical to compare models\n")
    
    ##number of observations
    if(is.null(nobs)) {
        nobs <- nrow(y1)
    }    

    ##extract LL
    LL.simple <- logLik(mod.simple)
    LL.complex <- logLik(mod.complex)

    ##extract number of estimated parameters
    K.simple <- length(coef(mod.simple))
    K.complex <- length(coef(mod.complex))

    ##residual df of complex model
    df.complex <- nobs - K.complex

    ##extract model formula
    ##simple model
    form.simpleLam <- formulaShort(mod.simple, unmarked.type = "lambda")
    form.simpleLam2 <- paste("lam(", form.simpleLam, ")", sep = "")

    form.simpleGam <- formulaShort(mod.simple, unmarked.type = "gamma")
    form.simpleGam2 <- paste("gam(", form.simpleGam, ")", sep = "")

    form.simpleOmega <- formulaShort(mod.simple, unmarked.type = "omega")
    form.simpleOmega2 <- paste("omega(", form.simpleOmega, ")", sep = "")

    form.simpleDet <- formulaShort(mod.simple, unmarked.type = "det")
    form.simpleDet2 <- paste("p(", form.simpleDet, ")", sep = "")

    if(exists("iota", mod.simple@estimates@estimates)){
        form.simpleIota <- formulaShort(mod.simple, unmarked.type = "iota")
        form.simpleIota2 <- paste("iota(", form.simpleIota, ")", sep = "")

        form.simple <- paste(form.simpleLam2, form.simpleGam2,
                             form.simpleOmega2, form.simpleIota2,
                             form.simpleDet2, sep = "")
    } else {
        form.simple <- paste(form.simpleLam2, form.simpleGam2,
                             form.simpleOmega2, 
                             form.simpleDet2, sep = "")
    }

    
    ##complex model
    form.complexLam <- formulaShort(mod.complex, unmarked.type = "lambda")
    form.complexLam2 <- paste("lam(", form.complexLam, ")", sep = "")

    form.complexGam <- formulaShort(mod.complex, unmarked.type = "gamma")
    form.complexGam2 <- paste("gam(", form.complexGam, ")", sep = "")

    form.complexOmega <- formulaShort(mod.complex, unmarked.type = "omega")
    form.complexOmega2 <- paste("omega(", form.complexOmega, ")", sep = "")

    form.complexDet <- formulaShort(mod.complex, unmarked.type = "det")
    form.complexDet2 <- paste("p(", form.complexDet, ")", sep = "")

    if(exists("iota", mod.complex@estimates@estimates)){
        form.complexIota <- formulaShort(mod.complex, unmarked.type = "iota")
        form.complexIota2 <- paste("iota(", form.complexIota, ")", sep = "")

        form.complex <- paste(form.complexLam2, form.complexGam2,
                             form.complexOmega2, form.complexIota2,
                             form.complexDet2, sep = "")
    } else {
        form.complex <- paste(form.complexLam2, form.complexGam2,
                             form.complexOmega2, 
                             form.complexDet2, sep = "")
    }
    
    
    ##- 2 * (logLik(simple) - logLik(complex))
    LR <- -2 * (LL.simple - LL.complex)
    ##difference in number of estimated parameters
    K.diff <- K.complex - K.simple

    ##use chi-square if no overdispersion
    if(c.hat == 1) {
        Chistat <- LR
        df <- K.diff
        pval <- 1 - pchisq(Chistat, df = df)

        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Chistat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "Chistat", "pval")
    
    } else {
        ##compute F statistic
        Fstat <- LR/((K.diff)*c.hat)

        ##compute df
        df.num <- K.diff
        df.denom <- df.complex

        ##P value
        pval <- 1 - pf(Fstat, df1 = df.num, df2 = df.denom)
        
        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Fstat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "F", "pval")
    }
    
    ##assemble in list
    outList <- list(form.simple = form.simple,
                    form.complex = form.complex,
                    c.hat = c.hat,
                    devMat = devMat)
    
    class(outList) <- c("anovaOD", "list")
    return(outList)
}



##distsampOpen
anovaOD.unmarkedFitDSO <- function(mod.simple, mod.complex, c.hat = 1, nobs = NULL, ...) {

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##check for distributions
    modFamily1 <- mod.simple@mixture
    modFamily2 <- mod.complex@mixture
    if(!identical(modFamily1, modFamily2)) stop("\nComparisons only appropriate for models using the same mixture distribution\n")

    if(!identical(modFamily1, "P") && !identical(modFamily1, "ZIP")) {
        if(c.hat > 1) stop("\ndistribution not appropriate for overdispersion correction\n")
    }

    ##extract response
    y1 <- mod.simple@data@y
    y2 <- mod.complex@data@y

    ##check that data are the same
    if(!identical(y1, y2)) stop("\nData set should be identical to compare models\n")
    
    ##number of observations
    if(is.null(nobs)) {
        nobs <- nrow(y1)
    }    

    ##extract LL
    LL.simple <- logLik(mod.simple)
    LL.complex <- logLik(mod.complex)

    ##extract number of estimated parameters
    K.simple <- length(coef(mod.simple))
    K.complex <- length(coef(mod.complex))

    ##residual df of complex model
    df.complex <- nobs - K.complex

    ##extract model formula
    ##simple model
    form.simpleLam <- formulaShort(mod.simple, unmarked.type = "lambda")
    form.simpleLam2 <- paste("lam(", form.simpleLam, ")", sep = "")

    form.simpleGam <- formulaShort(mod.simple, unmarked.type = "gamma")
    form.simpleGam2 <- paste("gam(", form.simpleGam, ")", sep = "")

    form.simpleOmega <- formulaShort(mod.simple, unmarked.type = "omega")
    form.simpleOmega2 <- paste("omega(", form.simpleOmega, ")", sep = "")

    form.simpleDet <- formulaShort(mod.simple, unmarked.type = "det")
    form.simpleDet2 <- paste("p(", form.simpleDet, ")", sep = "")

    if(exists("iota", mod.simple@estimates@estimates)){
        form.simpleIota <- formulaShort(mod.simple, unmarked.type = "iota")
        form.simpleIota2 <- paste("iota(", form.simpleIota, ")", sep = "")

        form.simple <- paste(form.simpleLam2, form.simpleGam2,
                             form.simpleOmega2, form.simpleIota2,
                             form.simpleDet2, sep = "")
    } else {
        form.simple <- paste(form.simpleLam2, form.simpleGam2,
                             form.simpleOmega2, 
                             form.simpleDet2, sep = "")
    }

    
    ##complex model
    form.complexLam <- formulaShort(mod.complex, unmarked.type = "lambda")
    form.complexLam2 <- paste("lam(", form.complexLam, ")", sep = "")

    form.complexGam <- formulaShort(mod.complex, unmarked.type = "gamma")
    form.complexGam2 <- paste("gam(", form.complexGam, ")", sep = "")

    form.complexOmega <- formulaShort(mod.complex, unmarked.type = "omega")
    form.complexOmega2 <- paste("omega(", form.complexOmega, ")", sep = "")

    form.complexDet <- formulaShort(mod.complex, unmarked.type = "det")
    form.complexDet2 <- paste("p(", form.complexDet, ")", sep = "")

    if(exists("iota", mod.complex@estimates@estimates)){
        form.complexIota <- formulaShort(mod.complex, unmarked.type = "iota")
        form.complexIota2 <- paste("iota(", form.complexIota, ")", sep = "")

        form.complex <- paste(form.complexLam2, form.complexGam2,
                             form.complexOmega2, form.complexIota2,
                             form.complexDet2, sep = "")
    } else {
        form.complex <- paste(form.complexLam2, form.complexGam2,
                             form.complexOmega2, 
                             form.complexDet2, sep = "")
    }
    
    
    ##- 2 * (logLik(simple) - logLik(complex))
    LR <- -2 * (LL.simple - LL.complex)
    ##difference in number of estimated parameters
    K.diff <- K.complex - K.simple

    ##use chi-square if no overdispersion
    if(c.hat == 1) {
        Chistat <- LR
        df <- K.diff
        pval <- 1 - pchisq(Chistat, df = df)

        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Chistat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "Chistat", "pval")
    
    } else {
        ##compute F statistic
        Fstat <- LR/((K.diff)*c.hat)

        ##compute df
        df.num <- K.diff
        df.denom <- df.complex

        ##P value
        pval <- 1 - pf(Fstat, df1 = df.num, df2 = df.denom)
        
        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Fstat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "F", "pval")
    }
    
    ##assemble in list
    outList <- list(form.simple = form.simple,
                    form.complex = form.complex,
                    c.hat = c.hat,
                    devMat = devMat)
    
    class(outList) <- c("anovaOD", "list")
    return(outList)
}



##goccu
anovaOD.unmarkedFitGOccu <- function(mod.simple, mod.complex, c.hat = 1, nobs = NULL, ...) {

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##check for distributions
    #modFamily1 <- mod.simple@mixture
    #modFamily2 <- mod.complex@mixture
    #if(!identical(modFamily1, modFamily2)) stop("\nComparisons only appropriate for models using the same mixture distribution\n")

    ##extract response
    y1 <- mod.simple@data@y
    y2 <- mod.complex@data@y

    ##check that data are the same
    if(!identical(y1, y2)) stop("\nData set should be identical to compare models\n")
    
    ##number of observations
    if(is.null(nobs)) {
        nobs <- nrow(y1)
    }    

    ##extract LL
    LL.simple <- logLik(mod.simple)
    LL.complex <- logLik(mod.complex)

    ##extract number of estimated parameters
    K.simple <- length(coef(mod.simple))
    K.complex <- length(coef(mod.complex))

    ##residual df of complex model
    df.complex <- nobs - K.complex

    ##extract model formula
    ##simple model
    form.simplePsi <- formulaShort(mod.simple, unmarked.type = "psi")
    form.simplePsi2 <- paste("psi(", form.simplePsi, ")", sep = "")

    form.simplePhi <- formulaShort(mod.simple, unmarked.type = "phi")
    form.simplePhi2 <- paste("phi(", form.simplePhi, ")", sep = "")

    form.simpleDet <- formulaShort(mod.simple, unmarked.type = "det")
    form.simpleDet2 <- paste("p(", form.simpleDet, ")", sep = "")
    form.simple <- paste(form.simplePsi2, form.simplePhi2,
                         form.simpleDet2, sep = "")
    
    ##complex model
    form.complexPsi <- formulaShort(mod.complex, unmarked.type = "psi")
    form.complexPsi2 <- paste("psi(", form.complexPsi, ")", sep = "")

    form.complexPhi <- formulaShort(mod.complex, unmarked.type = "phi")
    form.complexPhi2 <- paste("phi(", form.complexPhi, ")", sep = "")

    form.complexDet <- formulaShort(mod.complex, unmarked.type = "det")
    form.complexDet2 <- paste("p(", form.complexDet, ")", sep = "")
    form.complex <- paste(form.complexPsi2, form.complexPhi2,
                          form.complexDet2, sep = "")
    
    
    ##- 2 * (logLik(simple) - logLik(complex))
    LR <- -2 * (LL.simple - LL.complex)
    ##difference in number of estimated parameters
    K.diff <- K.complex - K.simple

    ##use chi-square if no overdispersion
    if(c.hat == 1) {
        Chistat <- LR
        df <- K.diff
        pval <- 1 - pchisq(Chistat, df = df)

        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Chistat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "Chistat", "pval")
    
    } else {
        ##compute F statistic
        Fstat <- LR/((K.diff)*c.hat)

        ##compute df
        df.num <- K.diff
        df.denom <- df.complex

        ##P value
        pval <- 1 - pf(Fstat, df1 = df.num, df2 = df.denom)
        
        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Fstat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "F", "pval")
    }
    
    ##assemble in list
    outList <- list(form.simple = form.simple,
                    form.complex = form.complex,
                    c.hat = c.hat,
                    devMat = devMat)
    
    class(outList) <- c("anovaOD", "list")
    return(outList)
}


##occuComm
anovaOD.unmarkedFitOccuComm <- function(mod.simple, mod.complex, c.hat = 1, nobs = NULL, ...) {

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##check for distributions
    ##modFamily1 <- mod.simple@mixture
    ##modFamily2 <- mod.complex@mixture
    ##if(!identical(modFamily1, modFamily2)) stop("\nComparisons only appropriate for models using the same mixture distribution\n")

    ##extract response
    y1 <- mod.simple@data@y
    y2 <- mod.complex@data@y

    ##check that data are the same
    if(!identical(y1, y2)) stop("\nData set should be identical to compare models\n")
    
    ##number of observations
    if(is.null(nobs)) {
        nobs <- nrow(y1)
    }    

    ##extract LL
    LL.simple <- logLik(mod.simple)
    LL.complex <- logLik(mod.complex)

    ##extract number of estimated parameters
    K.simple <- length(coef(mod.simple))
    K.complex <- length(coef(mod.complex))

    ##residual df of complex model
    df.complex <- nobs - K.complex

    ##extract model formula
    ##simple model
    ##form.simplePsi <- formulaShort(mod.simple, unmarked.type = "state")
    form.simplePsi <- names(coef(mod.simple, type = "state"))
    ##form.simplePsi2 <- paste("psi(", form.simplePsi, ")", sep = "")
    form.simplePsi2 <- paste(form.simplePsi, collapse = "+")

    ##form.simpleDet <- formulaShort(mod.simple, unmarked.type = "det")
    form.simpleDet <- names(coef(mod.simple, type = "det"))
    ##form.simpleDet2 <- paste("p(", form.simpleDet, ")", sep = "")
    form.simpleDet2 <- paste(form.simpleDet, collapse = "+")
    form.simple <- paste(form.simplePsi2, form.simpleDet2,
                         sep = "")
    
    ##complex model
    ##form.complexPsi <- formulaShort(mod.complex, unmarked.type = "psi")
    form.complexPsi <- names(coef(mod.complex, type = "state"))
    ##form.complexPsi2 <- paste("psi(", form.complexPsi, ")", sep = "")
    form.complexPsi2 <- paste(form.complexPsi, collapse = "+")

    ##form.complexDet <- formulaShort(mod.complex, unmarked.type = "det")
    form.complexDet <- names(coef(mod.complex, type = "det"))
    ##form.complexDet2 <- paste("p(", form.complexDet, ")", sep = "")
    form.complexDet2 <- paste(form.complexDet, collapse = "+")
    form.complex <- paste(form.complexPsi2, form.complexDet2,
                          sep = "")
    
    
    ##- 2 * (logLik(simple) - logLik(complex))
    LR <- -2 * (LL.simple - LL.complex)
    ##difference in number of estimated parameters
    K.diff <- K.complex - K.simple

    ##use chi-square if no overdispersion
    if(c.hat == 1) {
        Chistat <- LR
        df <- K.diff
        pval <- 1 - pchisq(Chistat, df = df)

        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Chistat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "Chistat", "pval")
    
    } else {
        ##compute F statistic
        Fstat <- LR/((K.diff)*c.hat)

        ##compute df
        df.num <- K.diff
        df.denom <- df.complex

        ##P value
        pval <- 1 - pf(Fstat, df1 = df.num, df2 = df.denom)
        
        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Fstat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "F", "pval")
    }
    
    ##assemble in list
    outList <- list(form.simple = form.simple,
                    form.complex = form.complex,
                    c.hat = c.hat,
                    devMat = devMat)
    
    class(outList) <- c("anovaOD", "list")
    return(outList)
}



##maxlike
anovaOD.maxlikeFit <- function(mod.simple, mod.complex, c.hat = 1, nobs = NULL, ...){

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##response variable
    y1 <- mod.simple$points.retained
    y2 <- mod.complex$points.retained

    ##number of observations
    if(is.null(nobs)) {
        nobs <- nrow(y1)
    }
    
    ##check that data are the same
    if(!identical(y1, y2)) stop("\nData set should be identical to compare models\n")
  
    ##extract log-likelihood
    LL1 <- logLik(mod.simple)
    LL2 <- logLik(mod.complex)
    LL.simple <- LL1[1]
    LL.complex <- LL2[1]

    ##extract number of estimated parameters
    K.simple <- length(coef(mod.simple))
    K.complex <- length(coef(mod.complex))
    
    ##residual df of complex model
    df.complex <- nobs - K.complex

    ##extract model formula
    simpleForm <- formula(mod.simple)
    form.simple <- paste("y ~", simpleForm[2])
    complexForm <- formula(mod.complex)
    form.complex <- paste("y ~", complexForm[2])
    
    ##- 2 * (logLik(simple) - logLik(complex))
    LR <- -2 * (LL.simple - LL.complex)
    ##difference in number of estimated parameters
    K.diff <- K.complex - K.simple

    ##use chi-square if no overdispersion
    if(c.hat == 1) {
        Chistat <- LR
        df <- K.diff
        pval <- 1 - pchisq(Chistat, df = df)

        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Chistat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "Chistat", "pval")
    
    } else {
        ##compute F statistic
        Fstat <- LR/((K.diff)*c.hat)

        ##compute df
        df.num <- K.diff
        df.denom <- df.complex

        ##P value
        pval <- 1 - pf(Fstat, df1 = df.num, df2 = df.denom)
        
        devMat <- matrix(data = c(K.simple, K.complex,
                                  LL.simple, LL.complex,
                                  NA, K.diff,
                                  NA, LR,
                                  NA, Fstat,
                                  NA, pval),
                         nrow = 2, ncol = 6)
        colnames(devMat) <- c("K", "logLik",
                              "Kdiff", "-2LL", "F", "pval")
    }
    
    ##assemble in list
    outList <- list(form.simple = form.simple,
                    form.complex = form.complex,
                    c.hat = c.hat,
                    devMat = devMat)
    class(outList) <- c("anovaOD", "list")
    return(outList)
}



##residual df of most complex model is required to compute F statistic
##this value is difficult to obtain for models of unmarked classes and glmm's
##potentially give option to user of specifying effective sample size (as in AICc)
##nobs



print.anovaOD <- function(x, digits = 4, ...) {

    ##extract information
    simple.text <- x$form.simple
    complex.text <- x$form.complex

    ##if text too long, truncate name
    if(nchar(simple.text) > 70) {
        simple.text <- paste(substr(x = simple.text,
                                    start = 1,
                                    stop = 70),
                             " ...", sep = "")
    }

    ##if text too long, truncate name
    if(nchar(complex.text) > 70) {
        complex.text <- paste(substr(x = complex.text,
                                     start = 1,
                                     stop = 70),
                              " ...", sep = "")
    }
        
    c.hat <- x$c.hat
    outMat <- x$devMat
    rownames(outMat) <- c("1", "2")

    ##check value of P and set to < 0.0001
    #pval <- ifelse(outMat[, "pval"] < 0.0001, "< 0.0001", outMat[, "pval"]) 
    
    ##replace NA's with spaces
    #outMat[is.na(outMat)] <- " "

    if(c.hat == 1) {
        ##add names
        colnames(outMat) <- c("K", "logLik",
                              "Kdiff", "-2logLik", "Chisq", "Pr(>Chisq)")

        cat("\nAnalysis of deviance table\n",
            sep = "")

    } else {
        ##add names
        colnames(outMat) <- c("K", "logLik",
                              "Kdiff", "-2logLik", "F", "Pr(>F)")
        
        cat("\nAnalysis of deviance table corrected for overdispersion\n",
            sep = "")
    }
    
    cat("\nSimple model  (1): ", simple.text,
        "\nComplex model (2): ", complex.text, "\n\n")
    
    printCoefmat(outMat, digits = digits, signif.stars = FALSE, cs.ind = 1,
                 na.print = "")
    
    if(c.hat > 1) { 
        cat("\n(c-hat = ", c.hat, ")", "\n", sep = "")
    }
    cat("\n")
}

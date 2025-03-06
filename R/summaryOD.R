##generic
summaryOD <- function(mod, c.hat = 1, conf.level = 0.95, 
                      out.type = "confint", ...){
    UseMethod("summaryOD", mod)
}



summaryOD.default <- function(mod, c.hat = 1, conf.level = 0.95, 
                              out.type = "confint", ...){
    stop("\nFunction not yet defined for this object class\n")
}



##summaryOD: summary with overdispersion to display CI or P-values
summaryOD.glm <- function(mod, c.hat = 1,
                          conf.level = 0.95,
                          out.type = "confint", ...){

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##check for distributions
    modFamily <- family(mod)$family
    if(!identical(modFamily, "poisson") && !identical(modFamily, "binomial")) {
        if(c.hat > 1) stop("\ndistribution not appropriate for overdispersion correction\n\n")
    }

    ##for binomial, check that number of trials > 1
    if(identical(modFamily, "binomial")) {
        if(!any(mod$prior.weights > 1)) stop("\nOverdispersion correction only appropriate for success/trials syntax\n\n")
    }
    
    ##extract coefs
    coefs <- coef(mod)
    ##extract SE's
    ses <- sqrt(diag(vcov(mod) * c.hat))

    ##number of estimated parameters
    nparms <- length(coefs)
    
    ##arrange in matrix
    outMat <- matrix(NA, nrow = nparms,
                     ncol = 4)
    outMat[, 1:2] <- cbind(coefs, ses)
    rownames(outMat) <- names(coefs)
    
    ##if interval
    if(identical(out.type, "confint")) {
        
        ##compute confidence intervals
        zstat <- qnorm(p = conf.level)
        outMat[, 3] <- outMat[, 1] - zstat * outMat[, 2]
        outMat[, 4] <- outMat[, 1] + zstat * outMat[, 2]

        ##add names
        colnames(outMat) <- c("estimate", "se", "lowlim", "upplim")
    }

    ##if htest
    if(identical(out.type, "nhst")) {

        ##compute P-value
        outMat[, 3] <- outMat[, 1]/outMat[, 2]
        outMat[, 4] <- 2 * pnorm(abs(outMat[, 3]), lower.tail = FALSE)

        ##add names
        colnames(outMat) <- c("estimate", "se", "z", "pvalue")
    }

    ##assemble in list
    outList <- list(out.type = out.type,
                    c.hat = c.hat,
                    conf.level = conf.level,
                    outMat = outMat)
    
    class(outList) <- c("summaryOD", "list")
    return(outList)
}



##occu
summaryOD.unmarkedFitOccu <- function(mod, c.hat = 1,
                                         conf.level = 0.95,
                                         out.type = "confint", ...){

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##extract coefs
    coefs <- coef(mod)
    ##extract SE's
    ses <- sqrt(diag(vcov(mod) * c.hat))

    ##number of estimated parameters
    nparms <- length(coefs)
    
    ##arrange in matrix
    outMat <- matrix(NA, nrow = nparms,
                     ncol = 4)
    outMat[, 1:2] <- cbind(coefs, ses)
    rownames(outMat) <- names(coefs)
    
    ##if interval
    if(identical(out.type, "confint")) {
        
        ##compute confidence intervals
        zstat <- qnorm(p = conf.level)
        outMat[, 3] <- outMat[, 1] - zstat * outMat[, 2]
        outMat[, 4] <- outMat[, 1] + zstat * outMat[, 2]

        ##add names
        colnames(outMat) <- c("estimate", "se", "lowlim", "upplim")
    }

    ##if htest
    if(identical(out.type, "nhst")) {

        ##compute P-value
        outMat[, 3] <- outMat[, 1]/outMat[, 2]
        outMat[, 4] <- 2 * pnorm(abs(outMat[, 3]), lower.tail = FALSE)

        ##add names
        colnames(outMat) <- c("estimate", "se", "z", "pvalue")
    }

    ##assemble in list
    outList <- list(out.type = out.type,
                    c.hat = c.hat,
                    conf.level = conf.level,
                    outMat = outMat)
    
    class(outList) <- c("summaryOD", "list")
    return(outList)
}



##colext
summaryOD.unmarkedFitColExt <- function(mod, c.hat = 1,
                                        conf.level = 0.95,
                                        out.type = "confint", ...){

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##extract coefs
    coefs <- coef(mod)
    ##extract SE's
    ses <- sqrt(diag(vcov(mod) * c.hat))

    ##number of estimated parameters
    nparms <- length(coefs)
    
    ##arrange in matrix
    outMat <- matrix(NA, nrow = nparms,
                     ncol = 4)
    outMat[, 1:2] <- cbind(coefs, ses)
    rownames(outMat) <- names(coefs)
    
    ##if interval
    if(identical(out.type, "confint")) {
        
        ##compute confidence intervals
        zstat <- qnorm(p = conf.level)
        outMat[, 3] <- outMat[, 1] - zstat * outMat[, 2]
        outMat[, 4] <- outMat[, 1] + zstat * outMat[, 2]

        ##add names
        colnames(outMat) <- c("estimate", "se", "lowlim", "upplim")
    }

    ##if htest
    if(identical(out.type, "nhst")) {

        ##compute P-value
        outMat[, 3] <- outMat[, 1]/outMat[, 2]
        outMat[, 4] <- 2 * pnorm(abs(outMat[, 3]), lower.tail = FALSE)

        ##add names
        colnames(outMat) <- c("estimate", "se", "z", "pvalue")
    }

    ##assemble in list
    outList <- list(out.type = out.type,
                    c.hat = c.hat,
                    conf.level = conf.level,
                    outMat = outMat)
    
    class(outList) <- c("summaryOD", "list")
    return(outList)
}



##occuRN
summaryOD.unmarkedFitOccuRN <- function(mod, c.hat = 1,
                                        conf.level = 0.95,
                                        out.type = "confint", ...){

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##extract coefs
    coefs <- coef(mod)
    ##extract SE's
    ses <- sqrt(diag(vcov(mod) * c.hat))

    ##number of estimated parameters
    nparms <- length(coefs)
    
    ##arrange in matrix
    outMat <- matrix(NA, nrow = nparms,
                     ncol = 4)
    outMat[, 1:2] <- cbind(coefs, ses)
    rownames(outMat) <- names(coefs)
    
    ##if interval
    if(identical(out.type, "confint")) {
        
        ##compute confidence intervals
        zstat <- qnorm(p = conf.level)
        outMat[, 3] <- outMat[, 1] - zstat * outMat[, 2]
        outMat[, 4] <- outMat[, 1] + zstat * outMat[, 2]

        ##add names
        colnames(outMat) <- c("estimate", "se", "lowlim", "upplim")
    }

    ##if htest
    if(identical(out.type, "nhst")) {

        ##compute P-value
        outMat[, 3] <- outMat[, 1]/outMat[, 2]
        outMat[, 4] <- 2 * pnorm(abs(outMat[, 3]), lower.tail = FALSE)

        ##add names
        colnames(outMat) <- c("estimate", "se", "z", "pvalue")
    }

    ##assemble in list
    outList <- list(out.type = out.type,
                    c.hat = c.hat,
                    conf.level = conf.level,
                    outMat = outMat)
    
    class(outList) <- c("summaryOD", "list")
    return(outList)
}



##pcount
summaryOD.unmarkedFitPCount <- function(mod, c.hat = 1,
                                        conf.level = 0.95,
                                        out.type = "confint", ...){

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##check for distributions
    modFamily <- mod@mixture
    if(!identical(modFamily, "P") && !identical(modFamily, "ZIP")) {
        if(c.hat > 1) stop("\ndistribution not appropriate for overdispersion correction\n\n")
    }

    ##extract coefs
    coefs <- coef(mod)
    ##extract SE's
    ses <- sqrt(diag(vcov(mod) * c.hat))

    ##number of estimated parameters
    nparms <- length(coefs)
    
    ##arrange in matrix
    outMat <- matrix(NA, nrow = nparms,
                     ncol = 4)
    outMat[, 1:2] <- cbind(coefs, ses)
    rownames(outMat) <- names(coefs)
    
    ##if interval
    if(identical(out.type, "confint")) {
        
        ##compute confidence intervals
        zstat <- qnorm(p = conf.level)
        outMat[, 3] <- outMat[, 1] - zstat * outMat[, 2]
        outMat[, 4] <- outMat[, 1] + zstat * outMat[, 2]

        ##add names
        colnames(outMat) <- c("estimate", "se", "lowlim", "upplim")
    }

    ##if htest
    if(identical(out.type, "nhst")) {

        ##compute P-value
        outMat[, 3] <- outMat[, 1]/outMat[, 2]
        outMat[, 4] <- 2 * pnorm(abs(outMat[, 3]), lower.tail = FALSE)

        ##add names
        colnames(outMat) <- c("estimate", "se", "z", "pvalue")
    }

    ##assemble in list
    outList <- list(out.type = out.type,
                    c.hat = c.hat,
                    conf.level = conf.level,
                    outMat = outMat)
    
    class(outList) <- c("summaryOD", "list")
    return(outList)
}



##pcountOpen
summaryOD.unmarkedFitPCO <- function(mod, c.hat = 1,
                                     conf.level = 0.95,
                                     out.type = "confint", ...){

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##check for distributions
    modFamily <- mod@mixture
    if(!identical(modFamily, "P") && !identical(modFamily, "ZIP")) {
        if(c.hat > 1) stop("\ndistribution not appropriate for overdispersion correction\n\n")
    }

    ##extract coefs
    coefs <- coef(mod)
    ##extract SE's
    ses <- sqrt(diag(vcov(mod) * c.hat))

    ##number of estimated parameters
    nparms <- length(coefs)
    
    ##arrange in matrix
    outMat <- matrix(NA, nrow = nparms,
                     ncol = 4)
    outMat[, 1:2] <- cbind(coefs, ses)
    rownames(outMat) <- names(coefs)
    
    ##if interval
    if(identical(out.type, "confint")) {
        
        ##compute confidence intervals
        zstat <- qnorm(p = conf.level)
        outMat[, 3] <- outMat[, 1] - zstat * outMat[, 2]
        outMat[, 4] <- outMat[, 1] + zstat * outMat[, 2]

        ##add names
        colnames(outMat) <- c("estimate", "se", "lowlim", "upplim")
    }

    ##if htest
    if(identical(out.type, "nhst")) {

        ##compute P-value
        outMat[, 3] <- outMat[, 1]/outMat[, 2]
        outMat[, 4] <- 2 * pnorm(abs(outMat[, 3]), lower.tail = FALSE)

        ##add names
        colnames(outMat) <- c("estimate", "se", "z", "pvalue")
    }

    ##assemble in list
    outList <- list(out.type = out.type,
                    c.hat = c.hat,
                    conf.level = conf.level,
                    outMat = outMat)
    
    class(outList) <- c("summaryOD", "list")
    return(outList)
}



##distsamp
summaryOD.unmarkedFitDS <- function(mod, c.hat = 1,
                                    conf.level = 0.95,
                                    out.type = "confint", ...){

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##extract coefs
    coefs <- coef(mod)
    ##extract SE's
    ses <- sqrt(diag(vcov(mod) * c.hat))

    ##number of estimated parameters
    nparms <- length(coefs)
    
    ##arrange in matrix
    outMat <- matrix(NA, nrow = nparms,
                     ncol = 4)
    outMat[, 1:2] <- cbind(coefs, ses)
    rownames(outMat) <- names(coefs)
    
    ##if interval
    if(identical(out.type, "confint")) {
        
        ##compute confidence intervals
        zstat <- qnorm(p = conf.level)
        outMat[, 3] <- outMat[, 1] - zstat * outMat[, 2]
        outMat[, 4] <- outMat[, 1] + zstat * outMat[, 2]

        ##add names
        colnames(outMat) <- c("estimate", "se", "lowlim", "upplim")
    }

    ##if htest
    if(identical(out.type, "nhst")) {

        ##compute P-value
        outMat[, 3] <- outMat[, 1]/outMat[, 2]
        outMat[, 4] <- 2 * pnorm(abs(outMat[, 3]), lower.tail = FALSE)

        ##add names
        colnames(outMat) <- c("estimate", "se", "z", "pvalue")
    }

    ##assemble in list
    outList <- list(out.type = out.type,
                    c.hat = c.hat,
                    conf.level = conf.level,
                    outMat = outMat)
    
    class(outList) <- c("summaryOD", "list")
    return(outList)
}



##gdistsamp
summaryOD.unmarkedFitGDS <- function(mod, c.hat = 1,
                                     conf.level = 0.95,
                                     out.type = "confint", ...){

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##check for distributions
    modFamily <- mod@mixture
    if(!identical(modFamily, "P")) {
        if(c.hat > 1) stop("\ndistribution not appropriate for overdispersion correction\n\n")
    }

    ##extract coefs
    coefs <- coef(mod)
    ##extract SE's
    ses <- sqrt(diag(vcov(mod) * c.hat))

    ##number of estimated parameters
    nparms <- length(coefs)
    
    ##arrange in matrix
    outMat <- matrix(NA, nrow = nparms,
                     ncol = 4)
    outMat[, 1:2] <- cbind(coefs, ses)
    rownames(outMat) <- names(coefs)
    
    ##if interval
    if(identical(out.type, "confint")) {
        
        ##compute confidence intervals
        zstat <- qnorm(p = conf.level)
        outMat[, 3] <- outMat[, 1] - zstat * outMat[, 2]
        outMat[, 4] <- outMat[, 1] + zstat * outMat[, 2]

        ##add names
        colnames(outMat) <- c("estimate", "se", "lowlim", "upplim")
    }

    ##if htest
    if(identical(out.type, "nhst")) {

        ##compute P-value
        outMat[, 3] <- outMat[, 1]/outMat[, 2]
        outMat[, 4] <- 2 * pnorm(abs(outMat[, 3]), lower.tail = FALSE)

        ##add names
        colnames(outMat) <- c("estimate", "se", "z", "pvalue")
    }

    ##assemble in list
    outList <- list(out.type = out.type,
                    c.hat = c.hat,
                    conf.level = conf.level,
                    outMat = outMat)
    
    class(outList) <- c("summaryOD", "list")
    return(outList)
}



##occuFP
summaryOD.unmarkedFitOccuFP <- function(mod, c.hat = 1,
                                        conf.level = 0.95,
                                        out.type = "confint", ...){

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##extract coefs
    coefs <- coef(mod)
    ##extract SE's
    ses <- sqrt(diag(vcov(mod) * c.hat))

    ##number of estimated parameters
    nparms <- length(coefs)
    
    ##arrange in matrix
    outMat <- matrix(NA, nrow = nparms,
                     ncol = 4)
    outMat[, 1:2] <- cbind(coefs, ses)
    rownames(outMat) <- names(coefs)
    
    ##if interval
    if(identical(out.type, "confint")) {
        
        ##compute confidence intervals
        zstat <- qnorm(p = conf.level)
        outMat[, 3] <- outMat[, 1] - zstat * outMat[, 2]
        outMat[, 4] <- outMat[, 1] + zstat * outMat[, 2]

        ##add names
        colnames(outMat) <- c("estimate", "se", "lowlim", "upplim")
    }

    ##if htest
    if(identical(out.type, "nhst")) {

        ##compute P-value
        outMat[, 3] <- outMat[, 1]/outMat[, 2]
        outMat[, 4] <- 2 * pnorm(abs(outMat[, 3]), lower.tail = FALSE)

        ##add names
        colnames(outMat) <- c("estimate", "se", "z", "pvalue")
    }

    ##assemble in list
    outList <- list(out.type = out.type,
                    c.hat = c.hat,
                    conf.level = conf.level,
                    outMat = outMat)
    
    class(outList) <- c("summaryOD", "list")
    return(outList)
}



##multinomPois
summaryOD.unmarkedFitMPois <- function(mod, c.hat = 1,
                                       conf.level = 0.95,
                                       out.type = "confint", ...){

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##extract coefs
    coefs <- coef(mod)
    ##extract SE's
    ses <- sqrt(diag(vcov(mod) * c.hat))

    ##number of estimated parameters
    nparms <- length(coefs)
    
    ##arrange in matrix
    outMat <- matrix(NA, nrow = nparms,
                     ncol = 4)
    outMat[, 1:2] <- cbind(coefs, ses)
    rownames(outMat) <- names(coefs)
    
    ##if interval
    if(identical(out.type, "confint")) {
        
        ##compute confidence intervals
        zstat <- qnorm(p = conf.level)
        outMat[, 3] <- outMat[, 1] - zstat * outMat[, 2]
        outMat[, 4] <- outMat[, 1] + zstat * outMat[, 2]

        ##add names
        colnames(outMat) <- c("estimate", "se", "lowlim", "upplim")
    }

    ##if htest
    if(identical(out.type, "nhst")) {

        ##compute P-value
        outMat[, 3] <- outMat[, 1]/outMat[, 2]
        outMat[, 4] <- 2 * pnorm(abs(outMat[, 3]), lower.tail = FALSE)

        ##add names
        colnames(outMat) <- c("estimate", "se", "z", "pvalue")
    }

    ##assemble in list
    outList <- list(out.type = out.type,
                    c.hat = c.hat,
                    conf.level = conf.level,
                    outMat = outMat)
    
    class(outList) <- c("summaryOD", "list")
    return(outList)
}



##gmultmix
summaryOD.unmarkedFitGMM <- function(mod, c.hat = 1,
                                     conf.level = 0.95,
                                     out.type = "confint", ...){

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##check for distributions
    modFamily <- mod@mixture
    if(!identical(modFamily, "P")) {
        if(c.hat > 1) stop("\ndistribution not appropriate for overdispersion correction\n\n")
    }

    ##extract coefs
    coefs <- coef(mod)
    ##extract SE's
    ses <- sqrt(diag(vcov(mod) * c.hat))

    ##number of estimated parameters
    nparms <- length(coefs)
    
    ##arrange in matrix
    outMat <- matrix(NA, nrow = nparms,
                     ncol = 4)
    outMat[, 1:2] <- cbind(coefs, ses)
    rownames(outMat) <- names(coefs)
    
    ##if interval
    if(identical(out.type, "confint")) {
        
        ##compute confidence intervals
        zstat <- qnorm(p = conf.level)
        outMat[, 3] <- outMat[, 1] - zstat * outMat[, 2]
        outMat[, 4] <- outMat[, 1] + zstat * outMat[, 2]

        ##add names
        colnames(outMat) <- c("estimate", "se", "lowlim", "upplim")
    }

    ##if htest
    if(identical(out.type, "nhst")) {

        ##compute P-value
        outMat[, 3] <- outMat[, 1]/outMat[, 2]
        outMat[, 4] <- 2 * pnorm(abs(outMat[, 3]), lower.tail = FALSE)

        ##add names
        colnames(outMat) <- c("estimate", "se", "z", "pvalue")
    }

    ##assemble in list
    outList <- list(out.type = out.type,
                    c.hat = c.hat,
                    conf.level = conf.level,
                    outMat = outMat)
    
    class(outList) <- c("summaryOD", "list")
    return(outList)
}



##gpcount
summaryOD.unmarkedFitGPC <- function(mod, c.hat = 1,
                                     conf.level = 0.95,
                                     out.type = "confint", ...){

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##check for distributions
    modFamily <- mod@mixture
    if(!identical(modFamily, "P")) {
        if(c.hat > 1) stop("\ndistribution not appropriate for overdispersion correction\n\n")
    }

    ##extract coefs
    coefs <- coef(mod)
    ##extract SE's
    ses <- sqrt(diag(vcov(mod) * c.hat))

    ##number of estimated parameters
    nparms <- length(coefs)
    
    ##arrange in matrix
    outMat <- matrix(NA, nrow = nparms,
                     ncol = 4)
    outMat[, 1:2] <- cbind(coefs, ses)
    rownames(outMat) <- names(coefs)
    
    ##if interval
    if(identical(out.type, "confint")) {
        
        ##compute confidence intervals
        zstat <- qnorm(p = conf.level)
        outMat[, 3] <- outMat[, 1] - zstat * outMat[, 2]
        outMat[, 4] <- outMat[, 1] + zstat * outMat[, 2]

        ##add names
        colnames(outMat) <- c("estimate", "se", "lowlim", "upplim")
    }

    ##if htest
    if(identical(out.type, "nhst")) {

        ##compute P-value
        outMat[, 3] <- outMat[, 1]/outMat[, 2]
        outMat[, 4] <- 2 * pnorm(abs(outMat[, 3]), lower.tail = FALSE)

        ##add names
        colnames(outMat) <- c("estimate", "se", "z", "pvalue")
    }

    ##assemble in list
    outList <- list(out.type = out.type,
                    c.hat = c.hat,
                    conf.level = conf.level,
                    outMat = outMat)
    
    class(outList) <- c("summaryOD", "list")
    return(outList)
}



##glmer
summaryOD.glmerMod <- function(mod, c.hat = 1,
                               conf.level = 0.95,
                               out.type = "confint", ...){

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##check for distributions
    modFamily <- fam.link.mer(mod)$family
    if(!identical(modFamily, "poisson") && !identical(modFamily, "binomial")) {
        if(c.hat > 1) stop("\ndistribution not appropriate for overdispersion correction\n\n")
    }

    ##for binomial, check that number of trials > 1
    if(identical(modFamily, "binomial")) {
        if(!any(mod@resp$weights > 1)) stop("\nOverdispersion correction only appropriate for success/trials syntax\n\n")
    }
    
    ##extract coefs
    coefs <- fixef(mod)
    ##extract SE's
    ses <- sqrt(diag(vcov(mod) * c.hat))

    ##number of estimated parameters
    nparms <- length(coefs)
    
    ##arrange in matrix
    outMat <- matrix(NA, nrow = nparms,
                     ncol = 4)
    outMat[, 1:2] <- cbind(coefs, ses)
    rownames(outMat) <- names(coefs)
    
    ##if interval
    if(identical(out.type, "confint")) {
        
        ##compute confidence intervals
        zstat <- qnorm(p = conf.level)
        outMat[, 3] <- outMat[, 1] - zstat * outMat[, 2]
        outMat[, 4] <- outMat[, 1] + zstat * outMat[, 2]

        ##add names
        colnames(outMat) <- c("estimate", "se", "lowlim", "upplim")
    }

    ##if htest
    if(identical(out.type, "nhst")) {

        ##compute P-value
        outMat[, 3] <- outMat[, 1]/outMat[, 2]
        outMat[, 4] <- 2 * pnorm(abs(outMat[, 3]), lower.tail = FALSE)

        ##add names
        colnames(outMat) <- c("estimate", "se", "z", "pvalue")
    }

    ##assemble in list
    outList <- list(out.type = out.type,
                    c.hat = c.hat,
                    conf.level = conf.level,
                    outMat = outMat)
    
    class(outList) <- c("summaryOD", "list")
    return(outList)
}



##occuMulti
summaryOD.unmarkedFitOccuMulti <- function(mod, c.hat = 1,
                                           conf.level = 0.95,
                                           out.type = "confint", ...){

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##extract coefs
    coefs <- coef(mod)
    ##extract SE's
    ses <- sqrt(diag(vcov(mod) * c.hat))

    ##number of estimated parameters
    nparms <- length(coefs)
    
    ##arrange in matrix
    outMat <- matrix(NA, nrow = nparms,
                     ncol = 4)
    outMat[, 1:2] <- cbind(coefs, ses)
    rownames(outMat) <- names(coefs)
    
    ##if interval
    if(identical(out.type, "confint")) {
        
        ##compute confidence intervals
        zstat <- qnorm(p = conf.level)
        outMat[, 3] <- outMat[, 1] - zstat * outMat[, 2]
        outMat[, 4] <- outMat[, 1] + zstat * outMat[, 2]

        ##add names
        colnames(outMat) <- c("estimate", "se", "lowlim", "upplim")
    }

    ##if htest
    if(identical(out.type, "nhst")) {

        ##compute P-value
        outMat[, 3] <- outMat[, 1]/outMat[, 2]
        outMat[, 4] <- 2 * pnorm(abs(outMat[, 3]), lower.tail = FALSE)

        ##add names
        colnames(outMat) <- c("estimate", "se", "z", "pvalue")
    }

    ##assemble in list
    outList <- list(out.type = out.type,
                    c.hat = c.hat,
                    conf.level = conf.level,
                    outMat = outMat)
    
    class(outList) <- c("summaryOD", "list")
    return(outList)
}



##occuMS
summaryOD.unmarkedFitOccuMS <- function(mod, c.hat = 1,
                                        conf.level = 0.95,
                                        out.type = "confint", ...){

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##extract coefs
    coefs <- coef(mod)
    ##extract SE's
    ses <- sqrt(diag(vcov(mod) * c.hat))

    ##number of estimated parameters
    nparms <- length(coefs)
    
    ##arrange in matrix
    outMat <- matrix(NA, nrow = nparms,
                     ncol = 4)
    outMat[, 1:2] <- cbind(coefs, ses)
    rownames(outMat) <- names(coefs)
    
    ##if interval
    if(identical(out.type, "confint")) {
        
        ##compute confidence intervals
        zstat <- qnorm(p = conf.level)
        outMat[, 3] <- outMat[, 1] - zstat * outMat[, 2]
        outMat[, 4] <- outMat[, 1] + zstat * outMat[, 2]

        ##add names
        colnames(outMat) <- c("estimate", "se", "lowlim", "upplim")
    }

    ##if htest
    if(identical(out.type, "nhst")) {

        ##compute P-value
        outMat[, 3] <- outMat[, 1]/outMat[, 2]
        outMat[, 4] <- 2 * pnorm(abs(outMat[, 3]), lower.tail = FALSE)

        ##add names
        colnames(outMat) <- c("estimate", "se", "z", "pvalue")
    }

    ##assemble in list
    outList <- list(out.type = out.type,
                    c.hat = c.hat,
                    conf.level = conf.level,
                    outMat = outMat)
    
    class(outList) <- c("summaryOD", "list")
    return(outList)
}



##occuTTD
summaryOD.unmarkedFitOccuTTD <- function(mod, c.hat = 1,
                                         conf.level = 0.95,
                                         out.type = "confint", ...){

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##extract coefs
    coefs <- coef(mod)
    ##extract SE's
    ses <- sqrt(diag(vcov(mod) * c.hat))

    ##number of estimated parameters
    nparms <- length(coefs)
    
    ##arrange in matrix
    outMat <- matrix(NA, nrow = nparms,
                     ncol = 4)
    outMat[, 1:2] <- cbind(coefs, ses)
    rownames(outMat) <- names(coefs)
    
    ##if interval
    if(identical(out.type, "confint")) {
        
        ##compute confidence intervals
        zstat <- qnorm(p = conf.level)
        outMat[, 3] <- outMat[, 1] - zstat * outMat[, 2]
        outMat[, 4] <- outMat[, 1] + zstat * outMat[, 2]

        ##add names
        colnames(outMat) <- c("estimate", "se", "lowlim", "upplim")
    }

    ##if htest
    if(identical(out.type, "nhst")) {

        ##compute P-value
        outMat[, 3] <- outMat[, 1]/outMat[, 2]
        outMat[, 4] <- 2 * pnorm(abs(outMat[, 3]), lower.tail = FALSE)

        ##add names
        colnames(outMat) <- c("estimate", "se", "z", "pvalue")
    }

    ##assemble in list
    outList <- list(out.type = out.type,
                    c.hat = c.hat,
                    conf.level = conf.level,
                    outMat = outMat)
    
    class(outList) <- c("summaryOD", "list")
    return(outList)
}



##multmixOpen
summaryOD.unmarkedFitMMO <- function(mod, c.hat = 1,
                                     conf.level = 0.95,
                                     out.type = "confint", ...){

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##extract coefs
    coefs <- coef(mod)
    ##extract SE's
    ses <- sqrt(diag(vcov(mod) * c.hat))

    ##number of estimated parameters
    nparms <- length(coefs)
    
    ##arrange in matrix
    outMat <- matrix(NA, nrow = nparms,
                     ncol = 4)
    outMat[, 1:2] <- cbind(coefs, ses)
    rownames(outMat) <- names(coefs)
    
    ##if interval
    if(identical(out.type, "confint")) {
        
        ##compute confidence intervals
        zstat <- qnorm(p = conf.level)
        outMat[, 3] <- outMat[, 1] - zstat * outMat[, 2]
        outMat[, 4] <- outMat[, 1] + zstat * outMat[, 2]

        ##add names
        colnames(outMat) <- c("estimate", "se", "lowlim", "upplim")
    }

    ##if htest
    if(identical(out.type, "nhst")) {

        ##compute P-value
        outMat[, 3] <- outMat[, 1]/outMat[, 2]
        outMat[, 4] <- 2 * pnorm(abs(outMat[, 3]), lower.tail = FALSE)

        ##add names
        colnames(outMat) <- c("estimate", "se", "z", "pvalue")
    }

    ##assemble in list
    outList <- list(out.type = out.type,
                    c.hat = c.hat,
                    conf.level = conf.level,
                    outMat = outMat)
    
    class(outList) <- c("summaryOD", "list")
    return(outList)
}



##distsampOpen
summaryOD.unmarkedFitDSO <- function(mod, c.hat = 1,
                                     conf.level = 0.95,
                                     out.type = "confint", ...){

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##extract coefs
    coefs <- coef(mod)
    ##extract SE's
    ses <- sqrt(diag(vcov(mod) * c.hat))

    ##number of estimated parameters
    nparms <- length(coefs)
    
    ##arrange in matrix
    outMat <- matrix(NA, nrow = nparms,
                     ncol = 4)
    outMat[, 1:2] <- cbind(coefs, ses)
    rownames(outMat) <- names(coefs)
    
    ##if interval
    if(identical(out.type, "confint")) {
        
        ##compute confidence intervals
        zstat <- qnorm(p = conf.level)
        outMat[, 3] <- outMat[, 1] - zstat * outMat[, 2]
        outMat[, 4] <- outMat[, 1] + zstat * outMat[, 2]

        ##add names
        colnames(outMat) <- c("estimate", "se", "lowlim", "upplim")
    }

    ##if htest
    if(identical(out.type, "nhst")) {

        ##compute P-value
        outMat[, 3] <- outMat[, 1]/outMat[, 2]
        outMat[, 4] <- 2 * pnorm(abs(outMat[, 3]), lower.tail = FALSE)

        ##add names
        colnames(outMat) <- c("estimate", "se", "z", "pvalue")
    }

    ##assemble in list
    outList <- list(out.type = out.type,
                    c.hat = c.hat,
                    conf.level = conf.level,
                    outMat = outMat)
    
    class(outList) <- c("summaryOD", "list")
    return(outList)
}



##goccu
summaryOD.unmarkedFitGOccu <- function(mod, c.hat = 1,
                                     conf.level = 0.95,
                                     out.type = "confint", ...){

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }


    ##extract coefs
    coefs <- coef(mod)
    ##extract SE's
    ses <- sqrt(diag(vcov(mod) * c.hat))

    ##number of estimated parameters
    nparms <- length(coefs)
    
    ##arrange in matrix
    outMat <- matrix(NA, nrow = nparms,
                     ncol = 4)
    outMat[, 1:2] <- cbind(coefs, ses)
    rownames(outMat) <- names(coefs)
    
    ##if interval
    if(identical(out.type, "confint")) {
        
        ##compute confidence intervals
        zstat <- qnorm(p = conf.level)
        outMat[, 3] <- outMat[, 1] - zstat * outMat[, 2]
        outMat[, 4] <- outMat[, 1] + zstat * outMat[, 2]

        ##add names
        colnames(outMat) <- c("estimate", "se", "lowlim", "upplim")
    }

    ##if htest
    if(identical(out.type, "nhst")) {

        ##compute P-value
        outMat[, 3] <- outMat[, 1]/outMat[, 2]
        outMat[, 4] <- 2 * pnorm(abs(outMat[, 3]), lower.tail = FALSE)

        ##add names
        colnames(outMat) <- c("estimate", "se", "z", "pvalue")
    }

    ##assemble in list
    outList <- list(out.type = out.type,
                    c.hat = c.hat,
                    conf.level = conf.level,
                    outMat = outMat)
    
    class(outList) <- c("summaryOD", "list")
    return(outList)
}



##occuComm
summaryOD.unmarkedFitOccuComm <- function(mod, c.hat = 1,
                                          conf.level = 0.95,
                                          out.type = "confint", ...){

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##extract coefs
    coefs <- coef(mod)
    ##extract SE's
    ses <- sqrt(diag(vcov(mod) * c.hat))

    ##number of estimated parameters
    nparms <- length(coefs)
    
    ##arrange in matrix
    outMat <- matrix(NA, nrow = nparms,
                     ncol = 4)
    outMat[, 1:2] <- cbind(coefs, ses)
    rownames(outMat) <- names(coefs)
    
    ##if interval
    if(identical(out.type, "confint")) {
        
        ##compute confidence intervals
        zstat <- qnorm(p = conf.level)
        outMat[, 3] <- outMat[, 1] - zstat * outMat[, 2]
        outMat[, 4] <- outMat[, 1] + zstat * outMat[, 2]

        ##add names
        colnames(outMat) <- c("estimate", "se", "lowlim", "upplim")
    }

    ##if htest
    if(identical(out.type, "nhst")) {

        ##compute P-value
        outMat[, 3] <- outMat[, 1]/outMat[, 2]
        outMat[, 4] <- 2 * pnorm(abs(outMat[, 3]), lower.tail = FALSE)

        ##add names
        colnames(outMat) <- c("estimate", "se", "z", "pvalue")
    }

    ##assemble in list
    outList <- list(out.type = out.type,
                    c.hat = c.hat,
                    conf.level = conf.level,
                    outMat = outMat)
    
    class(outList) <- c("summaryOD", "list")
    return(outList)
}



##glmmTMB
summaryOD.glmmTMB <- function(mod, c.hat = 1,
                              conf.level = 0.95,
                              out.type = "confint", ...){

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##check for distributions
    ##determine family of model
    fam <- family(mod)$family

    ##extract response
    
    ##if binomial, check if n > 1 for each case
    if(fam == "binomial") {
        resp <- mod$frame[, mod$modelInfo$respCol]
        if(!is.matrix(resp)) {
            if(!any(names(mod$frame) == "(weights)")) {
                stop("\nOverdispersion correction only appropriate for success/trials syntax\n\n")
            }
        }
    }

    ##Poisson or binomial
    if(!any(fam == c("poisson", "binomial"))) {
        stop("\ndistribution not appropriate for overdispersion correction\n")
    }


    ##extract coefs
    coefs <- fixef(mod)$cond
    ##extract SE's
    ses <- sqrt(diag(vcov(mod)$cond * c.hat))

    ##number of estimated parameters
    nparms <- length(coefs)
    
    ##arrange in matrix
    outMat <- matrix(NA, nrow = nparms,
                     ncol = 4)
    outMat[, 1:2] <- cbind(coefs, ses)
    rownames(outMat) <- names(coefs)
    
    ##if interval
    if(identical(out.type, "confint")) {
        
        ##compute confidence intervals
        zstat <- qnorm(p = conf.level)
        outMat[, 3] <- outMat[, 1] - zstat * outMat[, 2]
        outMat[, 4] <- outMat[, 1] + zstat * outMat[, 2]

        ##add names
        colnames(outMat) <- c("estimate", "se", "lowlim", "upplim")
    }

    ##if nhst
    if(identical(out.type, "nhst")) {

        ##compute P-value
        outMat[, 3] <- outMat[, 1]/outMat[, 2]
        outMat[, 4] <- 2 * pnorm(abs(outMat[, 3]), lower.tail = FALSE)

        ##add names
        colnames(outMat) <- c("estimate", "se", "z", "pvalue")
    }

    ##assemble in list
    outList <- list(out.type = out.type,
                    c.hat = c.hat,
                    conf.level = conf.level,
                    outMat = outMat)
    
    class(outList) <- c("summaryOD", "list")
    return(outList)
}



##maxlike
summaryOD.maxlikeFit <- function(mod, c.hat = 1,
                                 conf.level = 0.95,
                                 out.type = "confint", ...){

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##extract coefs
    coefs <- coef(mod)
    ##extract SE's
    ses <- sqrt(diag(vcov(mod) * c.hat))

    ##number of estimated parameters
    nparms <- length(coefs)
    
    ##arrange in matrix
    outMat <- matrix(NA, nrow = nparms,
                     ncol = 4)
    outMat[, 1:2] <- cbind(coefs, ses)
    rownames(outMat) <- names(coefs)
    
    ##if interval
    if(identical(out.type, "confint")) {
        
        ##compute confidence intervals
        zstat <- qnorm(p = conf.level)
        outMat[, 3] <- outMat[, 1] - zstat * outMat[, 2]
        outMat[, 4] <- outMat[, 1] + zstat * outMat[, 2]

        ##add names
        colnames(outMat) <- c("estimate", "se", "lowlim", "upplim")
    }

    ##if htest
    if(identical(out.type, "nhst")) {

        ##compute P-value
        outMat[, 3] <- outMat[, 1]/outMat[, 2]
        outMat[, 4] <- 2 * pnorm(abs(outMat[, 3]), lower.tail = FALSE)

        ##add names
        colnames(outMat) <- c("estimate", "se", "z", "pvalue")
    }

    ##assemble in list
    outList <- list(out.type = out.type,
                    c.hat = c.hat,
                    conf.level = conf.level,
                    outMat = outMat)
    
    class(outList) <- c("summaryOD", "list")
    return(outList)
}



##multinom
summaryOD.multinom <- function(mod, c.hat = 1,
                               conf.level = 0.95,
                               out.type = "confint", ...){

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    ##extract coefs
    coefs <- as.vector(coef(mod))
    ##extract SE's
    ses <- sqrt(diag(vcov(mod) * c.hat))
    ##extract names
    parm.names <- names(ses)
    
    ##number of estimated parameters
    nparms <- length(ses)

    ##arrange in matrix
    outMat <- matrix(NA, nrow = nparms,
                     ncol = 4)
    outMat[, 1:2] <- cbind(coefs, ses)
    rownames(outMat) <- parm.names
      
    ##if interval
    if(identical(out.type, "confint")) {
        
        ##compute confidence intervals
        zstat <- qnorm(p = conf.level)
        outMat[, 3] <- outMat[, 1] - zstat * outMat[, 2]
        outMat[, 4] <- outMat[, 1] + zstat * outMat[, 2]

        ##add names
        colnames(outMat) <- c("estimate", "se", "lowlim", "upplim")
    }

    ##if htest
    if(identical(out.type, "nhst")) {

        ##compute P-value
        outMat[, 3] <- outMat[, 1]/outMat[, 2]
        outMat[, 4] <- 2 * pnorm(abs(outMat[, 3]), lower.tail = FALSE)

        ##add names
        colnames(outMat) <- c("estimate", "se", "z", "pvalue")
    }

    ##assemble in list
    outList <- list(out.type = out.type,
                    c.hat = c.hat,
                    conf.level = conf.level,
                    outMat = outMat)
    
    class(outList) <- c("summaryOD", "list")
    return(outList)
}



##vglm
summaryOD.vglm <- function(mod, c.hat = 1,
                           conf.level = 0.95,
                           out.type = "confint", ...){

    if(c.hat > 4) warning("\nHigh overdispersion:  model fit is questionable\n")
    if(c.hat < 1) {
        warning("\nUnderdispersion: c-hat is fixed to 1\n")
        c.hat <- 1
    }

    
    ##check family of vglm to avoid problems
    fam.type <- mod@family@vfamily[1]
    if(!(fam.type == "poissonff" || fam.type == "binomialff" || fam.type == "multinomial")) stop("\nDistribution not supported by function\n")
    if(fam.type == "binomialff") {
        if(!any(mod@prior.weights > 1)) stop("\nOverdispersion correction only appropriate for success/trials syntax\n\n")
    }
        
    
    ##extract coefs
    coefs <- coef(mod)
    ##extract SE's
    ses <- sqrt(diag(vcov(mod) * c.hat))

    ##number of estimated parameters
    nparms <- length(coefs)
    
    ##arrange in matrix
    outMat <- matrix(NA, nrow = nparms,
                     ncol = 4)
    outMat[, 1:2] <- cbind(coefs, ses)
    rownames(outMat) <- names(coefs)
    
    ##if interval
    if(identical(out.type, "confint")) {
        
        ##compute confidence intervals
        zstat <- qnorm(p = conf.level)
        outMat[, 3] <- outMat[, 1] - zstat * outMat[, 2]
        outMat[, 4] <- outMat[, 1] + zstat * outMat[, 2]

        ##add names
        colnames(outMat) <- c("estimate", "se", "lowlim", "upplim")
    }

    ##if htest
    if(identical(out.type, "nhst")) {

        ##compute P-value
        outMat[, 3] <- outMat[, 1]/outMat[, 2]
        outMat[, 4] <- 2 * pnorm(abs(outMat[, 3]), lower.tail = FALSE)

        ##add names
        colnames(outMat) <- c("estimate", "se", "z", "pvalue")
    }

    ##assemble in list
    outList <- list(out.type = out.type,
                    c.hat = c.hat,
                    conf.level = conf.level,
                    outMat = outMat)
    
    class(outList) <- c("summaryOD", "list")
    return(outList)
}



print.summaryOD <- function(x, digits = 4, ...) {

    ##extract information
    conf.level <- x$conf.level
    c.hat <- x$c.hat
    out.type <- x$out.type
    outMat <- x$outMat

    if(identical(out.type, "confint")) {
        ##label for confidence limit
        lowLab <- paste("Lower ", conf.level * 100, "%", " CL", sep = "")
        uppLab <- paste("Upper ", conf.level * 100, "%", " CL", sep = "")
        
        ##add names
        colnames(outMat) <- c("Estimate", "Std. Error", lowLab, uppLab)

        if(c.hat <= 1) {
            cat("\nPrecision unadjusted for overdispersion:\n\n")
        } else {
            cat("\nPrecision adjusted for overdispersion:\n\n")
        }
        printCoefmat(outMat, digits = digits)
                    cat("\n(c-hat = ", c.hat, ")", "\n", sep = "")
        cat("\n")
    }


    
    if(identical(out.type, "nhst")) {

        ##label for P-values
        zLab <- "z value"
        pLab <- "Pr(>|z|)"
        
        ##add names
        colnames(outMat) <- c("Estimate", "Std. Error", zLab, pLab)

        if(c.hat <= 1) {
            cat("\nPrecision and hypothesis tests unadjusted for overdispersion:\n\n")
        } else {
            cat("\nPrecision and hypothesis tests adjusted for overdispersion:\n\n")
        }

        printCoefmat(outMat, digits = digits, signif.stars = FALSE)
        cat("\n(c-hat = ", c.hat, ")", "\n", sep = "")
        cat("\n")
    }
}

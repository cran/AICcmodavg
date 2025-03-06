##custom functions for user-supplied model input


##Custom AICc computation where user inputs logL, K, and nobs manually
##convenient when output imported from other software
AICcCustom <- function(logL, K, return.K = FALSE, second.ord = TRUE,
                       nobs = NULL, c.hat = 1) {
  if(is.null(nobs) && identical(second.ord, TRUE)) {
    stop("\nYou must supply a value for 'nobs' for the second-order AIC\n")
  } else {n <- nobs}
  
    LL <- logL
      
    if(c.hat == 1) {
      if(second.ord == TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))} else{AICc <- -2*LL+2*K}
    }
    if(c.hat > 1 && c.hat <= 4) {
      K <- K+1
      if(second.ord==TRUE) {
        AICc <- (-2*LL/c.hat)+2*K*(n/(n-K-1))
        ##adjust parameter count to include estimation of dispersion parameter
      } else{
        AICc <- (-2*LL/c.hat)+2*K}
#      cat("\nc-hat estimate was added to parameter count\n")
    }

    if(c.hat > 4) stop("High overdispersion and model fit is questionable\n")
    if(c.hat < 1) stop("You should set \'c.hat\' to 1 if < 1, but values << 1 might also indicate lack of fit\n")

  if(return.K == TRUE) AICc <- K
  AICc
}



##Custom model selection where user inputs logL, K, and nobs manually
##convenient when output imported from other software
aictabCustom <-
  function(logL, K, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1){
    
    ##check if modnames are not supplied
    if(is.null(modnames)) {
      modnames <- paste("Mod", 1:length(logL), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    }

    ##check that nobs is the same for all models
    if(length(unique(nobs)) > 1) stop("\nSample size must be identical for all models\n")
    
    ##check that logL, K, estimate, se are vectors of same length
    nlogL <- length(logL)
    nK <- length(K)

    if(!all(nlogL == c(nlogL, nK))) stop("\nArguments 'logL' and 'K' must be of equal length\n")

    ##create model selection table
    Results <- NULL
    Results <- data.frame(Modnames = modnames)                    #assign model names to first column
    Results$K <- AICcCustom(logL = logL, K = K, return.K = TRUE,
                            second.ord = second.ord, nobs = nobs,
                            c.hat = c.hat)  #extract number of parameters
    Results$AICc <- AICcCustom(logL = logL, K = K, second.ord = second.ord,
                               nobs = nobs, c.hat = c.hat)  #extract AICc                                      
    Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
    Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
    Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

    ##check if some models are redundant
    if(length(unique(Results$AICc)) != length(K)) warning("\nCheck model structure carefully as some models may be redundant\n")
         
    ##check if AICc and c.hat = 1
    if(second.ord == TRUE && c.hat == 1) {
      Results$LL <- logL
    }
    
    ##rename correctly to QAICc and add column for c-hat
    if(second.ord == TRUE && c.hat > 1) {
      colnames(Results) <- c("Modnames", "K", "QAICc", "Delta_QAICc", "ModelLik", "QAICcWt")
      LL <- logL
      Results$Quasi.LL <- LL/c.hat
      Results$c_hat <- c.hat
    }      

    ##rename correctly to AIC
    if(second.ord == FALSE && c.hat == 1) {
      colnames(Results)<-c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
      Results$LL <- logL
    }  

    ##rename correctly to QAIC and add column for c-hat
    if(second.ord == FALSE && c.hat > 1) {
      colnames(Results)<-c("Modnames", "K", "QAIC", "Delta_QAIC", "ModelLik", "QAICWt")
      LL <- logL
      Results$Quasi.LL <- LL/c.hat
      Results$c_hat<-c.hat
    }

    
    if(sort)  {
      Results <- Results[order(Results[, 4]),] 	  #if sort=TRUE, models are ranked based on Akaike weights
      Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
    } else {Results$Cum.Wt <- NULL}

    
    class(Results) <- c("aictab", "data.frame")
    return(Results)
  }



###################################################
###################################################
##BIC-related functions
##BIC
useBICCustom <-
    function(logL, K, return.K = FALSE, nobs = NULL, c.hat = 1){
    
        if(is.null(nobs)) {
            stop("\nYou must supply a value for 'nobs' for the BIC\n")
        } else {n <- nobs}
        
        LL <- logL
      
        if(c.hat == 1) {
            BIC <- -2*LL + K * log(n)
        }
        if(c.hat > 1 && c.hat <= 4) {
            K <- K+1
            BIC <- -2*LL/c.hat + K * log(n)
        }

        if(c.hat > 4) stop("High overdispersion and model fit is questionable\n")
        if(c.hat < 1) stop("You should set \'c.hat\' to 1 if < 1, but values << 1 might also indicate lack of fit\n")

        if(return.K == TRUE) BIC <- K #attributes the first element of BIC to K
        BIC
    }



##model selection with BIC
bictabCustom <-
    function(logL, K, modnames = NULL, nobs = NULL, sort = TRUE, c.hat = 1){
    
        ##check if modnames are not supplied
        if(is.null(modnames)) {
            modnames <- paste("Mod", 1:length(logL), sep = "")
            warning("\nModel names have been supplied automatically in the table\n")
        }
        
        ##check that nobs is the same for all models
        if(length(unique(nobs)) > 1) stop("\nSample size must be identical for all models\n")
    
        ##check that logL, K, estimate, se are vectors of same length
        nlogL <- length(logL)
        nK <- length(K)

        if(!all(nlogL == c(nlogL, nK))) stop("\nArguments 'logL' and 'K' must be of equal length\n")

        ##create model selection table
        Results <- NULL
        Results <- data.frame(Modnames = modnames)                    #assign model names to first column
        Results$K <- useBICCustom(logL = logL, K = K, return.K = TRUE,
                                  nobs = nobs, c.hat = c.hat)  #extract number of parameters
        Results$BIC <- useBICCustom(logL = logL, K = K, nobs = nobs,
                                 c.hat = c.hat)  #extract BIC                                      
        Results$Delta_BIC <- Results$BIC - min(Results$BIC)          #compute delta BIC
        Results$ModelLik <- exp(-0.5*Results$Delta_BIC)                #compute model likelihood required to compute Akaike weights
        Results$BICWt <- Results$ModelLik/sum(Results$ModelLik)        #compute BIC weights

        ##check if some models are redundant
        if(length(unique(Results$BIC)) != length(K)) warning("\nCheck model structure carefully as some models may be redundant\n")
         
        ##check if BIC and c.hat = 1
        if(c.hat == 1) {
            Results$LL <- logL
        }
    
        ##rename correctly to QBIC and add column for c-hat
        if(c.hat > 1) {
            colnames(Results) <- c("Modnames", "K", "QBIC", "Delta_QBIC", "ModelLik", "QBICWt")
            LL <- logL
            Results$Quasi.LL <- LL/c.hat
            Results$c_hat <- c.hat
        }      

    
        if(sort)  {
            Results <- Results[order(Results[, 4]),] 	  #if sort=TRUE, models are ranked based on Akaike weights
            Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
        } else {Results$Cum.Wt <- NULL}

        
        class(Results) <- c("bictab", "data.frame")
        return(Results)
    }



##Custom model averaging where user inputs estimates and SE's manually
##convenient when model type or SE's not available from predict methods
modavgCustom <-
  function(logL, K, modnames = NULL, estimate, se, second.ord = TRUE,
		 nobs = NULL, uncond.se = "revised", conf.level = 0.95, c.hat = 1,
		 useBIC = FALSE){

    ##check if modnames are not supplied
    if(is.null(modnames)) {
      modnames <- paste("Mod", 1:length(logL), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    }
    

    ##check that logL, K, estimate, se are vectors of same length
    nlogL <- length(logL)
    nK <- length(K)
    nestimate <- length(estimate)
    nse <- length(se)

    if(!all(nlogL == c(nlogL, nK, nestimate, nse))) stop("\nArguments 'logL', 'K', 'estimate', and 'se' must be of equal length\n")


    ##compute table
    if(!useBIC) {
        new_table <- aictabCustom(modnames = modnames, logL = logL, K = K,
                                  second.ord = second.ord, nobs = nobs, sort = FALSE,
                                  c.hat = c.hat)  #recompute AIC table and associated measures
    }

    if(useBIC) {
        new_table <- bictabCustom(modnames = modnames, logL = logL, K = K,
                                  nobs = nobs, sort = FALSE, c.hat = c.hat)  #recompute AIC table and associated measures
    }

    new_table$Estimate <- estimate
    new_table$SE <- se
        
    ##if c-hat is estimated adjust the SE's by multiplying with sqrt of c-hat
    if(c.hat > 1) {
        new_table$SE <-new_table$SE*sqrt(c.hat)
    } 


    if(!useBIC) {
        ##AICc
        ##compute model-averaged estimates, unconditional SE, and 95% CL
        if(c.hat == 1 && second.ord == TRUE) {
            Modavg_estimate <- sum(new_table$AICcWt*new_table$Estimate)
            
            ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
            if(identical(uncond.se, "old")) {
                Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Estimate- Modavg_estimate)^2))
            }
      
            ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
            if(identical(uncond.se, "revised")) {
                Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Estimate- Modavg_estimate)^2)))
            }
        }

    
    
        ##QAICc
        if(c.hat > 1 && second.ord == TRUE) {
            Modavg_estimate <- sum(new_table$QAICcWt*new_table$Estimate)
            
            ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
            if(identical(uncond.se, "old")) {      
                Uncond_SE <- sum(new_table$QAICcWt*sqrt(new_table$SE^2 + (new_table$Estimate- Modavg_estimate)^2))
            }
            
            ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
            if(identical(uncond.se, "revised")) {
                Uncond_SE <- sqrt(sum(new_table$QAICcWt*(new_table$SE^2 + (new_table$Estimate- Modavg_estimate)^2)))
            }
        }     
        


        ##AIC
        if(c.hat == 1 && second.ord == FALSE) {
            Modavg_estimate <- sum(new_table$AICWt*new_table$Estimate)
            
            ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
            if(identical(uncond.se, "old")) {
                Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Estimate- Modavg_estimate)^2))
            }
      
            ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
            if(identical(uncond.se, "revised")) {
                Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Estimate- Modavg_estimate)^2)))
            }
        }
        


        ##QAIC
        if(c.hat > 1 && second.ord == FALSE) {
            Modavg_estimate <- sum(new_table$QAICWt*new_table$Estimate)
            
            ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
            if(identical(uncond.se, "old")) {
                Uncond_SE <- sum(new_table$QAICWt*sqrt(new_table$SE^2 + (new_table$Estimate- Modavg_estimate)^2))
            }
            
            ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
            if(identical(uncond.se, "revised")) {
                Uncond_SE <- sqrt(sum(new_table$QAICWt*(new_table$SE^2 + (new_table$Estimate- Modavg_estimate)^2)))
            }  
        }     
    }

    if(useBIC){
        ##BIC
        if(c.hat == 1) {
            Modavg_estimate <- sum(new_table$BICWt*new_table$Estimate)
            
            ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
            if(identical(uncond.se, "old")) {
                Uncond_SE <- sum(new_table$BICWt*sqrt(new_table$SE^2 + (new_table$Estimate- Modavg_estimate)^2))
            }
      
            ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
            if(identical(uncond.se, "revised")) {
                Uncond_SE <- sqrt(sum(new_table$BICWt*(new_table$SE^2 + (new_table$Estimate- Modavg_estimate)^2)))
            }
        }
        


        ##QBIC
        if(c.hat > 1) {
            Modavg_estimate <- sum(new_table$QBICWt*new_table$Estimate)
            
            ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
            if(identical(uncond.se, "old")) {
                Uncond_SE <- sum(new_table$QBICWt*sqrt(new_table$SE^2 + (new_table$Estimate- Modavg_estimate)^2))
            }
            
            ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
            if(identical(uncond.se, "revised")) {
                Uncond_SE <- sqrt(sum(new_table$QBICWt*(new_table$SE^2 + (new_table$Estimate- Modavg_estimate)^2)))
            }  
        }     
        
    }
        
    zcrit <- qnorm(p = (1-conf.level)/2, lower.tail = FALSE)
    Lower_CL <- Modavg_estimate - zcrit*Uncond_SE
    Upper_CL <- Modavg_estimate + zcrit*Uncond_SE
    out.modavg <- list("Mod.avg.table" = new_table, "Mod.avg.est" = Modavg_estimate,
                       "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level,
                       "Lower.CL"= Lower_CL, "Upper.CL" = Upper_CL)
  
    class(out.modavg) <- c("modavgCustom", "list")
    return(out.modavg)

}



##print method
print.modavgCustom <- function(x, digits = 2, ...) {
ic <- colnames(x$Mod.avg.table)[3]
cat("\nMultimodel inference on manually-supplied parameter based on", ic, "\n")
cat("\n", ic, "table used to obtain model-averaged estimate:\n")
oldtab <- x$Mod.avg.table
if (any(names(oldtab)=="c_hat")) {cat("\t(c-hat estimate = ", oldtab$c_hat[1], ")\n")}
cat("\n")
if (any(names(oldtab)=="c_hat")) {
  nice.tab <- cbind(oldtab[,2], oldtab[,3], oldtab[,4], oldtab[,6],
                    oldtab[,9], oldtab[,10])
} else {nice.tab <- cbind(oldtab[,2], oldtab[,3], oldtab[,4], oldtab[,6],
                          oldtab[,8], oldtab[,9])
      }
    
colnames(nice.tab) <- c(colnames(oldtab)[c(2, 3, 4, 6)], "Estimate", "SE")
rownames(nice.tab) <- oldtab[, 1]
print(round(nice.tab, digits = digits))
cat("\nModel-averaged estimate:", eval(round(x$Mod.avg.est, digits = digits)), "\n")
cat("Unconditional SE:", eval(round(x$Uncond.SE, digits = digits)), "\n")
cat("",x$Conf.level*100, "% Unconditional confidence interval:", round(x$Lower.CL, digits = digits),
    ",", round(x$Upper.CL, digits = digits), "\n\n")
}



##function for generic information criteria
ictab <- function(ic, K, modnames = NULL, sort = TRUE, ic.name = NULL){

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
        modnames <- paste("Mod", 1:length(ic), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
    }
      
    ##check that logL, K, estimate, se are vectors of same length
    nic <- length(ic)
    nK <- length(K)

    if(!all(nic == c(nic, nK))) stop("\nArguments 'ic' and 'K' must be of equal length\n")

    Results <- NULL
    Results <- data.frame(Modnames = modnames)                    #assign model names to first column
    Results$K <- K
    Results$IC <- ic
    Results$Delta_IC <- Results$IC - min(Results$IC)            #compute delta IC
    Results$ModelLik <- exp(-0.5*Results$Delta_IC)                #compute model likelihood required to compute Akaike weights
    Results$ICWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights
    
    ##check if some models are redundant
    if(length(unique(Results$IC)) != length(K)) warning("\nCheck model structure carefully as some models may be redundant\n")

    ##check if name of IC is specified
    if(!is.null(ic.name)) {
        ##replace IC with ic.name
        names(Results) <- gsub(pattern = "IC", replacement = ic.name,
                               x = names(Results))
    }
    
    if(sort)  {
        Results <- Results[order(Results[, 4]),] 	  #if sort=TRUE, models are ranked based on delta IC
        Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
    } else {Results$Cum.Wt <- NULL}

    
    class(Results) <- c("ictab", "data.frame")
    return(Results)
}



print.ictab <- function(x, digits = 2, ...) {
    cat("\nModel selection based on ", colnames(x)[3], ":\n", sep = "")
    cat("\n")

    #check if Cum.Wt should be printed
    if(any(names(x) == "Cum.Wt")) {
        nice.tab <- cbind(x[, c(2:4, 6:7)])
        colnames(nice.tab) <- colnames(x)[c(2:4, 6:7)]
        rownames(nice.tab) <- x[, 1]
    } else {
        nice.tab <- cbind(x[, c(2:4, 6)])
        colnames(nice.tab) <- colnames(x)[c(2:4, 6)]
        rownames(nice.tab) <- x[, 1]
    }
    
    print(round(nice.tab, digits = digits)) #select rounding off with digits argument
    cat("\n")
  }



##model averaging for generic IC where user inputs estimates and SE's manually
modavgIC <- function(ic, K, modnames = NULL, estimate, se,
                     uncond.se = "revised", conf.level = 0.95, ic.name = NULL){

    ##check if modnames are not supplied
    if(is.null(modnames)) {
        modnames <- paste("Mod", 1:length(ic), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
    }
    

    ##check that logL, K, estimate, se are vectors of same length
    nic <- length(ic)
    nK <- length(K)
    nestimate <- length(estimate)
    nse <- length(se)

    if(!all(nic == c(nic, nK, nestimate, nse))) stop("\nArguments 'ic', 'K', 'estimate', and 'se' must be of equal length\n")

    ##compute table
    new_table <- ictab(ic = ic, K = K, modnames = modnames,
                       sort = FALSE, ic.name = ic.name)

    new_table$Estimate <- estimate
    new_table$SE <- se
    
    ##compute model-averaged estimates, unconditional SE, and 95% CL
    Modavg_estimate <- sum(new_table[, 6] * new_table$Estimate)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table[, 6] * sqrt(new_table$SE^2 + (new_table$Estimate- Modavg_estimate)^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table[, 6] * (new_table$SE^2 + (new_table$Estimate - Modavg_estimate)^2)))
    }
        
    zcrit <- qnorm(p = (1-conf.level)/2, lower.tail = FALSE)
    Lower_CL <- Modavg_estimate - zcrit*Uncond_SE
    Upper_CL <- Modavg_estimate + zcrit*Uncond_SE
    out.modavg <- list("Mod.avg.table" = new_table, "Mod.avg.est" = Modavg_estimate,
                       "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level,
                       "Lower.CL"= Lower_CL, "Upper.CL" = Upper_CL)
  
    class(out.modavg) <- c("modavgIC", "list")
    return(out.modavg)
}



##print method
print.modavgIC <- function(x, digits = 2, ...) {
    ic <- colnames(x$Mod.avg.table)[3]
    cat("\nMultimodel inference on manually-supplied parameter based on", ic, "\n")
    cat("\n", ic, "table used to obtain model-averaged estimate:\n")
    oldtab <- x$Mod.avg.table
    cat("\n")
    nice.tab <- cbind(oldtab[, c(2:4, 6:8)])
    
    colnames(nice.tab) <- colnames(oldtab)[c(2:4, 6:8)]
    rownames(nice.tab) <- oldtab[, 1]
    print(round(nice.tab, digits = digits))
    cat("\nModel-averaged estimate:", eval(round(x$Mod.avg.est, digits = digits)), "\n")
    cat("Unconditional SE:", eval(round(x$Uncond.SE, digits = digits)), "\n")
    cat("",x$Conf.level*100, "% Unconditional confidence interval:", round(x$Lower.CL, digits = digits),
        ",", round(x$Upper.CL, digits = digits), "\n\n")
}


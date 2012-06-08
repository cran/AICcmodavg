modavg.shrink.unmarked <- 
  function(cand.set, parm, modnames, c.hat = 1, conf.level = 0.95, second.ord = TRUE, nobs = NULL,
           uncond.se = "revised", parm.type = NULL){
    ##note that parameter is referenced differently from unmarked object - see labels( )

    ##check model types
    mod.class <- unlist(lapply(cand.set, FUN = function(i) class(i)[1]))
    mod.type <- unique(mod.class)

    if(length(mod.type) > 1) stop("\nThis function is not appropriate to model-average parameters from different model types\n")
    
    ##check for supported mod.type
    supp.class <- c("unmarkedFitOccu", "unmarkedFitColExt", "unmarkedFitOccuRN", "unmarkedFitPCount", "unmarkedFitPCO",
                    "unmarkedFitDS", "unmarkedFitGDS")
                  
    if(!any(supp.class == mod.type)) {stop("\nFunction not yet defined for this object class\n")}


    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavg for details\n")}
  
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"
   
    ##single-season occupancy model
    if(identical(mod.type, "unmarkedFitOccu")) {
      ##psi
      if(identical(parm.type, "psi")) {
        ##extract model formula for each model in cand.set
        mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$state)))
        parm.unmarked <- "psi"
        parm <- paste(parm.unmarked, "(", parm, ")", sep="")
      }
      ##detect
      if(identical(parm.type, "detect")) {
        mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$det)))
        parm.unmarked <- "p"
        parm <- paste(parm.unmarked, "(", parm, ")", sep="")
      }
    }

  
    ##multiseason occupancy model
    if(identical(mod.type, "unmarkedFitColExt")) {
      ##psi - initial occupancy
      if(identical(parm.type, "psi")) {
        ##extract model formula for each model in cand.set
        mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$psi)))
        ##create label for parm
        parm.unmarked <- "psi"
        parm <- paste(parm.unmarked, "(", parm, ")", sep="")
      }
      ##gamma - extinction
      if(identical(parm.type, "gamma")) {
        ##extract model formula for each model in cand.set
        mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$col)))
        ##create label for parm
        parm.unmarked <- "col"
        parm <- paste(parm.unmarked, "(", parm, ")", sep="")
      }
      ##epsilon - extinction
      if(identical(parm.type, "epsilon")) {
        ##extract model formula for each model in cand.set
        mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$ext)))
        ##create label for parm
        parm.unmarked <- "ext"
        parm <- paste(parm.unmarked, "(", parm, ")", sep="")
      }
      ##detect
      if(identical(parm.type, "detect")) {
        mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$det)))
        parm.unmarked <- "p"
        parm <- paste(parm.unmarked, "(", parm, ")", sep="")
      }
    }

  
    ##single season N-mixture model
    if(identical(mod.type, "unmarkedFitPCount")) {
      ##lambda - abundance
      if(identical(parm.type, "lambda")) {
        ##extract model formula for each model in cand.set
        mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$state)))
        ##create label for parm
        parm.unmarked <- "lam"
        parm <- paste(parm.unmarked, "(", parm, ")", sep="")
      }
      ##detect
      if(identical(parm.type, "detect")) {
        mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$det)))
        parm.unmarked <- "p"
        parm <- paste(parm.unmarked, "(", parm, ")", sep="")
      }
    }
    
    ##open version of N-mixture model
    if(identical(mod.type, "unmarkedFitPCO")) {
      ##lambda - abundance
      if(identical(parm.type, "lambda")) {
        ##extract model formula for each model in cand.set
        mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$lambda)))
        parm.unmarked <- "lam"
        parm <- paste(parm.unmarked, "(", parm, ")", sep="")
      }
      ##gamma - recruitment
      if(identical(parm.type, "gamma")) {
        ##extract model formula for each model in cand.set
        mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$gamma)))
        ##create label for parm
        parm.unmarked <- "gam"
        parm <- paste(parm.unmarked, "(", parm, ")", sep="")
      }
      ##omega - apparent survival
      if(identical(parm.type, "omega")) {
        ##extract model formula for each model in cand.set
        mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$omega)))
        ##create label for parm
        parm.unmarked <- "omega"
        parm <- paste(parm.unmarked, "(", parm, ")", sep="")
      }
      ##detect
      if(identical(parm.type, "detect")) {
        mod_formula<-lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$det)))
        parm.unmarked <- "p"
        parm <- paste(parm.unmarked, "(", parm, ")", sep="")
      }
    }

    ##Royle-Nichols heterogeneity model
    if(identical(mod.type, "unmarkedFitOccuRN")) {
      ##lambda - abundance
      if(identical(parm.type, "lambda")) {
        ##extract model formula for each model in cand.set
        mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$state)))
        parm.unmarked <- "lam"
        parm <- paste(parm.unmarked, "(", parm, ")", sep="")
      }
      ##detect
      if(identical(parm.type, "detect")) {
        mod_formula<-lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$det)))
        parm.unmarked <- "p"
        parm <- paste(parm.unmarked, "(", parm, ")", sep="")
      }
    }


      ##Distance sampling model
  if(identical(mod.type, "unmarkedFitDS")) {
    ##lambda - abundance
    if(identical(parm.type, "lambda")) {
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$state)))
      parm <- paste("lam", "(", parm, ")", sep="")
      if(!is.null(reversed.parm)) {reversed.parm <- paste("lam", "(", reversed.parm, ")", sep="")}
      not.include <- lapply(cand.set, FUN = function(i) i@formula[[3]])
    }
    ##detect
    if(identical(parm.type, "detect")) {
      if(identical(parm.type, "detect")) {
        stop("\nModel-averaging estimates of detection covariates not yet supported for unmarkedFitDS class\n")
      }
    }
  }



  ##Distance sampling model with availability
  if(identical(mod.type, "unmarkedFitDS")) {
    ##lambda - abundance
    if(identical(parm.type, "lambda")) {
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$state)))
      parm <- paste("lam", "(", parm, ")", sep="")
      if(!is.null(reversed.parm)) {reversed.parm <- paste("lam", "(", reversed.parm, ")", sep="")}
      not.include <- lapply(cand.set, FUN = function(i) i@formlist$lambdaformula)
    }
    ##detect
    if(identical(parm.type, "detect")) {
      if(identical(parm.type, "detect")) {
        stop("\nModel-averaging estimates of detection covariates not yet supported for unmarkedFitGDS class\n")
      }
    }
    ##availability
    if(identical(parm.type, "phi")) {
      stop("\nModel-averaging estimates of availability covariates not yet supported for unmarkedFitGDS class\n")
    }
  }


    ##NEED TO PASTE THE PARAMETER TYPE - INCLUDE THIS STEP ABOVE FOR EACH PARM.TYPE
    ##determine frequency of each term across models (except (Intercept) ) 
    pooled.terms <- unlist(mod_formula)
    ##remove intercept from vector
    no.int <- pooled.terms[which(pooled.terms != paste(parm.unmarked, "(Int)", sep = ""))]
    terms.freq <- table(no.int)
    if(length(unique(terms.freq)) > 1) stop("\n\nTo compute a shrinkage version of model-averaged estimate, each term must appear with the same frequency across models\n")

    ##check whether parm is involved in interaction
    parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
    inter.check <- ifelse(attr(regexpr(parm.inter[1], mod_formula, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                     mod_formula, fixed = TRUE), "match.length") == "-1", 0, 1)
    if(sum(inter.check) > 0) stop("\nParameter of interest should not be involved in interaction for shrinkage version of model-averaging to be appropriate\n")

    
    nmods <- length(cand.set)
  
    
    new_table <- aictab.unmarked(cand.set = cand.set, modnames = modnames, sort = FALSE, c.hat = c.hat,
                               second.ord = second.ord, nobs = nobs)  #recompute AIC table and associated measures
    new_table$Beta_est <- unlist(lapply(cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
    new_table$SE <- unlist(lapply(cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))

    ##replace NA's with 0
    new_table$Beta_est[is.na(new_table$Beta_est)] <- 0
    new_table$SE[is.na(new_table$SE)] <- 0

    ##add a check to determine if parameter occurs in any model
    if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("\nParameter not found in any of the candidate models\n") }

    
    ##if c-hat is estimated adjust the SE's by multiplying with sqrt of c-hat
    if(c.hat > 1) {
      new_table$SE <- new_table$SE*sqrt(c.hat)
    } 

    ##AICc
    ##compute model-averaged estimates, unconditional SE, and 95% CL
    if(c.hat == 1 && second.ord == TRUE) {
      Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }

    ##QAICc
    if(c.hat > 1 && second.ord == TRUE) {
      Modavg_beta <- sum(new_table$QAICcWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {      
        Uncond_SE <- sum(new_table$QAICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$QAICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }     

    
    ##AIC
    if(c.hat == 1 && second.ord == FALSE) {
      Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }


    ##QAIC
    if(c.hat > 1 && second.ord == FALSE) {
      Modavg_beta <- sum(new_table$QAICWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$QAICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$QAICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }  
    }     

  
    zcrit <- qnorm(p = (1 - conf.level)/2, lower.tail = FALSE)
    Lower_CL <- Modavg_beta - zcrit*Uncond_SE
    Upper_CL <- Modavg_beta + zcrit*Uncond_SE
    out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                       "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                       "Upper.CL" = Upper_CL)
    
    class(out.modavg) <- c("modavg.shrink", "list")
    return(out.modavg)
  }

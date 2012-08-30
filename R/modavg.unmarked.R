##added functionality for reversing parameters
##added additional argument parm.type = "psi", "gamma", "epsilon", "lambda", "omega", "detect"
##model type:  parameters labeled in unmarked - parameters labeled in AICcmodavg.unmarked
##single season: state, det - USE psi, detect
##multiseason model:  psi, col, ext, det - USE psi, gamma, epsilon, detect
##RN heterogeneity model: state, det - USE lambda, detect
##N-mixture: state, det - USE lambda, detect
##Open N-mixture: lambda, gamma, omega, det - USE lambda, gamma, omega, detect
##distsamp: state, det - USE lambda, detect
##gdistsamp: state, det, phi - USE lambda, detect, phi

modavg.unmarked <-
function(cand.set, parm, modnames, c.hat = 1, conf.level = 0.95, second.ord = TRUE, nobs = NULL, exclude = NULL,
         warn = TRUE, uncond.se = "revised", parm.type = NULL){
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
  
#####MODIFICATIONS BEGIN#######
  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)
  parm.strip <- parm #to use later
  
  ##if (Intercept) is chosen assign (Int) - for compatibility
  if(identical(parm, "(Intercept)")) parm <- "Int"

  ##reverse parm
  reversed.parm <- reverse.parm(parm)
  reversed.parm.strip <- reversed.parm #to use later
  exclude <- reverse.exclude(exclude = exclude)
#####MODIFICATIONS END######
  
  
  ##single-season occupancy model
  if(identical(mod.type, "unmarkedFitOccu")) {
    ##psi
    if(identical(parm.type, "psi")) {
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$state)))
      parm <- paste("psi", "(", parm, ")", sep="")
      if(!is.null(reversed.parm)) {reversed.parm <- paste("psi", "(", reversed.parm, ")", sep="")}
      not.include <- lapply(cand.set, FUN = function(i) i@formula[[3]])      
    }
    ##detect
    if(identical(parm.type, "detect")) {
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$det)))
      parm <- paste("p", "(", parm, ")", sep="")
      if(!is.null(reversed.parm)) {reversed.parm <- paste("p", "(", reversed.parm, ")", sep="")}
      not.include <- lapply(cand.set, FUN = function(i) i@formula[[2]])
    }
  }

  
  ##multiseason occupancy model
  if(identical(mod.type, "unmarkedFitColExt")) {
    ##psi - initial occupancy
    if(identical(parm.type, "psi")) {
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$psi)))
      ##create label for parm
      parm <- paste("psi", "(", parm, ")", sep="")
      if(!is.null(reversed.parm)) {reversed.parm <- paste("psi", "(", reversed.parm, ")", sep="")}
      not.include <- lapply(cand.set, FUN = function(i) i@psiformula)
    }
    ##gamma - extinction
    if(identical(parm.type, "gamma")) {
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$col)))
      ##create label for parm
      parm <- paste("col", "(", parm, ")", sep="")
      if(!is.null(reversed.parm)) {reversed.parm <- paste("col", "(", reversed.parm, ")", sep="")}
      not.include <- lapply(cand.set, FUN = function(i) i@gamformula)
    }
    ##epsilon - extinction
    if(identical(parm.type, "epsilon")) {
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$ext)))
      ##create label for parm
      parm <- paste("ext", "(", parm, ")", sep="")
      if(!is.null(reversed.parm)) {reversed.parm <- paste("ext", "(", reversed.parm, ")", sep="")}
      ##for epsilon
      not.include <- lapply(cand.set, FUN = function(i) i@epsformula) 
    }
    ##detect
    if(identical(parm.type, "detect")) {
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$det)))
      parm <- paste("p", "(", parm, ")", sep="")
      if(!is.null(reversed.parm)) {reversed.parm <- paste("p", "(", reversed.parm, ")", sep="")}
      not.include <- lapply(cand.set, FUN = function(i) i@detformula)
    }
  }

  
  ##single season N-mixture model
  if(identical(mod.type, "unmarkedFitPCount")) {
    ##lambda - abundance
    if(identical(parm.type, "lambda")) {
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$state)))
      ##create label for parm
      parm <- paste("lam", "(", parm, ")", sep="")
      if(!is.null(reversed.parm)) {reversed.parm <- paste("lam", "(", reversed.parm, ")", sep="")}
      not.include <- lapply(cand.set, FUN = function(i) i@formula[[3]])
    }
    ##detect
    if(identical(parm.type, "detect")) {
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$det)))
      parm <- paste("p", "(", parm, ")", sep="")
      if(!is.null(reversed.parm)) {reversed.parm <- paste("p", "(", reversed.parm, ")", sep="")}
      not.include <- lapply(cand.set, FUN = function(i) i@formula[[2]])
    }
  }
    
    ##open version of N-mixture model
    if(identical(mod.type, "unmarkedFitPCO")) {
      ##lambda - abundance
      if(identical(parm.type, "lambda")) {
        ##extract model formula for each model in cand.set
        mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$lambda)))
        parm <- paste("lam", "(", parm, ")", sep="")
        if(!is.null(reversed.parm)) {reversed.parm <- paste("lam", "(", reversed.parm, ")", sep="")}
        not.include <- lapply(cand.set, FUN = function(i) i@formlist$lambdaformula)
      }
      ##gamma - recruitment
      if(identical(parm.type, "gamma")) {
        ##extract model formula for each model in cand.set
        mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$gamma)))

        ##determine if same H0 on gamma (gamConst, gamAR, gamTrend)
        strip.gam <- sapply(mod_formula, FUN = function(i) unlist(strsplit(i, "\\("))[[1]])
        unique.gam <- unique(strip.gam)
        if(length(unique.gam) > 1) stop("\nDifferent formulations of gamma parameter occur among models:\n
beta estimates cannot be model-averaged\n")

        ##create label for parm
        parm <- paste(unique.gam, "(", parm, ")", sep="")
        if(!is.null(reversed.parm)) {reversed.parm <- paste(unique.gam, "(", reversed.parm, ")", sep="")}
        not.include <- lapply(cand.set, FUN = function(i) i@formlist$gammaformula)
      }
      ##omega - apparent survival
      if(identical(parm.type, "omega")) {
        ##extract model formula for each model in cand.set
        mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$omega)))   
        ##create label for parm
        parm <- paste("omega", "(", parm, ")", sep="")
        if(!is.null(reversed.parm)) {reversed.parm <- paste("omega", "(", reversed.parm, ")", sep="")}
        not.include <- lapply(cand.set, FUN = function(i) i@formlist$omegaformula)
      }
      ##detect
      if(identical(parm.type, "detect")) {
        mod_formula<-lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$det)))
        parm <- paste("p", "(", parm, ")", sep="")
        if(!is.null(reversed.parm)) {reversed.parm <- paste("p", "(", reversed.parm, ")", sep="")}
        not.include <- lapply(cand.set, FUN = function(i) i@formlist$pformula)
      }
    }



    ##Royle-Nichols heterogeneity model
    if(identical(mod.type, "unmarkedFitOccuRN")) {
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
        mod_formula<-lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$det)))
        parm <- paste("p", "(", parm, ")", sep="")
        if(!is.null(reversed.parm)) {reversed.parm <- paste("p", "(", reversed.parm, ")", sep="")}
        not.include <- lapply(cand.set, FUN = function(i) i@formula[[2]])
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
  if(identical(mod.type, "unmarkedFitGDS")) {
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
  

  nmods <- length(cand.set)
  
  ##setup matrix to indicate presence of parms in the model
  include <- matrix(NA, nrow=nmods, ncol=1)
  ##add a check for multiple instances of same variable in given model (i.e., interactions)
  include.check <- matrix(NA, nrow=nmods, ncol=1)

##################################
##################################
###ADDED A NEW OBJECT TO STRIP AWAY lam( ) from parm on line 35
###to enable search with regexpr( ) 
  
  
  ##iterate over each formula in mod_formula list
  for (i in 1:nmods) {
    idents <- NULL
    idents.check <- NULL
    form <- mod_formula[[i]]


######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
      ##iterate over each element of formula[[i]] in list
      if(is.null(reversed.parm)) {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j])
          ##added parm.strip here for regexpr( )
          idents.check[j] <- ifelse(attr(regexpr(parm.strip, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)  
        }
      } else {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm.strip, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm.strip, form[j],
                                                                  fixed=TRUE), "match.length")=="-1" , 0, 1)  
        }
      }
###MODIFICATIONS END
######################################################################################################
######################################################################################################

    
    include[i] <- ifelse(any(idents==1), 1, 0)
    include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
  }

  #####################################################
  #exclude == NULL; warn=TRUE:  warn that duplicates occur and stop
  if(is.null(exclude) && identical(warn, TRUE)) {
    #check for duplicates in same model
    if(any(include.check == "duplicates")) {
      stop("\nSome models include more than one instance of the parameter of interest. \n",
           "This may be due to the presence of interaction/polynomial terms, or variables\n",
             "with similar names:\n",
             "\tsee \"?modavg\" for details on variable specification and \"exclude\" argument\n")
      }

    }

    #exclude == NULL; warn=FALSE:  compute model-averaged beta estimate from models including variable of interest,
    #assuming that the variable is not involved in interaction or higher order polynomial (x^2, x^3, etc...),
    #warn that models were not excluded
    if(is.null(exclude) && identical(warn, FALSE)) {
      if(any(include.check == "duplicates")) {
        warning("\nMultiple instances of parameter of interest in given model is presumably\n",
                "not due to interaction or polynomial terms - these models will not be\n",
                "excluded from the computation of model-averaged estimate\n")
      }
      
    }

    #warn if exclude is neither a list nor NULL
    if(!is.null(exclude)) {
      if(!is.list(exclude)) {stop("\nItems in \"exclude\" must be specified as a list")}
    }


    #if exclude is list  
    if(is.list(exclude)) {

    #determine number of elements in exclude
      nexcl <- length(exclude)

      #check each formula for presence of exclude variable in not.include list
      #not.include <- lapply(cand.set, FUN = formula)

      #set up a new list with model formula
      forms <- list()
      for (i in 1:nmods) {
        form.tmp <- as.character(not.include[i]) #changed from other versions as formula returned is of different structure for unmarked objects
        if(attr(regexpr("\\+", form.tmp), "match.length")==-1) {
          forms[i] <- form.tmp
        } else {forms[i] <- strsplit(form.tmp, split=" \\+ ")}
      }

      #additional check to see whether some variable names include "+"
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\+", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nPlease avoid \"+\" in variable names\n")

      ##additional check to determine if intercept was removed from models
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\- 1", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nModels without intercept are not supported in this version, please use alternative parameterization\n")

 
      #search within formula for variables to exclude
      mod.exclude <- matrix(NA, nrow=nmods, ncol=nexcl)

      #iterate over each element in exclude list
      for (var in 1:nexcl) {

      #iterate over each formula in mod_formula list
        for (i in 1:nmods) {
          idents <- NULL
          form.excl <- forms[[i]]

          #iterate over each element of forms[[i]]
          for (j in 1:length(form.excl)) {
            idents[j] <- identical(exclude[var][[1]], form.excl[j])
          }
          mod.exclude[i,var] <- ifelse(any(idents == 1), 1, 0)
        }    
        
      }
  
      #determine outcome across all variables to exclude
      to.exclude <- rowSums(mod.exclude)
  
  
      #exclude models following models from model averaging  
      include[which(to.exclude >= 1)] <- 0
      
      
    }


 
  ##add a check to determine if include always == 0
  if (sum(include) == 0) {stop("\nParameter not found in any of the candidate models\n") }
  
  new.cand.set <- cand.set[which(include == 1)] #select models including a given parameter
  new.mod.name <- modnames[which(include == 1)]    #update model names

  new_table <- aictab.unmarked(cand.set = new.cand.set, modnames = new.mod.name, sort = FALSE, c.hat = c.hat,
                        second.ord = second.ord, nobs = nobs)  #recompute AIC table and associated measures
  new_table$Beta_est <- unlist(lapply(new.cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
  new_table$SE <- unlist(lapply(new.cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))

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
  
  class(out.modavg) <- c("modavg", "list")
  return(out.modavg)
}

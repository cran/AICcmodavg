importance <-
  function(cand.set, parm, modnames, c.hat = 1, second.ord = TRUE, nobs = NULL, parm.type = NULL){
    ##check if class is appropriate
    ##extract classes
    mod.class <- unlist(lapply(X = cand.set, FUN = class))
  
    ##check if all are identical
    check.class <- unique(mod.class)

    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##reverse parm
    reversed.parm <- reverse.parm(parm)

    ##add check for supported classes
    known <- rep(0, 9) #create an identifier of class type other than lm, glm, multinom, polr, lme, gls, mer, unmarked, or nls
    
    if(identical(check.class, "lm") || identical(check.class, c("glm", "lm")))  {
    
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN=function(i) rownames(summary(i)$coefficients))
      known[1] <- 1
    }
  
    if(identical(check.class, "lme")) {
      mod_formula <- lapply(cand.set, FUN=function(i) labels(summary(i)$coefficients$fixed))
      known[2] <- 1
    }

    if(identical(check.class, "gls")) {
      mod_formula <- lapply(cand.set, FUN=function(i) labels(summary(i)$coefficients))
      known[3] <- 1
    }

    if(identical(check.class, c("multinom", "nnet"))) {
      mod_formula <- lapply(cand.set, FUN=function(i) colnames(summary(i)$coefficients))
      known[4] <- 1
    }

    if(identical(check.class, "mer")) {
      mod_formula <- lapply(cand.set, FUN=function(i) labels(fixef(i)))
      known[5] <- 1
    }

    if(identical(check.class, "polr")) {
      mod_formula <- lapply(cand.set, FUN=function(i) rownames(summary(i)$coefficients))
      known[6] <- 1
    }

    
    ##determine if unmarked
    unmarked.class <- c("unmarkedFitOccu", "unmarkedFitColExt", "unmarkedFitOccuRN", "unmarkedFitPCount", "unmarkedFitPCO")
    if(any(sapply(unmarked.class, FUN = function(i) identical(i, check.class)))) {

      known[7] <- 1
      
      ##check for parm.type and stop if NULL
      if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavg for details\n")}

      ##if (Intercept) is chosen assign (Int) - for compatibility
      if(identical(parm, "(Intercept)")) parm <- "Int"
      
      ##single-season occupancy model
      if(identical(check.class, "unmarkedFitOccu")) {
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
      if(identical(check.class, "unmarkedFitColExt")) {
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
      if(identical(check.class, "unmarkedFitPCount")) {
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
      if(identical(check.class, "unmarkedFitPCO")) {
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
      if(identical(check.class, "unmarkedFitOccuRN")) {
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
    }
  

    ##warn if class is neither lm, glm, multinom, polr, lme, gls, nls, mer, nor unmarkedFit
    if(sum(known) < 1) {stop("Function not yet defined for this object class\n")}
    
    ##setup matrix to indicate presence of parm in the model
    include <- matrix(NA, nrow=length(cand.set), ncol=1)

    ##iterate over each formula in mod_formula list
    for (i in 1:length(cand.set)) {
      idents <- NULL
      form <- mod_formula[[i]]

      ##iterate over each element of formula[[i]] in list
      if(is.null(reversed.parm)) {
        for (j in 1:length(form)) {
          idents[j] <- identical(paste(parm), form[j])
        }
      } else {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
        }
      }
    

    include[i] <- ifelse(any(idents==1), 1, 0)
    }

    ##add a check to determine if include always == 0
    if (sum(include)==0) {stop("Parameter not found in any of the candidate models") }

    new_table <- aictab(cand.set = cand.set, modnames = modnames, sort = FALSE, c.hat = c.hat, second.ord = second.ord, nobs = nobs)  

    ##add a check to determine if the same number of models include and exlude parameter
    if (length(which(include == 1)) != length(which(include != 1)) ) {
      stop("Importance values are only meaningful when the number of models with and without parameter are equal")
    }

    w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
    w.minus <- 1 - w.plus
    imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
    class(imp) <- c("importance", "list")
    return(imp)
  }

##function for nicer printing of importance values
print.importance <- function(x, digits = 2, ...) {
  cat("\nImportance values of '", x$parm, "' :\n\n")
  cat("w+ (models including parameter): ", round(x$w.plus, digits = digits), "\n")
  cat("w- (models excluding parameter): ", round(x$w.minus, digits = digits), "\n")
  cat("\n")
}


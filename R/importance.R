##generic
importance <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, ...){
  cand.set <- formatCands(cand.set)
  UseMethod("importance", cand.set)
}



##default
importance.default <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, ...){
  stop("\nFunction not yet defined for this object class\n")
}



##aov
importance.AICaov.lm <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, ...){

  ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##reverse parm
    reversed.parm <- reverse.parm(parm)


    ##extract labels
    mod_formula <- lapply(cand.set, FUN=function(i) rownames(summary(i)$coefficients))


    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"


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
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

    new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                        sort = FALSE)  

    ##add a check to determine if the same number of models include and exlude parameter
    if (length(which(include == 1)) != length(which(include != 1)) ) {
      stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
    }

    w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
    w.minus <- 1 - w.plus
    imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
    class(imp) <- c("importance", "list")
    return(imp)
  }



##betareg
importance.AICbetareg <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, ...){

  ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##reverse parm
    reversed.parm <- reverse.parm(parm)


    ##extract labels
    ##determine if parameter is on mean or phi
    if(regexpr(pattern = "\\(phi\\)_", parm) == "-1") {
      mod_formula <- lapply(cand.set, FUN=function(i) rownames(summary(i)$coefficients$mean))
    } else {
      ##replace parm
      parm <- gsub(pattern = "\\(phi\\)_", "", parm)
      mod_formula <- lapply(cand.set, FUN=function(i) rownames(summary(i)$coefficients$precision))
    }

    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"


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
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

    new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                        sort = FALSE)  

    ##add a check to determine if the same number of models include and exlude parameter
    if (length(which(include == 1)) != length(which(include != 1)) ) {
      stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
    }

    w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
    w.minus <- 1 - w.plus
    imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
    class(imp) <- c("importance", "list")
    return(imp)
  }



##clm
importance.AICsclm.clm <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, ...){

  ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##reverse parm
    reversed.parm <- reverse.parm(parm)


    ##extract labels
    mod_formula <- lapply(cand.set, FUN = function(i) rownames(summary(i)$coefficients))

    
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"


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
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

    new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                        sort = FALSE)  

    ##add a check to determine if the same number of models include and exlude parameter
    if (length(which(include == 1)) != length(which(include != 1)) ) {
      stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
    }

    w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
    w.minus <- 1 - w.plus
    imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
    class(imp) <- c("importance", "list")
    return(imp)
  }


##clm
importance.AICclm <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, ...){

  ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##reverse parm
    reversed.parm <- reverse.parm(parm)


    ##extract labels
    mod_formula <- lapply(cand.set, FUN = function(i) rownames(summary(i)$coefficients))

    
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"


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
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

    new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                        sort = FALSE)  

    ##add a check to determine if the same number of models include and exlude parameter
    if (length(which(include == 1)) != length(which(include != 1)) ) {
      stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
    }

    w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
    w.minus <- 1 - w.plus
    imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
    class(imp) <- c("importance", "list")
    return(imp)
}


##clmm
importance.AICclmm <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, ...){

  ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##reverse parm
    reversed.parm <- reverse.parm(parm)


    ##extract labels
    mod_formula <- lapply(cand.set, FUN = function(i) rownames(summary(i)$coefficients))

    
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"


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
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

    new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                        sort = FALSE)  

    ##add a check to determine if the same number of models include and exlude parameter
    if (length(which(include == 1)) != length(which(include != 1)) ) {
      stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
    }

    w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
    w.minus <- 1 - w.plus
    imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
    class(imp) <- c("importance", "list")
    return(imp)
  }


##clogit
importance.AICclogit.coxph <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, ...){

  ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##reverse parm
    reversed.parm <- reverse.parm(parm)


    ##extract labels
    mod_formula <- lapply(cand.set, FUN=function(i) rownames(summary(i)$coefficients))

    
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"


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
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

    new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                        sort = FALSE)  

    ##add a check to determine if the same number of models include and exlude parameter
    if (length(which(include == 1)) != length(which(include != 1)) ) {
      stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
    }

    w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
    w.minus <- 1 - w.plus
    imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
    class(imp) <- c("importance", "list")
    return(imp)
  }



##coxme
importance.AICcoxme <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, ...){

  ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##reverse parm
    reversed.parm <- reverse.parm(parm)


    ##extract labels
    mod_formula <- lapply(cand.set, FUN=function(i) labels(coef(i)))

    
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"


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
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

    new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                        sort = FALSE)  

    ##add a check to determine if the same number of models include and exlude parameter
    if (length(which(include == 1)) != length(which(include != 1)) ) {
      stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
    }

    w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
    w.minus <- 1 - w.plus
    imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
    class(imp) <- c("importance", "list")
    return(imp)
  }



##coxph
importance.AICcoxph <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, ...){

  ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##reverse parm
    reversed.parm <- reverse.parm(parm)


    ##extract labels
    mod_formula <- lapply(cand.set, FUN=function(i) rownames(summary(i)$coefficients))

    
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"


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
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

    new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                        sort = FALSE)  

    ##add a check to determine if the same number of models include and exlude parameter
    if (length(which(include == 1)) != length(which(include != 1)) ) {
      stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
    }

    w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
    w.minus <- 1 - w.plus
    imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
    class(imp) <- c("importance", "list")
    return(imp)
  }



##glm
importance.AICglm.lm <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, c.hat = 1, ...){

  ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##reverse parm
    reversed.parm <- reverse.parm(parm)


    ##extract labels
    mod_formula <- lapply(cand.set, FUN=function(i) rownames(summary(i)$coefficients))


    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"


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
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

    new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                        sort = FALSE, c.hat = c.hat)  

    ##add a check to determine if the same number of models include and exlude parameter
    if (length(which(include == 1)) != length(which(include != 1)) ) {
      stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
    }

    w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
    w.minus <- 1 - w.plus
    imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
    class(imp) <- c("importance", "list")
    return(imp)
  }



##glmer
importance.AICglmerMod <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, ...){

  ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##reverse parm
    reversed.parm <- reverse.parm(parm)


    ##extract labels
    mod_formula <- lapply(cand.set, FUN = function(i) labels(fixef(i)))

    
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"


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
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

    new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                        sort = FALSE)  

    ##add a check to determine if the same number of models include and exlude parameter
    if (length(which(include == 1)) != length(which(include != 1)) ) {
      stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
    }

    w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
    w.minus <- 1 - w.plus
    imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
    class(imp) <- c("importance", "list")
    return(imp)
  }



##glmmTMB
importance.AICglmmTMB <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, c.hat = 1, ...){

  ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##reverse parm
    reversed.parm <- reverse.parm(parm)


    ##extract labels
    mod_formula <- lapply(cand.set, FUN = function(i) labels(fixef(i)$cond))

    
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"


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
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

    new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                        sort = FALSE, c.hat = c.hat)  

    ##add a check to determine if the same number of models include and exlude parameter
    if (length(which(include == 1)) != length(which(include != 1)) ) {
      stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
    }

    w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
    w.minus <- 1 - w.plus
    imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
    class(imp) <- c("importance", "list")
    return(imp)
}



##gls
importance.AICgls <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, ...){

  ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##reverse parm
    reversed.parm <- reverse.parm(parm)


    ##extract labels
    mod_formula <- lapply(cand.set, FUN=function(i) labels(summary(i)$coefficients))
    
    
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"


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
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

    new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                        sort = FALSE)  

    ##add a check to determine if the same number of models include and exlude parameter
    if (length(which(include == 1)) != length(which(include != 1)) ) {
      stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
    }

    w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
    w.minus <- 1 - w.plus
    imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
    class(imp) <- c("importance", "list")
    return(imp)
  }



##hurdle
importance.AIChurdle <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, ...){
  
  ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##reverse parm
    reversed.parm <- reverse.parm(parm)


    ##extract labels
    mod_formula <- lapply(cand.set, FUN=function(i) labels(coef(i)))


    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"


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
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

    new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                        sort = FALSE)  

    ##add a check to determine if the same number of models include and exlude parameter
    if (length(which(include == 1)) != length(which(include != 1)) ) {
      stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
    }

    w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
    w.minus <- 1 - w.plus
    imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
    class(imp) <- c("importance", "list")
    return(imp)
  }



##lm
importance.AIClm <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, ...){

  ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##reverse parm
    reversed.parm <- reverse.parm(parm)


    ##extract labels
    mod_formula <- lapply(cand.set, FUN=function(i) rownames(summary(i)$coefficients))


    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"


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
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

    new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                        sort = FALSE)  

    ##add a check to determine if the same number of models include and exlude parameter
    if (length(which(include == 1)) != length(which(include != 1)) ) {
      stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
    }

    w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
    w.minus <- 1 - w.plus
    imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
    class(imp) <- c("importance", "list")
    return(imp)
  }


##lme
importance.AIClme <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, ...){

  ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##reverse parm
    reversed.parm <- reverse.parm(parm)


    ##extract labels
    mod_formula <- lapply(cand.set, FUN=function(i) labels(summary(i)$coefficients$fixed))

    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"


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
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

    new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                        sort = FALSE)  

    ##add a check to determine if the same number of models include and exlude parameter
    if (length(which(include == 1)) != length(which(include != 1)) ) {
      stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
    }

    w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
    w.minus <- 1 - w.plus
    imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
    class(imp) <- c("importance", "list")
    return(imp)
  }
    


##lmekin
importance.AIClmekin <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, ...){

  ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##reverse parm
    reversed.parm <- reverse.parm(parm)


    ##extract labels
    mod_formula <- lapply(cand.set, FUN = function(i) labels(fixef(i)))

    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"


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
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

    new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                        sort = FALSE)  

    ##add a check to determine if the same number of models include and exlude parameter
    if (length(which(include == 1)) != length(which(include != 1)) ) {
      stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
    }

    w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
    w.minus <- 1 - w.plus
    imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
    class(imp) <- c("importance", "list")
    return(imp)
  }



##maxlike
importance.AICmaxlikeFit.list <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, c.hat = 1, ...){

  ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##reverse parm
    reversed.parm <- reverse.parm(parm)


    ##extract labels
    mod_formula <- lapply(cand.set, FUN = function(i) labels(fixef(i)))


    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"


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
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

    new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                        sort = FALSE)  

    ##add a check to determine if the same number of models include and exlude parameter
    if (length(which(include == 1)) != length(which(include != 1)) ) {
      stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
    }

    w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
    w.minus <- 1 - w.plus
    imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
    class(imp) <- c("importance", "list")
    return(imp)
  }




##mer
importance.AICmer <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, ...){

  ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##reverse parm
    reversed.parm <- reverse.parm(parm)


    ##extract labels
    mod_formula <- lapply(cand.set, FUN=function(i) labels(fixef(i)))

    
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"


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
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

    new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                        sort = FALSE)  

    ##add a check to determine if the same number of models include and exlude parameter
    if (length(which(include == 1)) != length(which(include != 1)) ) {
      stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
    }

    w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
    w.minus <- 1 - w.plus
    imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
    class(imp) <- c("importance", "list")
    return(imp)
  }
      


##lmerModLmerTest
importance.AIClmerModLmerTest <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, ...){
              
  ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##reverse parm
    reversed.parm <- reverse.parm(parm)


    ##extract labels
    mod_formula <- lapply(cand.set, FUN=function(i) labels(fixef(i)))

    
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"


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
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

    new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                        sort = FALSE)  

    ##add a check to determine if the same number of models include and exlude parameter
    if (length(which(include == 1)) != length(which(include != 1)) ) {
      stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
    }

    w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
    w.minus <- 1 - w.plus
    imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
    class(imp) <- c("importance", "list")
    return(imp)
}



##multinom
importance.AICmultinom.nnet <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, c.hat = 1, ...){

  ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##reverse parm
    reversed.parm <- reverse.parm(parm)


    ##extract labels
    mod_formula <- lapply(cand.set, FUN=function(i) colnames(summary(i)$coefficients))

    
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"


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
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

    new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                        sort = FALSE, c.hat = c.hat)  

    ##add a check to determine if the same number of models include and exlude parameter
    if (length(which(include == 1)) != length(which(include != 1)) ) {
      stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
    }

    w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
    w.minus <- 1 - w.plus
    imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
    class(imp) <- c("importance", "list")
    return(imp)
  }


   
##glm.nb
importance.AICnegbin.glm.lm <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, ...){

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##reverse parm
    reversed.parm <- reverse.parm(parm)


    ##extract labels
    mod_formula <- lapply(cand.set, FUN=function(i) rownames(summary(i)$coefficients))


    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"


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
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

    new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                        sort = FALSE)  

    ##add a check to determine if the same number of models include and exlude parameter
    if (length(which(include == 1)) != length(which(include != 1)) ) {
      stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
    }

    w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
    w.minus <- 1 - w.plus
    imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
    class(imp) <- c("importance", "list")
    return(imp)
  }



##nlmer
importance.AICnlmerMod <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, ...){

  ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##reverse parm
    reversed.parm <- reverse.parm(parm)


    ##extract labels
    mod_formula <- lapply(cand.set, FUN = function(i) labels(fixef(i)))

    
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"


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
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

    new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                        sort = FALSE)  

    ##add a check to determine if the same number of models include and exlude parameter
    if (length(which(include == 1)) != length(which(include != 1)) ) {
      stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
    }

    w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
    w.minus <- 1 - w.plus
    imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
    class(imp) <- c("importance", "list")
    return(imp)
  }



##polr
importance.AICpolr <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, ...){

  ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##reverse parm
    reversed.parm <- reverse.parm(parm)


    ##extract labels
    mod_formula <- lapply(cand.set, FUN=function(i) rownames(summary(i)$coefficients))

    
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"


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
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

    new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                        sort = FALSE)  

    ##add a check to determine if the same number of models include and exlude parameter
    if (length(which(include == 1)) != length(which(include != 1)) ) {
      stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
    }

    w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
    w.minus <- 1 - w.plus
    imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
    class(imp) <- c("importance", "list")
    return(imp)
  }
      
      

##rlm
importance.AICrlm.lm <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, ...){

  ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##reverse parm
    reversed.parm <- reverse.parm(parm)


    ##extract labels
    mod_formula <- lapply(cand.set, FUN=function(i) rownames(summary(i)$coefficients))

    
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"


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
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

    new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                        sort = FALSE)  

    ##add a check to determine if the same number of models include and exlude parameter
    if (length(which(include == 1)) != length(which(include != 1)) ) {
      stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
    }

    w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
    w.minus <- 1 - w.plus
    imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
    class(imp) <- c("importance", "list")
    return(imp)
  }



##survreg
importance.AICsurvreg <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, ...){

  ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##reverse parm
    reversed.parm <- reverse.parm(parm)


    ##extract labels
    mod_formula <- lapply(cand.set, FUN=function(i) names(summary(i)$coefficients))


    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"


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
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

    new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                        sort = FALSE)  

    ##add a check to determine if the same number of models include and exlude parameter
    if (length(which(include == 1)) != length(which(include != 1)) ) {
      stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
    }

    w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
    w.minus <- 1 - w.plus
    imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
    class(imp) <- c("importance", "list")
    return(imp)
  }



##colext
importance.AICunmarkedFitColExt <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, c.hat = 1,
                                               parm.type = NULL, ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)

  ##reverse parm
  reversed.parm <- reverse.parm(parm)

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?importance for details\n")}
    
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"
  
  
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
  if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }
  
  new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                      sort = FALSE, c.hat = c.hat)  

  ##add a check to determine if the same number of models include and exlude parameter
  if (length(which(include == 1)) != length(which(include != 1)) ) {
    stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
  }
  
  w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
  w.minus <- 1 - w.plus
  imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
  class(imp) <- c("importance", "list")
  return(imp)
}



##occu
importance.AICunmarkedFitOccu <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, c.hat = 1,
                                             parm.type = NULL, ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)

  ##reverse parm
  reversed.parm <- reverse.parm(parm)

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?importance for details\n")}
    
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"

    
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
  if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }
  
  new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                      sort = FALSE, c.hat = c.hat)  

  ##add a check to determine if the same number of models include and exlude parameter
  if (length(which(include == 1)) != length(which(include != 1)) ) {
    stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
  }
  
  w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
  w.minus <- 1 - w.plus
  imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
  class(imp) <- c("importance", "list")
  return(imp)
}


##occuFP
importance.AICunmarkedFitOccuFP <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, c.hat = 1,
                                            parm.type = NULL, ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)

  ##reverse parm
  reversed.parm <- reverse.parm(parm)

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?importance for details\n")}
    
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"
  
 
    
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
  
  ##false positives - fp
  if(identical(parm.type, "fp")) {
    mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$fp)))
    parm.unmarked <- "fp"
    parm <- paste(parm.unmarked, "(", parm, ")", sep="")
  }
  
   
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
  if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }
  
  new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                      sort = FALSE, c.hat = c.hat)  

  ##add a check to determine if the same number of models include and exlude parameter
  if (length(which(include == 1)) != length(which(include != 1)) ) {
    stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
  }
  
  w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
  w.minus <- 1 - w.plus
  imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
  class(imp) <- c("importance", "list")
  return(imp)
}






##occuRN
importance.AICunmarkedFitOccuRN <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, c.hat = 1,
                                               parm.type = NULL, ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)

  ##reverse parm
  reversed.parm <- reverse.parm(parm)

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?importance for details\n")}
    
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"
  
  
  
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
  if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }
  
  new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                      sort = FALSE, c.hat = c.hat)  

  ##add a check to determine if the same number of models include and exlude parameter
  if (length(which(include == 1)) != length(which(include != 1)) ) {
    stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
  }
  
  w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
  w.minus <- 1 - w.plus
  imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
  class(imp) <- c("importance", "list")
  return(imp)
}



##pcount
importance.AICunmarkedFitPCount <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, c.hat = 1,
                                               parm.type = NULL, ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)

  ##reverse parm
  reversed.parm <- reverse.parm(parm)

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?importance for details\n")}
    
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"
   
  
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
  if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }
  
  new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                      sort = FALSE, c.hat = c.hat)  

  ##add a check to determine if the same number of models include and exlude parameter
  if (length(which(include == 1)) != length(which(include != 1)) ) {
    stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
  }
  
  w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
  w.minus <- 1 - w.plus
  imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
  class(imp) <- c("importance", "list")
  return(imp)
}
    


##pcountOpen
importance.AICunmarkedFitPCO <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, c.hat = 1,
                                               parm.type = NULL, ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)

  ##reverse parm
  reversed.parm <- reverse.parm(parm)

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?importance for details\n")}
    
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"
 
  
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

      ##determine if same H0 on gamma (gamConst, gamAR, gamTrend)
      strip.gam <- sapply(mod_formula, FUN = function(i) unlist(strsplit(i, "\\("))[[1]])
      unique.gam <- unique(strip.gam)
      if(length(unique.gam) > 1) stop("\nDifferent formulations of gamma parameter occur among models:\n
beta estimates cannot be model-averaged\n")
        
      ##create label for parm
      parm.unmarked <- unique.gam
      parm <- paste(unique.gam, "(", parm, ")", sep="")
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
  if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }
  
  new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                      sort = FALSE, c.hat = c.hat)  

  ##add a check to determine if the same number of models include and exlude parameter
  if (length(which(include == 1)) != length(which(include != 1)) ) {
    stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
  }
  
  w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
  w.minus <- 1 - w.plus
  imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
  class(imp) <- c("importance", "list")
  return(imp)
}



##distsamp
importance.AICunmarkedFitDS <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, c.hat = 1,
                                           parm.type = NULL, ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)

  ##reverse parm
  reversed.parm <- reverse.parm(parm)

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?importance for details\n")}
    
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"
  
  
  ##lambda - abundance
  if(identical(parm.type, "lambda")) {
    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$state)))
    parm.unmarked <- "lam"
    parm <- paste(parm.unmarked, "(", parm, ")", sep="")
  }
  
  ##detect
  if(identical(parm.type, "detect")) {
         ##check for key function used
        keyid <- unique(sapply(cand.set, FUN = function(i) i@keyfun))
        if(length(keyid) > 1) stop("\nDifferent key functions used across models:\n",
                                   "cannot compute model-averaged estimate\n")
        if(identical(keyid, "uniform")) stop("\nDetection parameter not found in models\n")
        ##set key prefix used in coef( )
        if(identical(keyid, "halfnorm")) {
            parm.key <- "sigma"
        }
        if(identical(keyid, "hazard")) {
            parm.key <- "shape"
        }
        if(identical(keyid, "exp")) {
            parm.key <- "rate"
        }

      ##label for intercept - label different with this model type
      if(identical(parm, "Int")) {parm <- "(Intercept)"}

      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$det)))
      parm <- paste("p", "(", parm.key, "(", parm, "))", sep="")
      if(!is.null(reversed.parm)) {reversed.parm <- paste("p", "(", parm.key, "(", reversed.parm, "))", sep="")}
      ##stop("\nImportance values for detection covariates not yet supported for unmarkedFitDS class\n")
  }
  

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
  if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }
  
  new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                      sort = FALSE, c.hat = c.hat)  

  ##add a check to determine if the same number of models include and exlude parameter
  if (length(which(include == 1)) != length(which(include != 1)) ) {
    stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
  }
  
  w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
  w.minus <- 1 - w.plus
  imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
  class(imp) <- c("importance", "list")
  return(imp)
}

    

##gdistsamp
importance.AICunmarkedFitGDS <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, c.hat = 1,
                                            parm.type = NULL, ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)

  ##reverse parm
  reversed.parm <- reverse.parm(parm)

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?importance for details\n")}
    
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"
      
  ##lambda - abundance
  if(identical(parm.type, "lambda")) {
    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$state)))
    parm.unmarked <- "lam"
    parm <- paste(parm.unmarked, "(", parm, ")", sep="")
  }
  
  ##detect
  if(identical(parm.type, "detect")) {
    stop("\nImportance values for detection covariates not yet supported for unmarkedFitGDS class\n")
  }
  
  ##availability
  if(identical(parm.type, "phi")) {
    stop("\nImportance values for availability covariates not yet supported for unmarkedFitGDS class\n")
  }
  
  
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
  if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }
  
  new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                      sort = FALSE, c.hat = c.hat)  

  ##add a check to determine if the same number of models include and exlude parameter
  if (length(which(include == 1)) != length(which(include != 1)) ) {
    stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
  }
  
  w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
  w.minus <- 1 - w.plus
  imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
  class(imp) <- c("importance", "list")
  return(imp)
}


##multinomPois
importance.AICunmarkedFitMPois <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, c.hat = 1,
                                              parm.type = NULL, ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)

  ##reverse parm
  reversed.parm <- reverse.parm(parm)

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?importance for details\n")}
    
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"
  
  
  ##lambda - abundance
  if(identical(parm.type, "lambda")) {
    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$state)))
    parm.unmarked <- "lambda"
    parm <- paste(parm.unmarked, "(", parm, ")", sep="")
  }

  ##detect
  if(identical(parm.type, "detect")) {
    mod_formula<-lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$det)))
    parm.unmarked <- "p"
    parm <- paste(parm.unmarked, "(", parm, ")", sep="")
  }
  
    
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
  if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }
  
  new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                      sort = FALSE, c.hat = c.hat)  

  ##add a check to determine if the same number of models include and exlude parameter
  if (length(which(include == 1)) != length(which(include != 1)) ) {
    stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
  }
  
  w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
  w.minus <- 1 - w.plus
  imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
  class(imp) <- c("importance", "list")
  return(imp)
}



##gmultmix
importance.AICunmarkedFitGMM <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, c.hat = 1,
                                            parm.type = NULL, ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)

  ##reverse parm
  reversed.parm <- reverse.parm(parm)

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?importance for details\n")}
    
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"

    
  ##lambda - abundance
  if(identical(parm.type, "lambda")) {
    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$lambda)))
    parm.unmarked <- "lambda"
    parm <- paste(parm.unmarked, "(", parm, ")", sep="")
  }
  
  ##detect
  if(identical(parm.type, "detect")) {
    mod_formula<-lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$det)))
    parm.unmarked <- "p"
    parm <- paste(parm.unmarked, "(", parm, ")", sep="")
  }

  ##availability
  if(identical(parm.type, "phi")) {
    mod_formula<-lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$phi)))
    parm.unmarked <- "phi"
    parm <- paste(parm.unmarked, "(", parm, ")", sep="")
  }
   
  
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
  if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }
  
  new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                      sort = FALSE, c.hat = c.hat)  

  ##add a check to determine if the same number of models include and exlude parameter
  if (length(which(include == 1)) != length(which(include != 1)) ) {
    stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
  }
  
  w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
  w.minus <- 1 - w.plus
  imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
  class(imp) <- c("importance", "list")
  return(imp)
}


##gmultmix
importance.AICunmarkedFitGPC <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, c.hat = 1,
                                         parm.type = NULL, ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)

  ##reverse parm
  reversed.parm <- reverse.parm(parm)

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?importance for details\n")}
    
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"

  ##lambda - abundance
  if(identical(parm.type, "lambda")) {
    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$lambda)))
    parm.unmarked <- "lambda"
    parm <- paste(parm.unmarked, "(", parm, ")", sep="")
  }
  
  ##detect
  if(identical(parm.type, "detect")) {
    mod_formula<-lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$det)))
    parm.unmarked <- "p"
    parm <- paste(parm.unmarked, "(", parm, ")", sep="")
  }

  ##availability
  if(identical(parm.type, "phi")) {
    mod_formula<-lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$phi)))
    parm.unmarked <- "phi"
    parm <- paste(parm.unmarked, "(", parm, ")", sep="")
  }
  
  
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
  if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }
  
  new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                      sort = FALSE, c.hat = c.hat)  

  ##add a check to determine if the same number of models include and exlude parameter
  if (length(which(include == 1)) != length(which(include != 1)) ) {
    stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
  }
  
  w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
  w.minus <- 1 - w.plus
  imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
  class(imp) <- c("importance", "list")
  return(imp)
}


##occuMulti
importance.AICunmarkedFitOccuMulti <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, c.hat = 1,
                                               parm.type = NULL, ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  ##remove all leading and trailing white space and within parm
  ##parm <- gsub('[[:space:]]+', "", parm)

  ##reverse parm
  reversed.parm <- reverse.parm(parm)

    
    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?importance for details\n")}
    
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"
  
  
  ##psi
  if(identical(parm.type, "psi")) {
    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN = function(x) labels(coef(x@estimates@estimates$state)))
    parm.unmarked <- "psi"
    parm <- paste(parm.unmarked, "(", parm, ")", sep="")
  }

  ##detect
  if(identical(parm.type, "detect")) {
    mod_formula <- lapply(cand.set, FUN = function(x) labels(coef(x@estimates@estimates$det)))
    parm.unmarked <- "p"
    parm <- paste(parm.unmarked, "(", parm, ")", sep="")
  }
  
 
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
  if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }
  
  new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                      sort = FALSE, c.hat = c.hat)  

  ##add a check to determine if the same number of models include and exlude parameter
  if (length(which(include == 1)) != length(which(include != 1)) ) {
    stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
  }
  
  w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
  w.minus <- 1 - w.plus
  imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
  class(imp) <- c("importance", "list")
  return(imp)
}


##occuMS
importance.AICunmarkedFitOccuMS <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, c.hat = 1,
                                            parm.type = NULL, ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  ##remove all leading and trailing white space and within parm
  ##parm <- gsub('[[:space:]]+', "", parm)

  ##reverse parm
  reversed.parm <- reverse.parm(parm)

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?importance for details\n")}
    
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"
  
  
  ##psi
  if(identical(parm.type, "psi")) {
    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN = function(x) labels(coef(x@estimates@estimates$state)))
    parm.unmarked <- "psi"
    parm <- paste(parm.unmarked, "(", parm, ")", sep="")
  }

  ##detect
  if(identical(parm.type, "detect")) {
    mod_formula <- lapply(cand.set, FUN = function(x) labels(coef(x@estimates@estimates$det)))
    parm.unmarked <- "p"
    parm <- paste(parm.unmarked, "(", parm, ")", sep="")
  }

    ##transition
    if(identical(parm.type, "phi")) {
        ##check that parameter appears in all models
        nseasons <- unique(sapply(cand.set, FUN = function(i) i@data@numPrimary))
        if(nseasons == 1) {
            stop("\nParameter \'phi\' does not appear in single-season models\n")
        }
        mod_formula <- lapply(cand.set, FUN = function(x) labels(coef(x@estimates@estimates$transition)))
        parm.unmarked <- "phi"
        parm <- paste(parm.unmarked, "(", parm, ")", sep="")
    }

    
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
  if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }
  
  new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                      sort = FALSE, c.hat = c.hat)  

  ##add a check to determine if the same number of models include and exlude parameter
  if (length(which(include == 1)) != length(which(include != 1)) ) {
    stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
  }
  
  w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
  w.minus <- 1 - w.plus
  imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
  class(imp) <- c("importance", "list")
  return(imp)
}


##occuTTD
importance.AICunmarkedFitOccuTTD <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, c.hat = 1,
                                             parm.type = NULL, ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)

  ##reverse parm
  reversed.parm <- reverse.parm(parm)

   
    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?importance for details\n")}
    
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"

    
  ##psi - initial occupancy
  if(identical(parm.type, "psi")) {
    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$psi)))
    ##create label for parm
    parm.unmarked <- "psi"
    parm <- paste(parm.unmarked, "(", parm, ")", sep="")
  }

  ##gamma - colonization
  if(identical(parm.type, "gamma")) {
      nseasons <- unique(sapply(cand.set, FUN = function(i) i@data@numPrimary))
      if(nseasons == 1) {
          stop("\nParameter \'gamma\' does not appear in single-season models\n")
      }
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$col)))
      ##create label for parm
      parm.unmarked <- "col"
      parm <- paste(parm.unmarked, "(", parm, ")", sep="")
  }

  ##epsilon - extinction
  if(identical(parm.type, "epsilon")) {
      nseasons <- unique(sapply(cand.set, FUN = function(i) i@data@numPrimary))
      if(nseasons == 1) {
          stop("\nParameter \'epsilon\' does not appear in single-season models\n")
      }
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$ext)))
      ##create label for parm
      parm.unmarked <- "ext"
      parm <- paste(parm.unmarked, "(", parm, ")", sep="")
  }

  ##detect - lambda parameter is a rate of a species not detected in t to be detected at next time step
  if(identical(parm.type, "detect")) {
    mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$det)))
    parm.unmarked <- "lam"
    parm <- paste(parm.unmarked, "(", parm, ")", sep="")
  }
  
  
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
  if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }
  
  new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                      sort = FALSE, c.hat = c.hat)  

  ##add a check to determine if the same number of models include and exlude parameter
  if (length(which(include == 1)) != length(which(include != 1)) ) {
    stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
  }
  
  w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
  w.minus <- 1 - w.plus
  imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
  class(imp) <- c("importance", "list")
  return(imp)
}


##multmixOpen
importance.AICunmarkedFitMMO <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, c.hat = 1,
                                         parm.type = NULL, ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)

  ##reverse parm
  reversed.parm <- reverse.parm(parm)

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?importance for details\n")}
    
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"
      
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

        ##determine if same H0 on gamma (gamConst, gamAR, gamTrend)
        strip.gam <- sapply(mod_formula, FUN = function(i) unlist(strsplit(i, "\\("))[[1]])
        unique.gam <- unique(strip.gam)
        if(length(unique.gam) > 1) stop("\nDifferent formulations of gamma parameter occur among models:\n
beta estimates cannot be model-averaged\n")
        ##create label for parm
        parm.unmarked <- unique.gam
        parm <- paste(unique.gam, "(", parm, ")", sep="")
    }
  
  ##omega - apparent survival
  if(identical(parm.type, "omega")) {
    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$omega)))
    ##create label for parm
    parm.unmarked <- "omega"
    parm <- paste(parm.unmarked, "(", parm, ")", sep="")
  }
  
    ##iota (for immigration = TRUE with dynamics = "autoreg", "trend", "ricker", or "gompertz")
    if(identical(parm.type, "iota")) {
        not.include <- lapply(cand.set, FUN = function(i) i@formlist$iotaformula)
        ##check that parameter appears in all models
        parfreq <- sum(sapply(cand.set, FUN = function(i) any(names(i@estimates@estimates) == parm.type)))
        if(!identical(length(cand.set), parfreq)) {
            stop("\nParameter \'iota\' (parm.type = \"", parm.type, "\") does not appear in all models:",
                 "\ncannot compute model-averaged estimate across all models\n")
        }
        ##extract model formula for each model in cand.set
        mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$iota)))
        ##create label for parm
        parm.unmarked <- "iota"
        parm <- paste(parm.unmarked, "(", parm, ")", sep="")

    }

    ##detect
    if(identical(parm.type, "detect")) {
        mod_formula<-lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$det)))
        parm.unmarked <- "p"
        parm <- paste(parm.unmarked, "(", parm, ")", sep="")
    }
  

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
  if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }
  
  new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                      sort = FALSE, c.hat = c.hat)  

  ##add a check to determine if the same number of models include and exlude parameter
  if (length(which(include == 1)) != length(which(include != 1)) ) {
    stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
  }
  
  w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
  w.minus <- 1 - w.plus
  imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
  class(imp) <- c("importance", "list")
  return(imp)
}


##distsampOpen
importance.AICunmarkedFitDSO <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, c.hat = 1,
                                         parm.type = NULL, ...){

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
        if(is.null(names(cand.set))) {
            modnames <- paste("Mod", 1:length(cand.set), sep = "")
            warning("\nModel names have been supplied automatically in the table\n")
        } else {
            modnames <- names(cand.set)
        }
    }
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##reverse parm
    reversed.parm <- reverse.parm(parm)

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?importance for details\n")}
    
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"

  
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
        ##determine if same H0 on gamma (gamConst, gamAR, gamTrend)
        strip.gam <- sapply(mod_formula, FUN = function(i) unlist(strsplit(i, "\\("))[[1]])
        unique.gam <- unique(strip.gam)
        if(length(unique.gam) > 1) stop("\nDifferent formulations of gamma parameter occur among models:\n
beta estimates cannot be model-averaged\n")
        
        ##create label for parm
        parm.unmarked <- unique.gam
        parm <- paste(unique.gam, "(", parm, ")", sep="")
    }
  
  ##omega - apparent survival
  if(identical(parm.type, "omega")) {
    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$omega)))
    ##create label for parm
    parm.unmarked <- "omega"
    parm <- paste(parm.unmarked, "(", parm, ")", sep="")
  }
  
    ##iota (for immigration = TRUE with dynamics = "autoreg", "trend", "ricker", or "gompertz")
    if(identical(parm.type, "iota")) {
        not.include <- lapply(cand.set, FUN = function(i) i@formlist$iotaformula)
        ##check that parameter appears in all models
        parfreq <- sum(sapply(cand.set, FUN = function(i) any(names(i@estimates@estimates) == parm.type)))
        if(!identical(length(cand.set), parfreq)) {
            stop("\nParameter \'iota\' (parm.type = \"", parm.type, "\") does not appear in all models:",
                 "\ncannot compute model-averaged estimate across all models\n")
        }
        ##extract model formula for each model in cand.set
        mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$iota)))
        ##create label for parm
        parm.unmarked <- "iota"
        parm <- paste(parm.unmarked, "(", parm, ")", sep="")

    }
    
    ##detect
    if(identical(parm.type, "detect")) {
        mod_formula<-lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$det)))
        parm.unmarked <- "sigma"
        parm <- paste(parm.unmarked, "(", parm, ")", sep="")
    }
  
  
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
  if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }
  
  new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                      sort = FALSE, c.hat = c.hat)  

  ##add a check to determine if the same number of models include and exlude parameter
  if (length(which(include == 1)) != length(which(include != 1)) ) {
    stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
  }
  
  w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
  w.minus <- 1 - w.plus
  imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
  class(imp) <- c("importance", "list")
  return(imp)
}


##goccu
importance.AICunmarkedFitGOccu <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, c.hat = 1,
                                           parm.type = NULL, ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)

  ##reverse parm
  reversed.parm <- reverse.parm(parm)

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?importance for details\n")}
    
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"

  ##psi - occupancy
  if(identical(parm.type, "psi")) {
    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$psi)))
    parm.unmarked <- "psi"
    parm <- paste(parm.unmarked, "(", parm, ")", sep="")
  }
  
  ##detect
  if(identical(parm.type, "detect")) {
    mod_formula<-lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$det)))
    parm.unmarked <- "p"
    parm <- paste(parm.unmarked, "(", parm, ")", sep="")
  }

  ##availability
  if(identical(parm.type, "phi")) {
    mod_formula<-lapply(cand.set, FUN = function(i) labels(coef(i@estimates@estimates$phi)))
    parm.unmarked <- "phi"
    parm <- paste(parm.unmarked, "(", parm, ")", sep="")
  }
  
  
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
  if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }
  
  new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                      sort = FALSE, c.hat = c.hat)  

  ##add a check to determine if the same number of models include and exlude parameter
  if (length(which(include == 1)) != length(which(include != 1)) ) {
    stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
  }
  
  w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
  w.minus <- 1 - w.plus
  imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
  class(imp) <- c("importance", "list")
  return(imp)
}



##occuComm
importance.AICunmarkedFitOccuComm <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, c.hat = 1,
                                               parm.type = NULL, ...){

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
  ##remove all leading and trailing white space and within parm
  ##parm <- gsub('[[:space:]]+', "", parm)

  ##reverse parm
  reversed.parm <- reverse.parm(parm)

    
    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?importance for details\n")}
    
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"
  
  
  ##psi
  if(identical(parm.type, "psi")) {
    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN = function(x) labels(coef(x@estimates@estimates$state)))
    parm.unmarked <- "psi"
    parm <- paste(parm.unmarked, "(", parm, ")", sep="")
  }

  ##detect
  if(identical(parm.type, "detect")) {
    mod_formula <- lapply(cand.set, FUN = function(x) labels(coef(x@estimates@estimates$det)))
    parm.unmarked <- "p"
    parm <- paste(parm.unmarked, "(", parm, ")", sep="")
  }
  
 
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
  if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }
  
  new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                      sort = FALSE, c.hat = c.hat)  

  ##add a check to determine if the same number of models include and exlude parameter
  if (length(which(include == 1)) != length(which(include != 1)) ) {
    stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
  }
  
  w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
  w.minus <- 1 - w.plus
  imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
  class(imp) <- c("importance", "list")
  return(imp)
}



##vglm
importance.AICvglm <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, c.hat = 1, ...){

  ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##reverse parm
    reversed.parm <- reverse.parm(parm)


    ##extract labels
    mod_formula <- lapply(cand.set, FUN=function(i) labels(coefficients(i)))


    ##check whether parm is involved in interaction or if label changes for some models - e.g., ZIP models
    ##if : not already included
    if(regexpr(":", parm, fixed = TRUE) == -1) {

      ##if : not included
      parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
      inter.check <- ifelse(attr(regexpr(parm.inter[1], mod_formula, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                       mod_formula, fixed = TRUE), "match.length") == "-1", 0, 1)
      if(sum(inter.check) > 0) warning("\nLabel of parameter of interest seems to change across models:\n",
                                       "check model syntax for possible problems\n")
    } else {
      ##if : already included
      ##remove : from parm
      simple.parm <- unlist(strsplit(parm, split = ":"))[1]
      
      ##search for simple.parm and parm in model formulae
      no.colon <- sum(ifelse(attr(regexpr(simple.parm, mod_formula, fixed = TRUE), "match.length") != "-1", 1, 0))
      with.colon <- sum(ifelse(attr(regexpr(parm, mod_formula, fixed = TRUE), "match.length") != "-1", 0, 1))
        
      ##check if both are > 0
      if(no.colon > 0 && with.colon > 0) warning("\nLabel of parameter of interest seems to change across models:\n",
                                                 "check model syntax for possible problems\n")
    }
    
    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"


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
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

    new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                        sort = FALSE, c.hat = c.hat)  

    ##add a check to determine if the same number of models include and exlude parameter
    if (length(which(include == 1)) != length(which(include != 1)) ) {
      stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
    }

    w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
    w.minus <- 1 - w.plus
    imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
    class(imp) <- c("importance", "list")
    return(imp)
  }

      

##zeroinfl
importance.AICzeroinfl <- function(cand.set, parm, modnames = NULL, second.ord = TRUE, nobs = NULL, ...){
  
  ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }
    
    
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##reverse parm
    reversed.parm <- reverse.parm(parm)


    ##extract labels
    mod_formula <- lapply(cand.set, FUN=function(i) labels(coef(i)))


    ##if (Intercept) is chosen assign (Int) - for compatibility
    if(identical(parm, "(Intercept)")) parm <- "Int"


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
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

    new_table <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs,
                        sort = FALSE)  

    ##add a check to determine if the same number of models include and exlude parameter
    if (length(which(include == 1)) != length(which(include != 1)) ) {
      stop("\nImportance values are only meaningful when the number of models with and without parameter are equal\n")
    }

    w.plus <- sum(new_table[which(include == 1), 6]) #select models including a given parameter
    w.minus <- 1 - w.plus
    imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
    class(imp) <- c("importance", "list")
    return(imp)
  }

    

##function for nicer printing of importance values
print.importance <- function(x, digits = 2, ...) {
  cat("\nImportance values of '", x$parm, "':\n\n", sep = "")
  cat("w+ (models including parameter):", round(x$w.plus, digits = digits), "\n")
  cat("w- (models excluding parameter):", round(x$w.minus, digits = digits), "\n")
  cat("\n")
}

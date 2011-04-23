importance <-
  function(cand.set, parm, modnames, c.hat=1, second.ord=TRUE, nobs=NULL){
    ##check if class is appropriate
    ##extract classes
    mod.class <- unlist(lapply(X=cand.set, FUN=class))
  
    ##check if all are identical
    check.class <- unique(mod.class)
  
    if(identical(check.class, "lm") || identical(check.class, c("glm", "lm")))  {
    
      ##extract model formula for each model in cand.set
      mod_formula <- lapply(cand.set, FUN=function(i) rownames(summary(i)$coefficients))
    }
  
    if(identical(check.class, "lme")) {
      mod_formula <- lapply(cand.set, FUN=function(i) labels(summary(i)$coefficients$fixed))
    }

    if(identical(check.class, "gls")) {
      mod_formula <- lapply(cand.set, FUN=function(i) labels(summary(i)$coefficients))
    }

    if(identical(check.class, c("multinom", "nnet"))) {
      mod_formula <- lapply(cand.set, FUN=function(i) colnames(summary(i)$coefficients))
    }

    if(identical(check.class, "mer")) {
      mod_formula <- lapply(cand.set, FUN=function(i) rownames(summary(i)@coefs))
    }

    ##setup matrix to indicate presence of parm in the model
    include <- matrix(NA, nrow=length(cand.set), ncol=1)

    ##iterate over each formula in mod_formula list
    for (i in 1:length(cand.set)) {
      idents <- NULL
      form <- mod_formula[[i]]

      ##iterate over each element of formula[[i]] in list
      for (j in 1:length(form)) {
        idents[j] <- identical(paste(parm), form[j])
      }
      include[i] <- ifelse(any(idents==1), 1, 0)
    }

    ##add a check to determine if include always == 0
    if (sum(include)==0) {stop("Parameter not found in any of the candidate models") }

    new_table <- aictab(cand.set=cand.set, modnames=modnames, sort=FALSE, c.hat=c.hat, second.ord=second.ord, nobs=nobs)  

    ##add a check to determine if the same number of models include and exlude parameter
    if (length(which(include==1)) != length(which(include!=1)) ) {
      stop("Importance values are only meaningful when the number of models with and without parameter are equal")
    }

    w.plus <- sum(new_table[which(include==1), 6]) #select models including a given parameter
    w.minus <- 1 - w.plus
    imp <- list("parm" = parm, "w.plus" = w.plus, "w.minus" = w.minus)

  
    class(imp) <- c("importance", "list")
    return(imp)
  }

##function for nicer printing of importance values
print.importance <- function(x, digits = 2, ...) {
  cat("\nImportance values of '", x$parm, "' :\n\n")
  cat("w+ (models including parameter): ", round(x$w.plus, digits=digits), "\n")
  cat("w- (models excluding parameter): ", round(x$w.minus, digits=digits), "\n")
  cat("\n")
}


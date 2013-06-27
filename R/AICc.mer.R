AICc.mer <-
  function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL){
    ##check whether object is of appropriate class
    if(identical(paste(class(mod)), "mer")) {  #class(mod) yields different attributes
    
      if(is.null(nobs)) {
        n <- mod@dims["n"]
      } else {n <- nobs}
      
      LL <- logLik(mod)[1]
      K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model
      
      if(second.ord == TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))} else{AICc <- -2*LL+2*K}
                    
      if(return.K == TRUE) AICc[1] <- K #attributes the first element of AICc to K
      AICc

    } else {stop("This function is only appropriate with objects of \'mer\' class")}
  }



AICc.merMod <-
  function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL){
    ##check whether object is of appropriate class
    if(identical(paste(class(mod)), "lmerMod") || identical(paste(class(mod)), "glmerMod") || identical(paste(class(mod)), "nlmerMod")) {  #class(mod) yields different attributes
    
      if(is.null(nobs)) {
        n <- mod@devcomp$dims["n"]
        names(n) <- NULL
      } else {n <- nobs}
      
      LL <- logLik(mod)[1]
      K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model
      
      if(second.ord == TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))} else{AICc <- -2*LL+2*K}
                    
      if(return.K == TRUE) AICc[1] <- K #attributes the first element of AICc to K
      AICc

    } else {stop("This function is only appropriate with objects of \'mer\' class")}
  }

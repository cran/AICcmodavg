AICc.mer <-
  function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL){
    ##check whether object is of appropriate class
    if(identical(paste(class(mod)), "mer")) {  #class(mod) yields different attributes
    
      if(identical(nobs, NULL)) {
        ##add check for binomial family to sum weights
        if(identical(fam.link.mer(mod)$family, "binomial")) {
          success.trials <- NULL
          success.failure <- NULL

          ##for events/trials syntax
          if(!is.na(charmatch(x = "(weights)", table =  names(mod@frame)))){
            success.trials <- TRUE
            n <- sum(mod@frame$"(weights)")
          }

          ##for cbind(incidence, size - incidence) syntax 
          if(!is.na(charmatch(x = "cbind(incidence, size - incidence)", table = names(mod@frame)))){
            success.failure <- TRUE
            n <- sum(mod@frame$cbind)            
          }

          ##for binary case syntax
          if(is.null(success.trials) && is.null(success.failure)) {
            n <- mod@dims["n"]
          }

          ##for all other cases
        } else {          
          n <- mod@dims["n"]
        }
      } else {n <- nobs}
      
      LL <- logLik(mod)[1]
      K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model
      
      if(second.ord == TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))} else{AICc <- -2*LL+2*K}
                    
      if(return.K == TRUE) AICc[1] <- K #attributes the first element of AICc to K
      AICc

    } else {stop("This function is only appropriate with objects of \'mer\' class")}
  }

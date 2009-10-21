aictab.lme <-
  function(cand.set,modnames, sort=TRUE, second.ord=TRUE, nobs=NULL){  #specify whether table should be sorted or not by delta AICc
    Results <- NULL
#check if models were fit with same method (REML or ML)
    check_ML<-unlist(lapply(cand.set, FUN=function(i) i$method))
    
    if (any(check_ML!="ML")) {
      warning(paste("Model selection for fixed effects is only appropriate with method=ML:", "\n",
                    "REML (default) should only be used to select random effects for a constant set of fixed effects"))
    }
    
    check.method <- unique(check_ML)
    
    if(identical(check.method, c("ML", "REML"))) {
      stop("You should not have models fit with REML and ML in the same candidate model set")
    }
    
    Results<-data.frame(Modnames=modnames)                    #assign model names to first column
    Results$K<-unlist(lapply(cand.set, AICc.lme, return.K=TRUE, second.ord=second.ord, nobs=nobs))     #extract number of parameters
    Results$AICc<-unlist(lapply(cand.set, AICc.lme, return.K=FALSE, second.ord=second.ord, nobs=nobs))  #extract AICc                                      #
    Results$Delta_AICc<-Results$AICc-min(Results$AICc)            #compute delta AICc
    Results$ModelLik<-exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
    Results$AICcWt<-Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights
    if(sort)  {Results<-Results[rev(order(Results$AICcWt)),]}  	  #if sort=TRUE, models are ranked based on delta AICc
    Results$Cum.Wt<-cumsum(Results$AICcWt)                        #display cumulative sum Akaike weights

    #rename correctly to AIC
    if(second.ord==FALSE) {
      colnames(Results)<-c("Modnames", "K", "AIC", "Delta AIC", "ModelLik", "AICWt", "Cum.Wt")
    }
    
    class(Results) <- c("aictab", "data.frame")
    return(Results)
  }


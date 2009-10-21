aictab.polr <-
  function(cand.set, modnames, sort=TRUE, second.ord=TRUE, nobs=NULL){  #specify whether table should be sorted or not by delta AICc
    Results <- NULL
    Results<-data.frame(Modnames=modnames)                    #assign model names to first column
    Results$K<-unlist(lapply(X=cand.set, FUN=AICc.polr, return.K=TRUE, 
                             second.ord=second.ord, nobs=nobs))     #extract number of parameters
    Results$AICc<-unlist(lapply(X=cand.set, FUN=AICc.polr, return.K=FALSE, 
                                second.ord=second.ord, nobs=nobs))  #extract AICc                                      #
    Results$Delta_AICc<-Results$AICc-min(Results$AICc)            #compute delta AICc
    Results$ModelLik<-exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
    Results$AICcWt<-Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights
    if(sort)  Results<-Results[rev(order(Results$AICcWt)),] 	  #if sort=TRUE, models are ranked based on delta AICc
    Results$Cum.Wt<-cumsum(Results$AICcWt)                        #display cumulative sum Akaike weights

    #rename correctly to AIC
    if(second.ord==FALSE) {
      colnames(Results)<-c("Modnames", "K", "AIC", "Delta AIC", "ModelLik",
                           "AICWt", "Cum.Wt")
    }

    class(Results) <- c("aictab", "data.frame")
    return(Results)
  }


aictab.polr <-
  function(cand.set, modnames, sort=TRUE, second.ord=TRUE, nobs=NULL){  #specify whether table should be sorted or not by delta AICc

    ##add check to see whether response variable is the same for all models
    check.resp <- lapply(X = cand.set, FUN = function(b) formula(b)[2])
    if(length(unique(check.resp)) > 1) stop("\nYou must use the same response variable for all models\n")
    
    Results <- NULL
    Results<-data.frame(Modnames=modnames)                    #assign model names to first column
    Results$K<-unlist(lapply(X=cand.set, FUN=AICc.polr, return.K=TRUE, 
                             second.ord=second.ord, nobs=nobs))     #extract number of parameters
    Results$AICc<-unlist(lapply(X=cand.set, FUN=AICc.polr, return.K=FALSE, 
                                second.ord=second.ord, nobs=nobs))  #extract AICc                                      #
    Results$Delta_AICc<-Results$AICc-min(Results$AICc)            #compute delta AICc
    Results$ModelLik<-exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
    Results$AICcWt<-Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights
    Results$LL <- unlist(lapply(X=cand.set, FUN=function(i) logLik(i)[1]))      

    ##check if some models are redundant
    if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")
    
    #rename correctly to AIC
    if(second.ord==FALSE) {
      colnames(Results)[1:6]<-c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik",
                           "AICWt")
    }


    if(sort)  {
      Results<-Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
      Results$Cum.Wt<-cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
    } else {Results$Cum.Wt <- NULL}


    
    class(Results) <- c("aictab", "data.frame")
    return(Results)
  }


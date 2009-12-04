aictab.glm <-
  function(cand.set, modnames, sort=TRUE, c.hat=1, second.ord=TRUE, nobs=NULL){  #specify whether table should be sorted or not by delta AICc
    Results <- NULL
    Results<-data.frame(Modnames=modnames)                    #assign model names to first column
    Results$K<-unlist(lapply(X=cand.set, FUN=AICc.glm, return.K=TRUE, c.hat=c.hat,
                             second.ord=second.ord, nobs=nobs))     #extract number of parameters
    Results$AICc<-unlist(lapply(X=cand.set, FUN=AICc.glm, return.K=FALSE, c.hat=c.hat,
                                second.ord=second.ord, nobs=nobs))  #extract AICc                                      #
    Results$Delta_AICc<-Results$AICc-min(Results$AICc)            #compute delta AICc
    Results$ModelLik<-exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
    Results$AICcWt<-Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights
    
         
    ##check if AICc and c.hat = 1
    if(second.ord==TRUE && c.hat == 1) {
      Results$LL <- unlist(lapply(X=cand.set, FUN=function(i) logLik(i)[1]))
    }
    
    ##rename correctly to QAICc and add column for c-hat
    if(second.ord==TRUE && c.hat > 1) {
      colnames(Results) <- c("Modnames", "K", "QAICc", "Delta QAICc", "ModelLik", "QAICcWt")
      LL <- unlist(lapply(X=cand.set, FUN=function(i) logLik(i)[1]))
      Results$Quasi.LL <- LL/c.hat
      Results$c_hat <- c.hat
    }      

    ##rename correctly to AIC
    if(second.ord==FALSE && c.hat==1) {
      colnames(Results)<-c("Modnames", "K", "AIC", "Delta AIC", "ModelLik", "AICWt")
      Results$LL <- unlist(lapply(X=cand.set, FUN=function(i) logLik(i)[1]))      
    }  

    ##rename correctly to QAIC and add column for c-hat
    if(second.ord==FALSE && c.hat > 1) {
      colnames(Results)<-c("Modnames", "K", "QAIC", "Delta QAIC", "ModelLik", "QAICWt")
      LL <- unlist(lapply(X=cand.set, FUN=function(i) logLik(i)[1]))
      Results$Quasi.LL <- LL/c.hat
      Results$c_hat<-c.hat
    }      

    
    if(sort)  {
      Results<-Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
      Results$Cum.Wt<-cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
    } else {Results$Cum.Wt <- NULL}

    
    class(Results) <- c("aictab", "data.frame")
    return(Results)
  }

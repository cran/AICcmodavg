AICc.glm <-
  function(mod, return.K=FALSE, c.hat=1, second.ord=TRUE, nobs=NULL){
#check whether object is of appropriate class
    if(identical(class(mod), "lm") || identical(class(mod), c("glm", "lm"))) {
    
      if(identical(nobs, NULL)) {
        if(identical(family(mod)[[1]], "binomial") && any(mod$prior.weights!=1)) {
          n <- sum(mod$prior.weights) } else {n<-length(mod$fitted)}
      } else {n <- nobs}
      
      LL<-logLik(mod)[1]
      K<-attr(logLik(mod), "df")  #extract correct number of parameters included in model - this includes LM
      
      if(c.hat == 1) {
        if(second.ord==TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))} else{AICc <- -2*LL+2*K}
                    }
      if(c.hat > 1 && c.hat <= 4) {K<-K+1; if(second.ord==TRUE) {AICc <- (-2*LL/c.hat)+2*K*(n/(n-K-1))
      #adjust parameter count to include estimation of dispersion parameter
                                            } else{AICc <- (-2*LL/c.hat)+2*K}
                                        }
      if(c.hat > 4) stop("High overdispersion and model fit is questionable")

      #check if negative binomial and add 1 to K for estimation of theta if glm( ) was used
      if(!is.na(charmatch(x="Negative Binomial", table=family(mod)$family))) {
        if(is.null(mod$theta)) {K <- K+1;
                                if(second.ord==TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))
                                                    } else {AICc <- -2*LL+2*K}
                                                                                                    }
        if(c.hat!=1) stop("You should not use the c.hat argument with the negative binomial")
      }
#add 1 for theta parameter in negative binomial fit glm( )
#check if gamma and add 1 to K for estimation of shape parameter if glm( ) was used
#an extra condition must be added to avoid adding a parameter for theta with negative binomial when glm.nb( ) is fit which estimates the correct number of parameters
      if(return.K==TRUE) AICc[1]<-K #attributes the first element of AICc to K
      AICc



                   
    } else {stop("This function is only appropriate with either \'lm\' or \'glm\' classes")}
  
  }


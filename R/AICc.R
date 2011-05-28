AICc <-
  function(mod, return.K = FALSE, c.hat = 1, second.ord = TRUE, nobs = NULL) {
    aicc <- NULL

    ##determine if lm or glm
    if(identical(class(mod)[1], "lm") || identical(class(mod)[1], "glm")) {
      aicc <- AICc.glm(mod = mod, return.K = return.K, c.hat = c.hat, second.ord = second.ord, 
                       nobs = nobs)
    }   

    ##determine if multinom
    if(identical(class(mod)[1], "multinom")) {
      aicc <- AICc.mult(mod = mod, return.K = return.K, c.hat = c.hat, second.ord = second.ord, 
                        nobs = nobs)
    }

    ##determine if polr
    if(identical(class(mod), "polr")) {
      aicc <- AICc.polr(mod = mod, return.K = return.K, second.ord = second.ord,
                        nobs = nobs)
    }

    ##determine if lme
    if(identical(class(mod)[1], "lme")) {
      aicc <- AICc.lme(mod = mod, return.K = return.K, second.ord = second.ord,
                       nobs = nobs)
    }      

    ##determine if gls
    if(identical(class(mod)[1], "gls")) {
      aicc <- AICc.gls(mod = mod, return.K = return.K, second.ord = second.ord,
                       nobs = nobs)
    }

    ##determine if mer
    if(identical(class(mod)[1], "mer")) {
      aicc <- AICc.mer(mod = mod, return.K = return.K, second.ord = second.ord,
                       nobs = nobs)
    }


    ##determine if unmarked
    unmarked.class <- c("unmarkedFitOccu", "unmarkedFitColExt", "unmarkedFitOccuRN", "unmarkedFitPCount", "unmarkedFitPCO")
    if(any(sapply(unmarked.class, FUN = function(i) identical(i, class(mod))))) {
      aicc <- AICc.unmarked(mod = mod, return.K = return.K, c.hat = c.hat, second.ord = second.ord,
                            nobs = nobs)
    }
        
    

    ##if(class(mod)[1]=="nlm") {aicc <- AICc.nlm(mod, return.K)}      #determine if object from nlm optimizer
    return(aicc)
  }

##compute BIC
##generic
useBIC <- function(mod, return.K = FALSE, nobs = NULL, ...){
  UseMethod("useBIC", mod)
}

useBIC.default <- function(mod, return.K = FALSE, nobs = NULL, ...){
  stop("\nFunction not yet defined for this object class\n")
}



##aov objects
useBIC.aov <-
  function(mod, return.K = FALSE, nobs = NULL, ...){
    
    if(identical(nobs, NULL)) {n <- length(mod$fitted)} else {n <- nobs}
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model - this includes LM
    BIC <- -2*LL + K * log(n)
    if(return.K == TRUE) BIC[1] <- K #attributes the first element of BIC to K
    BIC
  }



##betareg objects
useBIC.betareg <-
  function(mod, return.K = FALSE, nobs = NULL, ...){
    
    if(identical(nobs, NULL)) {n <- length(mod$fitted)} else {n <- nobs}
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model - this includes LM
    BIC <- -2*LL + K * log(n)
    if(return.K == TRUE) BIC[1] <- K #attributes the first element of BIC to K
    BIC
  }



##clm objects
useBIC.clm <-
  function(mod, return.K = FALSE, nobs = NULL, ...){
    
    if(identical(nobs, NULL)) {n <- length(fitted(mod))} else {n <- nobs}
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model - this includes LM
    BIC <- -2*LL + K * log(n)
    if(return.K == TRUE) BIC[1] <- K #attributes the first element of BIC to K
    BIC       
  }


##clmm objects
useBIC.clmm <-
  function(mod, return.K = FALSE, nobs = NULL, ...){

    if(identical(nobs, NULL)) {n <- length(fitted(mod))} else {n <- nobs}
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model - this includes LM
    BIC <- -2*LL + K * log(n)
    if(return.K == TRUE) BIC[1] <- K #attributes the first element of BIC to K
    BIC
  }



##coxme objects
useBIC.coxme <- function(mod, return.K = FALSE, nobs = NULL, ...){
  
  if(identical(nobs, NULL)) {n <- length(mod$linear.predictor)} else {n <- nobs}
  LL <- extractLL(mod)[1]
  K <- attr(extractLL(mod), "df")  #extract correct number of parameters included in model
  BIC <- -2*LL + K * log(n)
  if(return.K==TRUE) BIC[1] <- K #attributes the first element of BIC to K
  BIC
}



##coxph objects
useBIC.coxph <- function(mod, return.K = FALSE, nobs = NULL, ...){
    
  if(identical(nobs, NULL)) {n <- length(residuals(mod))} else {n <- nobs}
  LL <- extractLL(mod)[1]
  K <- attr(extractLL(mod), "df")  #extract correct number of parameters included in model
  BIC <- -2*LL + K * log(n)
  if(return.K == TRUE) BIC[1] <- K #attributes the first element of BIC to K
  BIC
}



##fitdist (from fitdistrplus)
useBIC.fitdist <-
  function(mod, return.K = FALSE, nobs = NULL, ...){
    
    if(identical(nobs, NULL)) {n <- mod$n} else {n <- nobs}
    LL <- logLik(mod)
    K <- length(mod$estimate)
    BIC <- -2*LL + K * log(n)
    if(return.K == TRUE) BIC[1] <- K #attributes the first element of BIC to K
    BIC
  }



##fitdistr (from MASS)
useBIC.fitdistr <-
  function(mod, return.K = FALSE, nobs = NULL, ...){
    
    if(identical(nobs, NULL)) {n <- mod$n} else {n <- nobs}
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model - this includes LM
    BIC <- -2*LL + K * log(n)
    if(return.K == TRUE) BIC[1] <- K #attributes the first element of BIC to K
    BIC
  }



##glm and lm objects
useBIC.glm <-
    function(mod, return.K = FALSE, nobs = NULL, c.hat = 1, ...){
    
        if(is.null(nobs)) {
            n <- length(mod$fitted)
        } else {n <- nobs}
    
        LL <- logLik(mod)[1]
        K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model - this includes LM
      
        if(c.hat == 1) {
            BIC <- -2*LL + K * log(n)
        }
        if(c.hat > 1 && c.hat <= 4) {
            K <- K+1
            BIC <- -2*LL/c.hat + K * log(n)
        }

        if(c.hat > 4) stop("High overdispersion and model fit is questionable\n")
        if(c.hat < 1) stop("You should set \'c.hat\' to 1 if < 1, but values << 1 might also indicate lack of fit\n")

        ##check if negative binomial and add 1 to K for estimation of theta if glm( ) was used
        if(!is.na(charmatch(x="Negative Binomial", table=family(mod)$family))) {
            if(!identical(class(mod)[1], "negbin")) { #if not negbin, add + 1 because k of negbin was estimated glm.convert( ) screws up logLik
                K <- K+1
                BIC <- -2*LL + K * log(n)
            }
            if(c.hat != 1) stop("You should not use the c.hat argument with the negative binomial")
        }
        
        ##add 1 for theta parameter in negative binomial fit with glm( )

        ##check if gamma and add 1 to K for estimation of shape parameter if glm( ) was used
        if(identical(family(mod)$family, "Gamma") && c.hat > 1) stop("You should not use the c.hat argument with the gamma")
        
        ##an extra condition must be added to avoid adding a parameter for theta with negative binomial when glm.nb( ) is fit which estimates the correct number of parameters
        if(return.K == TRUE) BIC[1] <- K #attributes the first element of BIC to K
        BIC
    }



##glmmTMB
useBIC.glmmTMB <-
    function(mod, return.K = FALSE, nobs = NULL, c.hat = 1, ...){
    
        if(is.null(nobs)) {
            n <- nrow(mod$frame)
            names(n) <- NULL
        } else {n <- nobs}
        
        LL <- logLik(mod)[1]
        K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model - this includes LM
      
        if(c.hat == 1) {
            BIC <- -2*LL + K * log(n)
        }
        if(c.hat > 1 && c.hat <= 4) {
            K <- K+1
            BIC <- -2*LL/c.hat + K * log(n)
        }

        if(c.hat > 4) stop("High overdispersion and model fit is questionable\n")
        if(c.hat < 1) stop("You should set \'c.hat\' to 1 if < 1, but values << 1 might also indicate lack of fit\n")

        if(return.K == TRUE) BIC[1] <- K #attributes the first element of BIC to K
        BIC
    }



##gls objects
useBIC.gls <-
    function(mod, return.K = FALSE, nobs = NULL, ...){
    
        if(identical(nobs, NULL)) {n<-length(mod$fitted)} else {n <- nobs}
        LL <- logLik(mod)[1]
        K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model - this includes LM
        BIC <- -2*LL + K * log(n)
        if(return.K == TRUE) BIC[1] <- K #attributes the first element of BIC to K
        BIC
    }



##gnls objects
useBIC.gnls <-
  function(mod, return.K = FALSE, nobs = NULL, ...){

    if(identical(nobs, NULL)) {n <- length(fitted(mod))} else {n <- nobs}
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model
    BIC <- -2*LL + K * log(n)
    if(return.K == TRUE) BIC[1] <- K #attributes the first element of BIC to K
    BIC
  }



##hurdle objects
useBIC.hurdle <-
  function(mod, return.K = FALSE, nobs = NULL, ...){
    
    if(identical(nobs, NULL)) {n <- length(mod$fitted)} else {n <- nobs}
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model - this includes LM
    BIC <- -2*LL + K * log(n)
    if(return.K == TRUE) BIC[1] <- K #attributes the first element of BIC to K
    BIC
  }


##lavaan
useBIC.lavaan <-
  function(mod, return.K = FALSE, nobs = NULL, ...){
    
    if(identical(nobs, NULL)) {n <- mod@Data@nobs[[1]]} else {n <- nobs}
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model
    BIC <- -2*LL + K * log(n)
    if(return.K == TRUE) BIC[1] <- K #attributes the first element of BIC to K
    BIC
  }


##lm objects
useBIC.lm <-
  function(mod, return.K = FALSE, nobs = NULL, ...){
    
    if(identical(nobs, NULL)) {n <- length(mod$fitted)} else {n <- nobs}
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model - this includes LM
    BIC <- -2*LL + K * log(n)
    if(return.K == TRUE) BIC[1] <- K #attributes the first element of BIC to K
    BIC
  }



##lme objects
useBIC.lme <-
  function(mod, return.K = FALSE, nobs = NULL, ...){
    
    if(identical(nobs, NULL)) {n <- nrow(mod$fitted)} else {n <- nobs}
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model - this includes LM
    BIC <- -2*LL + K * log(n)
    if(return.K == TRUE) BIC[1] <- K #attributes the first element of BIC to K
    BIC
  }


##lmekin objects
useBIC.lmekin <-
  function(mod, return.K = FALSE, nobs = NULL, ...){

    if(identical(nobs, NULL)) {n <- length(mod$residuals)} else {n <- nobs}
    LL <- extractLL(mod)[1]
    K <- attr(extractLL(mod), "df")  #extract correct number of parameters included in model - this includes LM
    BIC <- -2*LL + K * log(n)
    if(return.K == TRUE) BIC[1] <- K #attributes the first element of BIC to K
    return(BIC)
  }


##maxlike objects
useBIC.maxlikeFit <- function(mod, return.K = FALSE, nobs = NULL, c.hat = 1, ...) {
  
  LL <- extractLL(mod)[1]
  K <- attr(extractLL(mod), "df")

  if(is.null(nobs)) {
    n <- nrow(mod$points.retained)
  } else {n <- nobs}
  
  BIC <- -2*LL + K * log(n)
  
  if(c.hat != 1) stop("\nThis function does not support overdispersion in \'maxlikeFit\' models\n")

  if(identical(return.K, TRUE)) {
    return(K)
  } else {return(BIC)}
}


##mer object
useBIC.mer <-
  function(mod, return.K = FALSE, nobs = NULL, ...){

    if(is.null(nobs)) {
      n <- mod@dims["n"]
    } else {n <- nobs}
      
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model
    
    BIC <- -2*LL + K * log(n)
    
    if(return.K == TRUE) BIC[1] <- K #attributes the first element of BIC to K
    BIC
  }


##merMod objects
useBIC.merMod <-
  function(mod, return.K = FALSE, nobs = NULL, ...){
 
    if(is.null(nobs)) {
      n <- mod@devcomp$dims["n"]
      names(n) <- NULL
    } else {n <- nobs}
    
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model
      
    BIC <- -2*LL + K * log(n)
                    
    if(return.K == TRUE) BIC[1] <- K #attributes the first element of BIC to K
    BIC
  }


##lmerModLmerTest objects
useBIC.lmerModLmerTest <-
  function(mod, return.K = FALSE, nobs = NULL, ...){
 
    if(is.null(nobs)) {
      n <- mod@devcomp$dims["n"]
      names(n) <- NULL
    } else {n <- nobs}
    
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model
      
    BIC <- -2*LL + K * log(n)
                    
    if(return.K == TRUE) BIC[1] <- K #attributes the first element of BIC to K
    BIC
  }


##mult objects
useBIC.multinom <-
  function(mod, return.K = FALSE, nobs = NULL, c.hat = 1, ...){

    if(identical(nobs, NULL)) {n<-length(mod$fitted)/length(mod$lev)} else {n <- nobs}
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model - this includes LM
      
    if(c.hat == 1) {
        BIC <- -2*LL + K * log(n)
    }
    if(c.hat > 1 && c.hat <= 4) {
      K <- K+1
      BIC <- -2*LL/c.hat + K * log(n)
  }
    if(c.hat > 4) stop("High overdispersion and model fit is questionable")
    
    if(return.K==TRUE) BIC[1]<-K #attributes the first element of BIC to K
    BIC
}


##nlme objects
useBIC.nlme <-
    function(mod, return.K = FALSE, nobs = NULL, ...){

        if(identical(nobs, NULL)) {n <- nrow(mod$fitted)} else {n <- nobs}
        LL <- logLik(mod)[1]
        K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model - this includes LM
        BIC <- -2*LL + K * log(n)
        if(return.K == TRUE) BIC[1] <- K #attributes the first element of BIC to K
        BIC
    }


##nls objects
useBIC.nls <-
  function(mod, return.K = FALSE, nobs = NULL, ...){

    if(identical(nobs, NULL)) {n <- length(fitted(mod))} else {n <- nobs}
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model
    BIC <- -2*LL + K * log(n)
    if(return.K == TRUE) BIC[1] <- K #attributes the first element of BIC to K
    BIC
  }


##polr objects
useBIC.polr <-
function(mod, return.K = FALSE, nobs = NULL, ...){

  if(identical(nobs, NULL)) {n<-length(mod$fitted)} else {n <- nobs}
  LL <- logLik(mod)[1]
  K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model - this includes LM
  BIC <- -2*LL + K * log(n)                                
  if(return.K == TRUE) BIC[1] <- K #attributes the first element of BIC to K
  BIC
}



##rlm objects
##only valid for M-estimation (Huber M-estimator)
##modified from Tharmaratnam and Claeskens 2013 (equation 8)
##useBIC.rlm <- function(mod, return.K = FALSE, nobs = NULL, ...)
##{

##  if(second.ord == TRUE) stop("\nOnly 'second.ord = FALSE' is supported for 'rlm' models\n")

##  ##extract design matrix
##  X <- model.matrix(mod)
  
##  ##extract scale
##  scale.m <- mod$s

##  ##extract threshold value
##  cval <- mod$k2
    
##  ##extract residuals
##  res <- residuals(mod)
##  res.scaled <- res/scale.m
##  n <- length(res)
  
##  ##partial derivatives based on Huber's loss function
##  dPsi <- ifelse(abs(res.scaled) <= cval, 2, 0)
##  Psi <- (ifelse(abs(res.scaled) <= cval, 2*res.scaled, 2*cval*sign(res.scaled)))^2
    
##  J <- (t(X) %*% diag(as.vector(dPsi)) %*% X * (1/(scale.m^2)))/n
##  inv.J <- solve(J)
  
##  ##variance
##  K.var <- (t(X) %*% diag(as.vector(Psi)) %*% X * (1/(scale.m^2)))/n
##  AIC <- 2*n*(log(scale.m)) + 2 * sum(diag(inv.J %*%(K.var)))
  
##  if(return.K) {AIC <- 2 * sum(diag(inv.J %*%(K.var)))}
##  return(AIC)
##}

##the estimator below extracts the estimates obtained from M- or MM-estimator
##and plugs them in the normal likelihood function
useBIC.rlm <-
  function(mod, return.K = FALSE, nobs = NULL, ...){

    if(identical(nobs, NULL)) {n <- length(mod$fitted)} else {n <- nobs}
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model - this includes LM
    BIC <- -2*LL + K * log(n)
    if(return.K == TRUE) BIC[1] <- K #attributes the first element of BIC to K
    BIC
  }



##survreg objects
useBIC.survreg <-
  function(mod, return.K = FALSE, nobs = NULL, ...){
    
    if(identical(nobs, NULL)) {n <- nrow(mod$y)} else {n <- nobs}
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model - this includes LM
    BIC <- -2*LL + K * log(n)
    if(return.K == TRUE) BIC[1] <- K #attributes the first element of BIC to K
    BIC
  }



##unmarkedFit objects
##create function to extract BIC from 'unmarkedFit'
useBIC.unmarkedFit <- function(mod, return.K = FALSE, nobs = NULL, c.hat = 1, ...) {
  
  LL <- extractLL(mod)[1]
  K <- attr(extractLL(mod), "df")
  
  if(is.null(nobs)) {
    n <- dim(mod@data@y)[1]
  } else {n <- nobs}
  
  if(c.hat == 1) {
    BIC <- -2*LL + K * log(n)
  }
  if(c.hat > 1 && c.hat <= 4) {
    ##adjust parameter count to include estimation of dispersion parameter
    K <- K + 1
    BIC <- -2*LL/c.hat + K * log(n)
  }

  if(c.hat > 4) stop("\nHigh overdispersion and model fit is questionable\n")
  if(c.hat < 1) stop("\nYou should set \'c.hat\' to 1 if < 1, but values << 1 might also indicate lack of fit\n")

  if(identical(return.K, TRUE)) {
    return(K)
  } else {return(BIC)}
}


##vglm objects
useBIC.vglm <- function(mod, return.K = FALSE, nobs = NULL, c.hat = 1, ...){
    
    if(is.null(nobs)) {
      n <- nrow(mod@fitted.values)
    } else {n <- nobs}
    
    LL <- extractLL(mod)[1]

    ##extract number of estimated parameters
    K <- attr(extractLL(mod), "df")
    
    if(c.hat !=1) {
      fam.name <- mod@family@vfamily
      if(fam.name != "poissonff" && fam.name != "binomialff") stop("\nOverdispersion correction only appropriate for Poisson or binomial models\n")
    }
    if(c.hat == 1) {
      BIC <- -2*LL + K * log(n)
    }
    if(c.hat > 1 && c.hat <= 4) {
      K <- K + 1
      BIC <- -2*LL/c.hat + K * log(n)
    }

    if(c.hat > 4) stop("High overdispersion and model fit is questionable\n")
    if(c.hat < 1) stop("You should set \'c.hat\' to 1 if < 1, but values << 1 might also indicate lack of fit\n")

    if(return.K == TRUE) BIC[1] <- K #attributes the first element of BIC to K
    BIC
  }



##zeroinfl objects
useBIC.zeroinfl <-
  function(mod, return.K = FALSE, nobs = NULL, ...){
    
    if(identical(nobs, NULL)) {n <- length(mod$fitted)} else {n <- nobs}
    LL <- logLik(mod)[1]
    K <- attr(logLik(mod), "df")  #extract correct number of parameters included in model - this includes LM
    BIC <- -2*LL + K * log(n)
    if(return.K == TRUE) BIC[1] <- K #attributes the first element of BIC to K
    BIC
}


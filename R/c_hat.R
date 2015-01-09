##create generic c_hat
c_hat <- function(mod, ...) {
  ##format list according to model class
  UseMethod("c_hat", mod)
}

##default to indicate when object class not supported
c_hat.default <- function(mod, ...) {
  stop("\nFunction not yet defined for this object class\n")
}



##function to compute c-hat from Poisson or binomial GLM with success/total syntax
c_hat.glm <- function(mod, ...){

  ##determine family of model
  fam <- family(mod)$family

  ##if binomial, check if n > 1 for each case
  if(fam == "binomial") {
    if(identical(unique(mod$prior.weights), 1)) {
      stop("\nWith a binomial distribution, the number of successes must be summarized for valid computation of c-hat\n")
    }
  }

  ##Poisson or binomial
  if(!any(fam == c("poisson", "binomial"))) {
    stop("\nEstimation of c-hat only valid for Poisson or binomial GLM's\n")
  }
  
  ##estimate Pearson chi-square
  chisq <- sum(residuals(mod, type = "pearson")^2)
  c_hat.est <- chisq/mod$df.residual
  return(c_hat.est)
}



##method for GLMM from lme4
c_hat.merMod <- function(mod, ...) {

  #determine family of model
  fam <- family(mod)$family

  ##if binomial, check if n > 1 for each case
  if(fam == "binomial") {
    if(identical(unique(mod@resp$n), 1)) {
      stop("\nWith a binomial distribution, the number of successes must be summarized for valid computation of c-hat\n")
    }
  }
      
  ##Poisson or binomial
  if(!any(fam == c("poisson", "binomial"))) {
    stop("\nEstimation of c-hat only valid for Poisson or binomial GLMM's\n")
  }
    
  ##number of parameters estimated
  n.parms <- attr(logLik(mod), "df")
  
  ##total number of observations
  n.obs <- nrow(model.frame(mod))

  ##residual df
  res.df <- n.obs - n.parms
  
  chisq <- sum(residuals(mod, type = "pearson")^2)
  c_hat.est <- chisq/res.df
  return(c_hat.est)
}

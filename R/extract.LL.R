##extract log-likelihood of model
extract.LL <- function(mod) {
  if(identical(class(mod), "coxph") || identical(class(mod), c("coxph.null", "coxph"))) {
    LL <- extract.LL.coxph(mod)
  } else {
  if(attr(regexpr(pattern = "unmarked", text = class(mod)), "match.length") == -1) {
    stop("\nObject must be of 'unmarkedFit' class\n")}
    LL <- extract.LL.unmarked(mod)
  }
  return(LL)
}
  
  
    
extract.LL.coxph <- function(mod) {
  if(identical(class(mod), "coxph") || identical(class(mod), c("coxph.null", "coxph"))) {
    coefs <- coef(mod)
    if(is.null(coefs)) {
      ncoefs <- 0
      LL <- mod$loglik[1] #when null model, only 1 log-likelihood value
    } else {
      ncoefs <- length(coefs)
      LL <- mod$loglik[2] #second value is the logLik at the solution
    }
    
    attr(LL, "df") <- ncoefs
  } else {stop("This function is only appropriate with the \'coxph\' class\n")}
  return(LL)
}


##extract logLik from 'unmarkedFit' objects
extract.LL.unmarked <- function(mod) {
  ##check if model is of 'unmarkedFit' class
  if(attr(regexpr(pattern = "unmarked", text = class(mod)), "match.length") == -1) stop("\nObject must be of 'unmarkedFit' class\n")
  LL <- -1*mod@negLogLike
  df <- length(mod@opt$par)

  attr(LL, "df") <- df
  return(LL)
}

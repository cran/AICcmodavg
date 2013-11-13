AICc.coxph <- function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL){
  ##check if correct class
  if(identical(class(mod), "coxph") || identical(class(mod), c("coxph.null", "coxph")) || identical(class(mod), c("clogit", "coxph"))) {
    
    if(identical(nobs, NULL)) {n <- length(residuals(mod))} else {n <- nobs}
    LL <- extract.LL.coxph(mod)[1]
    K <- attr(extract.LL.coxph(mod), "df")  #extract correct number of parameters included in model
    if(second.ord==TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))}  else{AICc <- -2*LL+2*K}
    if(return.K==TRUE) AICc[1] <- K #attributes the first element of AICc to K
    AICc
  } else {stop("This function is only appropriate with objects of \'coxph\' class")}
  
}

##extract logLik from 'unmarkedFit' objects
extract.LL.unmarked <- function(mod) {
  ##check if model is of 'unmarkedFit' class
  if(attr(regexpr(pattern = "unmarked", text = class(mod)), "match.length") == -1) stop("\nObject must be of 'unmarkedFit' class\n")
  LL <- -1*mod@negLogLike
  return(LL)
}


##create function to extract AICc from 'unmarkedFit'
AICc.unmarked <- function(mod, return.K = FALSE, c.hat = 1, second.ord = TRUE, nobs = NULL) {
  
  LL <- extract.LL.unmarked(mod)
  K <- length(mod@opt$par)
  
  if(is.null(nobs)) {
    n <- dim(mod@data@y)[1]
  }
  
  if(c.hat == 1) {
    if(second.ord == TRUE) {AICc <- -2 * LL + 2 * K * (n/(n - K - 1))} else {AICc <- -2*LL + 2*K}
  }
  if(c.hat > 1 && c.hat <= 4) {
    ##adjust parameter count to include estimation of dispersion parameter
    K <- K + 1
    if(second.ord == TRUE) {
      AICc <- (-2 * LL/c.hat) + 2 * K * (n/(n - K - 1))
    } else {
      AICc <- (-2 * LL/c.hat) + 2*K}
  }

  if(c.hat > 4) stop("High overdispersion and model fit is questionable\n")
  if(c.hat < 1) stop("You should set \'c.hat\' to 1 if < 1, but values << 1 might also indicate lack of fit\n")

  
  if(identical(return.K, TRUE)) {
    return(K)
  } else {return(AICc)}
}

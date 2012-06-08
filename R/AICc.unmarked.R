##create function to extract AICc from 'unmarkedFit'
AICc.unmarked <- function(mod, return.K = FALSE, c.hat = 1, second.ord = TRUE, nobs = NULL) {
  
  LL <- extract.LL.unmarked(mod)[1]
  K <- attr(extract.LL.unmarked(mod), "df")
  
  if(is.null(nobs)) {
    n <- dim(mod@data@y)[1]
  } else {n <- nobs}
  
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

  if(c.hat > 4) stop("\nHigh overdispersion and model fit is questionable\n")
  if(c.hat < 1) stop("\nYou should set \'c.hat\' to 1 if < 1, but values << 1 might also indicate lack of fit\n")

  
  if(identical(return.K, TRUE)) {
    return(K)
  } else {return(AICc)}
}

DIC <-
  function(mod, return.pD = FALSE) {
    dic <- NULL

    ##determine if bugs
    if(identical(class(mod), "bugs")) {
      dic <- DIC.bugs(mod = mod, return.pD = return.pD)
    }

    ##determine if rjags
    if(identical(class(mod), "rjags")) {
      dic <- DIC.rjags(mod = mod, return.pD = return.pD)
    }

    return(dic)
  }




##function to extract DIC of bugs model class
DIC.bugs <- function(mod, return.pD = FALSE){
#  DIC = posterior mean of the deviance + pD, where pD is the effective number of parameters
  if(identical(class(mod), "bugs")) {
    if(return.pD == FALSE){
      DIC <- mod$DIC
    } else {DIC <- mod$pD}
    return(DIC)
  } else {stop("This function is only appropriate with objects of \'bugs\' class")}
}


##function to extract DIC of bugs model class
DIC.rjags <- function(mod, return.pD = FALSE){
#  DIC = posterior mean of the deviance + pD, where pD is the effective number of parameters
  if(identical(class(mod), "rjags")) {
    if(return.pD == FALSE){
      DIC <- mod$BUGSoutput$DIC
    } else {DIC <- mod$BUGSoutput$pD}
    return(DIC)
  } else {stop("This function is only appropriate with objects of \'bugs\' class")}
}

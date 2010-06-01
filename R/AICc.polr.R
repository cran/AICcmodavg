AICc.polr <-
function(mod, return.K=FALSE, second.ord=TRUE, nobs=NULL){
#check whether object is of appropriate class
  if(identical(class(mod), "polr")) {
    
    if(identical(nobs, NULL)) {n<-length(mod$fitted)} else {n <- nobs}
    LL<-logLik(mod)[1]
    K<-attr(logLik(mod), "df")  #extract correct number of parameters included in model - this includes LM
      
    if(second.ord==TRUE) {
      AICc <- -2*LL+2*K*(n/(n-K-1))
    } else{
        AICc <- -2*LL+2*K
      }
                                 
    if(return.K==TRUE) AICc[1]<-K #attributes the first element of AICc to K
    AICc
     
                   
  } else {stop("This function is only appropriate with objects of \'polr\' class")}
  
}


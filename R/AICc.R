AICc <-
function(mod, return.K=FALSE, c.hat=1, second.ord=TRUE, nobs=NULL) {
  aicc <- NULL

  #determine if lm or glm
  if(class(mod)[1]=="lm" || class(mod)[1]=="glm") {
    aicc <- AICc.glm(mod, return.K, c.hat, second.ord, nobs)
  }   

  #determine if multinom
  if(class(mod)[1]=="multinom") {
    aicc <- AICc.mult(mod, return.K, c.hat, second.ord, nobs)
  }      

  #determine if lme
  if(class(mod)[1]=="lme") {
    aicc <- AICc.lme(mod, return.K, second.ord, nobs)
  }      


  #if(class(mod)[1]=="nlm") {aicc <- AICc.nlm(mod, return.K)}      #determine if object from nlm optimizer
return(aicc)
}

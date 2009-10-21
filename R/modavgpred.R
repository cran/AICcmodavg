modavgpred <-
function(cand.set, modnames, newdata, type="response", c.hat=1, gamdisp=NULL, second.ord=TRUE, nobs=NULL) {
#newdata is data frame with exact structure of the original data frame (same variable names and type)
if(type=="terms") {stop("The terms argument is not defined for this function")}
dispersion <- c.hat
if(c.hat>1) {dispersion <- c.hat }
if(!is.null(gamdisp)) {dispersion <- gamdisp}
if(c.hat>1 && !is.null(gamdisp)) {stop("You cannot specify values for both c.hat and gamdisp")}
#dispersion is the dispersion parameter - this influences the SE's (to specify dispersion parameter for either overdispersed Poisson or Gamma glm)
#type enables to specify either "response" (original scale = point estimate) or "link" (linear predictor)
  
#check if object is of "lm" or "glm" class
#extract classes
  mod.class <- unlist(lapply(X=cand.set, FUN=class))
#check if all are identical
  check.class <- unique(mod.class)
  
  if(identical(check.class, "lm") || identical(check.class, c("glm", "lm")))  {

#determine number of observations in data set
nobs <- dim(newdata)[1]

#store AICc table
AICctab<-aictab(cand.set=cand.set, modnames=modnames, c.hat=c.hat, second.ord=second.ord, nobs=nobs, sort=FALSE)

#create object to hold Model-averaged estimates and unconditional SE's
Mod.avg.out<-matrix(NA, nrow=nobs, ncol=2)
colnames(Mod.avg.out)<-c("Mod.avg.est", "Uncond.SE")

#begin loop - AICc
if(second.ord==TRUE && c.hat==1){
  for (obs in 1:nobs) {

#extract fitted value for observation obs
fit<-unlist(lapply(X=cand.set, FUN=function(i)predict(i, se.fit=TRUE, newdata=newdata[obs, ], type=type, dispersion=dispersion)$fit))
#extract SE for fitted value for observation obs
SE<-unlist(lapply(X=cand.set, FUN=function(i)predict(i, se.fit=TRUE, newdata=newdata[obs, ], type=type, dispersion=dispersion)$se.fit))

#create temporary data.frame to store fitted values and SE 
AICctmp<-AICctab
AICctmp$fit<-fit
AICctmp$SE<-SE

#compute model averaged prediction and store in output matrix
Mod.avg.out[obs, 1]<-sum(AICctmp$AICcWt*AICctmp$fit)
#compute unconditional SE and store in output matrix
Mod.avg.out[obs, 2]<-sum(AICctmp$AICcWt*sqrt(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2))
}
}


#create temporary data.frame to store fitted values and SE - QAICc
if(second.ord==TRUE && c.hat > 1) {
for (obs in 1:nobs) {

#extract fitted value for observation obs
fit<-unlist(lapply(X=cand.set, FUN=function(i)predict(i, se.fit=TRUE, newdata=newdata[obs, ], type=type, dispersion=dispersion)$fit))
#extract SE for fitted value for observation obs
SE<-unlist(lapply(X=cand.set, FUN=function(i)predict(i, se.fit=TRUE, newdata=newdata[obs, ], type=type, dispersion=dispersion)$se.fit))

QAICctmp<-AICctab
QAICctmp$fit<-fit
QAICctmp$SE<-SE

#compute model averaged prediction and store in output matrix
Mod.avg.out[obs, 1]<-sum(QAICctmp$QAICcWt*QAICctmp$fit)
#compute unconditional SE and store in output matrix
Mod.avg.out[obs, 2]<-sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2))
}
}


#create temporary data.frame to store fitted values and SE - AIC
if(second.ord==FALSE && c.hat==1) {
for (obs in 1:nobs) {

#extract fitted value for observation obs
fit<-unlist(lapply(X=cand.set, FUN=function(i)predict(i, se.fit=TRUE, newdata=newdata[obs, ], type=type, dispersion=dispersion)$fit))
#extract SE for fitted value for observation obs
SE<-unlist(lapply(X=cand.set, FUN=function(i)predict(i, se.fit=TRUE, newdata=newdata[obs, ], type=type, dispersion=dispersion)$se.fit))

AICtmp<-AICctab
AICtmp$fit<-fit
AICtmp$SE<-SE

#compute model averaged prediction and store in output matrix
Mod.avg.out[obs, 1]<-sum(AICtmp$AICWt*AICtmp$fit)
#compute unconditional SE and store in output matrix
Mod.avg.out[obs, 2]<-sum(AICtmp$AICWt*sqrt(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2))
}
}

#create temporary data.frame to store fitted values and SE - QAIC
if(second.ord==FALSE && c.hat > 1) {
for (obs in 1:nobs) {

#extract fitted value for observation obs
fit<-unlist(lapply(X=cand.set, FUN=function(i)predict(i, se.fit=TRUE, newdata=newdata[obs, ], type=type, dispersion=dispersion)$fit))
#extract SE for fitted value for observation obs
SE<-unlist(lapply(X=cand.set, FUN=function(i)predict(i, se.fit=TRUE, newdata=newdata[obs, ], type=type, dispersion=dispersion)$se.fit))

QAICtmp<-AICctab
QAICtmp$fit<-fit
QAICtmp$SE<-SE

#compute model averaged prediction and store in output matrix
Mod.avg.out[obs, 1]<-sum(QAICtmp$QAICWt*QAICtmp$fit)
#compute unconditional SE and store in output matrix
Mod.avg.out[obs, 2]<-sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2))
}
}

Mod.pred.list <- list("type" = type, "mod.avg.pred" = Mod.avg.out[,1], "uncond.se" = Mod.avg.out[,2])
class(Mod.pred.list) <- c("modavgpred", "list")
return(Mod.pred.list)

} else {stop("This function is only appropriate with either \'lm\' or \'glm\' classes\n")}
}

print.modavgpred <- function(x, digits = 4, ...) {
  cat("\nModel-averaged predictions on the", x$type, "scale based on entire model set:\n\n")
  nice.tab <- cbind(x$mod.avg.pred, x$uncond.se)
  colnames(nice.tab) <- c("mod.avg.pred", "uncond.se")
  nrows <- dim(nice.tab)[1]
  rownames(nice.tab) <- 1:nrows
  print(round(nice.tab, digits = digits))
  cat("\n")
}

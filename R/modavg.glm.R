modavg.glm <-
function(cand.set, parm, modnames, c.hat=1, gamdisp=NULL, conf.level=0.95, second.ord=TRUE, nobs=NULL){
#check if class is appropriate
#extract classes
  mod.class <- unlist(lapply(X=cand.set, FUN=class))
#check if all are identical
  check.class <- unique(mod.class)
  
  if(identical(check.class, "lm") || identical(check.class, c("glm", "lm")))  {

mod_formula<-lapply(cand.set, FUN=function(i) rownames(summary(i)$coefficients)) #extract model formula for each model in cand.set

#setup matrix to indicate presence of parm in the model
include <- matrix(NA, nrow=length(cand.set), ncol=1)
#add a check for multiple instances of same variable in given model (i.e., interactions)
include.check <- matrix(NA, nrow=length(cand.set), ncol=1)

#iterate over each formula in mod_formula list
for (i in 1:length(cand.set)) {
  idents <- NULL
  idents.check <- NULL
  form <- mod_formula[[i]]

  #iterate over each element of formula[[i]] in list
  for (j in 1:length(form)) {
    idents[j] <- identical(parm, form[j])
    idents.check[j] <- ifelse(attr(regexpr(parm, form[j]), "match.length")=="-1", 0, 1)  
  }
  include[i] <- ifelse(any(idents==1), 1, 0)
  include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
}

#check for duplicates in same model
if(any(include.check == "duplicates")) {
  warning("Some models include more than one instance of the parameter of interest:  these models were excluded during model averaging")
  #exclude models with duplicates from model averaging
  include[which(include.check=="duplicates")] <- 0}


#add a check to determine if include always == 0
if (sum(include)==0) {stop("Parameter not found in any of the candidate models") }

new.cand.set<-cand.set[which(include==1)] #select models including a given parameter
new.mod.name<-modnames[which(include==1)]    #update model names
##

new_table<-aictab.glm(cand.set=new.cand.set, modnames=new.mod.name, sort=FALSE, c.hat=c.hat, second.ord=second.ord, nobs=nobs)  #recompute AIC table and associated measures
new_table$Beta_est<-unlist(lapply(new.cand.set, FUN=function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
new_table$SE<-unlist(lapply(new.cand.set, FUN=function(i) summary(i)$coefficients[paste(parm),2]))
if(c.hat > 1) {new_table$SE<-new_table$SE*sqrt(c.hat)} #if c-hat is estimated adjust the SE's by multiplying with sqrt of c-hat
gam1<-unlist(lapply(new.cand.set, FUN=function(i) family(i)$family[1]=="Gamma")) #check for gamma regression models
if(any(gam1)==TRUE)  {new_table$SE<-unlist(lapply(new.cand.set, FUN=function(i) summary(i, dispersion=gamdisp)$coefficients[paste(parm),2]))} #correct SE's for estimates of gamma regressions
#compute model-averaged estimates, unconditional SE, and 95% CL
if(c.hat == 1 && second.ord == TRUE) {Modavg_beta<-sum(new_table$AICcWt*new_table$Beta_est)
                                      Uncond_SE<-sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))}
if(c.hat > 1 && second.ord == TRUE) {Modavg_beta<-sum(new_table$QAICcWt*new_table$Beta_est)
                                      Uncond_SE<-sum(new_table$QAICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))}     #if c-hat is estimated adjust table names
if(c.hat == 1 && second.ord == FALSE) {Modavg_beta<-sum(new_table$AICWt*new_table$Beta_est)
                                      Uncond_SE<-sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))}
if(c.hat > 1 && second.ord == FALSE) {Modavg_beta<-sum(new_table$QAICWt*new_table$Beta_est)
                                      Uncond_SE<-sum(new_table$QAICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))}     #if c-hat is estimated adjust table names

zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
Lower_CL<-Modavg_beta-zcrit*Uncond_SE
Upper_CL<-Modavg_beta+zcrit*Uncond_SE
out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta, "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL"= Lower_CL, "Upper.CL" = Upper_CL)
} else {stop("This function is only appropriate with either lm or glm classes")}

class(out.modavg) <- c("modavg", "list")
return(out.modavg)

}


modavg.lme <-
function(cand.set, parm, modnames, conf.level=0.95, second.ord=TRUE, nobs=NULL, exclude=NULL, warn=TRUE){

#check if class is appropriate
#extract classes
  mod.class <- unlist(lapply(X=cand.set, FUN=class))
#check if all are identical
  check.class <- unique(mod.class)
  
  if(!identical(check.class, "lme"))  {stop("This function is only appropriate with the \'lme\' class\n")}



mod_formula<-lapply(cand.set, FUN=function(i) labels(summary(i)$coefficients$fixed)) #extract model formula for each model in cand.set
#setup matrix to indicate presence of parms in the model
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
    idents.check[j] <- ifelse(is.na(match(parm, form[j])), 0, 1)  
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

new_table<-aictab.lme(cand.set=new.cand.set, modnames=new.mod.name, sort=FALSE, second.ord=second.ord, nobs=nobs)  #recompute AIC table and associated measures
new_table$Beta_est<-unlist(lapply(new.cand.set, FUN=function(i) fixef(i)[paste(parm)])) #extract beta estimate for parm
new_table$SE<-unlist(lapply(new.cand.set, FUN=function(i) sqrt(diag(vcov(i)))[paste(parm)]))
#compute model-averaged estimates, unconditional SE, and 95% CL
if(second.ord==TRUE) {Modavg_beta<-sum(new_table$AICcWt*new_table$Beta_est)
Uncond_SE<-sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
                    } else {Modavg_beta<-sum(new_table$AICWt*new_table$Beta_est)
Uncond_SE<-sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))}
  
zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
Lower_CL<-Modavg_beta-zcrit*Uncond_SE
Upper_CL<-Modavg_beta+zcrit*Uncond_SE
out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta, "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL, "Upper.CL" = Upper_CL)
  
  class(out.modavg) <- c("modavg", "list")
  return(out.modavg)
}

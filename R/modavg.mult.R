modavg.mult <-
function(cand.set, parm, modnames, c.hat=1, conf.level=0.95, second.ord=TRUE, nobs=NULL){
#check if class is appropriate
#extract classes
  mod.class <- unlist(lapply(X=cand.set, FUN=class))
#check if all are identical
  check.class <- unique(mod.class)
  
  if(identical(check.class, c("multinom", "nnet"))) {
    
#extract model formula for each model in cand.set    
    mod_formula<-lapply(cand.set, FUN=function(i) colnames(summary(i)$coefficients)) 

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
      include[which(include.check=="duplicates")] <- 0
    }


#add a check to determine if include always == 0
    if (sum(include)==0) {stop("Parameter not found in any of the candidate models") }

    new.cand.set<-cand.set[which(include==1)] #select models including a given parameter
    new.mod.name<-modnames[which(include==1)]    #update model names
##


#determine number of levels - 1
    mod.levels<-lapply(cand.set, FUN=function(i) rownames(summary(i)$coefficients)) #extract level of response variable 
    check.levels <- unlist(unique(mod.levels))


#recompute AIC table and associated measures
    new_table<-aictab.mult(cand.set=new.cand.set, modnames=new.mod.name, sort=FALSE, c.hat=c.hat,
                           second.ord=second.ord, nobs=nobs) 

#create object to store model-averaged estimate and SE's of k - 1 level of response
    out.est <- matrix(data=NA, nrow=length(check.levels), ncol=4)
    colnames(out.est) <- c("Mod.avg.est", "Uncond.SE", "Lower.CL", "Upper.CL")
    rownames(out.est) <- check.levels

#iterate over levels of response variable
    for (g in 1:length(check.levels)) {
  #extract beta estimate for parm
      new_table$Beta_est<-unlist(lapply(new.cand.set, FUN=function(i) coef(i)[check.levels[g], paste(parm)]))
  #extract SE of estimate for parm
      new_table$SE<-unlist(lapply(new.cand.set, FUN=function(i) summary(i)$standard.errors[check.levels[g], paste(parm)]))

#if c-hat is estimated adjust the SE's by multiplying with sqrt of c-hat
      if(c.hat > 1) {new_table$SE<-new_table$SE*sqrt(c.hat)} 

#compute model-averaged estimates, unconditional SE, and 95% CL
      if(c.hat == 1 && second.ord == TRUE) {
        Modavg_beta<-sum(new_table$AICcWt*new_table$Beta_est)
        Uncond_SE<-sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

#if c-hat is estimated compute values accordingly and adjust table names
      if(c.hat > 1 && second.ord == TRUE) {
        Modavg_beta<-sum(new_table$QAICcWt*new_table$Beta_est)
        Uncond_SE<-sum(new_table$QAICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      if(c.hat == 1 && second.ord == FALSE) {
        Modavg_beta<-sum(new_table$AICWt*new_table$Beta_est)
        Uncond_SE<-sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

#if c-hat is estimated compute values accordingly and adjust table names  
      if(c.hat > 1 && second.ord == FALSE) {
        Modavg_beta<-sum(new_table$QAICWt*new_table$Beta_est)
        Uncond_SE<-sum(new_table$QAICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }     

      out.est[g, 1] <- Modavg_beta
      out.est[g, 2] <- Uncond_SE
    }
     
    zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
    out.est[,3] <- out.est[,1] - zcrit*out.est[,2]
    out.est[,4] <- out.est[,1] + zcrit*out.est[,2]
    out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = out.est[,1], "Uncond.SE" = out.est[,2], "Conf.level" = conf.level, "Lower.CL"= out.est[,3], "Upper.CL" = out.est[,4])
} else {stop("This function is only appropriate with multinom classes")}

  class(out.modavg) <- c("modavg", "list")
  return(out.modavg)

}


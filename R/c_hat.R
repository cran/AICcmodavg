#function to compute c-hat from Poisson or binomial GLM with success/total syntax
c_hat <- function(mod){

#determine family of model
  fam <- family(mod)$family

#indicate whether family is of correct type
  known <- rep(0,2)
   
  if(fam == "poisson" ) {
    chisq <- sum(residuals(mod, type="pearson")^2)
    c_hat.est <- chisq/mod$df.residual
    known[1] <- 1
  }

  if(fam == "binomial") {
    if( any(mod$prior.weights!=1) ) {
      chisq <- sum(residuals(mod, type="pearson")^2)
      c_hat.est <- chisq/mod$df.residual
      known[2] <- 1
    } else {stop("\nWith a binomial GLM, the number of successes must be summarized for valid computation of c-hat\n")}

  }

#return an error if other than binomial or Poisson glm used
  if( sum(known) == 0 ) stop("\nModel needs to be of class glm with either Poisson or binomial distribution\n")

  return(c_hat.est)
}


#binomial glm
#set.seed(seed=10)
#resp<-rbinom(n=60, size=1, prob=0.5)
#set.seed(seed=10)
#treat<-as.factor(sample(c(rep(x="m", times=30), rep(x="f", times=30))))
#age <- as.factor(c(rep("young", 20), rep("med", 20), rep("old", 20)))
#each invidual has its own response (n = 1)
#mod1<-glm(resp~treat+age, family=binomial)
#c_hat(mod1)

#computing table to summarize successes
#table(resp, treat, age)
#dat2 <- as.data.frame(table(resp, treat, age)) #not quite what we want here
#data2 <- data.frame(success=c(9, 4, 2, 3, 5, 2), sex=c("f", "m", "f", "m", "f", "m"), age=c("med", "med", "old", "old", "young", "young"), total=c(13, 7, 10, 10, 7, 13))
#data2$prop <- data2$success/data2$total
#data2$fail <- data2$total-data2$success
#run model with success/total syntax using weights argument
#mod2 <- glm(prop~sex+age, family=binomial, weights=total, data=data2)
#c_hat(mod2)

#run model with other syntax cbind(success, fail)
#mod3 <- glm(cbind(success, fail)~sex+age, family=binomial, data=data2)
#c_hat(mod3)


modavg.shrink.mult <-
function(cand.set, parm, modnames, c.hat = 1, conf.level = 0.95, second.ord = TRUE, nobs = NULL,
         uncond.se = "revised"){

  ##check if class is appropriate
  ##extract classes
  mod.class <- unlist(lapply(X = cand.set, FUN = class))
  ##check if all are identical
  check.class <- unique(mod.class)
  
  if(!identical(check.class, c("multinom", "nnet"))) {stop("\nThis function is only appropriate with the \'multinom\' class\n")}


  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)
  
 
  ##extract model formula for each model in cand.set    
  mod_formula <- lapply(cand.set, FUN = function(i) colnames(summary(i)$coefficients)) 

  nmods <- length(cand.set)
  
  
  ##determine frequency of each term across models (except (Intercept) ) 
  pooled.terms <- unlist(mod_formula)
  ##remove intercept from vector
  no.int <- pooled.terms[which(pooled.terms != "(Intercept)")]
  terms.freq <- table(no.int)
  if(length(unique(terms.freq)) > 1) stop("\nTo compute a shrinkage version of model-averaged estimate, each term must appear with the same frequency across models\n")

  ##check whether parm is involved in interaction
  parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
  inter.check <- ifelse(attr(regexpr(parm.inter[1], pooled.terms, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                    pooled.terms, fixed = TRUE), "match.length") == "-1", 0, 1)
  if(sum(inter.check) > 0) stop("\nParameter of interest should not be involved in interaction for shrinkage version of model-averaging to be appropriate\n")


  ##determine number of levels - 1
  mod.levels <- lapply(cand.set, FUN=function(i) rownames(summary(i)$coefficients)) #extract level of response variable 
  check.levels <- unlist(unique(mod.levels))


  ##recompute AIC table and associated measures
  new_table <- aictab.mult(cand.set = cand.set, modnames = modnames, sort = FALSE,
                         c.hat = c.hat, second.ord = second.ord, nobs = nobs) 

  ##create object to store model-averaged estimate and SE's of k - 1 level of response
  out.est <- matrix(data = NA, nrow = length(check.levels), ncol = 4)
  colnames(out.est) <- c("Mod.avg.est", "Uncond.SE", "Lower.CL", "Upper.CL")
  rownames(out.est) <- check.levels

  ##iterate over levels of response variable
  for (g in 1:length(check.levels)) {

    ##extract coefficients from each model for given level
    coefs.levels <- lapply(cand.set, FUN = function(i) coef(i)[check.levels[g], ])
    ##extract coefficients from each model for all levels
    SE.all.levels <- lapply(cand.set, FUN = function(i) sqrt(diag(vcov(i))))
    id.coef <- lapply(coefs.levels, FUN = function(i) which(names(i) == paste(parm)))
    ##temporary matrix to hold estimates and SE's from models and set to 0 otherwise
    tmp.coef <- matrix(NA, ncol = 2, nrow = nmods)
    for(k in 1:nmods) {
          tmp.coef[k, 1] <- ifelse(length(id.coef[[k]]) != 0, coefs.levels[[k]][paste(parm)], 0)
          tmp.coef[k, 2] <- ifelse(length(id.coef[[k]]) != 0, SE.all.levels[[k]][paste(check.levels[g], ":",
                                    parm, sep="")], 0)
        }
          
    
    ##extract beta estimate for parm
    new_table$Beta_est <- tmp.coef[, 1]
    
    ##extract SE of estimate for parm
    new_table$SE <- tmp.coef[, 2]

    
    ##add a check to determine if parameter occurs in any model
    if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("\nParameter not found in any of the candidate models\n") }

      
    ##if c-hat is estimated adjust the SE's by multiplying with sqrt of c-hat
    if(c.hat > 1) {new_table$SE <- new_table$SE*sqrt(c.hat)} 

    ##compute model-averaged estimates, unconditional SE, and 95% CL
    #AICc
    if(c.hat == 1 && second.ord == TRUE) {
      Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

      #unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      #revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }


    #QAICc
    #if c-hat is estimated compute values accordingly and adjust table names
    if(c.hat > 1 && second.ord == TRUE) {
      Modavg_beta <- sum(new_table$QAICcWt*new_table$Beta_est)

      #unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$QAICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      #revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$QAICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      } 
    }

    
    #AIC
    if(c.hat == 1 && second.ord == FALSE) {
      Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)

      #unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }
      
      #revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }
    }
    

    #QAIC
    #if c-hat is estimated compute values accordingly and adjust table names  
    if(c.hat > 1 && second.ord == FALSE) {
      Modavg_beta <- sum(new_table$QAICWt*new_table$Beta_est)

      #unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$QAICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      #revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$QAICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      } 
    }
    

    out.est[g, 1] <- Modavg_beta 
    out.est[g, 2] <- Uncond_SE 
  }
     
  zcrit <- qnorm(p = (1-conf.level)/2, lower.tail = FALSE)
  out.est[,3] <- out.est[,1] - zcrit*out.est[,2]
  out.est[,4] <- out.est[,1] + zcrit*out.est[,2]
  out.modavg <- list("Parameter" = paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = out.est[,1],
                     "Uncond.SE" = out.est[,2], "Conf.level" = conf.level, "Lower.CL"= out.est[,3],
                     "Upper.CL" = out.est[,4])

  class(out.modavg) <- c("modavg.shrink", "list")
  return(out.modavg)

}

modavg.shrink.polr <-
function(cand.set, parm, modnames, conf.level = 0.95, second.ord = TRUE, nobs = NULL,
         uncond.se = "revised"){
  ##check if class is appropriate
  ##extract classes
  mod.class <- unlist(lapply(X = cand.set, FUN = class))
  ##check if all are identical
  check.class <- unique(mod.class)
  
  if(!identical(check.class, "polr")) {stop("\nThis function is only appropriate with the \'polr\' class\n")}

  
  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)

  
  ##extract model formula for each model in cand.set    
  mod_formula <- lapply(cand.set, FUN = function(i) rownames(summary(i)$coefficients)) 

  ##determine frequency of each term across models (except (Intercept) ) 
  pooled.terms <- unlist(mod_formula)
  ##remove intercept from vector
  no.int <- pooled.terms[attr(regexpr(pattern = "\\|", text = pooled.terms), "match.length") == -1 ]
  terms.freq <- table(no.int)
  if(length(unique(terms.freq)) > 1) stop("\nTo compute a shrinkage version of model-averaged estimate, each term must appear with the same frequency across models\n")


  ##check whether parm is involved in interaction
  parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
  inter.check <- ifelse(attr(regexpr(parm.inter[1], mod_formula, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                                   mod_formula, fixed = TRUE), "match.length") == "-1", 0, 1)
  if(sum(inter.check) > 0) stop("\nParameter of interest should not be involved in interaction for shrinkage version of model-averaging to be appropriate\n")


  ##compute table
  new_table <- aictab.polr(cand.set = cand.set, modnames = modnames, sort = FALSE,
                         second.ord = second.ord, nobs = nobs)  #recompute AIC table and associated measures
  
  ##add logical test to distinguish between intercepts and other coefs
  if(attr(regexpr(pattern = "\\|", text = parm), "match.length") == -1) {
    new_table$Beta_est <- unlist(lapply(cand.set, FUN = function(i) coef(i)[paste(parm)]))
  } else {new_table$Beta_est <- unlist(lapply(cand.set, FUN = function(i) (i)$zeta[paste(parm)])) }
        
  ##extract beta estimate for parm
  new_table$SE <- unlist(lapply(cand.set, FUN = function(i) sqrt(diag(vcov(i)))[paste(parm)]))

  ##replace NA's with 0
  new_table$Beta_est[is.na(new_table$Beta_est)] <- 0
  new_table$SE[is.na(new_table$SE)] <- 0

  ##add a check to determine if parameter occurs in any model
  if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("\nParameter not found in any of the candidate models\n") }

  
  
  ##AICc
  ##compute model-averaged estimates, unconditional SE, and 95% CL based on AICc
  if(second.ord == TRUE) {
    Modavg_beta<-sum(new_table$AICcWt*new_table$Beta_est)
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE<-sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE<-sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }   
  }

  
  ##AICc  
  ##compute model-averaged estimates, unconditional SE, and 95% CL based on AIC
  if(second.ord == FALSE) {
    Modavg_beta<-sum(new_table$AICWt*new_table$Beta_est)

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE<-sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE<-sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    }
  }

  
  zcrit <- qnorm(p = (1-conf.level)/2, lower.tail = FALSE)
  Lower_CL <- Modavg_beta - zcrit*Uncond_SE
  Upper_CL <- Modavg_beta + zcrit*Uncond_SE
  out.modavg <- list("Parameter" = paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                     "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL"= Lower_CL,
                     "Upper.CL" = Upper_CL)
  class(out.modavg) <- c("modavg.shrink", "list")
  return(out.modavg)
}

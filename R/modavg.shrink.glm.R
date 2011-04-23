modavg.shrink.glm <-
function(cand.set, parm, modnames, c.hat = 1, gamdisp = NULL, conf.level = 0.95, second.ord = TRUE, nobs = NULL,
         uncond.se = "revised"){
  
#check if class is appropriate
#extract classes
  mod.class <- unlist(lapply(X=cand.set, FUN=class))
#check if all are identical
  check.class <- unique(mod.class)

  #check that link function is the same for all models
  if(identical(check.class[1], "glm")) {
  check.link <- unlist(lapply(X = cand.set, FUN=function(i) i$family$link))
  unique.link <- unique(x=check.link)
  if(length(unique.link) > 1) stop(cat("\nIt is not appropriate to compute a model averaged beta estimate\n",
"from models using different link functions\n"))
}
  
  if(identical(check.class, "lm") || identical(check.class, c("glm", "lm")))  {

    ##check family of glm to avoid problems when requesting predictions with argument 'dispersion'
    fam.type <- unlist(lapply(cand.set, FUN = function(i) family(i)$family))
    fam.unique <- unique(fam.type)
    if(identical(fam.unique, "gaussian")) {disp <- NULL} else{disp <- 1}
    ##poisson, binomial, and negative binomial defaults to 1 (no separate parameter for variance)
    ##gamma is treated separately


    ##check for frequency of each terms    
    #extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN = function(i) rownames(summary(i)$coefficients)) #extract model formula for each model in cand.set

    ##determine frequency of each term across models (except (Intercept) ) 
    pooled.terms <- unlist(mod_formula)
    ##remove intercept from vector
    no.int <- pooled.terms[which(pooled.terms != "(Intercept)")]
    terms.freq <- table(no.int)
    if(length(unique(terms.freq)) > 1) stop("\n\nTo compute a shrinkage version of model-averaged estimate, each term must appear with the same frequency across models\n")


    ##check whether parm is involved in interaction
    parm.inter <- c(paste(parm, ":", sep = ""), paste(":", parm, sep = ""))
    inter.check <- ifelse(attr(regexpr(parm.inter[1], mod_formula, fixed = TRUE), "match.length") == "-1" & attr(regexpr(parm.inter[2],
                                                      mod_formula, fixed = TRUE), "match.length") == "-1", 0, 1)
    if(sum(inter.check) > 0) stop("\nParameter of interest should not be involved in interaction for shrinkage version of model-averaging to be appropriate\n")


    ##compute table
    new_table <- aictab.glm(cand.set = cand.set, modnames = modnames, sort = FALSE, c.hat = c.hat,
                            second.ord = second.ord, nobs = nobs)  #recompute AIC table and associated measures
    new_table$Beta_est <- unlist(lapply(cand.set, FUN = function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
    new_table$SE <- unlist(lapply(cand.set, FUN = function(i) sqrt(diag(vcov(i, dispersion = disp)))[paste(parm)]))

    ##replace NA's with 0
    new_table$Beta_est[is.na(new_table$Beta_est)] <- 0
    new_table$SE[is.na(new_table$SE)] <- 0

    ##add a check to determine if parameter occurs in any model
    if (isTRUE(all.equal(unique(new_table$Beta_est), 0))) {stop("Parameter not found in any of the candidate models") }

    
    #if c-hat is estimated adjust the SE's by multiplying with sqrt of c-hat
    if(c.hat > 1) {
      new_table$SE <-new_table$SE*sqrt(c.hat)
    } 

    gam1 <- unlist(lapply(cand.set, FUN=function(i) family(i)$family[1]=="Gamma")) #check for gamma regression models
    #correct SE's for estimates of gamma regressions
    if(any(gam1) == TRUE)  {
      ##check for specification of gamdisp argument
      if(is.null(gamdisp)) stop("You must specify a gamma dispersion parameter with gamma generalized linear models\n")
      new_table$SE <- unlist(lapply(cand.set,
                                  FUN = function(i) sqrt(diag(vcov(i, dispersion = gamdisp)))[paste(parm)]))
    } 



    #AICc
    #compute model-averaged estimates, unconditional SE, and 95% CL
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

    zcrit <- qnorm(p = (1-conf.level)/2, lower.tail = FALSE)
    Lower_CL <- Modavg_beta - zcrit*Uncond_SE
    Upper_CL <- Modavg_beta + zcrit*Uncond_SE
    out.modavg <- list("Parameter" = paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta, "Uncond.SE" = Uncond_SE,
                       "Conf.level" = conf.level, "Lower.CL"= Lower_CL, "Upper.CL" = Upper_CL)
  } else {stop("This function is only appropriate with either \'lm\' or \'glm\' classes\n")}

  class(out.modavg) <- c("modavg.shrink", "list")
  return(out.modavg)

}

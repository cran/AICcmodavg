modavg.effect.unmarked <- function(cand.set, modnames, newdata, type = "response", c.hat = 1, conf.level = 0.95,
                                   second.ord = TRUE, nobs = NULL, uncond.se = "revised", parm.type = NULL) {

  ##check model types
  mod.class <- unlist(lapply(cand.set, FUN = function(i) class(i)[1]))
  mod.type <- unique(mod.class)

  if(length(mod.type) > 1) stop("\nThis function is not appropriate to model-average effect sizes from different model types\n")

  ##check for supported mod.type
  supp.class <- c("unmarkedFitOccu", "unmarkedFitColExt", "unmarkedFitOccuRN", "unmarkedFitPCount", "unmarkedFitPCO",
                  "unmarkedFitDS", "unmarkedFitGDS")
                  
  if(!any(supp.class == mod.type)) {stop("\nFunction not yet defined for this object class\n")}
  
  ##rename values according to unmarked to extract from object
  ##psi
  if(identical(parm.type, "psi")) {
    if(identical(mod.type, "unmarkedFitOccu")) {parm.type1 <- "state"; parm.id <- "psi"}
    if(identical(mod.type, "unmarkedFitColExt")) {parm.type1 <- "psi"; parm.id <- "psi"}
  }

  ##gamma
  if(identical(parm.type, "gamma")) {
    if(identical(mod.type, "unmarkedFitColExt")) {parm.type1 <- "col"; parm.id <- "col"}
    if(identical(mod.type, "unmarkedFitPCO")) {parm.type1 <- "gamma"; parm.id <- "gam"}
  }

  ##epsilon
  if(identical(parm.type, "epsilon")) {
    if(identical(mod.type, "unmarkedFitColExt")) {parm.type1 <- "ext"; parm.id <- "ext"}
  }

  ##lambda
  if(identical(parm.type, "lambda")) {
    if(identical(mod.type, "unmarkedFitPCount")) {parm.type1 <- "state"; parm.id <- "lam"}
    if(identical(mod.type, "unmarkedFitPCO")) {parm.type1 <- "lambda"; parm.id <- "lam"}
    if(identical(mod.type, "unmarkedFitOccuRN")) {parm.type1 <- "state"; parm.id <- "lam"}
    if(identical(mod.type, "unmarkedFitDS")) {parm.type1 <- "state"; parm.id <- "lam"}
    if(identical(mod.type, "unmarkedFitGDS")) {parm.type1 <- "state"; parm.id <- "lam"}
  }

  ##omega
  if(identical(parm.type, "omega")) {
    if(identical(mod.type, "unmarkedFitPCO")) {parm.type1 <- "omega"; parm.id <- "omega"}
  }

  ##detect
  if(identical(parm.type, "detect")) parm.type1 <- "det"; parm.id <- "p"
  if(identical(mod.type, "unmarkedFitDS") && identical(parm.type, "detect")) stop("\nModel-averaging effect sizes on detection not yet supported for unmarkedFitDS class\n")
  if(identical(mod.type, "unmarkedFitGDS") && identical(parm.type, "detect")) stop("\nModel-averaging effect sizes on detection not yet supported for unmarkedFitGDS class\n")

  ##availability
  if(identical(parm.type, "phi")) parm.type1 <- "phi"
  if(identical(mod.type, "unmarkedFitGDS") && identical(parm.type, "phi")) stop("\nModel-averaging effect sizes on availability not yet supported for unmarkedFitGDS class\n")

    
     
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  ##check on newdata
  ##determine number of observations in new data set
  nobserv <- nrow(newdata)
  if(nobserv > 2) stop("\nCurrent maximum number of groups compared is 2:\nmodify newdata argument accordingly\n")

  ##determine number of columns in new data set
  ncolumns <- ncol(newdata)

  ##if only 1 column, add an additional column to avoid problems in computation
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
  
  ##determine which column varies
  uniques <- apply(X = newdata, MARGIN = 2, FUN = unique)
  lengths <- lapply(X = uniques, FUN = length)
  singles <- sapply(X = lengths, FUN = function(i) i > 1)
  if(sum(singles) != 1) stop("\nAll columns in 'newdata' should be constant, except the group variable\n")

  ##extract name of column
  var.id <- names(singles)[which(singles == TRUE)]
  
  ##determine name of groups compared
  group1 <- as.character(newdata[,paste(var.id)][1])
  group2 <- as.character(newdata[,paste(var.id)][2])
  
  
  ##number of models
  nmods <- length(modnames)


  ##compute predicted values
  ##point estimate
  if(identical(type, "response")) {
    ##extract fitted value for observation obs
    fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                         type = parm.type1)$Predicted)),
                  nrow = nmods, ncol = 2, byrow = TRUE)

    ##extract SE for fitted value for observation obs
    SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                        type = parm.type1)$SE)),
                 nrow = nmods, ncol = 2, byrow = TRUE)
  }

  ##link scale
  if(identical(type, "link")) {
    ##extract fitted value for observation obs
    fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                       type = parm.type1, backTransform = FALSE)$Predicted)),
                  nrow = nmods, ncol = 2, byrow = TRUE)

    ##extract SE for fitted value for observation obs
    SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                      type = parm.type1, backTransform = FALSE)$SE)),
                 nrow = nmods, ncol = 2, byrow = TRUE)
  }  

  
  ##difference between groups 
  diff <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.diff <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  ##store AICc table
  AICctab <- aictab.unmarked(cand.set = cand.set, modnames = modnames, c.hat = c.hat, second.ord = second.ord,
                             nobs = nobs, sort = FALSE)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")



  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
     
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$diff <- diff
    AICctmp$SE.diff <- SE.diff
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$diff)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.diff^2 + (AICctmp$diff - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
      QAICctmp$diff <- diff
      QAICctmp$SE.diff <- SE.diff * sqrt(c.hat)

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$diff)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.diff^2 + (QAICctmp$diff - Mod.avg.out[, 1])^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.diff^2 + (QAICctmp$diff - Mod.avg.out[, 1])^2)))  
      }
      ##store table
      AICc.out <- QAICctmp
    
    }
  
  
  
  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord == FALSE && c.hat == 1) {
    AICtmp <- AICctab
    AICtmp$diff <- diff
    AICtmp$SE.diff <- SE.diff

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$diff)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.diff^2 + (AICtmp$diff - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }

  
  ##begin loop - QAICc
  if(second.ord == FALSE && c.hat > 1){
      
    ##create temporary data.frame to store fitted values and SE 
    QAICtmp <- AICctab
    QAICtmp$diff <- diff
    QAICtmp$SE.diff <- SE.diff* sqrt(c.hat)
      
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(QAICtmp$QAICWt*QAICtmp$diff)
      
    ##compute unconditional SE and store in output matrix
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.diff^2 + (QAICtmp$diff - Mod.avg.out[, 1])^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.diff^2 + (QAICtmp$diff - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- QAICtmp
  }

  Group.variable <- paste(parm.id, "(", var.id, ")")
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]
  Mod.eff.list <- list("Group.variable" = Group.variable, "Group1" = group1,
                       "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out,
                       "Mod.avg.eff" = Mod.avg.out[,1], "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level,
                       "Lower.CL" = Lower.CL, "Upper.CL" = Upper.CL)
  class(Mod.eff.list) <- c("modavg.effect", "list")
  return(Mod.eff.list)  
}

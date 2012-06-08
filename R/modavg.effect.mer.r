modavg.effect.mer <- function(cand.set, modnames, newdata, type = "response", c.hat = 1, conf.level = 0.95,
                              second.ord = TRUE, nobs = NULL, uncond.se = "revised") {
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  if(c.hat != 1) {warning("\nThis function only allows \'c.hat = 1\' for \'mer\' class objects\n")}

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
  
    
  ##extract classes
  mod.class <- unlist(lapply(X=cand.set, FUN=class))
  ##check if all are identical
  check.class <- unique(mod.class)

  ##check that link function is the same for all models if linear predictor is used
  if(identical(type, "link")) {
    link.list <- unlist(lapply(X = cand.set, FUN = function(i) fam.link.mer(i)$link))
    check.link <- unique(link.list)
    if(length(check.link) > 1) stop("\nIt is not appropriate to compute a model-averaged beta estimate\n",
                                        "from models using different link functions\n")
  }

 
  if(!identical(check.class, "mer")) {stop("\nThis function is only appropriate with \'mer\' classes\n")}

       
  ##number of models
  nmods <- length(modnames)

    
  ##compute fitted values
  fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predictSE.mer(i, se.fit = TRUE, newdata = newdata, type = type)$fit)),
                nrow = nmods, ncol = 2, byrow = TRUE)
    
  ##compute SE's on fitted values
  SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predictSE.mer(i, se.fit = TRUE, newdata = newdata, type = type)$se.fit)),
               nrow = nmods, ncol = 2, byrow = TRUE)
  
  ##difference between groups 
  diff <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.diff <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, c.hat = c.hat, second.ord = second.ord,
                    nobs = nobs, sort = FALSE)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")


  ##begin loop - AICc
  if(second.ord==TRUE && c.hat==1){
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
  

  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord==FALSE && c.hat==1) {
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
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]
  Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                        "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1],
                        "Uncond.se" = Mod.avg.out[,2],  "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                        "Upper.CL" = Upper.CL)
  class(Mod.eff.list) <- c("modavg.effect", "list")
  return(Mod.eff.list)
}

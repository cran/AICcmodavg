##generic
modavgEffect <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                         nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                         ...){
  cand.set <- formatCands(cand.set)
  UseMethod("modavgEffect", cand.set)
}



##default
modavgEffect.default <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                 nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                 ...){
  stop("\nFunction not yet defined for this object class\n")
}



##aov
modavgEffect.AICaov.lm <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                   nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                   ...){
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
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
  varies <- sapply(X = lengths, FUN = function(i) i > 1)
##########################################
    ##CHANGES: add case when only a single variable appears in data frame
    if(ncol(newdata) == 1) {
        varies <- 1
    }

    ##add extractX to check that variables appearing in model also appear in data frame
    ##checkVariables <- extractX(cand.set)
    ##if(any(!checkVariables$predictors %in% names(newdata))) {
    ##    stop("\nAll predictors must appear in the 'newdata' data frame\n")
    ##}
##########################################
    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
    
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

  
    
  ##number of models
  nmods <- length(modnames)

  ##compute fitted values
  fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata)$fit)),
                nrow = nmods, ncol = 2, byrow = TRUE)
  
  ##compute SE's on fitted values
  SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata)$se.fit)),
               nrow = nmods, ncol = 2, byrow = TRUE)
    
  ##difference between groups 
  differ <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.differ <- sqrt(SE[, 1]^2 + SE[, 2]^2)

  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs, sort = FALSE)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  ##colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

  ##begin loop - AICc
  if(second.ord == TRUE){
                   
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$differ <- differ
    AICctmp$SE.differ <- SE.differ
      
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$differ)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2))
    }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2)))
      }
      ##store table
      AICc.out <- AICctmp
    }
  


    ##create temporary data.frame to store fitted values and SE - AIC
    if(second.ord==FALSE) {
      
      AICtmp <- AICctab
      AICtmp$differ <- differ
      AICtmp$SE.differ <- SE.differ

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$differ)

      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2)))
      }  
      ##store table
      AICc.out <- AICtmp
    }

    ##indicate scale of predictions
    type <- "response"

    zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
    Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
    Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]

    ##arrange in matrix
    predsOutMat <- matrix(data = c(Mod.avg.out[, 1], Mod.avg.out[, 2],
                                   Lower.CL, Upper.CL),
                          nrow = 1, ncol = 4)
    colnames(predsOutMat) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
    rownames(predsOutMat) <- "effect.size"

    Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                         "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1], 
                         "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                         "Upper.CL" = Upper.CL, "Matrix.output" = predsOutMat)
    class(Mod.eff.list) <- c("modavgEffect", "list")
    return(Mod.eff.list)
}



##glm
modavgEffect.AICglm.lm <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                   nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                   type = "response", c.hat = 1, gamdisp = NULL,
                                   ...){
  
  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

  
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
  varies <- sapply(X = lengths, FUN = function(i) i > 1)
##########################################
    ##CHANGES: add case when only a single variable appears in data frame
    if(ncol(newdata) == 1) {
        varies <- 1
    }

    ##add extractX to check that variables appearing in model also appear in data frame
    ##checkVariables <- extractX(cand.set)
    ##if(any(!checkVariables$predictors %in% names(newdata))) {
    ##    stop("\nAll predictors must appear in the 'newdata' data frame\n")
    ##}
##########################################

  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

  
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  if(type == "terms") {stop("\nThe terms argument is not defined for this function\n")}

  ##check family of glm to avoid problems when requesting predictions with argument 'dispersion'
  fam.type <- unlist(lapply(cand.set, FUN=function(i) family(i)$family))
  fam.unique <- unique(fam.type)
  if(identical(fam.unique, "gaussian")) {
    dispersion <- NULL  #set to NULL if gaussian is used
  } else{dispersion <- c.hat}
  ##poisson and binomial defaults to 1 (no separate parameter for variance)

  ##for negative binomial - reset to NULL
  if(any(regexpr("Negative Binomial", fam.type) != -1)) {
    dispersion <- NULL
    ##check for mixture of negative binomial and other
    ##number of models with negative binomial
    negbin.num <- sum(regexpr("Negative Binomial", fam.type) != -1)
    if(negbin.num < length(fam.type)) {
      stop("Function does not support mixture of negative binomial with other distribution")
    }
  }
  
  
    
###################CHANGES####
##############################
  if(c.hat > 1) {dispersion <- c.hat }
  if(!is.null(gamdisp)) {dispersion <- gamdisp}
  if(c.hat > 1 && !is.null(gamdisp)) {stop("\nYou cannot specify values for both \'c.hat\' and \'gamdisp\'\n")}
  ##dispersion is the dispersion parameter - this influences the SE's (to specify dispersion parameter for either overdispersed Poisson or Gamma glm)
  ##type enables to specify either "response" (original scale = point estimate) or "link" (linear predictor)
  
  ##check if object is of "lm" or "glm" class
  ##extract classes
  mod.class <- unlist(lapply(X = cand.set, FUN = class))
  ##check if all are identical
  check.class <- unique(mod.class)

  ##check that link function is the same for all models if linear predictor is used
  if(identical(type, "link")) {
    check.link <- unlist(lapply(X = cand.set, FUN = function(i) i$family$link))
    unique.link <- unique(x = check.link)
    if(length(unique.link) > 1) {stop("\nIt is not appropriate to compute a model averaged beta estimate\n",
                                          "with different link functions\n")}
  }

 
  ##check if model uses gamma distribution
  gam1 <- unlist(lapply(cand.set, FUN = function(i) family(i)$family[1] == "Gamma")) #check for gamma regression models
  ##correct SE's for estimates of gamma regressions when gamdisp is specified
  if(any(gam1) == TRUE)  {
    ##check for specification of gamdisp argument
    if(is.null(gamdisp)) stop("\nYou must specify a gamma dispersion parameter with gamma generalized linear models\n")
  }
  
  ##number of models
  nmods <- length(modnames)

  ##compute fitted values
  fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata, type = type,
                                                     dispersion = dispersion)$fit)), nrow = nmods, ncol = 2, byrow = TRUE)
    
  ##compute SE's on fitted values
  SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata, type = type,
                                                    dispersion = dispersion)$se.fit)), nrow = nmods, ncol = 2, byrow = TRUE)
  
  ##difference between groups 
  differ <- fit[, 1] - fit[, 2]
  
  ##SE on difference
  SE.differ <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  ##colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
             
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$differ <- differ
    AICctmp$SE.differ <- SE.differ
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$differ)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }
  
  ##create temporary data.frame to store fitted values and SE - QAICc
  if(second.ord==TRUE && c.hat > 1) {
    
    QAICctmp <- AICctab
    QAICctmp$differ <- differ
    QAICctmp$SE.differ <- SE.differ
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$differ)
      
    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.differ^2 + (QAICctmp$differ - Mod.avg.out[, 1])^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.differ^2 + (QAICctmp$differ - Mod.avg.out[, 1])^2)))  
    }
    ##store table
    AICc.out <- QAICctmp
    
  }

  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord == FALSE && c.hat == 1) {
    
    AICtmp <- AICctab
    AICtmp$differ <- differ
    AICtmp$SE.differ <- SE.differ

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$differ)
    
    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }
  
  ##create temporary data.frame to store fitted values and SE - QAIC
  if(second.ord == FALSE && c.hat > 1) {
          
    QAICtmp <- AICctab
    QAICtmp$differ <- differ
    QAICtmp$SE.differ <- SE.differ
      
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(QAICtmp$QAICWt*QAICtmp$differ)
      
    ##compute unconditional SE and store in output matrix
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.differ^2 + (QAICtmp$differ - Mod.avg.out[, 1])^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.differ^2 + (QAICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- QAICtmp
  }

    zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
    Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
    Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]

    ##arrange in matrix
    predsOutMat <- matrix(data = c(Mod.avg.out[, 1], Mod.avg.out[, 2],
                                   Lower.CL, Upper.CL),
                          nrow = 1, ncol = 4)
    colnames(predsOutMat) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
    rownames(predsOutMat) <- "effect.size"

    Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                         "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1], 
                         "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                         "Upper.CL" = Upper.CL, "Matrix.output" = predsOutMat)
    class(Mod.eff.list) <- c("modavgEffect", "list")
    return(Mod.eff.list)
}



##gls
modavgEffect.AICgls <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                ...){
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
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
  varies <- sapply(X = lengths, FUN = function(i) i > 1)
##########################################
    ##CHANGES: add case when only a single variable appears in data frame
    if(ncol(newdata) == 1) {
        varies <- 1
    }

    ##add extractX to check that variables appearing in model also appear in data frame
    ##checkVariables <- extractX(cand.set)
    ##if(any(!checkVariables$predictors %in% names(newdata))) {
    ##    stop("\nAll predictors must appear in the 'newdata' data frame\n")
    ##}
##########################################
    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

    
    ##number of models
    nmods <- length(modnames)

    ##compute fitted values
    fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE, newdata = newdata)$fit)),
                  nrow = nmods, ncol = 2, byrow = TRUE)
    
    ##compute SE's on fitted values
    SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE, newdata = newdata)$se.fit)),
                 nrow = nmods, ncol = 2, byrow = TRUE)
    
    ##difference between groups 
    differ <- fit[, 1] - fit[, 2]
    
    ##SE on difference
    SE.differ <- sqrt(SE[, 1]^2 + SE[, 2]^2)


    
    ##store AICc table
    AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs, sort = FALSE)

    #create object to hold Model-averaged estimates and unconditional SE's
    Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
    ##colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

    ##begin loop - AICc
    if(second.ord == TRUE){
                   
      ##create temporary data.frame to store fitted values and SE 
      AICctmp <- AICctab
      AICctmp$differ <- differ
      AICctmp$SE.differ <- SE.differ
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$differ)
      ##compute unconditional SE and store in output matrix
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2)))
      }
      ##store table
      AICc.out <- AICctmp
    }
  


    ##create temporary data.frame to store fitted values and SE - AIC
    if(second.ord==FALSE) {
      
      AICtmp <- AICctab
      AICtmp$differ <- differ
      AICtmp$SE.differ <- SE.differ

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$differ)

      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2)))
      }  
      ##store table
      AICc.out <- AICtmp
    }

    ##indicate scale of predictions
    type <- "response"

    zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
    Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
    Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]

    ##arrange in matrix
    predsOutMat <- matrix(data = c(Mod.avg.out[, 1], Mod.avg.out[, 2],
                                   Lower.CL, Upper.CL),
                          nrow = 1, ncol = 4)
    colnames(predsOutMat) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
    rownames(predsOutMat) <- "effect.size"

    Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                         "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1], 
                         "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                         "Upper.CL" = Upper.CL, "Matrix.output" = predsOutMat)
    class(Mod.eff.list) <- c("modavgEffect", "list")
    return(Mod.eff.list)
}



##lm
modavgEffect.AIClm <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                               nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                               ...){
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
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
  varies <- sapply(X = lengths, FUN = function(i) i > 1)
##########################################
    ##CHANGES: add case when only a single variable appears in data frame
    if(ncol(newdata) == 1) {
        varies <- 1
    }

    ##add extractX to check that variables appearing in model also appear in data frame
    ##checkVariables <- extractX(cand.set)
    ##if(any(!checkVariables$predictors %in% names(newdata))) {
    ##    stop("\nAll predictors must appear in the 'newdata' data frame\n")
    ##}
##########################################
    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

    
    ##number of models
    nmods <- length(modnames)

    ##compute fitted values
    fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata)$fit)),
                  nrow = nmods, ncol = 2, byrow = TRUE)
    
    ##compute SE's on fitted values
    SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata)$se.fit)),
                 nrow = nmods, ncol = 2, byrow = TRUE)
    
    ##difference between groups 
    differ <- fit[, 1] - fit[, 2]
    
    ##SE on difference
    SE.differ <- sqrt(SE[, 1]^2 + SE[, 2]^2)

    ##store AICc table
    AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs, sort = FALSE)

    #create object to hold Model-averaged estimates and unconditional SE's
    Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
    ##colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

    ##begin loop - AICc
    if(second.ord == TRUE){
                   
      ##create temporary data.frame to store fitted values and SE 
      AICctmp <- AICctab
      AICctmp$differ <- differ
      AICctmp$SE.differ <- SE.differ
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$differ)
      ##compute unconditional SE and store in output matrix
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2)))
      }
      ##store table
      AICc.out <- AICctmp
    }
  


    ##create temporary data.frame to store fitted values and SE - AIC
    if(second.ord==FALSE) {
      
      AICtmp <- AICctab
      AICtmp$differ <- differ
      AICtmp$SE.differ <- SE.differ

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$differ)

      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2)))
      }  
      ##store table
      AICc.out <- AICtmp
    }

    ##indicate scale of predictions
    type <- "response"

    zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
    Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
    Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]

    ##arrange in matrix
    predsOutMat <- matrix(data = c(Mod.avg.out[, 1], Mod.avg.out[, 2],
                                   Lower.CL, Upper.CL),
                          nrow = 1, ncol = 4)
    colnames(predsOutMat) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
    rownames(predsOutMat) <- "effect.size"

    Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                         "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1], 
                         "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                         "Upper.CL" = Upper.CL, "Matrix.output" = predsOutMat)
    class(Mod.eff.list) <- c("modavgEffect", "list")
    return(Mod.eff.list)
}



##lme
modavgEffect.AIClme <-
function(cand.set, modnames = NULL, newdata, second.ord = TRUE, nobs = NULL,
         uncond.se = "revised", conf.level = 0.95, ...){
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
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
  varies <- sapply(X = lengths, FUN = function(i) i > 1)
##########################################
    ##CHANGES: add case when only a single variable appears in data frame
    if(ncol(newdata) == 1) {
        varies <- 1
    }

    ##add extractX to check that variables appearing in model also appear in data frame
    ##checkVariables <- extractX(cand.set)
    ##if(any(!checkVariables$predictors %in% names(newdata))) {
    ##    stop("\nAll predictors must appear in the 'newdata' data frame\n")
    ##}
##########################################
    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }
  
    
  ##number of models
  nmods <- length(modnames)

  ##compute fitted values
  fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE, newdata = newdata)$fit)),
                nrow = nmods, ncol = 2, byrow = TRUE)
    
  ##compute SE's on fitted values
  SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE, newdata = newdata)$se.fit)),
               nrow = nmods, ncol = 2, byrow = TRUE)
    
  ##difference between groups 
  differ <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.differ <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs, sort = FALSE)
  
  
  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  ##colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")
  
  ##begin loop - AICc
  if(second.ord == TRUE){
    
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$differ <- differ
    AICctmp$SE.differ <- SE.differ
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$differ)
    ##compute unconditional SE and store in output matrix
    
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }
  
  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord==FALSE) {
    
    AICtmp <- AICctab
    AICtmp$differ <- differ
    AICtmp$SE.differ <- SE.differ

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$differ)
    
    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }
  
  ##indicate scale of predictions
  type <- "response"
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]

    ##arrange in matrix
    predsOutMat <- matrix(data = c(Mod.avg.out[, 1], Mod.avg.out[, 2],
                                   Lower.CL, Upper.CL),
                          nrow = 1, ncol = 4)
    colnames(predsOutMat) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
    rownames(predsOutMat) <- "effect.size"

    Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                         "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1], 
                         "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                         "Upper.CL" = Upper.CL, "Matrix.output" = predsOutMat)  
    class(Mod.eff.list) <- c("modavgEffect", "list")
    return(Mod.eff.list)
}



##mer  - lme4 version < 1
modavgEffect.AICmer <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                type = "response", ...) {

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
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
  varies <- sapply(X = lengths, FUN = function(i) i > 1)
##########################################
    ##CHANGES: add case when only a single variable appears in data frame
    if(ncol(newdata) == 1) {
        varies <- 1
    }

    ##add extractX to check that variables appearing in model also appear in data frame
    ##checkVariables <- extractX(cand.set)
    ##if(any(!checkVariables$predictors %in% names(newdata))) {
    ##    stop("\nAll predictors must appear in the 'newdata' data frame\n")
    ##}
##########################################
   
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

    
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

 
  ##number of models
  nmods <- length(modnames)

    
  ##compute fitted values
  fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE, newdata = newdata, type = type)$fit)),
                nrow = nmods, ncol = 2, byrow = TRUE)
    
  ##compute SE's on fitted values
  SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE, newdata = newdata, type = type)$se.fit)),
               nrow = nmods, ncol = 2, byrow = TRUE)
  
  ##difference between groups 
  differ <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.differ <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  ##colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")


  ##begin loop - AICc
  if(second.ord==TRUE){
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$differ <- differ
    AICctmp$SE.differ <- SE.differ
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$differ)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }
  

  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord==FALSE) {
    AICtmp <- AICctab
    AICtmp$differ <- differ
    AICtmp$SE.differ <- SE.differ

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$differ)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]

    ##arrange in matrix
    predsOutMat <- matrix(data = c(Mod.avg.out[, 1], Mod.avg.out[, 2],
                                   Lower.CL, Upper.CL),
                          nrow = 1, ncol = 4)
    colnames(predsOutMat) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
    rownames(predsOutMat) <- "effect.size"

    Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                         "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1], 
                         "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                         "Upper.CL" = Upper.CL, "Matrix.output" = predsOutMat)  
    class(Mod.eff.list) <- c("modavgEffect", "list")
    return(Mod.eff.list)
}



##glmerMod
modavgEffect.AICglmerMod <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                     nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                     type = "response", ...) {

  ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
        if(is.null(names(cand.set))) {
            modnames <- paste("Mod", 1:length(cand.set), sep = "")
            warning("\nModel names have been supplied automatically in the table\n")
        } else {
            modnames <- names(cand.set)
        }
    }
    
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
  varies <- sapply(X = lengths, FUN = function(i) i > 1)
##########################################
    ##CHANGES: add case when only a single variable appears in data frame
    if(ncol(newdata) == 1) {
        varies <- 1
    }

    ##add extractX to check that variables appearing in model also appear in data frame
    ##checkVariables <- extractX(cand.set)
    ##if(any(!checkVariables$predictors %in% names(newdata))) {
    ##    stop("\nAll predictors must appear in the 'newdata' data frame\n")
    ##}
##########################################
    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

    
    
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

 
  ##number of models
  nmods <- length(modnames)

    
  ##compute fitted values
  fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE, newdata = newdata, type = type)$fit)),
                nrow = nmods, ncol = 2, byrow = TRUE)
    
  ##compute SE's on fitted values
  SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE, newdata = newdata, type = type)$se.fit)),
               nrow = nmods, ncol = 2, byrow = TRUE)
  
  ##difference between groups 
  differ <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.differ <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  ##colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")


  ##begin loop - AICc
  if(second.ord==TRUE){
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$differ <- differ
    AICctmp$SE.differ <- SE.differ
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$differ)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }
  

  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord==FALSE) {
    AICtmp <- AICctab
    AICtmp$differ <- differ
    AICtmp$SE.differ <- SE.differ

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$differ)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]

    ##arrange in matrix
    predsOutMat <- matrix(data = c(Mod.avg.out[, 1], Mod.avg.out[, 2],
                                   Lower.CL, Upper.CL),
                          nrow = 1, ncol = 4)
    colnames(predsOutMat) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
    rownames(predsOutMat) <- "effect.size"

    Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                         "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1], 
                         "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                         "Upper.CL" = Upper.CL, "Matrix.output" = predsOutMat)
    class(Mod.eff.list) <- c("modavgEffect", "list")
    return(Mod.eff.list)
}



##lmerMod
modavgEffect.AIClmerMod <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                    nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                    ...) {

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
        if(is.null(names(cand.set))) {
            modnames <- paste("Mod", 1:length(cand.set), sep = "")
            warning("\nModel names have been supplied automatically in the table\n")
        } else {
            modnames <- names(cand.set)
        }
  }
    
    
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
  varies <- sapply(X = lengths, FUN = function(i) i > 1)
##########################################
    ##CHANGES: add case when only a single variable appears in data frame
    if(ncol(newdata) == 1) {
        varies <- 1
    }

    ##add extractX to check that variables appearing in model also appear in data frame
    ##checkVariables <- extractX(cand.set)
    ##if(any(!checkVariables$predictors %in% names(newdata))) {
    ##    stop("\nAll predictors must appear in the 'newdata' data frame\n")
    ##}
##########################################
    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

    
  ##extract classes
  mod.class <- unlist(lapply(X=cand.set, FUN=class))
  ##check if all are identical
  check.class <- unique(mod.class)

  ##number of models
  nmods <- length(modnames)
  
    
  ##compute fitted values
  fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE, newdata = newdata)$fit)),
                nrow = nmods, ncol = 2, byrow = TRUE)
    
  ##compute SE's on fitted values
  SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE, newdata = newdata)$se.fit)),
               nrow = nmods, ncol = 2, byrow = TRUE)
  
  ##difference between groups 
  differ <- fit[, 1] - fit[, 2]
  
  ##SE on difference
  SE.differ <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  ##colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

  ##begin loop - AICc
  if(second.ord==TRUE){
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$differ <- differ
    AICctmp$SE.differ <- SE.differ
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$differ)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }
  
  
  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord==FALSE) {
    AICtmp <- AICctab
    AICtmp$differ <- differ
    AICtmp$SE.differ <- SE.differ

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$differ)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }

  ##scale of predictions
  type <- "response"
    
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]

    ##arrange in matrix
    predsOutMat <- matrix(data = c(Mod.avg.out[, 1], Mod.avg.out[, 2],
                                   Lower.CL, Upper.CL),
                          nrow = 1, ncol = 4)
    colnames(predsOutMat) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
    rownames(predsOutMat) <- "effect.size"

    Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                         "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1], 
                         "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                         "Upper.CL" = Upper.CL, "Matrix.output" = predsOutMat)
    class(Mod.eff.list) <- c("modavgEffect", "list")
    return(Mod.eff.list)
}



##lmerModLmerTest
modavgEffect.AIClmerModLmerTest <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                            nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                            ...) {

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
        if(is.null(names(cand.set))) {
            modnames <- paste("Mod", 1:length(cand.set), sep = "")
            warning("\nModel names have been supplied automatically in the table\n")
        } else {
            modnames <- names(cand.set)
        }
  }
    
    
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
  varies <- sapply(X = lengths, FUN = function(i) i > 1)
##########################################
    ##CHANGES: add case when only a single variable appears in data frame
    if(ncol(newdata) == 1) {
        varies <- 1
    }

    ##add extractX to check that variables appearing in model also appear in data frame
    ##checkVariables <- extractX(cand.set)
    ##if(any(!checkVariables$predictors %in% names(newdata))) {
    ##    stop("\nAll predictors must appear in the 'newdata' data frame\n")
    ##}
##########################################
    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

    
  ##extract classes
  mod.class <- unlist(lapply(X=cand.set, FUN=class))
  ##check if all are identical
  check.class <- unique(mod.class)

  ##number of models
  nmods <- length(modnames)
  
    
  ##compute fitted values
  fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE, newdata = newdata)$fit)),
                nrow = nmods, ncol = 2, byrow = TRUE)
    
  ##compute SE's on fitted values
  SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE, newdata = newdata)$se.fit)),
               nrow = nmods, ncol = 2, byrow = TRUE)
  
  ##difference between groups 
  differ <- fit[, 1] - fit[, 2]
  
  ##SE on difference
  SE.differ <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  ##colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

  ##begin loop - AICc
  if(second.ord==TRUE){
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$differ <- differ
    AICctmp$SE.differ <- SE.differ
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$differ)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }
  
  
  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord==FALSE) {
    AICtmp <- AICctab
    AICtmp$differ <- differ
    AICtmp$SE.differ <- SE.differ

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$differ)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }

  ##scale of predictions
  type <- "response"
    
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]
    
    ##arrange in matrix
    predsOutMat <- matrix(data = c(Mod.avg.out[, 1], Mod.avg.out[, 2],
                                   Lower.CL, Upper.CL),
                          nrow = 1, ncol = 4)
    colnames(predsOutMat) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
    rownames(predsOutMat) <- "effect.size"

    Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                         "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1], 
                         "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                         "Upper.CL" = Upper.CL, "Matrix.output" = predsOutMat)
    class(Mod.eff.list) <- c("modavgEffect", "list")
    return(Mod.eff.list)
}



##glm.nb
modavgEffect.AICnegbin.glm.lm <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                          nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                          type = "response", ...){
  
    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
        if(is.null(names(cand.set))) {
            modnames <- paste("Mod", 1:length(cand.set), sep = "")
            warning("\nModel names have been supplied automatically in the table\n")
        } else {
            modnames <- names(cand.set)
        }
    }


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
    varies <- sapply(X = lengths, FUN = function(i) i > 1)
##########################################
    ##CHANGES: add case when only a single variable appears in data frame
    if(ncol(newdata) == 1) {
        varies <- 1
    }

    ##add extractX to check that variables appearing in model also appear in data frame
    ##checkVariables <- extractX(cand.set)
    ##if(any(!checkVariables$predictors %in% names(newdata))) {
    ##    stop("\nAll predictors must appear in the 'newdata' data frame\n")
    ##}
##########################################

    ##extract name of column
    if(sum(varies) == 1) {
        var.id <- names(varies)[which(varies == TRUE)]
                
        ##determine name of groups compared
        group1 <- as.character(newdata[,paste(var.id)][1])
        group2 <- as.character(newdata[,paste(var.id)][2])

    } else {
        ##warn that no single variable defines groups
        warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
        ##use generic names
        var.id <- "Groups"
        group1 <- "group 1"
        group2 <- "group 2"
    }

  
    ##newdata is data frame with exact structure of the original data frame (same variable names and type)
    if(type == "terms") {stop("\nThe terms argument is not defined for this function\n")}
  
    
###################CHANGES####
##############################
    ##type enables to specify either "response" (original scale = point estimate) or "link" (linear predictor)
  
    ##check if object is of "lm" or "glm" class
    ##extract classes
    mod.class <- unlist(lapply(X = cand.set, FUN = class))
    ##check if all are identical
    check.class <- unique(mod.class)

    ##check that link function is the same for all models if linear predictor is used
    if(identical(type, "link")) {
        check.link <- unlist(lapply(X = cand.set, FUN = function(i) i$family$link))
        unique.link <- unique(x = check.link)
        if(length(unique.link) > 1) {stop("\nIt is not appropriate to compute a model averaged beta estimate\n",
                                          "with different link functions\n")}
    }

 
    ##number of models
    nmods <- length(modnames)
    
    ##compute fitted values
    fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                                              type = type)$fit)),
                  nrow = nmods, ncol = 2, byrow = TRUE)
    
    ##compute SE's on fitted values
    SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                                             type = type)$se.fit)),
                 nrow = nmods, ncol = 2, byrow = TRUE)
    
    ##difference between groups 
    differ <- fit[, 1] - fit[, 2]
  
    ##SE on difference
    SE.differ <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
    ##store AICc table
    AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                      nobs = nobs, sort = FALSE)

    ##create object to hold Model-averaged estimates and unconditional SE's
    Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
    ##colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

    ##begin loop - AICc
    if(second.ord == TRUE){
             
        ##create temporary data.frame to store fitted values and SE 
        AICctmp <- AICctab
        AICctmp$differ <- differ
        AICctmp$SE.differ <- SE.differ
    
        ##compute model averaged prediction and store in output matrix
        Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$differ)
        ##compute unconditional SE and store in output matrix
      
        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
        if(identical(uncond.se, "old")) {
            Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2))
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
            Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2)))
        }
        ##store table
        AICc.out <- AICctmp
    }
    

    ##create temporary data.frame to store fitted values and SE - AIC
    if(second.ord == FALSE) {
    
        AICtmp <- AICctab
        AICtmp$differ <- differ
        AICtmp$SE.differ <- SE.differ

        ##compute model averaged prediction and store in output matrix
        Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$differ)
    
        ##compute unconditional SE and store in output matrix
        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
        if(identical(uncond.se, "old")) {
            Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2))
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
            Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2)))
        }  
        ##store table
        AICc.out <- AICtmp
    }
  
    zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
    Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
    Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]

    ##arrange in matrix
    predsOutMat <- matrix(data = c(Mod.avg.out[, 1], Mod.avg.out[, 2],
                                   Lower.CL, Upper.CL),
                          nrow = 1, ncol = 4)
    colnames(predsOutMat) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
    rownames(predsOutMat) <- "effect.size"

    Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                         "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1], 
                         "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                         "Upper.CL" = Upper.CL, "Matrix.output" = predsOutMat)
    class(Mod.eff.list) <- c("modavgEffect", "list")
    return(Mod.eff.list)
}



##rlm
modavgEffect.AICrlm.lm <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                   nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                   ...){
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

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
  varies <- sapply(X = lengths, FUN = function(i) i > 1)
##########################################
    ##CHANGES: add case when only a single variable appears in data frame
    if(ncol(newdata) == 1) {
        varies <- 1
    }

    ##add extractX to check that variables appearing in model also appear in data frame
    ##checkVariables <- extractX(cand.set)
    ##if(any(!checkVariables$predictors %in% names(newdata))) {
    ##    stop("\nAll predictors must appear in the 'newdata' data frame\n")
    ##}
##########################################
    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

  
  ##number of models
  nmods <- length(modnames)
    
  ##compute fitted values
  fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE, newdata = newdata)$fit)),
                nrow = nmods, ncol = 2, byrow = TRUE)
    
  ##compute SE's on fitted values
  SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE, newdata = newdata)$se.fit)),
               nrow = nmods, ncol = 2, byrow = TRUE)
    
  ##difference between groups 
  differ <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.differ <- sqrt(SE[, 1]^2 + SE[, 2]^2)
    
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs, sort = FALSE)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  ##colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")
  
  ##begin loop - AICc
  if(second.ord == TRUE){
                   
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$differ <- differ
    AICctmp$SE.differ <- SE.differ
      
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$differ)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2))
    }
    
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }
  
  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord==FALSE) {
      
    AICtmp <- AICctab
    AICtmp$differ <- differ
    AICtmp$SE.differ <- SE.differ

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$differ)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2))
    }
    
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }
  
  ##indicate scale of predictions
  type <- "response"

  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]

    ##arrange in matrix
    predsOutMat <- matrix(data = c(Mod.avg.out[, 1], Mod.avg.out[, 2],
                                   Lower.CL, Upper.CL),
                          nrow = 1, ncol = 4)
    colnames(predsOutMat) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
    rownames(predsOutMat) <- "effect.size"

    Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                         "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1], 
                         "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                         "Upper.CL" = Upper.CL, "Matrix.output" = predsOutMat)
    class(Mod.eff.list) <- c("modavgEffect", "list")
    return(Mod.eff.list)
}



##survreg
modavgEffect.AICsurvreg <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                    nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                    type = "response", ...){
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

  ##check that distribution is the same for all models
  if(identical(type, "link")) {
    check.dist <- sapply(X = cand.set, FUN = function(i) i$dist)
    unique.dist <- unique(x = check.dist)
    if(length(unique.dist) > 1) stop("\nFunction does not support model-averaging effect size on link scale using different distributions\n")
  }
  
  
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
  varies <- sapply(X = lengths, FUN = function(i) i > 1)
##########################################
    ##CHANGES: add case when only a single variable appears in data frame
    if(ncol(newdata) == 1) {
        varies <- 1
    }

    ##add extractX to check that variables appearing in model also appear in data frame
    ##checkVariables <- extractX(cand.set)
    ##if(any(!checkVariables$predictors %in% names(newdata))) {
    ##    stop("\nAll predictors must appear in the 'newdata' data frame\n")
    ##}
##########################################

    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }
    
    
    ##number of models
    nmods <- length(modnames)

    ##compute fitted values
    fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata, type = type)$fit)),
                  nrow = nmods, ncol = 2, byrow = TRUE)
    
    ##compute SE's on fitted values
    SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata, type = type)$se.fit)),
                 nrow = nmods, ncol = 2, byrow = TRUE)
    
    ##difference between groups 
    differ <- fit[, 1] - fit[, 2]
    
    ##SE on difference
    SE.differ <- sqrt(SE[, 1]^2 + SE[, 2]^2)

    ##store AICc table
    AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs, sort = FALSE)

    #create object to hold Model-averaged estimates and unconditional SE's
    Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
    ##colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

    ##begin loop - AICc
    if(second.ord == TRUE){
                   
      ##create temporary data.frame to store fitted values and SE 
      AICctmp <- AICctab
      AICctmp$differ <- differ
      AICctmp$SE.differ <- SE.differ
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$differ)
      ##compute unconditional SE and store in output matrix
      
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2)))
      }
      ##store table
      AICc.out <- AICctmp
    }
  


    ##create temporary data.frame to store fitted values and SE - AIC
    if(second.ord==FALSE) {
      
      AICtmp <- AICctab
      AICtmp$differ <- differ
      AICtmp$SE.differ <- SE.differ

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$differ)

      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2)))
      }  
      ##store table
      AICc.out <- AICtmp
    }

    zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
    Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
    Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]

    ##arrange in matrix
    predsOutMat <- matrix(data = c(Mod.avg.out[, 1], Mod.avg.out[, 2],
                                   Lower.CL, Upper.CL),
                          nrow = 1, ncol = 4)
    colnames(predsOutMat) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
    rownames(predsOutMat) <- "effect.size"

    Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                         "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1], 
                         "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                         "Upper.CL" = Upper.CL, "Matrix.output" = predsOutMat)
    
    class(Mod.eff.list) <- c("modavgEffect", "list")
    return(Mod.eff.list)
}



##occu
modavgEffect.AICunmarkedFitOccu <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                            nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                            type = "response", c.hat = 1, parm.type = NULL,
                                            ...) {

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgEffect for details\n")}

    
  ##rename values according to unmarked to extract from object
  ##psi
  if(identical(parm.type, "psi")) {
    parm.type1 <- "state"; parm.id <- "psi"
    
  }

    ##detect
    if(identical(parm.type, "detect")) {parm.type1 <- "det"; parm.id <- "p"}


##################
    ##extract link function
    check.link <- sapply(X = cand.set, FUN = function(i) eval(parse(text = paste("i@estimates@estimates$",
                                                                                 parm.type1, "@invlink",
                                                                                 sep = ""))))
    unique.link <- unique(check.link)
    select.link <- unique.link[1]
    
    if(identical(type, "link")) {
        if(length(unique.link) > 1) {stop("\nIt is not appropriate to compute a model averaged linear predictor\n",
                                          "with different link functions\n")}
    }
##################    
    
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
  varies <- sapply(X = lengths, FUN = function(i) i > 1)
##########################################
    ##CHANGES: add case when only a single variable appears in data frame
    if(ncol(newdata) == 1) {
        varies <- 1
    }

    ##add extractX to check that variables appearing in model also appear in data frame
    ##checkVariables <- extractX(cand.set, parm.type = parm.type)
    ##if(any(!checkVariables$predictors %in% names(newdata))) {
    ##    stop("\nAll predictors must appear in the 'newdata' data frame\n")
    ##}
##########################################
    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

  
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
  differ <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.differ <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  ##colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
     
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$differ <- differ
    AICctmp$SE.differ <- SE.differ
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$differ)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
      QAICctmp$differ <- differ
      QAICctmp$SE.differ <- SE.differ * sqrt(c.hat)

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$differ)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.differ^2 + (QAICctmp$differ - Mod.avg.out[, 1])^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.differ^2 + (QAICctmp$differ - Mod.avg.out[, 1])^2)))  
      }
      ##store table
      AICc.out <- QAICctmp
    
    }
  
  
  
  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord == FALSE && c.hat == 1) {
    AICtmp <- AICctab
    AICtmp$differ <- differ
    AICtmp$SE.differ <- SE.differ

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$differ)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }

  
  ##begin loop - QAICc
  if(second.ord == FALSE && c.hat > 1){
      
    ##create temporary data.frame to store fitted values and SE 
    QAICtmp <- AICctab
    QAICtmp$differ <- differ
    QAICtmp$SE.differ <- SE.differ* sqrt(c.hat)
      
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(QAICtmp$QAICWt*QAICtmp$differ)
      
    ##compute unconditional SE and store in output matrix
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.differ^2 + (QAICtmp$differ - Mod.avg.out[, 1])^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.differ^2 + (QAICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- QAICtmp
  }

  Group.variable <- paste(parm.id, "(", var.id, ")")
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]

    ##arrange in matrix
    predsOutMat <- matrix(data = c(Mod.avg.out[, 1], Mod.avg.out[, 2],
                                   Lower.CL, Upper.CL),
                          nrow = 1, ncol = 4)
    colnames(predsOutMat) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
    rownames(predsOutMat) <- "effect.size"

    Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                         "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1], 
                         "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                         "Upper.CL" = Upper.CL, "Matrix.output" = predsOutMat)

    class(Mod.eff.list) <- c("modavgEffect", "list")
    return(Mod.eff.list)  
}



##colext
modavgEffect.AICunmarkedFitColExt <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                              nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                              type = "response", c.hat = 1, parm.type = NULL,
                                              ...) {

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgEffect for details\n")}


  ##rename values according to unmarked to extract from object
  ##psi
  if(identical(parm.type, "psi")) {
    parm.type1 <- "psi"; parm.id <- "psi"
  }

  ##gamma
  if(identical(parm.type, "gamma")) {
    parm.type1 <- "col"; parm.id <- "col"
  }

  ##epsilon
  if(identical(parm.type, "epsilon")) {
    parm.type1 <- "ext"; parm.id <- "ext"
  }

  ##detect
  if(identical(parm.type, "detect")) {parm.type1 <- "det"; parm.id <- "p"}


##################
    ##extract link function
    check.link <- sapply(X = cand.set, FUN = function(i) eval(parse(text = paste("i@estimates@estimates$",
                                                                                 parm.type1, "@invlink",
                                                                                 sep = ""))))
    unique.link <- unique(check.link)
    select.link <- unique.link[1]
    
    if(identical(type, "link")) {
        if(length(unique.link) > 1) {stop("\nIt is not appropriate to compute a model averaged linear predictor\n",
                                          "with different link functions\n")}
    }
##################
    
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
  varies <- sapply(X = lengths, FUN = function(i) i > 1)
##########################################
    ##CHANGES: add case when only a single variable appears in data frame
    if(ncol(newdata) == 1) {
        varies <- 1
    }

    ##add extractX to check that variables appearing in model also appear in data frame
    ##checkVariables <- extractX(cand.set, parm.type = parm.type)
    ##if(any(!checkVariables$predictors %in% names(newdata))) {
    ##    stop("\nAll predictors must appear in the 'newdata' data frame\n")
    ##}
##########################################
    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

  
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
  differ <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.differ <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  ##colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
     
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$differ <- differ
    AICctmp$SE.differ <- SE.differ
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$differ)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
      QAICctmp$differ <- differ
      QAICctmp$SE.differ <- SE.differ * sqrt(c.hat)

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$differ)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.differ^2 + (QAICctmp$differ - Mod.avg.out[, 1])^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.differ^2 + (QAICctmp$differ - Mod.avg.out[, 1])^2)))  
      }
      ##store table
      AICc.out <- QAICctmp
    
    }
  
  
  
  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord == FALSE && c.hat == 1) {
    AICtmp <- AICctab
    AICtmp$differ <- differ
    AICtmp$SE.differ <- SE.differ

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$differ)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }

  
  ##begin loop - QAICc
  if(second.ord == FALSE && c.hat > 1){
      
    ##create temporary data.frame to store fitted values and SE 
    QAICtmp <- AICctab
    QAICtmp$differ <- differ
    QAICtmp$SE.differ <- SE.differ* sqrt(c.hat)
      
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(QAICtmp$QAICWt*QAICtmp$differ)
      
    ##compute unconditional SE and store in output matrix
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.differ^2 + (QAICtmp$differ - Mod.avg.out[, 1])^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.differ^2 + (QAICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- QAICtmp
  }

  Group.variable <- paste(parm.id, "(", var.id, ")")
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]

    ##arrange in matrix
    predsOutMat <- matrix(data = c(Mod.avg.out[, 1], Mod.avg.out[, 2],
                                   Lower.CL, Upper.CL),
                          nrow = 1, ncol = 4)
    colnames(predsOutMat) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
    rownames(predsOutMat) <- "effect.size"

    Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                         "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1], 
                         "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                         "Upper.CL" = Upper.CL, "Matrix.output" = predsOutMat)
    class(Mod.eff.list) <- c("modavgEffect", "list")
    return(Mod.eff.list)  
}



##occuRN
modavgEffect.AICunmarkedFitOccuRN <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                              nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                              type = "response", c.hat = 1, parm.type = NULL,
                                              ...) {

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgEffect for details\n")}

    
  ##rename values according to unmarked to extract from object
  ##lambda
  if(identical(parm.type, "lambda")) {
    parm.type1 <- "state"; parm.id <- "lam"
  }

  ##detect
  if(identical(parm.type, "detect")) {parm.type1 <- "det"; parm.id <- "p"}

##################
    ##extract link function
    check.link <- sapply(X = cand.set, FUN = function(i) eval(parse(text = paste("i@estimates@estimates$",
                                                                                 parm.type1, "@invlink",
                                                                                 sep = ""))))
    unique.link <- unique(check.link)
    select.link <- unique.link[1]
    
    if(identical(type, "link")) {
        if(length(unique.link) > 1) {stop("\nIt is not appropriate to compute a model averaged linear predictor\n",
                                          "with different link functions\n")}
    }
##################
    
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
  varies <- sapply(X = lengths, FUN = function(i) i > 1)
##########################################
    ##CHANGES: add case when only a single variable appears in data frame
    if(ncol(newdata) == 1) {
        varies <- 1
    }

    ##add extractX to check that variables appearing in model also appear in data frame
    ##checkVariables <- extractX(cand.set, parm.type = parm.type)
    ##if(any(!checkVariables$predictors %in% names(newdata))) {
    ##    stop("\nAll predictors must appear in the 'newdata' data frame\n")
    ##}
##########################################
    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

  
  ##number of models
  nmods <- length(modnames)


  ##compute predicted values
  ##point estimate
  if(identical(type, "response")) {
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
  differ <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.differ <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  ##colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
     
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$differ <- differ
    AICctmp$SE.differ <- SE.differ
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$differ)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
      QAICctmp$differ <- differ
      QAICctmp$SE.differ <- SE.differ * sqrt(c.hat)

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$differ)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.differ^2 + (QAICctmp$differ - Mod.avg.out[, 1])^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.differ^2 + (QAICctmp$differ - Mod.avg.out[, 1])^2)))  
      }
      ##store table
      AICc.out <- QAICctmp
    
    }
  
  
  
  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord == FALSE && c.hat == 1) {
    AICtmp <- AICctab
    AICtmp$differ <- differ
    AICtmp$SE.differ <- SE.differ

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$differ)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }

  
  ##begin loop - QAICc
  if(second.ord == FALSE && c.hat > 1){
      
    ##create temporary data.frame to store fitted values and SE 
    QAICtmp <- AICctab
    QAICtmp$differ <- differ
    QAICtmp$SE.differ <- SE.differ* sqrt(c.hat)
      
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(QAICtmp$QAICWt*QAICtmp$differ)
      
    ##compute unconditional SE and store in output matrix
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.differ^2 + (QAICtmp$differ - Mod.avg.out[, 1])^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.differ^2 + (QAICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- QAICtmp
  }

  Group.variable <- paste(parm.id, "(", var.id, ")")
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]

    ##arrange in matrix
    predsOutMat <- matrix(data = c(Mod.avg.out[, 1], Mod.avg.out[, 2],
                                   Lower.CL, Upper.CL),
                          nrow = 1, ncol = 4)
    colnames(predsOutMat) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
    rownames(predsOutMat) <- "effect.size"

    Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                         "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1], 
                         "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                         "Upper.CL" = Upper.CL, "Matrix.output" = predsOutMat)  
    class(Mod.eff.list) <- c("modavgEffect", "list")
    return(Mod.eff.list)  
}



##pcount
modavgEffect.AICunmarkedFitPCount <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                              nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                              type = "response", c.hat = 1, parm.type = NULL,
                                              ...) {

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgEffect for details\n")}

  ##rename values according to unmarked to extract from object
  ##lambda
  if(identical(parm.type, "lambda")) {
    parm.type1 <- "state"; parm.id <- "lam"
    
    ##check mixture type for mixture models
    mixture.type <- sapply(X = cand.set, FUN = function(i) i@mixture)
    unique.mixture <- unique(mixture.type)
  }
   

  ##detect
  if(identical(parm.type, "detect")) {parm.type1 <- "det"; parm.id <- "p"}
    
##################
    ##extract link function
    check.link <- sapply(X = cand.set, FUN = function(i) eval(parse(text = paste("i@estimates@estimates$",
                                                                                 parm.type1, "@invlink",
                                                                                 sep = ""))))
    unique.link <- unique(check.link)
    select.link <- unique.link[1]
    
    if(identical(type, "link")) {
        if(length(unique.link) > 1) {stop("\nIt is not appropriate to compute a model averaged linear predictor\n",
                                          "with different link functions\n")}
    }
##################
    
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
  varies <- sapply(X = lengths, FUN = function(i) i > 1)
##########################################
    ##CHANGES: add case when only a single variable appears in data frame
    if(ncol(newdata) == 1) {
        varies <- 1
    }

    ##add extractX to check that variables appearing in model also appear in data frame
    ##checkVariables <- extractX(cand.set, parm.type = parm.type)
    ##if(any(!checkVariables$predictors %in% names(newdata))) {
    ##    stop("\nAll predictors must appear in the 'newdata' data frame\n")
    ##}
##########################################
    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

  
  
  ##number of models
  nmods <- length(modnames)


  ##compute predicted values
  ##point estimate
  if(identical(type, "response")) {
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
  differ <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.differ <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  ##colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
     
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$differ <- differ
    AICctmp$SE.differ <- SE.differ
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$differ)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
      QAICctmp$differ <- differ
      QAICctmp$SE.differ <- SE.differ * sqrt(c.hat)

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$differ)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.differ^2 + (QAICctmp$differ - Mod.avg.out[, 1])^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.differ^2 + (QAICctmp$differ - Mod.avg.out[, 1])^2)))  
      }
      ##store table
      AICc.out <- QAICctmp
    
    }
  
  
  
  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord == FALSE && c.hat == 1) {
    AICtmp <- AICctab
    AICtmp$differ <- differ
    AICtmp$SE.differ <- SE.differ

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$differ)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }

  
  ##begin loop - QAICc
  if(second.ord == FALSE && c.hat > 1){
      
    ##create temporary data.frame to store fitted values and SE 
    QAICtmp <- AICctab
    QAICtmp$differ <- differ
    QAICtmp$SE.differ <- SE.differ* sqrt(c.hat)
      
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(QAICtmp$QAICWt*QAICtmp$differ)
      
    ##compute unconditional SE and store in output matrix
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.differ^2 + (QAICtmp$differ - Mod.avg.out[, 1])^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.differ^2 + (QAICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- QAICtmp
  }

  Group.variable <- paste(parm.id, "(", var.id, ")")
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]

    ##arrange in matrix
    predsOutMat <- matrix(data = c(Mod.avg.out[, 1], Mod.avg.out[, 2],
                                   Lower.CL, Upper.CL),
                          nrow = 1, ncol = 4)
    colnames(predsOutMat) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
    rownames(predsOutMat) <- "effect.size"

    Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                         "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1], 
                         "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                         "Upper.CL" = Upper.CL, "Matrix.output" = predsOutMat)
    class(Mod.eff.list) <- c("modavgEffect", "list")
    return(Mod.eff.list)  
}



##unmarkedFitPCO
modavgEffect.AICunmarkedFitPCO <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                           type = "response", c.hat = 1, parm.type = NULL,
                                           ...) {
  
  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgEffect for details\n")}

  ##rename values according to unmarked to extract from object
  ##gamma
  if(identical(parm.type, "gamma")) {
    parm.type1 <- "gamma"; parm.id <- "gam"
  }

  ##lambda
  if(identical(parm.type, "lambda")) {
    parm.type1 <- "lambda"; parm.id <- "lam"
    
    ##check mixture type for mixture models
    mixture.type <- sapply(X = cand.set, FUN = function(i) i@mixture)
    unique.mixture <- unique(mixture.type)
    if(length(unique.mixture) > 1) {
      if(any(unique.mixture == "ZIP")) stop("\nThis function does not yet support mixing ZIP with other distributions\n")
    } else {
      mixture.id <- unique(mixture.type)
      if(identical(unique.mixture, "ZIP")) {
        if(identical(type, "link")) stop("\nLink scale not yet supported for ZIP mixtures\n")
      }
    }
  }
  

  ##omega
  if(identical(parm.type, "omega")) {
    parm.type1 <- "omega"; parm.id <- "omega"
  }

  ##iota (for immigration = TRUE with dynamics = "autoreg", "trend", "ricker", or "gompertz")
  if(identical(parm.type, "iota")) {
      parm.type1 <- "iota"; parm.id <- "iota"
      ##check that parameter appears in all models
      parfreq <- sum(sapply(cand.set, FUN = function(i) any(names(i@estimates@estimates) == parm.type1)))
      if(!identical(length(cand.set), parfreq)) {
          stop("\nParameter \'", parm.type1, "\' (parm.type = \"", parm.type, "\") does not appear in all models:",
               "\ncannot compute model-averaged effect size across all models\n")
      }
  }

  ##detect
  if(identical(parm.type, "detect")) {parm.type1 <- "det"; parm.id <- "p"}

##################
    ##extract link function
    check.link <- sapply(X = cand.set, FUN = function(i) eval(parse(text = paste("i@estimates@estimates$",
                                                                                 parm.type1, "@invlink",
                                                                                 sep = ""))))
    unique.link <- unique(check.link)
    select.link <- unique.link[1]
    
    if(identical(type, "link")) {
        if(length(unique.link) > 1) {stop("\nIt is not appropriate to compute a model averaged linear predictor\n",
                                          "with different link functions\n")}
    }
##################
    
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
  varies <- sapply(X = lengths, FUN = function(i) i > 1)
##########################################
    ##CHANGES: add case when only a single variable appears in data frame
    if(ncol(newdata) == 1) {
        varies <- 1
    }

    ##add extractX to check that variables appearing in model also appear in data frame
    ##checkVariables <- extractX(cand.set, parm.type = parm.type)
    ##if(any(!checkVariables$predictors %in% names(newdata))) {
    ##    stop("\nAll predictors must appear in the 'newdata' data frame\n")
    ##}
##########################################
    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

  
  
  ##number of models
  nmods <- length(modnames)


  ##compute predicted values
  ##point estimate
  if(identical(type, "response")) {
    ##extract fitted value for observation obs
    if(identical(parm.type, "lambda") && identical(mixture.id, "ZIP")) {
      fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE,
                                           newdata = newdata)$fit)),
                    nrow = nmods, ncol = 2, byrow = TRUE)
      
      SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE,
                                          newdata = newdata)$se.fit)),
                   nrow = nmods, ncol = 2, byrow = TRUE)
      
    } else {
      fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                         type = parm.type1)$Predicted)),
                    nrow = nmods, ncol = 2, byrow = TRUE)

    ##extract SE for fitted value for observation obs
    SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                        type = parm.type1)$SE)),
                 nrow = nmods, ncol = 2, byrow = TRUE)
    }
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
  differ <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.differ <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  ##colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
     
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$differ <- differ
    AICctmp$SE.differ <- SE.differ
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$differ)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
      QAICctmp$differ <- differ
      QAICctmp$SE.differ <- SE.differ * sqrt(c.hat)

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$differ)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.differ^2 + (QAICctmp$differ - Mod.avg.out[, 1])^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.differ^2 + (QAICctmp$differ - Mod.avg.out[, 1])^2)))  
      }
      ##store table
      AICc.out <- QAICctmp
    
    }
  
  
  
  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord == FALSE && c.hat == 1) {
    AICtmp <- AICctab
    AICtmp$differ <- differ
    AICtmp$SE.differ <- SE.differ

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$differ)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }

  
  ##begin loop - QAICc
  if(second.ord == FALSE && c.hat > 1){
      
    ##create temporary data.frame to store fitted values and SE 
    QAICtmp <- AICctab
    QAICtmp$differ <- differ
    QAICtmp$SE.differ <- SE.differ* sqrt(c.hat)
      
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(QAICtmp$QAICWt*QAICtmp$differ)
      
    ##compute unconditional SE and store in output matrix
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.differ^2 + (QAICtmp$differ - Mod.avg.out[, 1])^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.differ^2 + (QAICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- QAICtmp
  }

  Group.variable <- paste(parm.id, "(", var.id, ")")
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]

    ##arrange in matrix
    predsOutMat <- matrix(data = c(Mod.avg.out[, 1], Mod.avg.out[, 2],
                                   Lower.CL, Upper.CL),
                          nrow = 1, ncol = 4)
    colnames(predsOutMat) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
    rownames(predsOutMat) <- "effect.size"

    Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                         "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1], 
                         "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                         "Upper.CL" = Upper.CL, "Matrix.output" = predsOutMat)
    class(Mod.eff.list) <- c("modavgEffect", "list")
    return(Mod.eff.list)  
}



##DS
modavgEffect.AICunmarkedFitDS <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                          nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                          type = "response", c.hat = 1, parm.type = NULL,
                                          ...) {

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgEffect for details\n")}

  ##rename values according to unmarked to extract from object

  ##lambda
  if(identical(parm.type, "lambda")) {
    parm.type1 <- "state"; parm.id <- "lam"
    ##check mixture type for mixture models
  }

  ##detect
  if(identical(parm.type, "detect")) {
      parm.type1 <- "det"; parm.id <- "p"
      ##check for key function used
      keyid <- unique(sapply(cand.set, FUN = function(i) i@keyfun))
      if(any(keyid == "uniform")) stop("\nDetection parameter not found in some models\n")
  }

##################
    ##extract link function
    check.link <- sapply(X = cand.set, FUN = function(i) eval(parse(text = paste("i@estimates@estimates$",
                                                                                 parm.type1, "@invlink",
                                                                                 sep = ""))))
    unique.link <- unique(check.link)
    select.link <- unique.link[1]
    
    if(identical(type, "link")) {
        if(length(unique.link) > 1) {stop("\nIt is not appropriate to compute a model averaged linear predictor\n",
                                          "with different link functions\n")}
    }
##################
    
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
  varies <- sapply(X = lengths, FUN = function(i) i > 1)
##########################################
    ##CHANGES: add case when only a single variable appears in data frame
    if(ncol(newdata) == 1) {
        varies <- 1
    }

    ##add extractX to check that variables appearing in model also appear in data frame
    ##checkVariables <- extractX(cand.set, parm.type = parm.type)
    ##if(any(!checkVariables$predictors %in% names(newdata))) {
    ##    stop("\nAll predictors must appear in the 'newdata' data frame\n")
    ##}
##########################################
    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

  
  ##number of models
  nmods <- length(modnames)


  ##compute predicted values
  ##point estimate
  if(identical(type, "response")) {
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
  differ <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.differ <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  ##colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
     
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$differ <- differ
    AICctmp$SE.differ <- SE.differ
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$differ)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
      QAICctmp$differ <- differ
      QAICctmp$SE.differ <- SE.differ * sqrt(c.hat)

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$differ)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.differ^2 + (QAICctmp$differ - Mod.avg.out[, 1])^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.differ^2 + (QAICctmp$differ - Mod.avg.out[, 1])^2)))  
      }
      ##store table
      AICc.out <- QAICctmp
    
    }
  
  
  
  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord == FALSE && c.hat == 1) {
    AICtmp <- AICctab
    AICtmp$differ <- differ
    AICtmp$SE.differ <- SE.differ

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$differ)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }

  
  ##begin loop - QAICc
  if(second.ord == FALSE && c.hat > 1){
      
    ##create temporary data.frame to store fitted values and SE 
    QAICtmp <- AICctab
    QAICtmp$differ <- differ
    QAICtmp$SE.differ <- SE.differ* sqrt(c.hat)
      
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(QAICtmp$QAICWt*QAICtmp$differ)
      
    ##compute unconditional SE and store in output matrix
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.differ^2 + (QAICtmp$differ - Mod.avg.out[, 1])^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.differ^2 + (QAICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- QAICtmp
  }

  Group.variable <- paste(parm.id, "(", var.id, ")")
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]

    ##arrange in matrix
    predsOutMat <- matrix(data = c(Mod.avg.out[, 1], Mod.avg.out[, 2],
                                   Lower.CL, Upper.CL),
                          nrow = 1, ncol = 4)
    colnames(predsOutMat) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
    rownames(predsOutMat) <- "effect.size"
    
    Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                         "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1], 
                         "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                         "Upper.CL" = Upper.CL, "Matrix.output" = predsOutMat)  
    class(Mod.eff.list) <- c("modavgEffect", "list")
    return(Mod.eff.list)  
}



##gdistsamp
modavgEffect.AICunmarkedFitGDS <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                           type = "response", c.hat = 1, parm.type = NULL,
                                           ...) {
  
  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgEffect for details\n")}

  ##rename values according to unmarked to extract from object

  ##lambda
  if(identical(parm.type, "lambda")) {
    parm.type1 <- "lambda"; parm.id <- "lam"
  }

  ##detect
  if(identical(parm.type, "detect")) {
      parm.type1 <- "det"; parm.id <- "p"
      ##check for key function used
      keyid <- unique(sapply(cand.set, FUN = function(i) i@keyfun))
      if(any(keyid == "uniform")) stop("\nDetection parameter not found in some models\n")
  }
    
  ##availability
  if(identical(parm.type, "phi")) {parm.type1 <- "phi"; parm.id <- "phi"}

##################
    ##extract link function
    check.link <- sapply(X = cand.set, FUN = function(i) eval(parse(text = paste("i@estimates@estimates$",
                                                                                 parm.type1, "@invlink",
                                                                                 sep = ""))))
    unique.link <- unique(check.link)
    select.link <- unique.link[1]
    
    if(identical(type, "link")) {
        if(length(unique.link) > 1) {stop("\nIt is not appropriate to compute a model averaged linear predictor\n",
                                          "with different link functions\n")}
    }
##################
    
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
  varies <- sapply(X = lengths, FUN = function(i) i > 1)
##########################################
    ##CHANGES: add case when only a single variable appears in data frame
    if(ncol(newdata) == 1) {
        varies <- 1
    }

    ##add extractX to check that variables appearing in model also appear in data frame
    ##checkVariables <- extractX(cand.set, parm.type = parm.type)
    ##if(any(!checkVariables$predictors %in% names(newdata))) {
    ##    stop("\nAll predictors must appear in the 'newdata' data frame\n")
    ##}
##########################################
    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

  
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
  differ <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.differ <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  ##colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
     
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$differ <- differ
    AICctmp$SE.differ <- SE.differ
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$differ)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
      QAICctmp$differ <- differ
      QAICctmp$SE.differ <- SE.differ * sqrt(c.hat)

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$differ)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.differ^2 + (QAICctmp$differ - Mod.avg.out[, 1])^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.differ^2 + (QAICctmp$differ - Mod.avg.out[, 1])^2)))  
      }
      ##store table
      AICc.out <- QAICctmp
    
    }
  
  
  
  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord == FALSE && c.hat == 1) {
    AICtmp <- AICctab
    AICtmp$differ <- differ
    AICtmp$SE.differ <- SE.differ

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$differ)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }

  
  ##begin loop - QAICc
  if(second.ord == FALSE && c.hat > 1){
      
    ##create temporary data.frame to store fitted values and SE 
    QAICtmp <- AICctab
    QAICtmp$differ <- differ
    QAICtmp$SE.differ <- SE.differ* sqrt(c.hat)
      
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(QAICtmp$QAICWt*QAICtmp$differ)
      
    ##compute unconditional SE and store in output matrix
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.differ^2 + (QAICtmp$differ - Mod.avg.out[, 1])^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.differ^2 + (QAICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- QAICtmp
  }

  Group.variable <- paste(parm.id, "(", var.id, ")")
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]

    ##arrange in matrix
    predsOutMat <- matrix(data = c(Mod.avg.out[, 1], Mod.avg.out[, 2],
                                   Lower.CL, Upper.CL),
                          nrow = 1, ncol = 4)
    colnames(predsOutMat) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
    rownames(predsOutMat) <- "effect.size"

    Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                         "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1], 
                         "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                         "Upper.CL" = Upper.CL, "Matrix.output" = predsOutMat)
    class(Mod.eff.list) <- c("modavgEffect", "list")
    return(Mod.eff.list)  
}



##occuFP
modavgEffect.AICunmarkedFitOccuFP <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                              nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                              type = "response", c.hat = 1, parm.type = NULL,
                                              ...) {

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgEffect for details\n")}

  ##rename values according to unmarked to extract from object
  ##psi
  if(identical(parm.type, "psi")) {
    parm.type1 <- "state"; parm.id <- "psi"
  }


  ##detect
  if(identical(parm.type, "detect")) {parm.type1 <- "det"; parm.id <- "p"}

  ##false positives
  if(identical(parm.type, "falsepos") || identical(parm.type, "fp")) {parm.type1 <- "fp"; parm.id <- "fp"}

  ##certain detections
  if(identical(parm.type, "certain")) {
      parm.type1 <- "b"; parm.id <- "b"
      ##check that parameter appears in all models
      parfreq <- sum(sapply(cand.set, FUN = function(i) any(names(i@estimates@estimates) == parm.type1)))
      if(!identical(length(cand.set), parfreq)) {
          stop("\nParameter \'", parm.type1, "\' (parm.type = \"", parm.type, "\") does not appear in all models:",
               "\ncannot compute model-averaged effect size across all models\n")
      }
  }


##################
    ##extract link function
    check.link <- sapply(X = cand.set, FUN = function(i) eval(parse(text = paste("i@estimates@estimates$",
                                                                                 parm.type1, "@invlink",
                                                                                 sep = ""))))
    unique.link <- unique(check.link)
    select.link <- unique.link[1]
    
    if(identical(type, "link")) {
        if(length(unique.link) > 1) {stop("\nIt is not appropriate to compute a model averaged linear predictor\n",
                                          "with different link functions\n")}
    }
##################
    
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
  varies <- sapply(X = lengths, FUN = function(i) i > 1)
##########################################
    ##CHANGES: add case when only a single variable appears in data frame
    if(ncol(newdata) == 1) {
        varies <- 1
    }

    ##add extractX to check that variables appearing in model also appear in data frame
    ##checkVariables <- extractX(cand.set, parm.type = parm.type)
    ##if(any(!checkVariables$predictors %in% names(newdata))) {
    ##    stop("\nAll predictors must appear in the 'newdata' data frame\n")
    ##}
##########################################
    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

  
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
  differ <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.differ <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  ##colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
     
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$differ <- differ
    AICctmp$SE.differ <- SE.differ
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$differ)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
      QAICctmp$differ <- differ
      QAICctmp$SE.differ <- SE.differ * sqrt(c.hat)

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$differ)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.differ^2 + (QAICctmp$differ - Mod.avg.out[, 1])^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.differ^2 + (QAICctmp$differ - Mod.avg.out[, 1])^2)))  
      }
      ##store table
      AICc.out <- QAICctmp
    
    }
  
  
  
  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord == FALSE && c.hat == 1) {
    AICtmp <- AICctab
    AICtmp$differ <- differ
    AICtmp$SE.differ <- SE.differ

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$differ)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }

  
  ##begin loop - QAICc
  if(second.ord == FALSE && c.hat > 1){
      
    ##create temporary data.frame to store fitted values and SE 
    QAICtmp <- AICctab
    QAICtmp$differ <- differ
    QAICtmp$SE.differ <- SE.differ* sqrt(c.hat)
      
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(QAICtmp$QAICWt*QAICtmp$differ)
      
    ##compute unconditional SE and store in output matrix
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.differ^2 + (QAICtmp$differ - Mod.avg.out[, 1])^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.differ^2 + (QAICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- QAICtmp
  }

  Group.variable <- paste(parm.id, "(", var.id, ")")
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]

    ##arrange in matrix
    predsOutMat <- matrix(data = c(Mod.avg.out[, 1], Mod.avg.out[, 2],
                                   Lower.CL, Upper.CL),
                          nrow = 1, ncol = 4)
    colnames(predsOutMat) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
    rownames(predsOutMat) <- "effect.size"

    Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                         "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1], 
                         "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                         "Upper.CL" = Upper.CL, "Matrix.output" = predsOutMat)
    class(Mod.eff.list) <- c("modavgEffect", "list")
    return(Mod.eff.list)  
}



##multinomPois
modavgEffect.AICunmarkedFitMPois <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                             nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                             type = "response", c.hat = 1, parm.type = NULL,
                                             ...) {

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgEffect for details\n")}

  ##rename values according to unmarked to extract from object
  ##lambda
  if(identical(parm.type, "lambda")) {
    parm.type1 <- "state"; parm.id <- "lam"
    ##set check to NULL for other models
    mixture.id <- NULL
  }
  

  ##detect
  if(identical(parm.type, "detect")) {parm.type1 <- "det"; parm.id <- "p"}

##################
    ##extract link function
    check.link <- sapply(X = cand.set, FUN = function(i) eval(parse(text = paste("i@estimates@estimates$",
                                                                                 parm.type1, "@invlink",
                                                                                 sep = ""))))
    unique.link <- unique(check.link)
    select.link <- unique.link[1]
    
    if(identical(type, "link")) {
        if(length(unique.link) > 1) {stop("\nIt is not appropriate to compute a model averaged linear predictor\n",
                                          "with different link functions\n")}
    }
##################
    
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
  varies <- sapply(X = lengths, FUN = function(i) i > 1)
##########################################
    ##CHANGES: add case when only a single variable appears in data frame
    if(ncol(newdata) == 1) {
        varies <- 1
    }

    ##add extractX to check that variables appearing in model also appear in data frame
    ##checkVariables <- extractX(cand.set, parm.type = parm.type)
    ##if(any(!checkVariables$predictors %in% names(newdata))) {
    ##    stop("\nAll predictors must appear in the 'newdata' data frame\n")
    ##}
##########################################
    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

  
  ##number of models
  nmods <- length(modnames)


  ##compute predicted values
  ##point estimate
  if(identical(type, "response")) {
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
  differ <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.differ <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  ##colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
     
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$differ <- differ
    AICctmp$SE.differ <- SE.differ
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$differ)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
      QAICctmp$differ <- differ
      QAICctmp$SE.differ <- SE.differ * sqrt(c.hat)

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$differ)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.differ^2 + (QAICctmp$differ - Mod.avg.out[, 1])^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.differ^2 + (QAICctmp$differ - Mod.avg.out[, 1])^2)))  
      }
      ##store table
      AICc.out <- QAICctmp
    
    }
  
  
  
  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord == FALSE && c.hat == 1) {
    AICtmp <- AICctab
    AICtmp$differ <- differ
    AICtmp$SE.differ <- SE.differ

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$differ)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }

  
  ##begin loop - QAICc
  if(second.ord == FALSE && c.hat > 1){
      
    ##create temporary data.frame to store fitted values and SE 
    QAICtmp <- AICctab
    QAICtmp$differ <- differ
    QAICtmp$SE.differ <- SE.differ* sqrt(c.hat)
      
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(QAICtmp$QAICWt*QAICtmp$differ)
      
    ##compute unconditional SE and store in output matrix
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.differ^2 + (QAICtmp$differ - Mod.avg.out[, 1])^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.differ^2 + (QAICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- QAICtmp
  }

  Group.variable <- paste(parm.id, "(", var.id, ")")
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]

    ##arrange in matrix
    predsOutMat <- matrix(data = c(Mod.avg.out[, 1], Mod.avg.out[, 2],
                                   Lower.CL, Upper.CL),
                          nrow = 1, ncol = 4)
    colnames(predsOutMat) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
    rownames(predsOutMat) <- "effect.size"

    Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                         "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1], 
                         "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                         "Upper.CL" = Upper.CL, "Matrix.output" = predsOutMat)
    class(Mod.eff.list) <- c("modavgEffect", "list")
    return(Mod.eff.list)  
}



##gmultmix
modavgEffect.AICunmarkedFitGMM <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                           type = "response", c.hat = 1, parm.type = NULL,
                                           ...) {

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgEffect for details\n")}

  ##rename values according to unmarked to extract from object
  ##lambda
  if(identical(parm.type, "lambda")) {
    parm.type1 <- "lambda"; parm.id <- "lam"
  }

  ##detect
  if(identical(parm.type, "detect")) {parm.type1 <- "det"; parm.id <- "p"}
  
  ##availability
  if(identical(parm.type, "phi")) {parm.type1 <- "phi"; parm.id <- "phi"}

##################
    ##extract link function
    check.link <- sapply(X = cand.set, FUN = function(i) eval(parse(text = paste("i@estimates@estimates$",
                                                                                 parm.type1, "@invlink",
                                                                                 sep = ""))))
    unique.link <- unique(check.link)
    select.link <- unique.link[1]
    
    if(identical(type, "link")) {
        if(length(unique.link) > 1) {stop("\nIt is not appropriate to compute a model averaged linear predictor\n",
                                          "with different link functions\n")}
    }
##################
    
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
  varies <- sapply(X = lengths, FUN = function(i) i > 1)
##########################################
    ##CHANGES: add case when only a single variable appears in data frame
    if(ncol(newdata) == 1) {
        varies <- 1
    }

    ##add extractX to check that variables appearing in model also appear in data frame
    ##checkVariables <- extractX(cand.set, parm.type = parm.type)
    ##if(any(!checkVariables$predictors %in% names(newdata))) {
    ##    stop("\nAll predictors must appear in the 'newdata' data frame\n")
    ##}
##########################################
    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

  
  
  ##number of models
  nmods <- length(modnames)


  ##compute predicted values
  ##point estimate
  if(identical(type, "response")) {
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
  differ <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.differ <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  ##colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
     
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$differ <- differ
    AICctmp$SE.differ <- SE.differ
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$differ)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
      QAICctmp$differ <- differ
      QAICctmp$SE.differ <- SE.differ * sqrt(c.hat)

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$differ)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.differ^2 + (QAICctmp$differ - Mod.avg.out[, 1])^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.differ^2 + (QAICctmp$differ - Mod.avg.out[, 1])^2)))  
      }
      ##store table
      AICc.out <- QAICctmp
    
    }
  
  
  
  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord == FALSE && c.hat == 1) {
    AICtmp <- AICctab
    AICtmp$differ <- differ
    AICtmp$SE.differ <- SE.differ

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$differ)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }

  
  ##begin loop - QAICc
  if(second.ord == FALSE && c.hat > 1){
      
    ##create temporary data.frame to store fitted values and SE 
    QAICtmp <- AICctab
    QAICtmp$differ <- differ
    QAICtmp$SE.differ <- SE.differ* sqrt(c.hat)
      
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(QAICtmp$QAICWt*QAICtmp$differ)
      
    ##compute unconditional SE and store in output matrix
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.differ^2 + (QAICtmp$differ - Mod.avg.out[, 1])^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.differ^2 + (QAICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- QAICtmp
  }

  Group.variable <- paste(parm.id, "(", var.id, ")")
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]

    ##arrange in matrix
    predsOutMat <- matrix(data = c(Mod.avg.out[, 1], Mod.avg.out[, 2],
                                   Lower.CL, Upper.CL),
                          nrow = 1, ncol = 4)
    colnames(predsOutMat) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
    rownames(predsOutMat) <- "effect.size"

    Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                         "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1], 
                         "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                         "Upper.CL" = Upper.CL, "Matrix.output" = predsOutMat)
    class(Mod.eff.list) <- c("modavgEffect", "list")
    return(Mod.eff.list)  
}



##gpcount
modavgEffect.AICunmarkedFitGPC <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                           type = "response", c.hat = 1, parm.type = NULL,
                                           ...) {

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgEffect for details\n")}

  ##rename values according to unmarked to extract from object

  ##lambda
  if(identical(parm.type, "lambda")) {
    parm.type1 <- "lambda"; parm.id <- "lam"
  }

  ##detect
  if(identical(parm.type, "detect")) {parm.type1 <- "det"; parm.id <- "p"}
  
  ##availability
  if(identical(parm.type, "phi")) {parm.type1 <- "phi"; parm.id <- "phi"}

##################
    ##extract link function
    check.link <- sapply(X = cand.set, FUN = function(i) eval(parse(text = paste("i@estimates@estimates$",
                                                                                 parm.type1, "@invlink",
                                                                                 sep = ""))))
    unique.link <- unique(check.link)
    select.link <- unique.link[1]
    
    if(identical(type, "link")) {
        if(length(unique.link) > 1) {stop("\nIt is not appropriate to compute a model averaged linear predictor\n",
                                          "with different link functions\n")}
    }
##################
    
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
  varies <- sapply(X = lengths, FUN = function(i) i > 1)
##########################################
    ##CHANGES: add case when only a single variable appears in data frame
    if(ncol(newdata) == 1) {
        varies <- 1
    }

    ##add extractX to check that variables appearing in model also appear in data frame
    ##checkVariables <- extractX(cand.set, parm.type = parm.type)
    ##if(any(!checkVariables$predictors %in% names(newdata))) {
    ##    stop("\nAll predictors must appear in the 'newdata' data frame\n")
    ##}
##########################################
    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

    
  ##number of models
  nmods <- length(modnames)


  ##compute predicted values
  ##point estimate
  if(identical(type, "response")) {
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
  differ <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.differ <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  ##colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
     
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$differ <- differ
    AICctmp$SE.differ <- SE.differ
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$differ)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
      QAICctmp$differ <- differ
      QAICctmp$SE.differ <- SE.differ * sqrt(c.hat)

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$differ)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.differ^2 + (QAICctmp$differ - Mod.avg.out[, 1])^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.differ^2 + (QAICctmp$differ - Mod.avg.out[, 1])^2)))  
      }
      ##store table
      AICc.out <- QAICctmp
    
    }
  
  
  
  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord == FALSE && c.hat == 1) {
    AICtmp <- AICctab
    AICtmp$differ <- differ
    AICtmp$SE.differ <- SE.differ

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$differ)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }

  
  ##begin loop - QAICc
  if(second.ord == FALSE && c.hat > 1){
      
    ##create temporary data.frame to store fitted values and SE 
    QAICtmp <- AICctab
    QAICtmp$differ <- differ
    QAICtmp$SE.differ <- SE.differ* sqrt(c.hat)
      
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(QAICtmp$QAICWt*QAICtmp$differ)
      
    ##compute unconditional SE and store in output matrix
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.differ^2 + (QAICtmp$differ - Mod.avg.out[, 1])^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.differ^2 + (QAICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- QAICtmp
  }

  Group.variable <- paste(parm.id, "(", var.id, ")")
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]

    ##arrange in matrix
    predsOutMat <- matrix(data = c(Mod.avg.out[, 1], Mod.avg.out[, 2],
                                   Lower.CL, Upper.CL),
                          nrow = 1, ncol = 4)
    colnames(predsOutMat) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
    rownames(predsOutMat) <- "effect.size"

    Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                         "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1], 
                         "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                         "Upper.CL" = Upper.CL, "Matrix.output" = predsOutMat)
    class(Mod.eff.list) <- c("modavgEffect", "list")
    return(Mod.eff.list)  
}



##unmarkedFitMMO
modavgEffect.AICunmarkedFitMMO <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                           type = "response", c.hat = 1, parm.type = NULL,
                                           ...) {
  
  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgEffect for details\n")}


  ##rename values according to unmarked to extract from object
  ##gamma
  if(identical(parm.type, "gamma")) {
    parm.type1 <- "gamma"; parm.id <- "gam"
  }

  ##lambda
  if(identical(parm.type, "lambda")) {
    parm.type1 <- "lambda"; parm.id <- "lam"
    
    ##check mixture type for mixture models
    mixture.type <- sapply(X = cand.set, FUN = function(i) i@mixture)
    unique.mixture <- unique(mixture.type)
    if(length(unique.mixture) > 1) {
      if(any(unique.mixture == "ZIP")) stop("\nThis function does not yet support mixing ZIP with other distributions\n")
    } else {
      mixture.id <- unique(mixture.type)
      if(identical(unique.mixture, "ZIP")) {
        if(identical(type, "link")) stop("\nLink scale not yet supported for ZIP mixtures\n")
      }
    }
  }
  

  ##omega
  if(identical(parm.type, "omega")) {
      parm.type1 <- "omega"; parm.id <- "omega"
  }

  ##iota (for immigration = TRUE with dynamics = "autoreg", "trend", "ricker", or "gompertz")
  if(identical(parm.type, "iota")) {
      parm.type1 <- "iota"; parm.id <- "iota"
      ##check that parameter appears in all models
      parfreq <- sum(sapply(cand.set, FUN = function(i) any(names(i@estimates@estimates) == parm.type1)))
      if(!identical(length(cand.set), parfreq)) {
          stop("\nParameter \'", parm.type1, "\' (parm.type = \"", parm.type, "\") does not appear in all models:",
               "\ncannot compute model-averaged effect size across all models\n")
      }
  }

  ##detect
  if(identical(parm.type, "detect")) {parm.type1 <- "det"; parm.id <- "p"}

##################
    ##extract link function
    check.link <- sapply(X = cand.set, FUN = function(i) eval(parse(text = paste("i@estimates@estimates$",
                                                                                 parm.type1, "@invlink",
                                                                                 sep = ""))))
    unique.link <- unique(check.link)
    select.link <- unique.link[1]
    
    if(identical(type, "link")) {
        if(length(unique.link) > 1) {stop("\nIt is not appropriate to compute a model averaged linear predictor\n",
                                          "with different link functions\n")}
    }
##################
    
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
  varies <- sapply(X = lengths, FUN = function(i) i > 1)
##########################################
    ##CHANGES: add case when only a single variable appears in data frame
    if(ncol(newdata) == 1) {
        varies <- 1
    }

    ##add extractX to check that variables appearing in model also appear in data frame
    ##checkVariables <- extractX(cand.set, parm.type = parm.type)
    ##if(any(!checkVariables$predictors %in% names(newdata))) {
    ##    stop("\nAll predictors must appear in the 'newdata' data frame\n")
    ##}
##########################################
    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

  
  
  ##number of models
  nmods <- length(modnames)


  ##compute predicted values
  ##point estimate
  if(identical(type, "response")) {
    ##extract fitted value for observation obs
    if(identical(parm.type, "lambda") && identical(mixture.id, "ZIP")) {
      fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE,
                                           newdata = newdata)$fit)),
                    nrow = nmods, ncol = 2, byrow = TRUE)
      
      SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE,
                                          newdata = newdata)$se.fit)),
                   nrow = nmods, ncol = 2, byrow = TRUE)
      
    } else {
      fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                         type = parm.type1)$Predicted)),
                    nrow = nmods, ncol = 2, byrow = TRUE)

    ##extract SE for fitted value for observation obs
    SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                        type = parm.type1)$SE)),
                 nrow = nmods, ncol = 2, byrow = TRUE)
    }
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
  differ <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.differ <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  ##colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
     
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$differ <- differ
    AICctmp$SE.differ <- SE.differ
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$differ)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
      QAICctmp$differ <- differ
      QAICctmp$SE.differ <- SE.differ * sqrt(c.hat)

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$differ)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.differ^2 + (QAICctmp$differ - Mod.avg.out[, 1])^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.differ^2 + (QAICctmp$differ - Mod.avg.out[, 1])^2)))  
      }
      ##store table
      AICc.out <- QAICctmp
    
    }
  
  
  
  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord == FALSE && c.hat == 1) {
    AICtmp <- AICctab
    AICtmp$differ <- differ
    AICtmp$SE.differ <- SE.differ

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$differ)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }

  
  ##begin loop - QAICc
  if(second.ord == FALSE && c.hat > 1){
      
    ##create temporary data.frame to store fitted values and SE 
    QAICtmp <- AICctab
    QAICtmp$differ <- differ
    QAICtmp$SE.differ <- SE.differ* sqrt(c.hat)
      
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(QAICtmp$QAICWt*QAICtmp$differ)
      
    ##compute unconditional SE and store in output matrix
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.differ^2 + (QAICtmp$differ - Mod.avg.out[, 1])^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.differ^2 + (QAICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- QAICtmp
  }

  Group.variable <- paste(parm.id, "(", var.id, ")")
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]

    ##arrange in matrix
    predsOutMat <- matrix(data = c(Mod.avg.out[, 1], Mod.avg.out[, 2],
                                   Lower.CL, Upper.CL),
                          nrow = 1, ncol = 4)
    colnames(predsOutMat) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
    rownames(predsOutMat) <- "effect.size"

    Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                         "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1], 
                         "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                         "Upper.CL" = Upper.CL, "Matrix.output" = predsOutMat)
    class(Mod.eff.list) <- c("modavgEffect", "list")
    return(Mod.eff.list)  
}



##unmarkedFitDSO
modavgEffect.AICunmarkedFitDSO <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                           type = "response", c.hat = 1, parm.type = NULL,
                                           ...) {
  
  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgEffect for details\n")}

  ##rename values according to unmarked to extract from object
  ##gamma
  if(identical(parm.type, "gamma")) {
    parm.type1 <- "gamma"; parm.id <- "gam"
  }

  ##lambda
  if(identical(parm.type, "lambda")) {
    parm.type1 <- "lambda"; parm.id <- "lam"
    
    ##check mixture type for mixture models
    mixture.type <- sapply(X = cand.set, FUN = function(i) i@mixture)
    unique.mixture <- unique(mixture.type)
    if(length(unique.mixture) > 1) {
      if(any(unique.mixture == "ZIP")) stop("\nThis function does not yet support mixing ZIP with other distributions\n")
    } else {
      mixture.id <- unique(mixture.type)
      if(identical(unique.mixture, "ZIP")) {
        if(identical(type, "link")) stop("\nLink scale not yet supported for ZIP mixtures\n")
      }
    }
  }
  

  ##omega
  if(identical(parm.type, "omega")) {
    parm.type1 <- "omega"; parm.id <- "omega"
  }

  ##iota (for immigration = TRUE with dynamics = "autoreg", "trend", "ricker", or "gompertz")
  if(identical(parm.type, "iota")) {
      parm.type1 <- "iota"; parm.id <- "iota"
      ##check that parameter appears in all models
      parfreq <- sum(sapply(cand.set, FUN = function(i) any(names(i@estimates@estimates) == parm.type1)))
      if(!identical(length(cand.set), parfreq)) {
          stop("\nParameter \'", parm.type1, "\' (parm.type = \"", parm.type, "\") does not appear in all models:",
               "\ncannot compute model-averaged effect size across all models\n")
      }
  }

    ##detect
    if(identical(parm.type, "detect")) {
        parm.type1 <- "det"
        parm.id <- "p"
        ##check for key function used
        keyid <- unique(sapply(cand.set, FUN = function(i) i@keyfun))
        if(any(keyid == "uniform")) stop("\nDetection parameter not found in some models\n")
    }

##################
    ##extract link function
    check.link <- sapply(X = cand.set, FUN = function(i) eval(parse(text = paste("i@estimates@estimates$",
                                                                                 parm.type1, "@invlink",
                                                                                 sep = ""))))
    unique.link <- unique(check.link)
    select.link <- unique.link[1]
    
    if(identical(type, "link")) {
        if(length(unique.link) > 1) {stop("\nIt is not appropriate to compute a model averaged linear predictor\n",
                                          "with different link functions\n")}
    }
##################
    
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
  varies <- sapply(X = lengths, FUN = function(i) i > 1)
##########################################
    ##CHANGES: add case when only a single variable appears in data frame
    if(ncol(newdata) == 1) {
        varies <- 1
    }

    ##add extractX to check that variables appearing in model also appear in data frame
    ##checkVariables <- extractX(cand.set, parm.type = parm.type)
    ##if(any(!checkVariables$predictors %in% names(newdata))) {
    ##    stop("\nAll predictors must appear in the 'newdata' data frame\n")
    ##}
##########################################
    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

  
  
  ##number of models
  nmods <- length(modnames)


  ##compute predicted values
  ##point estimate
  if(identical(type, "response")) {
    ##extract fitted value for observation obs
    if(identical(parm.type, "lambda") && identical(mixture.id, "ZIP")) {
      fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE,
                                           newdata = newdata)$fit)),
                    nrow = nmods, ncol = 2, byrow = TRUE)
      
      SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE,
                                          newdata = newdata)$se.fit)),
                   nrow = nmods, ncol = 2, byrow = TRUE)
      
    } else {
      fit <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                         type = parm.type1)$Predicted)),
                    nrow = nmods, ncol = 2, byrow = TRUE)

    ##extract SE for fitted value for observation obs
    SE <- matrix(data = unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                        type = parm.type1)$SE)),
                 nrow = nmods, ncol = 2, byrow = TRUE)
    }
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
  differ <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.differ <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  ##colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
     
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$differ <- differ
    AICctmp$SE.differ <- SE.differ
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$differ)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
      QAICctmp$differ <- differ
      QAICctmp$SE.differ <- SE.differ * sqrt(c.hat)

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$differ)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.differ^2 + (QAICctmp$differ - Mod.avg.out[, 1])^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.differ^2 + (QAICctmp$differ - Mod.avg.out[, 1])^2)))  
      }
      ##store table
      AICc.out <- QAICctmp
    
    }
  
  
  
  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord == FALSE && c.hat == 1) {
    AICtmp <- AICctab
    AICtmp$differ <- differ
    AICtmp$SE.differ <- SE.differ

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$differ)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }

  
  ##begin loop - QAICc
  if(second.ord == FALSE && c.hat > 1){
      
    ##create temporary data.frame to store fitted values and SE 
    QAICtmp <- AICctab
    QAICtmp$differ <- differ
    QAICtmp$SE.differ <- SE.differ* sqrt(c.hat)
      
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(QAICtmp$QAICWt*QAICtmp$differ)
      
    ##compute unconditional SE and store in output matrix
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.differ^2 + (QAICtmp$differ - Mod.avg.out[, 1])^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.differ^2 + (QAICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- QAICtmp
  }

  Group.variable <- paste(parm.id, "(", var.id, ")")
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]

    ##arrange in matrix
    predsOutMat <- matrix(data = c(Mod.avg.out[, 1], Mod.avg.out[, 2],
                                   Lower.CL, Upper.CL),
                          nrow = 1, ncol = 4)
    colnames(predsOutMat) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
    rownames(predsOutMat) <- "effect.size"

    Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                         "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1], 
                         "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                         "Upper.CL" = Upper.CL, "Matrix.output" = predsOutMat)
    class(Mod.eff.list) <- c("modavgEffect", "list")
    return(Mod.eff.list)  
}



##occuTTD
modavgEffect.AICunmarkedFitOccuTTD <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                               nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                               type = "response", c.hat = 1, parm.type = NULL,
                                               ...) {

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgEffect for details\n")}

  ##rename values according to unmarked to extract from object
  ##psi
  if(identical(parm.type, "psi")) {
    parm.type1 <- "psi"; parm.id <- "psi"
  }
    
  ##gamma
  if(identical(parm.type, "gamma")) {
      nseasons <- unique(sapply(cand.set, FUN = function(i) i@data@numPrimary))
      if(nseasons == 1) {
          stop("\nParameter \'gamma\' does not appear in single-season models\n")
      }

      parm.type1 <- "col"; parm.id <- "col"
  }

  ##epsilon
  if(identical(parm.type, "epsilon")) {
      nseasons <- unique(sapply(cand.set, FUN = function(i) i@data@numPrimary))
      if(nseasons == 1) {
          stop("\nParameter \'epsilon\' does not appear in single-season models\n")
      }

      parm.type1 <- "ext"; parm.id <- "ext"
  }

    ##detect
    if(identical(parm.type, "detect")) {parm.type1 <- "det"; parm.id <- "p"}

##################
    ##extract link function
    check.link <- sapply(X = cand.set, FUN = function(i) eval(parse(text = paste("i@estimates@estimates$",
                                                                                 parm.type1, "@invlink",
                                                                                 sep = ""))))
    unique.link <- unique(check.link)
    select.link <- unique.link[1]
    
    if(identical(type, "link")) {
        if(length(unique.link) > 1) {stop("\nIt is not appropriate to compute a model averaged linear predictor\n",
                                          "with different link functions\n")}
    }
##################
    
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
  varies <- sapply(X = lengths, FUN = function(i) i > 1)
##########################################
    ##CHANGES: add case when only a single variable appears in data frame
    if(ncol(newdata) == 1) {
        varies <- 1
    }

    ##add extractX to check that variables appearing in model also appear in data frame
    ##checkVariables <- extractX(cand.set, parm.type = parm.type)
    ##if(any(!checkVariables$predictors %in% names(newdata))) {
    ##    stop("\nAll predictors must appear in the 'newdata' data frame\n")
    ##}
##########################################
    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

  
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
  differ <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.differ <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  ##colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
     
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$differ <- differ
    AICctmp$SE.differ <- SE.differ
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$differ)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
      QAICctmp$differ <- differ
      QAICctmp$SE.differ <- SE.differ * sqrt(c.hat)

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$differ)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.differ^2 + (QAICctmp$differ - Mod.avg.out[, 1])^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.differ^2 + (QAICctmp$differ - Mod.avg.out[, 1])^2)))  
      }
      ##store table
      AICc.out <- QAICctmp
    
    }
  
  
  
  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord == FALSE && c.hat == 1) {
    AICtmp <- AICctab
    AICtmp$differ <- differ
    AICtmp$SE.differ <- SE.differ

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$differ)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }

  
  ##begin loop - QAICc
  if(second.ord == FALSE && c.hat > 1){
      
    ##create temporary data.frame to store fitted values and SE 
    QAICtmp <- AICctab
    QAICtmp$differ <- differ
    QAICtmp$SE.differ <- SE.differ* sqrt(c.hat)
      
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(QAICtmp$QAICWt*QAICtmp$differ)
      
    ##compute unconditional SE and store in output matrix
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.differ^2 + (QAICtmp$differ - Mod.avg.out[, 1])^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.differ^2 + (QAICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- QAICtmp
  }

  Group.variable <- paste(parm.id, "(", var.id, ")")
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]

    ##arrange in matrix
    predsOutMat <- matrix(data = c(Mod.avg.out[, 1], Mod.avg.out[, 2],
                                   Lower.CL, Upper.CL),
                          nrow = 1, ncol = 4)
    colnames(predsOutMat) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
    rownames(predsOutMat) <- "effect.size"

    Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                         "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1], 
                         "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                         "Upper.CL" = Upper.CL, "Matrix.output" = predsOutMat)
    class(Mod.eff.list) <- c("modavgEffect", "list")
    return(Mod.eff.list)  
}



##occuMS
modavgEffect.AICunmarkedFitOccuMS <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                              nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                              type = "response", c.hat = 1, parm.type = NULL,
                                              ...) {

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
        if(is.null(names(cand.set))) {
            modnames <- paste("Mod", 1:length(cand.set), sep = "")
            warning("\nModel names have been supplied automatically in the table\n")
        } else {
            modnames <- names(cand.set)
        }
    }


    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgEffect for details\n")}

    ##check if type = "link"
    if(identical(type, "link")) stop("\nLink scale predictions not yet supported for this model type\n")
   
    
    
    ##rename values according to unmarked to extract from object
    ##psi
    if(identical(parm.type, "psi")) {
        parm.type1 <- "state"
        parm.id <- "psi"
        ##because the same elements have different labels in different parts of the results of this object type
        parm.type.alt <- parm.type
    }
    
    ##transition
    if(identical(parm.type, "phi")) {
        ##check that parameter appears in all models
        nseasons <- unique(sapply(cand.set, FUN = function(i) i@data@numPrimary))
        if(nseasons == 1) {
            stop("\nParameter \'phi\' does not appear in single-season models\n")
        }
        parm.id <- "phi"
        parm.type1 <- "transition"
        ##because the same elements have different labels in different parts of the results of this object type
        parm.type.alt <- parm.type
    }

    ##detect
    if(identical(parm.type, "detect")) {
        parm.type1 <- "det"
        parm.id <- "p"
        ##because the same elements have different labels in different parts of the results of this object type
        parm.type.alt <- parm.type1
    }


##################
    ##extract link function
    check.link <- sapply(X = cand.set, FUN = function(i) eval(parse(text = paste("i@estimates@estimates$",
                                                                                 parm.type1, "@invlink",
                                                                                 sep = ""))))
    unique.link <- unique(check.link)
    select.link <- unique.link[1]
    
    if(identical(type, "link")) {
        if(length(unique.link) > 1) {stop("\nIt is not appropriate to compute a model averaged linear predictor\n",
                                          "with different link functions\n")}
    }
##################    
    
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
    varies <- sapply(X = lengths, FUN = function(i) i > 1)
##########################################
    ##CHANGES: add case when only a single variable appears in data frame
    if(ncol(newdata) == 1) {
        varies <- 1
    }

    ##add extractX to check that variables appearing in model also appear in data frame
    ##checkVariables <- extractX(cand.set, parm.type = parm.type)
    ##if(any(!checkVariables$predictors %in% names(newdata))) {
    ##    stop("\nAll predictors must appear in the 'newdata' data frame\n")
    ##}
##########################################
    
    ##extract name of column
    if(sum(varies) == 1) {
        var.id <- names(varies)[which(varies == TRUE)]
                
        ##determine name of groups compared
        group1 <- as.character(newdata[,paste(var.id)][1])
        group2 <- as.character(newdata[,paste(var.id)][2])
        
    } else {
        ##warn that no single variable defines groups
        warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
        ##use generic names
        var.id <- "Groups"
        group1 <- "group 1"
        group2 <- "group 2"
    }

  
    ##number of models
    nmods <- length(modnames)


    ##compute predicted values
    ##point estimate
    ##extract predicted values
    predsList <- lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                               type = parm.type.alt, ...))
    
    ##determine number of parameters
    parmFirst <- predsList[[1]]
    parmNames <- names(parmFirst)
    nparms <- length(parmNames)


    ##lists to store predictions and SE's
    predsEstList <- vector("list", nparms)
    names(predsEstList) <- parmNames
    predsSEList <- vector("list", nparms)
    names(predsSEList) <- parmNames

    ##iterate over each parm
    for(k in 1:nparms) {
        predsEstList[[k]] <- lapply(predsList, FUN = function(i) i[[k]]$Predicted)
        predsSEList[[k]] <- lapply(predsList, FUN = function(i) i[[k]]$SE)
    }
    
    
    ##organize in an nobs x nmodels x nparms array
    predsEst <- array(unlist(predsEstList), dim = c(nobserv, length(cand.set), nparms),
                      dimnames = list(1:nobserv, modnames, parmNames))
    predsSE <- array(unlist(predsSEList), dim = c(nobserv, length(cand.set), nparms),
                     dimnames = list(1:nobserv, modnames, parmNames))
    
    
    ##adjust for overdispersion if c-hat > 1
    if(c.hat > 1) {predsSE <- predsSE * sqrt(c.hat)}
    
      
    ##difference between groups 
    differList <- vector("list", nparms)

    for(k in 1:nparms) {
        differList[[k]] <- predsEst[1, , k] - predsEst[2, , k]
    }
    
    ##SE on difference
    SE.differList <- vector("list", nparms)
    for(k in 1:nparms) {
        SE.differList[[k]] <- sqrt(predsSE[1, , k]^2 + predsSE[2, , k]^2)
    }


    
    ##store AICc table
    AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                      nobs = nobs, sort = FALSE, c.hat = c.hat)

    ##prepare list for model-averaged predictions and SE's
    predsOut <-  array(data = NA, dim = c(1, 4, nparms),
                       dimnames = list(1,
                                       c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL"),
                                       parmNames))

    ##begin loop - AICc
    if(second.ord == TRUE && c.hat == 1){
    
        ##create temporary data.frame to store fitted values and SE 
        AICctmp <- AICctab

        ##compute unconditional SE and store in output matrix
        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    
        if(identical(uncond.se, "old")) {
            for(j in 1:nparms) { 
                predsOut[, 1, j] <- sum(AICctmp$AICcWt * differList[[j]])
                predsOut[, 2, j] <- sum(AICctmp$AICcWt * sqrt(SE.differList[[j]]^2 + (differList[[j]] - predsOut[, 1, j])^2))
            }
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
            for(j in 1:nparms) { 
                predsOut[, 1, j] <- sum(AICctmp$AICcWt * differList[[j]])
                predsOut[, 2, j] <- sqrt(sum(AICctmp$AICcWt * (SE.differList[[j]]^2 + (differList[[j]] - predsOut[, 1, j])^2)))
            }
        }
    }


    ##begin loop - QAICc
    if(second.ord == TRUE && c.hat > 1){

        ##create temporary data.frame to store fitted values and SE 
        QAICctmp <- AICctab
        
        ##compute unconditional SE and store in output matrix
        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    
        if(identical(uncond.se, "old")) {
            for(j in 1:nparms) { 
                predsOut[, 1, j] <- sum(QAICctmp$QAICcWt * differList[[j]])
                predsOut[, 2, j] <- sum(QAICctmp$QAICcWt * sqrt(SE.differList[[j]]^2 + (differList[[j]] - predsOut[, 1, j])^2))
            }
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
            for(j in 1:nparms) { 
                predsOut[, 1, j] <- sum(QAICctmp$QAICcWt * differList[[j]])
                predsOut[, 2, j] <- sqrt(sum(QAICctmp$QAICcWt * (SE.differList[[j]]^2 + (differList[[j]] - predsOut[, 1, j])^2)))
            }
        }
    }


    ##begin loop - AIC
    if(second.ord == FALSE && c.hat == 1) {

        
        ##create temporary data.frame to store fitted values and SE 
        AICtmp <- AICctab

        ##compute unconditional SE and store in output matrix
        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    
        if(identical(uncond.se, "old")) {
            for(j in 1:nparms) { 
                predsOut[, 1, j] <- sum(AICtmp$AICWt * differList[[j]])
                predsOut[, 2, j] <- sum(AICtmp$AICWt * sqrt(SE.differList[[j]]^2 + (differList[[j]] - predsOut[, 1, j])^2))
            }
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
            for(j in 1:nparms) { 
                predsOut[, 1, j] <- sum(AICtmp$AICWt * differList[[j]])
                predsOut[, 2, j] <- sqrt(sum(AICtmp$AICWt * (SE.differList[[j]]^2 + (differList[[j]] - predsOut[, 1, j])^2)))
            }
        }
    }


    
    ##begin loop - QAICc
    if(second.ord == FALSE && c.hat > 1){
        
        ##create temporary data.frame to store fitted values and SE 
        QAICtmp <- AICctab

        ##compute unconditional SE and store in output matrix
        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    
        if(identical(uncond.se, "old")) {
            for(j in 1:nparms) { 
                predsOut[, 1, j] <- sum(QAICtmp$QAICWt * differList[[j]])
                predsOut[, 2, j] <- sum(QAICtmp$QAICWt * sqrt(SE.differList[[j]]^2 + (differList[[j]] - predsOut[, 1, j])^2))
            }
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
            for(j in 1:nparms) { 
                predsOut[, 1, j] <- sum(QAICtmp$QAICWt * differList[[j]])
                predsOut[, 2, j] <- sqrt(sum(QAICtmp$QAICWt * (SE.differList[[j]]^2 + (differList[[j]] - predsOut[, 1, j])^2)))
            }
        }
    }

    
    ##compute confidence interval
    zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)

    for(j in 1:nparms){
        predsOut[, 3, j] <- predsOut[, 1, j] - zcrit * predsOut[, 2, j]
        predsOut[, 4, j] <- predsOut[, 1, j] + zcrit * predsOut[, 2, j]
    }

    ##format to matrix
    arrayToMat <- apply(predsOut, 2L, c)

    #if(is.vector(arrayToMat)) {
    #    predsOutMat <- matrix(arrayToMat, nrow = 1)
    #    colnames(predsOutMat) <- names(arrayToMat)
    #    AICctab$differ <- unlist(differList)
    #    AICctab$SE.differ <- unlist(SE.differList)
    #} else {

    predsOutMat <- arrayToMat
    #}
    
    ##create label for rows
    rownames(predsOutMat) <- parmNames

    ##convert array to list
    ##predsOutList <- lapply(seq(dim(predsOut)[3]), function(i) predsOut[ , , i])
    ##names(predsOutList) <- parmNames

    ##store table
    AICc.out <- AICctab
    
    Group.variable <- paste(parm.id, "(", var.id, ")")
    
    ##organize as list
    Mod.eff.list <- list("Group.variable" = Group.variable, "Group1" = group1,
                         "Group2" = group2, "Type" = type,
                         "Mod.avg.table" = AICc.out,
                         "Mod.avg.eff" = predsOut[, 1, ],
                         "Uncond.se" = predsOut[, 2, ],
                         "Conf.level" = conf.level,
                         "Lower.CL" = predsOut[, 3, ],
                         "Upper.CL" = predsOut[, 4, ],
                         "Matrix.output" = predsOutMat)
    class(Mod.eff.list) <- c("modavgEffect", "list")
    return(Mod.eff.list)  
}



##occuMulti
modavgEffect.AICunmarkedFitOccuMulti <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                                 nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                                 type = "response", c.hat = 1, parm.type = NULL,
                                                 ...) {

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
        if(is.null(names(cand.set))) {
            modnames <- paste("Mod", 1:length(cand.set), sep = "")
            warning("\nModel names have been supplied automatically in the table\n")
        } else {
            modnames <- names(cand.set)
        }
    }


    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgEffect for details\n")}

    ##check if type = "link"
    if(identical(type, "link")) stop("\nLink scale predictions not yet supported for this model type\n")
   
    
    
    ##rename values according to unmarked to extract from object
    ##psi
    if(identical(parm.type, "psi")) {
        parm.type1 <- "state"
        parm.id <- "psi"
    }
    
    
    ##detect
    if(identical(parm.type, "detect")) {
        parm.type1 <- "det"
        parm.id <- "p"
    }


##################
    ##extract link function
    check.link <- sapply(X = cand.set, FUN = function(i) eval(parse(text = paste("i@estimates@estimates$",
                                                                                 parm.type1, "@invlink",
                                                                                 sep = ""))))
    unique.link <- unique(check.link)
    select.link <- unique.link[1]
    
    if(identical(type, "link")) {
        if(length(unique.link) > 1) {stop("\nIt is not appropriate to compute a model averaged linear predictor\n",
                                          "with different link functions\n")}
    }
##################    
    
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
    varies <- sapply(X = lengths, FUN = function(i) i > 1)
##########################################
    ##CHANGES: add case when only a single variable appears in data frame
    if(ncol(newdata) == 1) {
        varies <- 1
    }

    ##add extractX to check that variables appearing in model also appear in data frame
    ##checkVariables <- extractX(cand.set, parm.type = parm.type)
    ##if(any(!checkVariables$predictors %in% names(newdata))) {
    ##    stop("\nAll predictors must appear in the 'newdata' data frame\n")
    ##}
##########################################
    
    ##extract name of column
    if(sum(varies) == 1) {
        var.id <- names(varies)[which(varies == TRUE)]
                
        ##determine name of groups compared
        group1 <- as.character(newdata[,paste(var.id)][1])
        group2 <- as.character(newdata[,paste(var.id)][2])
        
    } else {
        ##warn that no single variable defines groups
        warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
        ##use generic names
        var.id <- "Groups"
        group1 <- "group 1"
        group2 <- "group 2"
    }

  
    ##number of models
    nmods <- length(modnames)


    ##compute predicted values
    ##point estimate
    ##extract predicted values
    predsList <- lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata,
                                                               type = parm.type1, ...))
    
    ##check structure of predsList
    ##if predictions for all psi parameters
    if(is.matrix(predsList[[1]]$Predicted) && identical(parm.type, "psi")) {

        ##determine number of parameters
        ##check if species argument was provided in call to function
        parmFirst <- predsList[[1]]$Predicted
        parmNames <- colnames(parmFirst)
        nparms <- length(parmNames)
  
        ##lists to store predictions and SE's
        predsEstList <- vector("list", nparms)
        names(predsEstList) <- parmNames
        predsSEList <- vector("list", nparms)
        names(predsSEList) <- parmNames

        ##iterate over each parm
        predsEstList <- lapply(predsList, FUN = function(i) i$Predicted)
        predsSEList <- lapply(predsList, FUN = function(i) i$SE)

        for(k in 1:nparms) {
            predsEstList[[k]] <- lapply(predsList, FUN = function(i) i$Predicted[, k])
            predsSEList[[k]] <- lapply(predsList, FUN = function(i) i$SE[, k])
        }

        ##organize in an nobs x nmodels x nparms array
        predsEst <- array(unlist(predsEstList), dim = c(nobserv, length(cand.set), nparms),
                          dimnames = list(1:nobserv, modnames, parmNames))
        predsSE <- array(unlist(predsSEList), dim = c(nobserv, length(cand.set), nparms),
                         dimnames = list(1:nobserv, modnames, parmNames))
    }
    
    ##if predictions for single species
    if(!is.matrix(predsList[[1]]$Predicted) && identical(parm.type, "psi")) {
    
        parmNames <- parm.type
        nparms <- length(parmNames)

        predsEstMat <- sapply(predsList, FUN = function(i) i$Predicted)
        predsSEMat <- sapply(predsList, FUN = function(i) i$SE)
        
        predsEst <- array(predsEstMat, dim = c(nobserv, length(cand.set), nparms),
                          dimnames = list(1:nobserv, modnames, parmNames))
        predsSE <- array(predsSEMat, dim = c(nobserv, length(cand.set), nparms),
                         dimnames = list(1:nobserv, modnames, parmNames))
    }

    ##if predictions for detection
    if(identical(parm.type, "detect")) {
        parmFirst <- predsList[[1]]

        if(!is.data.frame(parmFirst)) {
            orig.parmNames <- names(parmFirst)
            parmNames <- paste("p", orig.parmNames, sep = "-")
            nparms <- length(parmNames)
            
            ##iterate over species
            predsEst <- array(NA, dim = c(nobserv, length(cand.set), nparms),
                              dimnames = list(1:nobserv, modnames, parmNames))
            predsSE <- array(NA, dim = c(nobserv, length(cand.set), nparms),
                             dimnames = list(1:nobserv, modnames, parmNames))

            
            ##iterate over each parm
            for(k in 1:nparms) {
                predsEst[, , k] <- sapply(predsList, FUN = function(i) i[[k]]$Predicted)
                predsSE[, , k] <- sapply(predsList, FUN = function(i) i[[k]]$SE)
            }
        } else {
            ##single parameter p
            parmNames <- "p"
            nparms <- 1

            ##iterate over species
            predsEst <- array(NA, dim = c(nobserv, length(cand.set), nparms),
                              dimnames = list(1:nobserv, modnames, parmNames))
            predsSE <- array(NA, dim = c(nobserv, length(cand.set), nparms),
                             dimnames = list(1:nobserv, modnames, parmNames))

            predsEst[, , "p"] <- sapply(predsList, FUN = function(i) i[, "Predicted"])
            predsSE[, , "p"] <- sapply(predsList, FUN = function(i) i[, "SE"])
        }
    }

    
    ##adjust for overdispersion if c-hat > 1
    if(c.hat > 1) {predsSE <- predsSE * sqrt(c.hat)}
    
      
    ##difference between groups 
    differList <- vector("list", nparms)

    for(k in 1:nparms) {
        differList[[k]] <- predsEst[1, , k] - predsEst[2, , k]
    }
    
    ##SE on difference
    SE.differList <- vector("list", nparms)
    for(k in 1:nparms) {
        SE.differList[[k]] <- sqrt(predsSE[1, , k]^2 + predsSE[2, , k]^2)
    }


    
    ##store AICc table
    AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                      nobs = nobs, sort = FALSE, c.hat = c.hat)

    ##prepare list for model-averaged predictions and SE's
    predsOut <-  array(data = NA, dim = c(1, 4, nparms),
                       dimnames = list(1,
                                       c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL"),
                                       parmNames))

    ##begin loop - AICc
    if(second.ord == TRUE && c.hat == 1){
    
        ##create temporary data.frame to store fitted values and SE 
        AICctmp <- AICctab

        ##compute unconditional SE and store in output matrix
        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    
        if(identical(uncond.se, "old")) {
            for(j in 1:nparms) { 
                predsOut[, 1, j] <- sum(AICctmp$AICcWt * differList[[j]])
                predsOut[, 2, j] <- sum(AICctmp$AICcWt * sqrt(SE.differList[[j]]^2 + (differList[[j]] - predsOut[, 1, j])^2))
            }
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
            for(j in 1:nparms) { 
                predsOut[, 1, j] <- sum(AICctmp$AICcWt * differList[[j]])
                predsOut[, 2, j] <- sqrt(sum(AICctmp$AICcWt * (SE.differList[[j]]^2 + (differList[[j]] - predsOut[, 1, j])^2)))
            }
        }
    }


    ##begin loop - QAICc
    if(second.ord == TRUE && c.hat > 1){

        ##create temporary data.frame to store fitted values and SE 
        QAICctmp <- AICctab
        
        ##compute unconditional SE and store in output matrix
        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    
        if(identical(uncond.se, "old")) {
            for(j in 1:nparms) { 
                predsOut[, 1, j] <- sum(QAICctmp$QAICcWt * differList[[j]])
                predsOut[, 2, j] <- sum(QAICctmp$QAICcWt * sqrt(SE.differList[[j]]^2 + (differList[[j]] - predsOut[, 1, j])^2))
            }
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
            for(j in 1:nparms) { 
                predsOut[, 1, j] <- sum(QAICctmp$QAICcWt * differList[[j]])
                predsOut[, 2, j] <- sqrt(sum(QAICctmp$QAICcWt * (SE.differList[[j]]^2 + (differList[[j]] - predsOut[, 1, j])^2)))
            }
        }
    }


    ##begin loop - AIC
    if(second.ord == FALSE && c.hat == 1) {

        
        ##create temporary data.frame to store fitted values and SE 
        AICtmp <- AICctab

        ##compute unconditional SE and store in output matrix
        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    
        if(identical(uncond.se, "old")) {
            for(j in 1:nparms) { 
                predsOut[, 1, j] <- sum(AICtmp$AICWt * differList[[j]])
                predsOut[, 2, j] <- sum(AICtmp$AICWt * sqrt(SE.differList[[j]]^2 + (differList[[j]] - predsOut[, 1, j])^2))
            }
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
            for(j in 1:nparms) { 
                predsOut[, 1, j] <- sum(AICtmp$AICWt * differList[[j]])
                predsOut[, 2, j] <- sqrt(sum(AICtmp$AICWt * (SE.differList[[j]]^2 + (differList[[j]] - predsOut[, 1, j])^2)))
            }
        }
    }


    
    ##begin loop - QAICc
    if(second.ord == FALSE && c.hat > 1){
        
        ##create temporary data.frame to store fitted values and SE 
        QAICtmp <- AICctab

        ##compute unconditional SE and store in output matrix
        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    
        if(identical(uncond.se, "old")) {
            for(j in 1:nparms) { 
                predsOut[, 1, j] <- sum(QAICtmp$QAICWt * differList[[j]])
                predsOut[, 2, j] <- sum(QAICtmp$QAICWt * sqrt(SE.differList[[j]]^2 + (differList[[j]] - predsOut[, 1, j])^2))
            }
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
            for(j in 1:nparms) { 
                predsOut[, 1, j] <- sum(QAICtmp$QAICWt * differList[[j]])
                predsOut[, 2, j] <- sqrt(sum(QAICtmp$QAICWt * (SE.differList[[j]]^2 + (differList[[j]] - predsOut[, 1, j])^2)))
            }
        }
    }

    
    ##compute confidence interval
    zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)

    for(j in 1:nparms){
        predsOut[, 3, j] <- predsOut[, 1, j] - zcrit * predsOut[, 2, j]
        predsOut[, 4, j] <- predsOut[, 1, j] + zcrit * predsOut[, 2, j]
    }

    ##format to matrix
    arrayToMat <- apply(predsOut, 2L, c)

    if(is.vector(arrayToMat)) {
        predsOutMat <- matrix(arrayToMat, nrow = 1)
        colnames(predsOutMat) <- names(arrayToMat)
        AICctab$differ <- unlist(differList)
        AICctab$SE.differ <- unlist(SE.differList)
    } else {

        predsOutMat <- arrayToMat
    }
    
    ##create label for rows
    rownames(predsOutMat) <- parmNames
    
    ##store table
    AICc.out <- AICctab
    
    Group.variable <- paste(parm.id, "(", var.id, ")")
    
    ##organize as list
    Mod.eff.list <- list("Group.variable" = Group.variable, "Group1" = group1,
                         "Group2" = group2, "Type" = type,
                         "Mod.avg.table" = AICc.out,
                         "Mod.avg.eff" = predsOut[, 1, ],
                         "Uncond.se" = predsOut[, 2, ],
                         "Conf.level" = conf.level,
                         "Lower.CL" = predsOut[, 3, ],
                         "Upper.CL" = predsOut[, 4, ],
                         "Matrix.output" = predsOutMat)
    class(Mod.eff.list) <- c("modavgEffect", "list")
    return(Mod.eff.list)  
}



##goccu
modavgEffect.AICunmarkedFitGOccu <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                             nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                             type = "response", c.hat = 1, parm.type = NULL,
                                             ...) {

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgEffect for details\n")}

  ##rename values according to unmarked to extract from object

  ##psi
  if(identical(parm.type, "psi")) {
    parm.type1 <- "psi"; parm.id <- "psi"
  }

  ##detect
  if(identical(parm.type, "detect")) {parm.type1 <- "det"; parm.id <- "p"}
  
  ##availability
  if(identical(parm.type, "phi")) {parm.type1 <- "phi"; parm.id <- "phi"}

##################
    ##extract link function
    check.link <- sapply(X = cand.set, FUN = function(i) eval(parse(text = paste("i@estimates@estimates$",
                                                                                 parm.type1, "@invlink",
                                                                                 sep = ""))))
    unique.link <- unique(check.link)
    select.link <- unique.link[1]
    
    if(identical(type, "link")) {
        if(length(unique.link) > 1) {stop("\nIt is not appropriate to compute a model averaged linear predictor\n",
                                          "with different link functions\n")}
    }
##################
    
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
  varies <- sapply(X = lengths, FUN = function(i) i > 1)
##########################################
    ##CHANGES: add case when only a single variable appears in data frame
    if(ncol(newdata) == 1) {
        varies <- 1
    }

    ##add extractX to check that variables appearing in model also appear in data frame
    ##checkVariables <- extractX(cand.set, parm.type = parm.type)
    ##if(any(!checkVariables$predictors %in% names(newdata))) {
    ##    stop("\nAll predictors must appear in the 'newdata' data frame\n")
    ##}
##########################################
    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

    
  ##number of models
  nmods <- length(modnames)


  ##compute predicted values
  ##point estimate
  if(identical(type, "response")) {
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
  differ <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.differ <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  ##colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
     
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$differ <- differ
    AICctmp$SE.differ <- SE.differ
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$differ)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
      QAICctmp$differ <- differ
      QAICctmp$SE.differ <- SE.differ * sqrt(c.hat)

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$differ)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.differ^2 + (QAICctmp$differ - Mod.avg.out[, 1])^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.differ^2 + (QAICctmp$differ - Mod.avg.out[, 1])^2)))  
      }
      ##store table
      AICc.out <- QAICctmp
    
    }
  
  
  
  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord == FALSE && c.hat == 1) {
    AICtmp <- AICctab
    AICtmp$differ <- differ
    AICtmp$SE.differ <- SE.differ

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$differ)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }

  
  ##begin loop - QAICc
  if(second.ord == FALSE && c.hat > 1){
      
    ##create temporary data.frame to store fitted values and SE 
    QAICtmp <- AICctab
    QAICtmp$differ <- differ
    QAICtmp$SE.differ <- SE.differ* sqrt(c.hat)
      
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(QAICtmp$QAICWt*QAICtmp$differ)
      
    ##compute unconditional SE and store in output matrix
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.differ^2 + (QAICtmp$differ - Mod.avg.out[, 1])^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.differ^2 + (QAICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- QAICtmp
  }

  Group.variable <- paste(parm.id, "(", var.id, ")")
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]

    ##arrange in matrix
    predsOutMat <- matrix(data = c(Mod.avg.out[, 1], Mod.avg.out[, 2],
                                   Lower.CL, Upper.CL),
                          nrow = 1, ncol = 4)
    colnames(predsOutMat) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
    rownames(predsOutMat) <- "effect.size"

    Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                         "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1], 
                         "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                         "Upper.CL" = Upper.CL, "Matrix.output" = predsOutMat)
    class(Mod.eff.list) <- c("modavgEffect", "list")
    return(Mod.eff.list)  
}



##occuComm
modavgEffect.AICunmarkedFitOccuComm <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                                nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                                type = "response", c.hat = 1, parm.type = NULL,
                                                ...) {

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgEffect for details\n")}

    
    ##rename values according to unmarked to extract from object
    ##psi
    if(identical(parm.type, "psi")) {
        parm.type1 <- "state"; parm.id <- "psi"
    
    }

    ##detect
    if(identical(parm.type, "detect")) {parm.type1 <- "det"; parm.id <- "p"}


##################
    ##extract link function
    check.link <- sapply(X = cand.set, FUN = function(i) eval(parse(text = paste("i@estimates@estimates$",
                                                                                 parm.type1, "@invlink",
                                                                                 sep = ""))))
    unique.link <- unique(check.link)
    select.link <- unique.link[1]
    
    if(identical(type, "link")) {
        if(length(unique.link) > 1) {stop("\nIt is not appropriate to compute a model averaged linear predictor\n",
                                          "with different link functions\n")}
    }
##################    
    
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
  varies <- sapply(X = lengths, FUN = function(i) i > 1)
##########################################
    ##CHANGES: add case when only a single variable appears in data frame
    if(ncol(newdata) == 1) {
        varies <- 1
    }

    ##add extractX to check that variables appearing in model also appear in data frame
    ##checkVariables <- extractX(cand.set, parm.type = parm.type)
    ##if(any(!checkVariables$predictors %in% names(newdata))) {
    ##    stop("\nAll predictors must appear in the 'newdata' data frame\n")
    ##}
##########################################
    
  ##extract name of column
  if(sum(varies) == 1) {
    var.id <- names(varies)[which(varies == TRUE)]
                
    ##determine name of groups compared
    group1 <- as.character(newdata[,paste(var.id)][1])
    group2 <- as.character(newdata[,paste(var.id)][2])

  } else {
    ##warn that no single variable defines groups
    warning("\nGroups do not seem to be defined by a single variable.\n Function proceeding with generic group names\n")
    ##use generic names
    var.id <- "Groups"
    group1 <- "group 1"
    group2 <- "group 2"
  }

  
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
  differ <- fit[, 1] - fit[, 2]
    
  ##SE on difference
  SE.differ <- sqrt(SE[, 1]^2 + SE[, 2]^2)
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = 1, ncol = 2)
  ##colnames(Mod.avg.out) <- c("Mod.avg.diff", "Uncond.SE")

  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
     
    ##create temporary data.frame to store fitted values and SE 
    AICctmp <- AICctab
    AICctmp$differ <- differ
    AICctmp$SE.differ <- SE.differ
    
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICctmp$AICcWt*AICctmp$differ)
    ##compute unconditional SE and store in output matrix
      
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.differ^2 + (AICctmp$differ - Mod.avg.out[, 1])^2)))
    }
    ##store table
    AICc.out <- AICctmp
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
      QAICctmp$differ <- differ
      QAICctmp$SE.differ <- SE.differ * sqrt(c.hat)

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$differ)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.differ^2 + (QAICctmp$differ - Mod.avg.out[, 1])^2))
      }
      
      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.differ^2 + (QAICctmp$differ - Mod.avg.out[, 1])^2)))  
      }
      ##store table
      AICc.out <- QAICctmp
    
    }
  
  
  
  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord == FALSE && c.hat == 1) {
    AICtmp <- AICctab
    AICtmp$differ <- differ
    AICtmp$SE.differ <- SE.differ

    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(AICtmp$AICWt*AICtmp$differ)

    ##compute unconditional SE and store in output matrix
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.differ^2 + (AICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- AICtmp
  }

  
  ##begin loop - QAICc
  if(second.ord == FALSE && c.hat > 1){
      
    ##create temporary data.frame to store fitted values and SE 
    QAICtmp <- AICctab
    QAICtmp$differ <- differ
    QAICtmp$SE.differ <- SE.differ* sqrt(c.hat)
      
    ##compute model averaged prediction and store in output matrix
    Mod.avg.out[, 1] <- sum(QAICtmp$QAICWt*QAICtmp$differ)
      
    ##compute unconditional SE and store in output matrix
    if(identical(uncond.se, "old")) {
      Mod.avg.out[, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.differ^2 + (QAICtmp$differ - Mod.avg.out[, 1])^2))
    }
      
    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.differ^2 + (QAICtmp$differ - Mod.avg.out[, 1])^2)))
    }  
    ##store table
    AICc.out <- QAICtmp
  }

  Group.variable <- paste(parm.id, "(", var.id, ")")
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower.CL <- Mod.avg.out[, 1] - zcrit * Mod.avg.out[, 2]
  Upper.CL <- Mod.avg.out[, 1] + zcrit * Mod.avg.out[, 2]

    ##arrange in matrix
    predsOutMat <- matrix(data = c(Mod.avg.out[, 1], Mod.avg.out[, 2],
                                   Lower.CL, Upper.CL),
                          nrow = 1, ncol = 4)
    colnames(predsOutMat) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
    rownames(predsOutMat) <- "effect.size"

    Mod.eff.list <- list("Group.variable" = var.id, "Group1" = group1,
                         "Group2" = group2, "Type" = type, "Mod.avg.table" = AICc.out, "Mod.avg.eff" = Mod.avg.out[,1], 
                         "Uncond.se" = Mod.avg.out[,2], "Conf.level" = conf.level, "Lower.CL" = Lower.CL,
                         "Upper.CL" = Upper.CL, "Matrix.output" = predsOutMat)

    class(Mod.eff.list) <- c("modavgEffect", "list")
    return(Mod.eff.list)  
}



print.modavgEffect <- function(x, digits = 2, ...) {

    ##rework Group.variable labels
    old.type <- x$Group.variable
    stripped.type <- unlist(strsplit(old.type, split = "\\("))

    ic <- colnames(x$Mod.avg.table)[3]
  
    cat("\nModel-averaged effect size on the", x$Type, "scale based on entire model set:\n\n")
  
    ##extract elements
    if(length(stripped.type) == 1) {
        cat("\nMultimodel inference on \"", paste(x$Group.variable, x$Group1, sep = ""), " - ",
            paste(x$Group.variable, x$Group2, sep = ""), "\" based on ", ic, "\n", sep = "")
        
        ##if unmarkedFit model, then print differently
    } else {
        ##extract parameter name
        parm.type <- gsub("(^ +)|( +$)", "", stripped.type[1])
    
        ##extract Group.variable name
        var.id <- gsub("(^ +)|( +$)", "", unlist(strsplit(stripped.type[2], "\\)"))[1])

        cat("\nMultimodel inference on \"", paste(parm.type, "(", var.id, x$Group1, ")", sep = ""), " - ",
            paste(parm.type, "(", var.id, x$Group2, ")", sep = ""), "\" based on ", ic, "\n", sep = "")
    }

  
    cat("\n", ic, " table used to obtain model-averaged effect size:\n", sep = "")
    oldtab <- x$Mod.avg.table
    if (any(names(oldtab)=="c_hat")) {cat("\t(c-hat estimate = ", oldtab$c_hat[1], ")\n", sep = "")}
    cat("\n")

    ##check if result is a scalar or vector
    if(length(x$Mod.avg.eff) == 1) {
        if (any(names(oldtab)=="c_hat")) {
            nice.tab <- cbind(oldtab[,2], oldtab[,3], oldtab[,4], oldtab[,6],
                              oldtab[,9], oldtab[,10])
        } else {
            nice.tab <- cbind(oldtab[,2], oldtab[,3], oldtab[,4], oldtab[,6],
                              oldtab[,8], oldtab[,9])
        }

        colnames(nice.tab) <- c(colnames(oldtab)[c(2,3,4,6)], paste("Effect(", x$Group1, " - ", x$Group2, ")", sep = ""), "SE")
        rownames(nice.tab) <- oldtab[,1]
        print(round(nice.tab, digits=digits))
        cat("\nModel-averaged effect size:", eval(round(x$Mod.avg.eff, digits=digits)), "\n")
        cat("Unconditional SE:", eval(round(x$Uncond.se, digits=digits)), "\n")
        cat("",x$Conf.level * 100, "% Unconditional confidence interval: ", round(x$Lower.CL, digits=digits),
            ", ", round(x$Upper.CL, digits=digits), "\n\n", sep = "")
    } else { ##if result from occuMulti or occuMS

        ##extract parameter names
        parmNames <- names(x$Mod.avg.eff)
        nparms <- length(parmNames)
        nice.tab <- cbind(oldtab[,2], oldtab[,3], oldtab[,4], oldtab[,6])

        colnames(nice.tab) <- colnames(oldtab)[c(2,3,4,6)]
        rownames(nice.tab) <- oldtab[,1]
        print(round(nice.tab, digits=digits))

        if(nparms <= 3) {
            ##iterate over each parameter
            for(k in 1:nparms) {
                cat("\nModel-averaged effect size for ", parmNames[k], ": ", round(x$Mod.avg.eff[k], digits=digits), "\n", sep = "")
                cat("Unconditional SE for ", parmNames[k], ": ", eval(round(x$Uncond.se[k], digits=digits)), "\n", sep = "")
                cat("",x$Conf.level * 100, "% Unconditional confidence interval for ", parmNames[k], ": ",
                    round(x$Lower.CL[k], digits=digits), ", ", round(x$Upper.CL[k], digits=digits), "\n",
                    "---", sep = "")
            }
        } else {
            cat("\n")
            cat("Model-averaged effect sizes:\n\n")
            nice.mat <- x$Matrix.output
            colnames(nice.mat) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
            print(round(nice.mat, digits = digits))
        }
        cat("\n")
    }
}


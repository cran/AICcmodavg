##generic
modavgPred <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                       nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                       ...) {
  cand.set <- formatCands(cand.set)
  UseMethod("modavgPred", cand.set)
}



##default
modavgPred.default <- function(cand.set, modnames = NULL, newdata,
                               second.ord = TRUE, nobs = NULL,
                               uncond.se = "revised", conf.level = 0.95, ...) {
stop("\nFunction not yet defined for this object class\n")
}



##aov
##lm
modavgPred.AICaov.lm <-
  function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95, ...) {

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
    
###################CHANGES####
##############################

    ##determine number of observations in new data set
    nobserv <- dim(newdata)[1]

    ##determine number of columns in new data set
    ncolumns <- dim(newdata)[2]

    ##if only 1 column, add an additional column to avoid problems in computation with predict( )
    if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
 
    ##store AICc table
    AICctab <- aictab(cand.set = cand.set, modnames = modnames,
                      second.ord = second.ord, nobs = nobs, sort = FALSE)

    ##create object to hold Model-averaged estimates and unconditional SE's
    Mod.avg.out <- matrix(NA, nrow = nobserv, ncol = 2)
    colnames(Mod.avg.out) <- c("Mod.avg.est", "Uncond.SE")


    ##begin loop - AICc
    if(second.ord == TRUE){
      for (obs in 1:nobserv) {

        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE,
                                             newdata = newdata[obs, ])$fit))

        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE,
                                            newdata = newdata[obs, ])$se.fit))


        ##create temporary data.frame to store fitted values and SE 
        AICctmp <- AICctab
        AICctmp$fit <- fit
        AICctmp$SE <- SE

        ##compute model averaged prediction and store in output matrix
        Mod.avg.out[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit)
        ##compute unconditional SE and store in output matrix

        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
        if(identical(uncond.se, "old")) {
          Mod.avg.out[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE^2 + (AICctmp$fit - Mod.avg.out[obs, 1])^2))
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
          Mod.avg.out[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE^2 + (AICctmp$fit - Mod.avg.out[obs, 1])^2)))
        }
      }
    }


    ##create temporary data.frame to store fitted values and SE - AIC
    if(second.ord == FALSE) {
      for (obs in 1:nobserv) {

        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE,
                                             newdata = newdata[obs, ])$fit))
        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE,
                                            newdata = newdata[obs, ])$se.fit))

        AICtmp <- AICctab
        AICtmp$fit <- fit
        AICtmp$SE <- SE

        ##compute model averaged prediction and store in output matrix
        Mod.avg.out[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit)

        ##compute unconditional SE and store in output matrix
        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
        if(identical(uncond.se, "old")) {
          Mod.avg.out[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE^2 + (AICtmp$fit - Mod.avg.out[obs, 1])^2))
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
          Mod.avg.out[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE^2 + (AICtmp$fit - Mod.avg.out[obs, 1])^2)))
        }  
        
      }
    }

  type <- "response"

  
  ##############changes start here
  mod.avg.pred <- Mod.avg.out[, 1]
  uncond.se <- Mod.avg.out[, 2]
  
  ##compute confidence interval
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  lower.CL <- mod.avg.pred - zcrit * uncond.se
  upper.CL <- mod.avg.pred + zcrit * uncond.se

  ##create matrix
  matrix.output <- matrix(data = c(mod.avg.pred, uncond.se, lower.CL, upper.CL),
                          nrow = nrow(newdata), ncol = 4)
  colnames(matrix.output) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
  
  ##organize as list
  Mod.pred.list <- list("type" = type, "mod.avg.pred" = mod.avg.pred, "uncond.se" = uncond.se,
                        "conf.level" = conf.level, "lower.CL" = lower.CL, "upper.CL" = upper.CL,
                        "matrix.output" = matrix.output)
  class(Mod.pred.list) <- c("modavgPred", "list")
  return(Mod.pred.list)
}



##glm
modavgPred.AICglm.lm <-
  function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           type = "response", c.hat = 1, gamdisp = NULL, ...) {

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
  if(type == "terms") {stop("\nThe terms argument is not defined for this function\n")}

  ##check family of glm to avoid problems when requesting predictions with argument 'dispersion'
  fam.type <- unlist(lapply(cand.set, FUN=function(i) family(i)$family))
  fam.unique <- unique(fam.type)
  if(identical(fam.unique, "gaussian")) {
    dispersion <- NULL  #set to NULL if gaussian is used
  } else{dispersion <- c.hat}
  ##poisson, binomial, and negative binomial defaults to 1 (no separate parameter for variance)

  ##for negative binomial - reset to NULL
  if(any(regexpr("Negative Binomial", fam.type) != -1)) {
    dispersion <- NULL
    ##check for mixture of negative binomial and other
    ##number of models with negative binomial
    negbin.num <- sum(regexpr("Negative Binomial", fam.type) != -1)
    if(negbin.num < length(fam.type)) {
      stop("Function does not support mixture of negative binomial with other distributions in model set")
    }
  }
  
     
###################CHANGES####
##############################
  if(c.hat > 1) {dispersion <- c.hat }
  if(!is.null(gamdisp)) {dispersion <- gamdisp}
  if(c.hat > 1 && !is.null(gamdisp)) {stop("\nYou cannot specify values for both \'c.hat\' and \'gamdisp\'\n")}
  ##dispersion is the dispersion parameter - this influences the SE's (to specify dispersion parameter for either overdispersed Poisson or Gamma glm)
  ##type enables to specify either "response" (original scale = point estimate) or "link" (linear predictor)
  

  ##check that link function is the same for all models if linear predictor is used
  check.link <- unlist(lapply(X = cand.set, FUN = function(i) i$family$link))
  unique.link <- unique(x = check.link)
  if(identical(type, "link")) {
    if(length(unique.link) > 1) {stop("\nIt is not appropriate to compute a model averaged linear predictor\n",
                                      "with different link functions\n")}
  }
  ##extract inverse link function
  link.inv <- unlist(lapply(X = cand.set, FUN = function(i) i$family$linkinv))[[1]]

    ##check if model uses gamma distribution
    gam1 <- unlist(lapply(cand.set, FUN = function(i) family(i)$family[1] == "Gamma")) #check for gamma regression models
    ##correct SE's for estimates of gamma regressions when gamdisp is specified
    if(any(gam1) == TRUE)  {
      ##check for specification of gamdisp argument
      if(is.null(gamdisp)) stop("\nYou must specify a gamma dispersion parameter with gamma generalized linear models\n")
    }
  
    
    ##determine number of observations in new data set
    nobserv <- dim(newdata)[1]

    ##determine number of columns in new data set
    ncolumns <- dim(newdata)[2]

    ##if only 1 column, add an additional column to avoid problems in computation with predict( )
    if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
 
    ##store AICc table
    AICctab <- aictab(cand.set = cand.set, modnames = modnames,
                      second.ord = second.ord, nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = nobserv, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.est", "Uncond.SE")

  ##required for CI
  if(identical(type, "response")) {
    Mod.avg.out.link <- matrix(NA, nrow = nobserv, ncol = 2)
    colnames(Mod.avg.out.link) <- c("Mod.avg.est", "Uncond.SE")
  }
  
    ##begin loop - AICc
    if(second.ord == TRUE && c.hat == 1){
      for (obs in 1:nobserv) {

        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = type, dispersion = dispersion)$fit))
                
        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                            type = type, dispersion = dispersion)$se.fit))

        ##create temporary data.frame to store fitted values and SE 
        AICctmp <- AICctab
        AICctmp$fit <- fit
        AICctmp$SE <- SE

        
        ##required for CI
        if(identical(type, "response")) {
          if(length(unique.link) == 1) {
            AICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                              type = "link", dispersion = dispersion)$fit))
            AICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                             type = "link", dispersion = dispersion)$se.fit))
          } else {
            warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
            AICctmp$fit.link <- NA
            AICctmp$SE.link <- NA
          }
        }
        
        
        
        ##compute model averaged prediction and store in output matrix
        Mod.avg.out[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit)

        ##compute unconditional SE and store in output matrix
        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
        if(identical(uncond.se, "old")) {
          Mod.avg.out[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE^2 + (AICctmp$fit - Mod.avg.out[obs, 1])^2))
          if(identical(type, "response")) {
            Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
            Mod.avg.out.link[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
          }
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
          Mod.avg.out[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE^2 + (AICctmp$fit - Mod.avg.out[obs, 1])^2)))
          if(identical(type, "response")) {
            Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
            Mod.avg.out.link[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
          }
        }
      }
    }




    ##create temporary data.frame to store fitted values and SE - QAICc
    if(second.ord==TRUE && c.hat > 1) {
      for (obs in 1:nobserv) {

        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN=function(i) predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = type, dispersion = dispersion)$fit))
        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                            type = type, dispersion = dispersion)$se.fit))

        QAICctmp <- AICctab
        QAICctmp$fit <- fit
        QAICctmp$SE <- SE

        ##required for CI
        if(identical(type, "response")) {
          if(length(unique.link) == 1) {
            QAICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                               type = "link", dispersion = dispersion)$fit))
            QAICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                              type = "link", dispersion = dispersion)$se.fit))
          } else {
            warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
            QAICctmp$fit.link <- NA
            QAICctmp$SE.link <- NA
          }
        }


        
        ##compute model averaged prediction and store in output matrix
        Mod.avg.out[obs, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$fit)

        ##compute unconditional SE and store in output matrix
        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
        if(identical(uncond.se, "old")) {
          Mod.avg.out[obs, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE^2 + (QAICctmp$fit - Mod.avg.out[obs, 1])^2))
          if(identical(type, "response")) {
            Mod.avg.out.link[obs, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$fit.link)
            Mod.avg.out.link[obs, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.link^2 + (QAICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
          }
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
          Mod.avg.out[obs, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE^2 + (QAICctmp$fit - Mod.avg.out[obs, 1])^2)))
          if(identical(type, "response")) {
            Mod.avg.out.link[obs, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$fit.link)
            Mod.avg.out.link[obs, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.link^2 + (QAICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
          }
        }
      }
    }




    ##create temporary data.frame to store fitted values and SE - AIC
    if(second.ord == FALSE && c.hat == 1) {
      for (obs in 1:nobserv) {

        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = type, dispersion = dispersion)$fit))
        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                            type = type, dispersion = dispersion)$se.fit))

        AICtmp <- AICctab
        AICtmp$fit <- fit
        AICtmp$SE <- SE

        ##required for CI
        if(identical(type, "response")) {
          if(length(unique.link) == 1) {
            AICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                              type = "link", dispersion = dispersion)$fit))
            AICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                            type = "link", dispersion = dispersion)$se.fit))

          } else {
            warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
            AICtmp$fit.link <- NA
            AICtmp$SE.link <- NA
          }
        }
        
        ##compute model averaged prediction and store in output matrix
        Mod.avg.out[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit)

        ##compute unconditional SE and store in output matrix
        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
        if(identical(uncond.se, "old")) {
          Mod.avg.out[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE^2 + (AICtmp$fit - Mod.avg.out[obs, 1])^2))
          if(identical(type, "response")) {
            Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
            Mod.avg.out.link[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
          }
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
          Mod.avg.out[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE^2 + (AICtmp$fit - Mod.avg.out[obs, 1])^2)))
          if(identical(type, "response")) {
            Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
            Mod.avg.out.link[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
          }
        }  
        
      }
    }




    ##create temporary data.frame to store fitted values and SE - QAIC
    if(second.ord == FALSE && c.hat > 1) {
      for (obs in 1:nobserv) {

        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                         type = type, dispersion = dispersion)$fit))
        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                        type = type, dispersion = dispersion)$se.fit))

        QAICtmp <- AICctab
        QAICtmp$fit <- fit
        QAICtmp$SE <- SE

        ##required for CI
        if(identical(type, "response")) {
          if(length(unique.link) == 1) {
            QAICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                              type = "link", dispersion = dispersion)$fit))
            QAICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                             type = "link", dispersion = dispersion)$se.fit))
          } else {
            warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
            QAICtmp$fit.link <- NA
            QAICtmp$SE.link <- NA
          }
        }
        
        ##compute model averaged prediction and store in output matrix
        Mod.avg.out[obs, 1] <- sum(QAICtmp$QAICWt*QAICtmp$fit)

        ##compute unconditional SE and store in output matrix
        if(identical(uncond.se, "old")) {
          Mod.avg.out[obs, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2))
          if(identical(type, "response")) {
            Mod.avg.out.link[obs, 1] <- sum(QAICtmp$QAICWt*QAICtmp$fit.link)
            Mod.avg.out.link[obs, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.link^2 + (QAICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
          }
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
          Mod.avg.out[obs, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2)))
          if(identical(type, "response")) {
            Mod.avg.out.link[obs, 1] <- sum(QAICtmp$QAICWt*QAICtmp$fit.link)
            Mod.avg.out.link[obs, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.link^2 + (QAICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
          }
        }  
      }
    }

##############changes start here
  mod.avg.pred <- Mod.avg.out[, 1]
  uncond.se <- Mod.avg.out[, 2]
  
  ##compute confidence interval
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)

  if(identical(type, "link")) {
    lower.CL <- mod.avg.pred - zcrit * uncond.se
    upper.CL <- mod.avg.pred + zcrit * uncond.se
  } else {
    lower.CL <- link.inv(Mod.avg.out.link[, 1] - zcrit * Mod.avg.out.link[, 2])
    upper.CL <- link.inv(Mod.avg.out.link[, 1] + zcrit * Mod.avg.out.link[, 2])
  }
  
  ##create matrix
  matrix.output <- matrix(data = c(mod.avg.pred, uncond.se, lower.CL, upper.CL),
                          nrow = nrow(newdata), ncol = 4)
  colnames(matrix.output) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
  
  ##organize as list
  Mod.pred.list <- list("type" = type, "mod.avg.pred" = mod.avg.pred, "uncond.se" = uncond.se,
                        "conf.level" = conf.level, "lower.CL" = lower.CL, "upper.CL" = upper.CL,
                        "matrix.output" = matrix.output)
  class(Mod.pred.list) <- c("modavgPred", "list")
  return(Mod.pred.list)
}



##lm
modavgPred.AIClm <-
  function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95, ...) {

  ##check if named list if modnames are not supplied
  if(is.null(modnames)) {
    if(is.null(names(cand.set))) {
      modnames <- paste("Mod", 1:length(cand.set), sep = "")
      warning("\nModel names have been supplied automatically in the table\n")
    } else {
      modnames <- names(cand.set)
    }
  }
    
    
###################CHANGES####
##############################

    ##determine number of observations in new data set
    nobserv <- dim(newdata)[1]

    ##determine number of columns in new data set
    ncolumns <- dim(newdata)[2]

    ##if only 1 column, add an additional column to avoid problems in computation with predict( )
    if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
 
    ##store AICc table
    AICctab <- aictab(cand.set = cand.set, modnames = modnames,
                      second.ord = second.ord, nobs = nobs, sort = FALSE)

    ##create object to hold Model-averaged estimates and unconditional SE's
    Mod.avg.out <- matrix(NA, nrow = nobserv, ncol = 2)
    colnames(Mod.avg.out) <- c("Mod.avg.est", "Uncond.SE")


    ##begin loop - AICc
    if(second.ord == TRUE){
      for (obs in 1:nobserv) {

        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE,
                                             newdata = newdata[obs, ])$fit))

        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE,
                                            newdata = newdata[obs, ])$se.fit))


        ##create temporary data.frame to store fitted values and SE 
        AICctmp <- AICctab
        AICctmp$fit <- fit
        AICctmp$SE <- SE

        ##compute model averaged prediction and store in output matrix
        Mod.avg.out[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit)
        ##compute unconditional SE and store in output matrix

        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
        if(identical(uncond.se, "old")) {
          Mod.avg.out[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE^2 + (AICctmp$fit - Mod.avg.out[obs, 1])^2))
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
          Mod.avg.out[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE^2 + (AICctmp$fit - Mod.avg.out[obs, 1])^2)))
        }
      }
    }


    ##create temporary data.frame to store fitted values and SE - AIC
    if(second.ord == FALSE) {
      for (obs in 1:nobserv) {

        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE,
                                             newdata = newdata[obs, ])$fit))
        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE,
                                            newdata = newdata[obs, ])$se.fit))

        AICtmp <- AICctab
        AICtmp$fit <- fit
        AICtmp$SE <- SE

        ##compute model averaged prediction and store in output matrix
        Mod.avg.out[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit)

        ##compute unconditional SE and store in output matrix
        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
        if(identical(uncond.se, "old")) {
          Mod.avg.out[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE^2 + (AICtmp$fit - Mod.avg.out[obs, 1])^2))
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
          Mod.avg.out[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE^2 + (AICtmp$fit - Mod.avg.out[obs, 1])^2)))
        }  
        
      }
    }

  type <- "response"

  ##############changes start here
  mod.avg.pred <- Mod.avg.out[, 1]
  uncond.se <- Mod.avg.out[, 2]
  
  ##compute confidence interval
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  lower.CL <- mod.avg.pred - zcrit * uncond.se
  upper.CL <- mod.avg.pred + zcrit * uncond.se

  ##create matrix
  matrix.output <- matrix(data = c(mod.avg.pred, uncond.se, lower.CL, upper.CL),
                          nrow = nrow(newdata), ncol = 4)
  colnames(matrix.output) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
  
  ##organize as list
  Mod.pred.list <- list("type" = type, "mod.avg.pred" = mod.avg.pred, "uncond.se" = uncond.se,
                        "conf.level" = conf.level, "lower.CL" = lower.CL, "upper.CL" = upper.CL,
                        "matrix.output" = matrix.output)
  class(Mod.pred.list) <- c("modavgPred", "list")
  return(Mod.pred.list)
}



##gls
modavgPred.AICgls <-
function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
         nobs = NULL, uncond.se = "revised", conf.level = 0.95, ...) {

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

  ##determine number of observations in new data set
  nobserv <- dim(newdata)[1]

  ##determine number of columns in new data set
  ncolumns <- dim(newdata)[2]

  ##if only 1 column, add an additional column to avoid problems in computation with predictSE.mer( )
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
 
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs, sort = FALSE)

    ##create object to hold Model-averaged estimates and unconditional SE's
    Mod.avg.out <- matrix(NA, nrow = nobserv, ncol = 2)
    colnames(Mod.avg.out) <- c("Mod.avg.est", "Uncond.SE")


    ##begin loop - AICc
    if(second.ord==TRUE){
      for (obs in 1:nobserv) {

        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i) predictSE(i, se.fit = TRUE, newdata = newdata[obs, ])$fit))

        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i) predictSE(i, se.fit = TRUE, newdata = newdata[obs, ])$se.fit))


        ##create temporary data.frame to store fitted values and SE 
        AICctmp <- AICctab
        AICctmp$fit <- fit
        AICctmp$SE <- SE

        ##compute model averaged prediction and store in output matrix
        Mod.avg.out[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit)
        ##compute unconditional SE and store in output matrix

        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
        if(identical(uncond.se, "old")) {
          Mod.avg.out[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE^2 + (AICctmp$fit - Mod.avg.out[obs, 1])^2))
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
          Mod.avg.out[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2)))
        }
      }
    }



    ##create temporary data.frame to store fitted values and SE - AIC
    if(second.ord == FALSE) {
      for (obs in 1:nobserv) {

        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i) predictSE(i, se.fit = TRUE, newdata = newdata[obs, ])$fit))
        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i) predictSE(i, se.fit = TRUE, newdata = newdata[obs, ])$se.fit))

        AICtmp <- AICctab
        AICtmp$fit <- fit
        AICtmp$SE <- SE

        ##compute model averaged prediction and store in output matrix
        Mod.avg.out[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit)

        ##compute unconditional SE and store in output matrix
        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
        if(identical(uncond.se, "old")) {
          Mod.avg.out[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE^2 + (AICtmp$fit - Mod.avg.out[obs, 1])^2))
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
          Mod.avg.out[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE^2 + (AICtmp$fit - Mod.avg.out[obs, 1])^2)))
        }  
        
      }
    }

  type <- "response"

  
  ##############changes start here
  mod.avg.pred <- Mod.avg.out[, 1]
  uncond.se <- Mod.avg.out[, 2]
  
  ##compute confidence interval
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  lower.CL <- mod.avg.pred - zcrit * uncond.se
  upper.CL <- mod.avg.pred + zcrit * uncond.se

  ##create matrix
  matrix.output <- matrix(data = c(mod.avg.pred, uncond.se, lower.CL, upper.CL),
                          nrow = nrow(newdata), ncol = 4)
  colnames(matrix.output) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
  
  ##organize as list
  Mod.pred.list <- list("type" = type, "mod.avg.pred" = mod.avg.pred, "uncond.se" = uncond.se,
                        "conf.level" = conf.level, "lower.CL" = lower.CL, "upper.CL" = upper.CL,
                        "matrix.output" = matrix.output)
  class(Mod.pred.list) <- c("modavgPred", "list")
  return(Mod.pred.list)
}



##lme
modavgPred.AIClme <-
  function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95, ...) {

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

  ##determine number of observations in new data set
  nobserv <- dim(newdata)[1]

  ##determine number of columns in new data set
  ncolumns <- dim(newdata)[2]

  ##if only 1 column, add an additional column to avoid problems in computation with predictSE( )
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord, nobs = nobs, sort = FALSE)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = nobserv, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.est", "Uncond.SE")

  ##begin loop - AICc
  if(second.ord==TRUE){
    for (obs in 1:nobserv) {

      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i) predictSE(i, se.fit = TRUE, newdata = newdata[obs, ])$fit))

      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i) predictSE(i, se.fit = TRUE, newdata = newdata[obs, ])$se.fit))


      ##create temporary data.frame to store fitted values and SE 
      AICctmp <- AICctab
      AICctmp$fit <- fit
      AICctmp$SE <- SE

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE^2 + (AICctmp$fit - Mod.avg.out[obs, 1])^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE^2 + (AICctmp$fit - Mod.avg.out[obs, 1])^2)))
      }
    }
  }



  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord == FALSE) {
    for (obs in 1:nobserv) {

      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i) predictSE(i, se.fit = TRUE, newdata = newdata[obs, ])$fit))
      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i) predictSE(i, se.fit = TRUE, newdata = newdata[obs, ])$se.fit))

      AICtmp <- AICctab
      AICtmp$fit <- fit
      AICtmp$SE <- SE

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE^2 + (AICtmp$fit - Mod.avg.out[obs, 1])^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE^2 + (AICtmp$fit - Mod.avg.out[obs, 1])^2)))
      }  
        
    }
  }

  type <- "response"

  
  ##############changes start here
  mod.avg.pred <- Mod.avg.out[, 1]
  uncond.se <- Mod.avg.out[, 2]
  
  ##compute confidence interval
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  lower.CL <- mod.avg.pred - zcrit * uncond.se
  upper.CL <- mod.avg.pred + zcrit * uncond.se

  ##create matrix
  matrix.output <- matrix(data = c(mod.avg.pred, uncond.se, lower.CL, upper.CL),
                          nrow = nrow(newdata), ncol = 4)
  colnames(matrix.output) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
  
  ##organize as list
  Mod.pred.list <- list("type" = type, "mod.avg.pred" = mod.avg.pred, "uncond.se" = uncond.se,
                        "conf.level" = conf.level, "lower.CL" = lower.CL, "upper.CL" = upper.CL,
                        "matrix.output" = matrix.output)
  class(Mod.pred.list) <- c("modavgPred", "list")
  return(Mod.pred.list)
}


##mer
modavgPred.AICmer <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                              nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                              type = "response", c.hat = 1, ...) {

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
  if(c.hat != 1) {warning("\nThis function only allows \'c.hat = 1\' for \'mer\' class objects\n")}

  ##check that link function is the same for all models if linear predictor is used
  check.link <- unlist(lapply(X = cand.set, FUN = function(i) fam.link.mer(i)$link))
  unique.link <- unique(check.link)
  if(identical(type, "link")) {
    if(length(unique.link) > 1) stop("\nIt is not appropriate to compute a model-averaged beta estimate\n",
                                    "from models using different link functions\n")
  }
  ##extract inverse link function
  link.inv <- unlist(lapply(X = cand.set, FUN = function(i) i@resp$family$linkinv))[[1]]
  
  ##determine number of observations in data set
  nobserv <- dim(newdata)[1]

  ##determine number of columns in data set
  ncolumns <- dim(newdata)[2]

  ##if only 1 column, add an additional column to avoid problems in computation with predictSE.mer( )
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = nobserv, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.est", "Uncond.SE")

  ##required for CI
  if(identical(type, "response")) {
    Mod.avg.out.link <- matrix(NA, nrow = nobserv, ncol = 2)
    colnames(Mod.avg.out.link) <- c("Mod.avg.est", "Uncond.SE")
  }

  
  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
    for (obs in 1:nobserv) {

      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i) predictSE(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = type, level = 0)$fit))

      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i) predictSE(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = type, level = 0)$se.fit))


      ##create temporary data.frame to store fitted values and SE 
      AICctmp <- AICctab
      AICctmp$fit <- fit
      AICctmp$SE <- SE

      
      ##required for CI
      if(identical(type, "response")) {
        if(length(unique.link) == 1) {
          AICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i) predictSE(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                            type = "link", level = 0)$fit))
          AICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i) predictSE(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                           type = "link", level = 0)$se.fit))
        } else {
          warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
          AICctmp$fit.link <- NA
          AICctmp$SE.link <- NA
        }
      }


      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE^2 + (AICctmp$fit - Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE^2 + (AICctmp$fit - Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }
    }
  }




  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord == FALSE && c.hat == 1) {
    for (obs in 1:nobserv) {
      
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i) predictSE(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = type, level = 0)$fit))
      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i) predictSE(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = type, level = 0)$se.fit))

      AICtmp <- AICctab
      AICtmp$fit <- fit
      AICtmp$SE <- SE


      ##required for CI
      if(identical(type, "response")) {
        if(length(unique.link) == 1) {
          AICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i) predictSE(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                           type = "link", level = 0)$fit))
          AICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i) predictSE(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                          type = "link", level = 0)$se.fit))

          } else {
            warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
            AICtmp$fit.link <- NA
            AICtmp$SE.link <- NA
          }
        }

      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit)

      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE^2 + (AICtmp$fit - Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE^2 + (AICtmp$fit - Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }  
    }
  }

  ##############changes start here
  mod.avg.pred <- Mod.avg.out[, 1]
  uncond.se <- Mod.avg.out[, 2]
  
  ##compute confidence interval
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)

  if(identical(type, "link")) {
    lower.CL <- mod.avg.pred - zcrit * uncond.se
    upper.CL <- mod.avg.pred + zcrit * uncond.se
  } else {
    lower.CL <- link.inv(Mod.avg.out.link[, 1] - zcrit * Mod.avg.out.link[, 2])
    upper.CL <- link.inv(Mod.avg.out.link[, 1] + zcrit * Mod.avg.out.link[, 2])
  }
  
  ##create matrix
  matrix.output <- matrix(data = c(mod.avg.pred, uncond.se, lower.CL, upper.CL),
                          nrow = nrow(newdata), ncol = 4)
  colnames(matrix.output) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
  
  ##organize as list
  Mod.pred.list <- list("type" = type, "mod.avg.pred" = mod.avg.pred, "uncond.se" = uncond.se,
                        "conf.level" = conf.level, "lower.CL" = lower.CL, "upper.CL" = upper.CL,
                        "matrix.output" = matrix.output)
  class(Mod.pred.list) <- c("modavgPred", "list")
  return(Mod.pred.list)
}



##glmerMod
modavgPred.AICglmerMod <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                   nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                   type = "response", c.hat = 1, ...) {

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
    

  if(c.hat != 1) {warning("\nThis function only allows \'c.hat = 1\' for \'glmerMod\' class objects\n")}

  ##check that link function is the same for all models if linear predictor is used
  check.link <- unlist(lapply(X = cand.set, FUN = function(i) fam.link.mer(i)$link))
  unique.link <- unique(check.link)
  if(identical(type, "link")) {
    if(length(unique.link) > 1) stop("\nIt is not appropriate to compute a model-averaged beta estimate\n",
                                    "from models using different link functions\n")
  }
  ##extract inverse link function
  link.inv <- unlist(lapply(X = cand.set, FUN = function(i) i@resp$family$linkinv))[[1]]
 
  ##determine number of observations in data set
  nobserv <- dim(newdata)[1]

  ##determine number of columns in data set
  ncolumns <- dim(newdata)[2]

  ##if only 1 column, add an additional column to avoid problems in computation with predictSE.mer( )
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = nobserv, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.est", "Uncond.SE")

  ##required for CI
  if(identical(type, "response")) {
    Mod.avg.out.link <- matrix(NA, nrow = nobserv, ncol = 2)
    colnames(Mod.avg.out.link) <- c("Mod.avg.est", "Uncond.SE")
  }

  
  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
    for (obs in 1:nobserv) {

      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i) predictSE(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = type, level = 0)$fit))

      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i) predictSE(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = type, level = 0)$se.fit))


      ##create temporary data.frame to store fitted values and SE 
      AICctmp <- AICctab
      AICctmp$fit <- fit
      AICctmp$SE <- SE

           
      ##required for CI
      if(identical(type, "response")) {
        if(length(unique.link) == 1) {
          AICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i) predictSE(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                            type = "link", level = 0)$fit))
          AICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i) predictSE(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                           type = "link", level = 0)$se.fit))
        } else {
          warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
          AICctmp$fit.link <- NA
          AICctmp$SE.link <- NA
        }
      }

      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE^2 + (AICctmp$fit - Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE^2 + (AICctmp$fit - Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }
    }
  }




  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord == FALSE && c.hat == 1) {
    for (obs in 1:nobserv) {
      
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i) predictSE(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = type, level = 0)$fit))
      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i) predictSE(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = type, level = 0)$se.fit))

      AICtmp <- AICctab
      AICtmp$fit <- fit
      AICtmp$SE <- SE
      
      
      ##required for CI
      if(identical(type, "response")) {
        if(length(unique.link) == 1) {
          AICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i) predictSE(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                           type = "link", level = 0)$fit))
          AICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i) predictSE(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                          type = "link", level = 0)$se.fit))
          
        } else {
          warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
          AICtmp$fit.link <- NA
          AICtmp$SE.link <- NA
        }
      }

      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit)

      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE^2 + (AICtmp$fit - Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE^2 + (AICtmp$fit - Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }  
    }
  }


##############changes start here
  mod.avg.pred <- Mod.avg.out[, 1]
  uncond.se <- Mod.avg.out[, 2]
  
  ##compute confidence interval
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)

  if(identical(type, "link")) {
    lower.CL <- mod.avg.pred - zcrit * uncond.se
    upper.CL <- mod.avg.pred + zcrit * uncond.se
  } else {
    lower.CL <- link.inv(Mod.avg.out.link[, 1] - zcrit * Mod.avg.out.link[, 2])
    upper.CL <- link.inv(Mod.avg.out.link[, 1] + zcrit * Mod.avg.out.link[, 2])
  }
  
  ##create matrix
  matrix.output <- matrix(data = c(mod.avg.pred, uncond.se, lower.CL, upper.CL),
                          nrow = nrow(newdata), ncol = 4)
  colnames(matrix.output) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
  
  ##organize as list
  Mod.pred.list <- list("type" = type, "mod.avg.pred" = mod.avg.pred, "uncond.se" = uncond.se,
                        "conf.level" = conf.level, "lower.CL" = lower.CL, "upper.CL" = upper.CL,
                        "matrix.output" = matrix.output)
  class(Mod.pred.list) <- c("modavgPred", "list")
  return(Mod.pred.list)
}



##lmerMod
modavgPred.AIClmerMod <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                     nobs = NULL, uncond.se = "revised", conf.level = 0.95, ...) {
  
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
    

  ##determine number of observations in data set
  nobserv <- dim(newdata)[1]

  ##determine number of columns in data set
  ncolumns <- dim(newdata)[2]

  ##if only 1 column, add an additional column to avoid problems in computation with predictSE.mer( )
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = nobserv, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.est", "Uncond.SE")


  ##begin loop - AICc
  if(second.ord == TRUE){
    for (obs in 1:nobserv) {

      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i) predictSE(i, se.fit = TRUE, newdata = newdata[obs, ],
                                           level = 0)$fit))
      
      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i) predictSE(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          level = 0)$se.fit))


      ##create temporary data.frame to store fitted values and SE 
      AICctmp <- AICctab
      AICctmp$fit <- fit
      AICctmp$SE <- SE

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE^2 + (AICctmp$fit - Mod.avg.out[obs, 1])^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE^2 + (AICctmp$fit - Mod.avg.out[obs, 1])^2)))
      }
    }
  }




  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord == FALSE) {
    for (obs in 1:nobserv) {
      
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i) predictSE(i, se.fit = TRUE, newdata = newdata[obs, ],
                                           level = 0)$fit))
      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i) predictSE(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          level = 0)$se.fit))

      AICtmp <- AICctab
      AICtmp$fit <- fit
      AICtmp$SE <- SE

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit)

      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE^2 + (AICtmp$fit - Mod.avg.out[obs, 1])^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE^2 + (AICtmp$fit - Mod.avg.out[obs, 1])^2)))
      }  
    }
  }

  type <- "response"

  
  ##############changes start here
  mod.avg.pred <- Mod.avg.out[, 1]
  uncond.se <- Mod.avg.out[, 2]
  
  ##compute confidence interval
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  lower.CL <- mod.avg.pred - zcrit * uncond.se
  upper.CL <- mod.avg.pred + zcrit * uncond.se

  ##create matrix
  matrix.output <- matrix(data = c(mod.avg.pred, uncond.se, lower.CL, upper.CL),
                          nrow = nrow(newdata), ncol = 4)
  colnames(matrix.output) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
  
  ##organize as list
  Mod.pred.list <- list("type" = type, "mod.avg.pred" = mod.avg.pred, "uncond.se" = uncond.se,
                        "conf.level" = conf.level, "lower.CL" = lower.CL, "upper.CL" = upper.CL,
                        "matrix.output" = matrix.output)
  class(Mod.pred.list) <- c("modavgPred", "list")
  return(Mod.pred.list)
}



##lmerModLmerTest
modavgPred.AIClmerModLmerTest <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                          nobs = NULL, uncond.se = "revised", conf.level = 0.95, ...) {
  
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
    

  ##determine number of observations in data set
  nobserv <- dim(newdata)[1]

  ##determine number of columns in data set
  ncolumns <- dim(newdata)[2]

  ##if only 1 column, add an additional column to avoid problems in computation with predictSE.mer( )
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
  
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                    nobs = nobs, sort = FALSE)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = nobserv, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.est", "Uncond.SE")


  ##begin loop - AICc
  if(second.ord == TRUE){
    for (obs in 1:nobserv) {

      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i) predictSE(i, se.fit = TRUE, newdata = newdata[obs, ],
                                           level = 0)$fit))
      
      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i) predictSE(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          level = 0)$se.fit))


      ##create temporary data.frame to store fitted values and SE 
      AICctmp <- AICctab
      AICctmp$fit <- fit
      AICctmp$SE <- SE

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE^2 + (AICctmp$fit - Mod.avg.out[obs, 1])^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE^2 + (AICctmp$fit - Mod.avg.out[obs, 1])^2)))
      }
    }
  }


  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord == FALSE) {
    for (obs in 1:nobserv) {
      
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i) predictSE(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                     level = 0)$fit))
      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i) predictSE(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                    level = 0)$se.fit))

      AICtmp <- AICctab
      AICtmp$fit <- fit
      AICtmp$SE <- SE

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit)

      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE^2 + (AICtmp$fit - Mod.avg.out[obs, 1])^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE^2 + (AICtmp$fit - Mod.avg.out[obs, 1])^2)))
      }  
    }
  }

  type <- "response"

  
  ##############changes start here
  mod.avg.pred <- Mod.avg.out[, 1]
  uncond.se <- Mod.avg.out[, 2]
  
  ##compute confidence interval
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  lower.CL <- mod.avg.pred - zcrit * uncond.se
  upper.CL <- mod.avg.pred + zcrit * uncond.se

  ##create matrix
  matrix.output <- matrix(data = c(mod.avg.pred, uncond.se, lower.CL, upper.CL),
                          nrow = nrow(newdata), ncol = 4)
  colnames(matrix.output) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
  
  ##organize as list
  Mod.pred.list <- list("type" = type, "mod.avg.pred" = mod.avg.pred, "uncond.se" = uncond.se,
                        "conf.level" = conf.level, "lower.CL" = lower.CL, "upper.CL" = upper.CL,
                        "matrix.output" = matrix.output)
  class(Mod.pred.list) <- c("modavgPred", "list")
  return(Mod.pred.list)
}



##glm.nb
modavgPred.AICnegbin.glm.lm <-
  function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
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
  if(type == "terms") {stop("\nThe terms argument is not defined for this function\n")}
  
     
###################CHANGES####
##############################
      ##type enables to specify either "response" (original scale = point estimate) or "link" (linear predictor)
  

      ##check that link function is the same for all models if linear predictor is used
      check.link <- unlist(lapply(X = cand.set, FUN = function(i) i$family$link))
      unique.link <- unique(x = check.link)
      if(identical(type, "link")) {
          if(length(unique.link) > 1) {stop("\nIt is not appropriate to compute a model averaged linear predictor\n",
                                            "with different link functions\n")}
      }
      ##extract inverse link function
      link.inv <- unlist(lapply(X = cand.set, FUN = function(i) i$family$linkinv))[[1]]

      ##determine number of observations in new data set
      nobserv <- dim(newdata)[1]

      ##determine number of columns in new data set
      ncolumns <- dim(newdata)[2]

      ##if only 1 column, add an additional column to avoid problems in computation with predict( )
      if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
 
      ##store AICc table
      AICctab <- aictab(cand.set = cand.set, modnames = modnames,
                        second.ord = second.ord, nobs = nobs, sort = FALSE)

      ##create object to hold Model-averaged estimates and unconditional SE's
      Mod.avg.out <- matrix(NA, nrow = nobserv, ncol = 2)
      colnames(Mod.avg.out) <- c("Mod.avg.est", "Uncond.SE")

      ##required for CI
      if(identical(type, "response")) {
          Mod.avg.out.link <- matrix(NA, nrow = nobserv, ncol = 2)
          colnames(Mod.avg.out.link) <- c("Mod.avg.est", "Uncond.SE")
      }
  
      ##begin loop - AICc
      if(second.ord == TRUE){
          for (obs in 1:nobserv) {

              ##extract fitted value for observation obs
              fit <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                           type = type)$fit))
                
              ##extract SE for fitted value for observation obs
              SE <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                          type = type)$se.fit))

              ##create temporary data.frame to store fitted values and SE 
              AICctmp <- AICctab
              AICctmp$fit <- fit
              AICctmp$SE <- SE

        
              ##required for CI
              if(identical(type, "response")) {
                  if(length(unique.link) == 1) {
                      AICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE,
                                                                                                newdata = newdata[obs, ],
                                                                                                type = "link")$fit))
                      AICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE,
                                                                                               newdata = newdata[obs, ],
                                                                                               type = "link")$se.fit))
                  } else {
                      warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
                      AICctmp$fit.link <- NA
                      AICctmp$SE.link <- NA
                  }
              }
        
              
              ##compute model averaged prediction and store in output matrix
              Mod.avg.out[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit)

              ##compute unconditional SE and store in output matrix
              ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
              if(identical(uncond.se, "old")) {
                  Mod.avg.out[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE^2 + (AICctmp$fit - Mod.avg.out[obs, 1])^2))
                  if(identical(type, "response")) {
                      Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
                      Mod.avg.out.link[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
                  }
              }

              ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
              if(identical(uncond.se, "revised")) {
                  Mod.avg.out[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE^2 + (AICctmp$fit - Mod.avg.out[obs, 1])^2)))
                  if(identical(type, "response")) {
                      Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
                      Mod.avg.out.link[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
                  }
              }
          }
      }


    ##create temporary data.frame to store fitted values and SE - AIC
      if(second.ord == FALSE) {
          for (obs in 1:nobserv) {

              ##extract fitted value for observation obs
              fit <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                           type = type)$fit))
              ##extract SE for fitted value for observation obs
              SE <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                          type = type)$se.fit))

              AICtmp <- AICctab
              AICtmp$fit <- fit
              AICtmp$SE <- SE

              ##required for CI
              if(identical(type, "response")) {
                  if(length(unique.link) == 1) {
                      AICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE,
                                                                                               newdata = newdata[obs, ],
                                                                                               type = "link")$fit))
                      AICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE,
                                                                                              newdata = newdata[obs, ],
                                                                                              type = "link")$se.fit))

                  } else {
                      warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
                      AICtmp$fit.link <- NA
                      AICtmp$SE.link <- NA
                  }
              }
        
              ##compute model averaged prediction and store in output matrix
              Mod.avg.out[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit)

              ##compute unconditional SE and store in output matrix
              ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
              if(identical(uncond.se, "old")) {
                  Mod.avg.out[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE^2 + (AICtmp$fit - Mod.avg.out[obs, 1])^2))
                  if(identical(type, "response")) {
                      Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
                      Mod.avg.out.link[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
                  }
              }

              ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
              if(identical(uncond.se, "revised")) {
                  Mod.avg.out[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE^2 + (AICtmp$fit - Mod.avg.out[obs, 1])^2)))
                  if(identical(type, "response")) {
                      Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
                      Mod.avg.out.link[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
                  }
              }  
              
          }
      }




##############changes start here
      mod.avg.pred <- Mod.avg.out[, 1]
      uncond.se <- Mod.avg.out[, 2]
  
      ##compute confidence interval
      zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)

      if(identical(type, "link")) {
          lower.CL <- mod.avg.pred - zcrit * uncond.se
          upper.CL <- mod.avg.pred + zcrit * uncond.se
      } else {
          lower.CL <- link.inv(Mod.avg.out.link[, 1] - zcrit * Mod.avg.out.link[, 2])
          upper.CL <- link.inv(Mod.avg.out.link[, 1] + zcrit * Mod.avg.out.link[, 2])
      }
  
      ##create matrix
      matrix.output <- matrix(data = c(mod.avg.pred, uncond.se, lower.CL, upper.CL),
                              nrow = nrow(newdata), ncol = 4)
      colnames(matrix.output) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
  
      ##organize as list
      Mod.pred.list <- list("type" = type, "mod.avg.pred" = mod.avg.pred, "uncond.se" = uncond.se,
                            "conf.level" = conf.level, "lower.CL" = lower.CL, "upper.CL" = upper.CL,
                            "matrix.output" = matrix.output)
      class(Mod.pred.list) <- c("modavgPred", "list")
      return(Mod.pred.list)
  }



##rlm
modavgPred.AICrlm.lm <-
  function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95, ...) {
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
    nobserv <- dim(newdata)[1]

    ##determine number of columns in new data set
    ncolumns <- dim(newdata)[2]

    ##if only 1 column, add an additional column to avoid problems in computation with predictSE.mer( )
    if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
 
    ##store AICc table
    AICctab <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                      nobs = nobs, sort = FALSE)

    ##create object to hold Model-averaged estimates and unconditional SE's
    Mod.avg.out <- matrix(NA, nrow = nobserv, ncol = 2)
    colnames(Mod.avg.out) <- c("Mod.avg.est", "Uncond.SE")

    ##begin loop - AICc
    if(second.ord == TRUE){
      for (obs in 1:nobserv) {

        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE, newdata = newdata[obs, ])$fit))

        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE, newdata = newdata[obs, ])$se.fit))

        ##create temporary data.frame to store fitted values and SE 
        AICctmp <- AICctab
        AICctmp$fit <- fit
        AICctmp$SE <- SE

        ##compute model averaged prediction and store in output matrix
        Mod.avg.out[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit)
        ##compute unconditional SE and store in output matrix

        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
        if(identical(uncond.se, "old")) {
          Mod.avg.out[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE^2 + (AICctmp$fit - Mod.avg.out[obs, 1])^2))
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
          Mod.avg.out[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE^2 + (AICctmp$fit - Mod.avg.out[obs, 1])^2)))
        }
      }
    }


    ##create temporary data.frame to store fitted values and SE - AIC
    if(second.ord == FALSE) {
      for (obs in 1:nobserv) {

        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE, newdata = newdata[obs, ])$fit))
        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE, newdata = newdata[obs, ])$se.fit))

        AICtmp <- AICctab
        AICtmp$fit <- fit
        AICtmp$SE <- SE

        ##compute model averaged prediction and store in output matrix
        Mod.avg.out[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit)

        ##compute unconditional SE and store in output matrix
        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
        if(identical(uncond.se, "old")) {
          Mod.avg.out[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE^2 + (AICtmp$fit - Mod.avg.out[obs, 1])^2))
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
          Mod.avg.out[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE^2 + (AICtmp$fit - Mod.avg.out[obs, 1])^2)))
        }  
        
      }
    }

    type <- "response"

  
##############changes start here
    mod.avg.pred <- Mod.avg.out[, 1]
    uncond.se <- Mod.avg.out[, 2]
  
    ##compute confidence interval
    zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
    lower.CL <- mod.avg.pred - zcrit * uncond.se
    upper.CL <- mod.avg.pred + zcrit * uncond.se
    
    ##create matrix
    matrix.output <- matrix(data = c(mod.avg.pred, uncond.se, lower.CL, upper.CL),
                            nrow = nrow(newdata), ncol = 4)
    colnames(matrix.output) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
    
    ##organize as list
    Mod.pred.list <- list("type" = type, "mod.avg.pred" = mod.avg.pred, "uncond.se" = uncond.se,
                          "conf.level" = conf.level, "lower.CL" = lower.CL, "upper.CL" = upper.CL,
                          "matrix.output" = matrix.output)
    
    class(Mod.pred.list) <- c("modavgPred", "list")
    return(Mod.pred.list)
}



##survreg
modavgPred.AICsurvreg <-
  function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95, type = "response", ...) {

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
    if(type == "terms") {stop("\nThe terms argument is not defined for this function\n")}

  
###################CHANGES####
##############################
    ##type enables to specify either "response" (original scale = point estimate) or "link" (linear predictor)
    

    ##check that distribution is the same for all models
    check.dist <- sapply(X = cand.set, FUN = function(i) i$dist)
    unique.dist <- unique(x = check.dist)
    if(identical(type, "link")) {
      if(length(unique.dist) > 1) stop("\nFunction does not support model-averaging linear predictors using different distributions\n")
    }
    ##extract inverse link function
    link.inv <- sapply(X = cand.set, FUN = function(i) survreg.distributions[[unique.dist]]$itrans)[[1]]
    
  
    ##determine number of observations in new data set
    nobserv <- dim(newdata)[1]

    ##determine number of columns in new data set
    ncolumns <- dim(newdata)[2]

    ##if only 1 column, add an additional column to avoid problems in computation with predictSE.mer( )
    if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
 
    ##store AICc table
    AICctab <- aictab(cand.set = cand.set, modnames = modnames,
                      second.ord = second.ord, nobs = nobs, sort = FALSE)

    ##create object to hold Model-averaged estimates and unconditional SE's
    Mod.avg.out <- matrix(NA, nrow = nobserv, ncol = 2)
    colnames(Mod.avg.out) <- c("Mod.avg.est", "Uncond.SE")

    ##required for CI
    if(identical(type, "response")) {
      Mod.avg.out.link <- matrix(NA, nrow = nobserv, ncol = 2)
      colnames(Mod.avg.out.link) <- c("Mod.avg.est", "Uncond.SE")
    }

    
    ##begin loop - AICc
    if(second.ord == TRUE){
      for (obs in 1:nobserv) {

        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = type)$fit))

        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                            type = type)$se.fit))


        ##create temporary data.frame to store fitted values and SE 
        AICctmp <- AICctab
        AICctmp$fit <- fit
        AICctmp$SE <- SE


        ##required for CI
        if(identical(type, "response")) {
          if(length(unique.dist) == 1) {
            AICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                              type = "link")$fit))
            AICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                             type = "link")$se.fit))
          } else {
            warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
            AICctmp$fit.link <- NA
            AICctmp$SE.link <- NA
          }
        }

        
        ##compute model averaged prediction and store in output matrix
        Mod.avg.out[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit)
        ##compute unconditional SE and store in output matrix

        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
        if(identical(uncond.se, "old")) {
          Mod.avg.out[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE^2 + (AICctmp$fit - Mod.avg.out[obs, 1])^2))
          if(identical(type, "response")) {
            Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
            Mod.avg.out.link[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
          }
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
          Mod.avg.out[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE^2 + (AICctmp$fit - Mod.avg.out[obs, 1])^2)))
          if(identical(type, "response")) {
            Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
            Mod.avg.out.link[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
          }
        }
      }
    }



    ##create temporary data.frame to store fitted values and SE - AIC
    if(second.ord == FALSE) {
      for (obs in 1:nobserv) {

        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = type)$fit))
        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                            type = type)$se.fit))

        AICtmp <- AICctab
        AICtmp$fit <- fit
        AICtmp$SE <- SE


        ##required for CI
        if(identical(type, "response")) {
          if(length(unique.dist) == 1) {
            AICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                             type = "link")$fit))
            AICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i) predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                            type = "link")$se.fit))
            
          } else {
            warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
            AICtmp$fit.link <- NA
            AICtmp$SE.link <- NA
          }
        }

        ##compute model averaged prediction and store in output matrix
        Mod.avg.out[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit)

        ##compute unconditional SE and store in output matrix
        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
        if(identical(uncond.se, "old")) {
          Mod.avg.out[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE^2 + (AICtmp$fit - Mod.avg.out[obs, 1])^2))
          if(identical(type, "response")) {
            Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
            Mod.avg.out.link[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
          }
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
          Mod.avg.out[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE^2 + (AICtmp$fit - Mod.avg.out[obs, 1])^2)))
          if(identical(type, "response")) {
            Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
            Mod.avg.out.link[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
          }
        }  
      }
    }


##############changes start here
    mod.avg.pred <- Mod.avg.out[, 1]
    uncond.se <- Mod.avg.out[, 2]
  
    ##compute confidence interval
    zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)

    if(identical(type, "link")) {
      lower.CL <- mod.avg.pred - zcrit * uncond.se
      upper.CL <- mod.avg.pred + zcrit * uncond.se
    } else {
      lower.CL <- link.inv(Mod.avg.out.link[, 1] - zcrit * Mod.avg.out.link[, 2])
      upper.CL <- link.inv(Mod.avg.out.link[, 1] + zcrit * Mod.avg.out.link[, 2])
    }
    
    ##create matrix
    matrix.output <- matrix(data = c(mod.avg.pred, uncond.se, lower.CL, upper.CL),
                            nrow = nrow(newdata), ncol = 4)
    colnames(matrix.output) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
  
    ##organize as list
    Mod.pred.list <- list("type" = type, "mod.avg.pred" = mod.avg.pred, "uncond.se" = uncond.se,
                          "conf.level" = conf.level, "lower.CL" = lower.CL, "upper.CL" = upper.CL,
                          "matrix.output" = matrix.output)
    class(Mod.pred.list) <- c("modavgPred", "list")
    return(Mod.pred.list)
  }






##occu
modavgPred.AICunmarkedFitOccu <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                          nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                          type = "response", c.hat = 1, parm.type = NULL, ...) {

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
  if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgPred for details\n")}

  
  ##rename values according to unmarked to extract from object
  ##psi
  if(identical(parm.type, "psi")) {
    parm.type1 <- "state"
    
  }

  
    ##detect
    if(identical(parm.type, "detect")) {
        parm.type1 <- "det"
        
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
    
    ##extract inverse link function
    if(identical(select.link, "logistic")) {
        link.inv <- plogis
    }

    ##extract inverse link function
    if(identical(select.link, "cloglog")) {
        link.inv <- function(x) 1 - exp(-exp(x))
    }

    
    ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  
    ##determine number of observations in new data set
    nobserv <- dim(newdata)[1]
  
  ##determine number of columns in new data set
  ncolumns <- dim(newdata)[2]

  ##if only 1 column, add an additional column to avoid problems in computation with predict( )
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
 
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames,
                    second.ord = second.ord, nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = nobserv, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.est", "Uncond.SE")
  

  ##required for CI
  if(identical(type, "response")) {
    Mod.avg.out.link <- matrix(NA, nrow = nobserv, ncol = 2)
    colnames(Mod.avg.out.link) <- c("Mod.avg.est", "Uncond.SE")
  }
  
  
  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
    
    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      AICctmp <- AICctab
      
      if(identical(type, "response")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1)$Predicted))
        
        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                            type = parm.type1)$SE))
        if(length(unique.link) == 1) {
          AICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))
          AICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                             type = parm.type1, backTransform = FALSE)$SE))
        } else {
          warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
          AICctmp$fit.link <- NA
          AICctmp$SE.link <- NA
        }
      }
    
      
      if(identical(type, "link")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))

        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          type = parm.type1, backTransform = FALSE)$SE))
      }
      

      AICctmp$fit <- fit
      AICctmp$SE <- SE
    
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }        
      }
    }
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
    
    for (obs in 1:nobserv) {
      
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
  
      if(identical(type, "response")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                           type = parm.type1)$Predicted))
      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          type = parm.type1)$SE))
      if(length(unique.link) == 1) {
        QAICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                           type = parm.type1, backTransform = FALSE)$Predicted))
        QAICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                          type = parm.type1, backTransform = FALSE)$SE)) * sqrt(c.hat)
      } else {
        warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
        QAICctmp$fit.link <- NA
        QAICctmp$SE.link <- NA
      }
    }
      
      if(identical(type, "link")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                           type = parm.type1, backTransform = FALSE)$Predicted))
      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          type = parm.type1, backTransform = FALSE)$SE))
    }

      
      QAICctmp$fit <- fit
      QAICctmp$SE <- SE * sqrt(c.hat)
      
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(QAICctmp$QAICcWt * QAICctmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(QAICctmp$QAICcWt * sqrt(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.link^2 + (QAICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(QAICctmp$QAICcWt *(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.link^2 + (QAICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }
    }
  }
  
  
  
  ##begin loop - AIC
  if(second.ord == FALSE && c.hat == 1) {

    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      AICtmp <- AICctab
      
      if(identical(type, "response")) {
          fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                               type = parm.type1)$Predicted))

          ##extract SE for fitted value for observation obs
          SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                              type = parm.type1)$SE))
          if(length(unique.link) == 1) {
            AICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                             type = parm.type1, backTransform = FALSE)$Predicted))
            AICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                            type = parm.type1, backTransform = FALSE)$SE))
            
          } else {
            warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
            AICtmp$fit.link <- NA
            AICtmp$SE.link <- NA
          }
        }
      

      if(identical(type, "link")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))
        
        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                            type = parm.type1, backTransform = FALSE)$SE))
      }
        
      
      AICtmp$fit <- fit
      AICtmp$SE <- SE

      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }  
    }
  }

  
  ##begin loop - QAICc
  if(second.ord == FALSE && c.hat > 1){
     
    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      QAICtmp <- AICctab
      
      if(identical(type, "response")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = parm.type1)$Predicted))

      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = parm.type1)$SE))
      
      if(length(unique.link) == 1) {
        QAICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                          type = parm.type1, backTransform = FALSE)$Predicted))
        QAICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                         type = parm.type1, backTransform = FALSE)$SE)) * sqrt(c.hat)
      } else {
        warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
        QAICtmp$fit.link <- NA
        QAICtmp$SE.link <- NA
      }
    }
      
      if(identical(type, "link")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = parm.type1, backTransform = FALSE)$Predicted))

      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = parm.type1, backTransform = FALSE)$SE))
    }
      
      
      QAICtmp$fit <- fit
      QAICtmp$SE <- SE * sqrt(c.hat)
      
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(QAICtmp$QAICWt * QAICtmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(QAICtmp$QAICWt * sqrt(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICtmp$QAICWt*QAICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.link^2 + (QAICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(QAICtmp$QAICWt *(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICtmp$QAICWt*QAICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.link^2 + (QAICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }
    }
  }

##############changes start here
  mod.avg.pred <- Mod.avg.out[, 1]
  uncond.se <- Mod.avg.out[, 2]
  
  ##compute confidence interval
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)

  if(identical(type, "link")) {
    lower.CL <- mod.avg.pred - zcrit * uncond.se
    upper.CL <- mod.avg.pred + zcrit * uncond.se
  } else {
    lower.CL <- link.inv(Mod.avg.out.link[, 1] - zcrit * Mod.avg.out.link[, 2])
    upper.CL <- link.inv(Mod.avg.out.link[, 1] + zcrit * Mod.avg.out.link[, 2])
  }
  
  ##create matrix
  matrix.output <- matrix(data = c(mod.avg.pred, uncond.se, lower.CL, upper.CL),
                          nrow = nrow(newdata), ncol = 4)
  colnames(matrix.output) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
  
  ##organize as list
  Mod.pred.list <- list("type" = type, "mod.avg.pred" = mod.avg.pred, "uncond.se" = uncond.se,
                        "conf.level" = conf.level, "lower.CL" = lower.CL, "upper.CL" = upper.CL,
                        "matrix.output" = matrix.output)
  class(Mod.pred.list) <- c("modavgPred", "list")
  return(Mod.pred.list)
}



##colext
modavgPred.AICunmarkedFitColExt <-
  function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           type = "response", c.hat = 1, parm.type = NULL, ...) {

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
  if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgPred for details\n")}
    
  ##rename values according to unmarked to extract from object
  ##psi
  if(identical(parm.type, "psi")) {
      parm.type1 <- "psi"
      
  }

      ##gamma
      if(identical(parm.type, "gamma")) {
          parm.type1 <- "col"
      
  }

      ##epsilon
      if(identical(parm.type, "epsilon")) {
          parm.type1 <- "ext"

      }

      ##detect
      if(identical(parm.type, "detect")) {
          parm.type1 <- "det"

      }

      
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
      
      ##extract inverse link function
      if(identical(select.link, "logistic")) {
          link.inv <- plogis
      }
  
      ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  
      ##determine number of observations in new data set
      nobserv <- dim(newdata)[1]
  
  ##determine number of columns in new data set
  ncolumns <- dim(newdata)[2]

  ##if only 1 column, add an additional column to avoid problems in computation with predict( )
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
 
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames,
                    second.ord = second.ord, nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = nobserv, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.est", "Uncond.SE")
  

  ##required for CI
  if(identical(type, "response")) {
    Mod.avg.out.link <- matrix(NA, nrow = nobserv, ncol = 2)
    colnames(Mod.avg.out.link) <- c("Mod.avg.est", "Uncond.SE")
  }
  
  
  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
    
    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      AICctmp <- AICctab
      
      if(identical(type, "response")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1)$Predicted))
        
        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                            type = parm.type1)$SE))
        if(length(unique.link) == 1) {
          AICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))
          AICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                           type = parm.type1, backTransform = FALSE)$SE))
        } else {
          warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
          AICctmp$fit.link <- NA
          AICctmp$SE.link <- NA
        }
      }
    
      
      if(identical(type, "link")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))

        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          type = parm.type1, backTransform = FALSE)$SE))
      }
      

      AICctmp$fit <- fit
      AICctmp$SE <- SE
    
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }        
      }
    }
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
    
    for (obs in 1:nobserv) {
      
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
  
      if(identical(type, "response")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                           type = parm.type1)$Predicted))
      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          type = parm.type1)$SE))
      if(length(unique.link) == 1) {
        QAICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                           type = parm.type1, backTransform = FALSE)$Predicted))
        QAICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                          type = parm.type1, backTransform = FALSE)$SE)) * sqrt(c.hat)
      } else {
        warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
        QAICctmp$fit.link <- NA
        QAICctmp$SE.link <- NA
      }
    }
      
      if(identical(type, "link")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                           type = parm.type1, backTransform = FALSE)$Predicted))
      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          type = parm.type1, backTransform = FALSE)$SE))
    }

      
      QAICctmp$fit <- fit
      QAICctmp$SE <- SE * sqrt(c.hat)
      
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(QAICctmp$QAICcWt * QAICctmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(QAICctmp$QAICcWt * sqrt(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.link^2 + (QAICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(QAICctmp$QAICcWt *(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.link^2 + (QAICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }
    }
  }
  
  
  
  ##begin loop - AIC
  if(second.ord == FALSE && c.hat == 1) {

    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      AICtmp <- AICctab
      
      if(identical(type, "response")) {
          fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                               type = parm.type1)$Predicted))

          ##extract SE for fitted value for observation obs
          SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                              type = parm.type1)$SE))
          if(length(unique.link) == 1) {
            AICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                             type = parm.type1, backTransform = FALSE)$Predicted))
            AICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                            type = parm.type1, backTransform = FALSE)$SE))
            
          } else {
            warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
            AICtmp$fit.link <- NA
            AICtmp$SE.link <- NA
          }
        }
      

      if(identical(type, "link")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))
        
        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                            type = parm.type1, backTransform = FALSE)$SE))
      }
        
      
      AICtmp$fit <- fit
      AICtmp$SE <- SE

      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }  
    }
  }

  
  ##begin loop - QAIC
  if(second.ord == FALSE && c.hat > 1){
     
    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      QAICtmp <- AICctab
      
      if(identical(type, "response")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = parm.type1)$Predicted))

      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = parm.type1)$SE))
      
      if(length(unique.link) == 1) {
        QAICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                          type = parm.type1, backTransform = FALSE)$Predicted))
        QAICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                         type = parm.type1, backTransform = FALSE)$SE)) * sqrt(c.hat)
      } else {
        warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
        QAICtmp$fit.link <- NA
        QAICtmp$SE.link <- NA
      }
    }
      
      if(identical(type, "link")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = parm.type1, backTransform = FALSE)$Predicted))

      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = parm.type1, backTransform = FALSE)$SE))
    }
      
      
      QAICtmp$fit <- fit
      QAICtmp$SE <- SE * sqrt(c.hat)
      
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(QAICtmp$QAICWt * QAICtmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(QAICtmp$QAICWt * sqrt(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICtmp$QAICWt*QAICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.link^2 + (QAICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(QAICtmp$QAICWt *(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICtmp$QAICWt*QAICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.link^2 + (QAICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }
    }
  }

##############changes start here
  mod.avg.pred <- Mod.avg.out[, 1]
  uncond.se <- Mod.avg.out[, 2]
  
  ##compute confidence interval
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)

  if(identical(type, "link")) {
    lower.CL <- mod.avg.pred - zcrit * uncond.se
    upper.CL <- mod.avg.pred + zcrit * uncond.se
  } else {
    lower.CL <- link.inv(Mod.avg.out.link[, 1] - zcrit * Mod.avg.out.link[, 2])
    upper.CL <- link.inv(Mod.avg.out.link[, 1] + zcrit * Mod.avg.out.link[, 2])
  }
  
  ##create matrix
  matrix.output <- matrix(data = c(mod.avg.pred, uncond.se, lower.CL, upper.CL),
                          nrow = nrow(newdata), ncol = 4)
  colnames(matrix.output) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
  
  ##organize as list
  Mod.pred.list <- list("type" = type, "mod.avg.pred" = mod.avg.pred, "uncond.se" = uncond.se,
                        "conf.level" = conf.level, "lower.CL" = lower.CL, "upper.CL" = upper.CL,
                        "matrix.output" = matrix.output)
  class(Mod.pred.list) <- c("modavgPred", "list")
  return(Mod.pred.list)
}



##occuRN
modavgPred.AICunmarkedFitOccuRN <-
  function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           type = "response", c.hat = 1, parm.type = NULL, ...) {

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
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgPred for details\n")}

    
    ##rename values according to unmarked to extract from object

    ##lambda
    if(identical(parm.type, "lambda")) {
        parm.type1 <- "state"
        
    }
  
    ##detect
      if(identical(parm.type, "detect")) {
          parm.type1 <- "det"
          
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

      
      ##extract inverse link function
      if(identical(select.link, "logistic")) {
          link.inv <- plogis
      }
      if(identical(select.link, "exp")) {
          link.inv <- exp
      }
  
    ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  
    ##determine number of observations in new data set
    nobserv <- dim(newdata)[1]
  
    ##determine number of columns in new data set
    ncolumns <- dim(newdata)[2]

    ##if only 1 column, add an additional column to avoid problems in computation with predict( )
    if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
 
    ##store AICc table
    AICctab <- aictab(cand.set = cand.set, modnames = modnames,
                      second.ord = second.ord, nobs = nobs, sort = FALSE, c.hat = c.hat)

    ##create object to hold Model-averaged estimates and unconditional SE's
    Mod.avg.out <- matrix(NA, nrow = nobserv, ncol = 2)
    colnames(Mod.avg.out) <- c("Mod.avg.est", "Uncond.SE")
  

    ##required for CI
    if(identical(type, "response")) {
      Mod.avg.out.link <- matrix(NA, nrow = nobserv, ncol = 2)
      colnames(Mod.avg.out.link) <- c("Mod.avg.est", "Uncond.SE")
    }
  
  
    ##begin loop - AICc
    if(second.ord == TRUE && c.hat == 1){
    
      for (obs in 1:nobserv) {
        
        ##create temporary data.frame to store fitted values and SE 
        AICctmp <- AICctab
      
        if(identical(type, "response")) {
          ##extract fitted value for observation obs
          fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                               type = parm.type1)$Predicted))
          
          ##extract SE for fitted value for observation obs
          SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                              type = parm.type1)$SE))
          if(length(unique.link) == 1) {
            AICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                              type = parm.type1, backTransform = FALSE)$Predicted))
            AICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                             type = parm.type1, backTransform = FALSE)$SE))
          } else {
            warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
            AICctmp$fit.link <- NA
            AICctmp$SE.link <- NA
          }
        }
    
      
        if(identical(type, "link")) {
          ##extract fitted value for observation obs
          fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                               type = parm.type1, backTransform = FALSE)$Predicted))
          
          ##extract SE for fitted value for observation obs
          SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                              type = parm.type1, backTransform = FALSE)$SE))
        }
      

        AICctmp$fit <- fit
        AICctmp$SE <- SE
    
      
        ##compute model averaged prediction and store in output matrix
        Mod.avg.out[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit)
        ##compute unconditional SE and store in output matrix

        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
        if(identical(uncond.se, "old")) {
          Mod.avg.out[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2))
          if(identical(type, "response")) {
            Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
            Mod.avg.out.link[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
          }
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
          Mod.avg.out[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2)))
          if(identical(type, "response")) {
            Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
            Mod.avg.out.link[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
          }        
        }
      }
    }
      

    ##begin loop - QAICc
    if(second.ord == TRUE && c.hat > 1){
      
      for (obs in 1:nobserv) {
        
        ##create temporary data.frame to store fitted values and SE 
        QAICctmp <- AICctab
  
        if(identical(type, "response")) {
          ##extract fitted value for observation obs
          fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                               type = parm.type1)$Predicted))
          ##extract SE for fitted value for observation obs
          SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                              type = parm.type1)$SE))
          if(length(unique.link) == 1) {
            QAICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                               type = parm.type1, backTransform = FALSE)$Predicted))
            QAICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                              type = parm.type1, backTransform = FALSE)$SE)) * sqrt(c.hat)
          } else {
            warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
            QAICctmp$fit.link <- NA
            QAICctmp$SE.link <- NA
          }
        }
      
        if(identical(type, "link")) {
          ##extract fitted value for observation obs
          fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                      type = parm.type1, backTransform = FALSE)$Predicted))
          ##extract SE for fitted value for observation obs
          SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                     type = parm.type1, backTransform = FALSE)$SE))
        }

        
        QAICctmp$fit <- fit
        QAICctmp$SE <- SE * sqrt(c.hat)
        
        
        ##compute model averaged prediction and store in output matrix
        Mod.avg.out[obs, 1] <- sum(QAICctmp$QAICcWt * QAICctmp$fit)
        ##compute unconditional SE and store in output matrix
        
        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
        if(identical(uncond.se, "old")) {
          Mod.avg.out[obs, 2] <- sum(QAICctmp$QAICcWt * sqrt(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2))
          if(identical(type, "response")) {
            Mod.avg.out.link[obs, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$fit.link)
            Mod.avg.out.link[obs, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.link^2 + (QAICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
          }
        }
        
        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
          Mod.avg.out[obs, 2] <- sqrt(sum(QAICctmp$QAICcWt *(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2)))
          if(identical(type, "response")) {
            Mod.avg.out.link[obs, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$fit.link)
            Mod.avg.out.link[obs, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.link^2 + (QAICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
          }
        }
      }
    }
    
  
  
    ##begin loop - AIC
    if(second.ord == FALSE && c.hat == 1) {
      
      for (obs in 1:nobserv) {
        
        ##create temporary data.frame to store fitted values and SE 
        AICtmp <- AICctab
        
        if(identical(type, "response")) {
          fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                               type = parm.type1)$Predicted))
          
          ##extract SE for fitted value for observation obs
          SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                              type = parm.type1)$SE))
          if(length(unique.link) == 1) {
            AICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                             type = parm.type1, backTransform = FALSE)$Predicted))
            AICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                            type = parm.type1, backTransform = FALSE)$SE))
            
          } else {
            warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
            AICtmp$fit.link <- NA
            AICtmp$SE.link <- NA
          }
        }
      

      if(identical(type, "link")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))
        
        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                            type = parm.type1, backTransform = FALSE)$SE))
      }
        
      
      AICtmp$fit <- fit
      AICtmp$SE <- SE

      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }  
    }
  }

  
  ##begin loop - QAIC
  if(second.ord == FALSE && c.hat > 1){
     
    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      QAICtmp <- AICctab
      
      if(identical(type, "response")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = parm.type1)$Predicted))

      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = parm.type1)$SE))
      
      if(length(unique.link) == 1) {
        QAICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                          type = parm.type1, backTransform = FALSE)$Predicted))
        QAICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                         type = parm.type1, backTransform = FALSE)$SE)) * sqrt(c.hat)
      } else {
        warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
        QAICtmp$fit.link <- NA
        QAICtmp$SE.link <- NA
      }
    }
      
      if(identical(type, "link")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = parm.type1, backTransform = FALSE)$Predicted))

      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = parm.type1, backTransform = FALSE)$SE))
    }
      
      
      QAICtmp$fit <- fit
      QAICtmp$SE <- SE * sqrt(c.hat)
      
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(QAICtmp$QAICWt * QAICtmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(QAICtmp$QAICWt * sqrt(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICtmp$QAICWt*QAICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.link^2 + (QAICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(QAICtmp$QAICWt *(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICtmp$QAICWt*QAICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.link^2 + (QAICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }
    }
  }

##############changes start here
  mod.avg.pred <- Mod.avg.out[, 1]
  uncond.se <- Mod.avg.out[, 2]
  
  ##compute confidence interval
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)

  if(identical(type, "link")) {
    lower.CL <- mod.avg.pred - zcrit * uncond.se
    upper.CL <- mod.avg.pred + zcrit * uncond.se
  } else {
    lower.CL <- link.inv(Mod.avg.out.link[, 1] - zcrit * Mod.avg.out.link[, 2])
    upper.CL <- link.inv(Mod.avg.out.link[, 1] + zcrit * Mod.avg.out.link[, 2])
  }
  
  ##create matrix
  matrix.output <- matrix(data = c(mod.avg.pred, uncond.se, lower.CL, upper.CL),
                          nrow = nrow(newdata), ncol = 4)
  colnames(matrix.output) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
  
  ##organize as list
  Mod.pred.list <- list("type" = type, "mod.avg.pred" = mod.avg.pred, "uncond.se" = uncond.se,
                        "conf.level" = conf.level, "lower.CL" = lower.CL, "upper.CL" = upper.CL,
                        "matrix.output" = matrix.output)
  class(Mod.pred.list) <- c("modavgPred", "list")
  return(Mod.pred.list)
}

    


##pcount
modavgPred.AICunmarkedFitPCount <-
  function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           type = "response", c.hat = 1, parm.type = NULL, ...) {

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
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgPred for details\n")}
    

    ##rename values according to unmarked to extract from object
    
    ##lambda
    if(identical(parm.type, "lambda")) {
      parm.type1 <- "state"
      ##check mixture type for mixture models
      mixture.type <- sapply(X = cand.set, FUN = function(i) i@mixture)
      unique.mixture <- unique(mixture.type)

    }
    
  
    ##detect
      if(identical(parm.type, "detect")) {
          parm.type1 <- "det"

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
    
      ##extract inverse link function
      if(identical(select.link, "logistic")) {
          link.inv <- plogis
      }
      if(identical(select.link, "exp")) {
          link.inv <- exp
      }

      
    ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  
    ##determine number of observations in new data set
    nobserv <- dim(newdata)[1]
  
    ##determine number of columns in new data set
    ncolumns <- dim(newdata)[2]

    ##if only 1 column, add an additional column to avoid problems in computation with predict( )
    if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
 
    ##store AICc table
    AICctab <- aictab(cand.set = cand.set, modnames = modnames,
                      second.ord = second.ord, nobs = nobs, sort = FALSE, c.hat = c.hat)

    ##create object to hold Model-averaged estimates and unconditional SE's
    Mod.avg.out <- matrix(NA, nrow = nobserv, ncol = 2)
    colnames(Mod.avg.out) <- c("Mod.avg.est", "Uncond.SE")

    
    ##required for CI
    if(identical(type, "response")) {
      Mod.avg.out.link <- matrix(NA, nrow = nobserv, ncol = 2)
      colnames(Mod.avg.out.link) <- c("Mod.avg.est", "Uncond.SE")
    }

    
    ##begin loop - AICc
    if(second.ord == TRUE && c.hat == 1){

        for (obs in 1:nobserv) {

        ##create temporary data.frame to store fitted values and SE 
        AICctmp <- AICctab
        
        if(identical(type, "response")) {
          
            fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                        type = parm.type1)$Predicted))

            ##extract SE for fitted value for observation obs
            SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                       type = parm.type1)$SE))

            
            if(length(unique.link) == 1) {
                AICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                                         type = parm.type1, backTransform = FALSE)$Predicted))
                AICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                                        type = parm.type1, backTransform = FALSE)$SE)) * sqrt(c.hat)
            } else {
                warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
                AICctmp$fit.link <- NA
                AICctmp$SE.link <- NA
            } 
        }
      
        
            if(identical(type, "link")) {
                ##extract fitted value for observation obs
                fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                            type = parm.type1, backTransform = FALSE)$Predicted))
                
                ##extract SE for fitted value for observation obs
                SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                           type = parm.type1, backTransform = FALSE)$SE))
            }
      
            AICctmp$fit <- fit
            AICctmp$SE <- SE
      
            ##compute model averaged prediction and store in output matrix
            Mod.avg.out[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit)
            ##compute unconditional SE and store in output matrix

            ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
            if(identical(uncond.se, "old")) {
                Mod.avg.out[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2))
                if(identical(type, "response")) {
                    Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
                    Mod.avg.out.link[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
                }
            }

            ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
            if(identical(uncond.se, "revised")) {
                Mod.avg.out[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2)))
                if(identical(type, "response")) {
                    Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
                    Mod.avg.out.link[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
                }        
            }
        }
    }
  

      ##begin loop - QAICc
      if(second.ord == TRUE && c.hat > 1){
          for (obs in 1:nobserv) {

              ##create temporary data.frame to store fitted values and SE
              QAICctmp <- AICctab

              if(identical(type, "response")) {
          
                  ##extract fitted value for observation obs
                  fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                              type = parm.type1)$Predicted))
            
                  ##extract SE for fitted value for observation obs
                  SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                             type = parm.type1)$SE))


                  if(length(unique.link) == 1) {
                      QAICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                                                type = parm.type1, backTransform = FALSE)$Predicted))
                      QAICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                                               type = parm.type1, backTransform = FALSE)$SE)) * sqrt(c.hat)
                  } else {
                      warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
                      QAICctmp$fit.link <- NA
                      QAICctmp$SE.link <- NA
                  }
              }
              
        
              if(identical(type, "link")) {
                  ##extract fitted value for observation obs
                  fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                              type = parm.type1, backTransform = FALSE)$Predicted))

                  ##extract SE for fitted value for observation obs
                  SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                             type = parm.type1, backTransform = FALSE)$SE))
              }
      
              ##create temporary data.frame to store fitted values and SE 
              QAICctmp$fit <- fit
              QAICctmp$SE <- SE * sqrt(c.hat)
      
              ##compute model averaged prediction and store in output matrix
              Mod.avg.out[obs, 1] <- sum(QAICctmp$QAICcWt * QAICctmp$fit)
              ##compute unconditional SE and store in output matrix

              ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
              if(identical(uncond.se, "old")) {
                  Mod.avg.out[obs, 2] <- sum(QAICctmp$QAICcWt * sqrt(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2))
                  if(identical(type, "response")) {
                      Mod.avg.out.link[obs, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$fit.link)
                      Mod.avg.out.link[obs, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.link^2 + (QAICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
                  }
              }

              ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
              if(identical(uncond.se, "revised")) {
                  Mod.avg.out[obs, 2] <- sqrt(sum(QAICctmp$QAICcWt *(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2)))
                  if(identical(type, "response")) {
                      Mod.avg.out.link[obs, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$fit.link)
                      Mod.avg.out.link[obs, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.link^2 + (QAICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
                  }
              }
          }
      }
  
  
  
      ##create temporary data.frame to store fitted values and SE - AIC
      if(second.ord == FALSE && c.hat == 1) {
          for (obs in 1:nobserv) {

              ##create temporary data.frame to store fitted values and SE 
              AICtmp <- AICctab
        
              if(identical(type, "response")) {
        
                  fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                              type = parm.type1)$Predicted))
                  
                  ##extract SE for fitted value for observation obs
                  SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                             type = parm.type1)$SE))

                  if(length(unique.link) == 1) {
                      AICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                                              type = parm.type1, backTransform = FALSE)$Predicted))
                      AICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                                             type = parm.type1, backTransform = FALSE)$SE))
                  } else {
                      warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
                      AICtmp$fit.link <- NA
                      AICtmp$SE.link <- NA
                  }
              }
   

              if(identical(type, "link")) {
                  ##extract fitted value for observation obs
                  fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                              type = parm.type1, backTransform = FALSE)$Predicted))
                  
                  ##extract SE for fitted value for observation obs
                  SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                             type = parm.type1, backTransform = FALSE)$SE))
              }
        

              AICtmp$fit <- fit
              AICtmp$SE <- SE

              ##compute model averaged prediction and store in output matrix
              Mod.avg.out[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit)
        
              ##compute unconditional SE and store in output matrix
              ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
              if(identical(uncond.se, "old")) {
                  Mod.avg.out[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2))
                  if(identical(type, "response")) {
                      Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
                      Mod.avg.out.link[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
                  }
              }

              ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
              if(identical(uncond.se, "revised")) {
                  Mod.avg.out[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2)))
                  if(identical(type, "response")) {
                      Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
                      Mod.avg.out.link[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
                  }        
              }  
          }
      }

  
      ##begin loop - QAIC
      if(second.ord == FALSE && c.hat > 1){
          for (obs in 1:nobserv) {

              ##create temporary data.frame to store fitted values and SE 
              QAICtmp <- AICctab
        
              if(identical(type, "response")) {
          
                  ##extract fitted value for observation obs
                  fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                              type = parm.type1)$Predicted))
                  
                  ##extract SE for fitted value for observation obs
                  SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                             type = parm.type1)$SE))

                  if(length(unique.link) == 1) {
                      QAICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                                               type = parm.type1, backTransform = FALSE)$Predicted))
                      QAICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                                              type = parm.type1, backTransform = FALSE)$SE)) * sqrt(c.hat)
                  } else {
                      warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
                      QAICtmp$fit.link <- NA
                      QAICtmp$SE.link <- NA
                  }
              }
                    
              if(identical(type, "link")) {
                  ##extract fitted value for observation obs
                  fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                          type = parm.type1, backTransform = FALSE)$Predicted))
          
                  ##extract SE for fitted value for observation obs
                  SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                             type = parm.type1, backTransform = FALSE)$SE))
              }
      
              ##create temporary data.frame to store fitted values and SE 
              QAICtmp$fit <- fit
              QAICtmp$SE <- SE * sqrt(c.hat)
              
              ##compute model averaged prediction and store in output matrix
              Mod.avg.out[obs, 1] <- sum(QAICtmp$QAICWt * QAICtmp$fit)
              ##compute unconditional SE and store in output matrix
              
              ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
              if(identical(uncond.se, "old")) {
                  Mod.avg.out[obs, 2] <- sum(QAICtmp$QAICWt * sqrt(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2))
                  if(identical(type, "response")) {
                      Mod.avg.out.link[obs, 1] <- sum(QAICtmp$QAICWt*QAICtmp$fit.link)
                      Mod.avg.out.link[obs, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.link^2 + (QAICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
                  }
              }

              ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
              if(identical(uncond.se, "revised")) {
                  Mod.avg.out[obs, 2] <- sqrt(sum(QAICtmp$QAICWt *(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2)))
                  if(identical(type, "response")) {
                      Mod.avg.out.link[obs, 1] <- sum(QAICtmp$QAICWt*QAICtmp$fit.link)
                      Mod.avg.out.link[obs, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.link^2 + (QAICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
                  }        
              }
          }
      }
    


##############changes start here
      mod.avg.pred <- Mod.avg.out[, 1]
      uncond.se <- Mod.avg.out[, 2]
  
      ##compute confidence interval
      zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)

      if(identical(type, "link")) {
          lower.CL <- mod.avg.pred - zcrit * uncond.se
          upper.CL <- mod.avg.pred + zcrit * uncond.se
      } else {
          lower.CL <- link.inv(Mod.avg.out.link[, 1] - zcrit * Mod.avg.out.link[, 2])
          upper.CL <- link.inv(Mod.avg.out.link[, 1] + zcrit * Mod.avg.out.link[, 2])
      }
  
      ##create matrix
      matrix.output <- matrix(data = c(mod.avg.pred, uncond.se, lower.CL, upper.CL),
                              nrow = nrow(newdata), ncol = 4)
      colnames(matrix.output) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
      
      ##organize as list
      Mod.pred.list <- list("type" = type, "mod.avg.pred" = mod.avg.pred, "uncond.se" = uncond.se,
                            "conf.level" = conf.level, "lower.CL" = lower.CL, "upper.CL" = upper.CL,
                            "matrix.output" = matrix.output)
      class(Mod.pred.list) <- c("modavgPred", "list")
      return(Mod.pred.list)
  }



##pcountOpen
modavgPred.AICunmarkedFitPCO <-
  function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           type = "response", c.hat = 1, parm.type = NULL, ...) {

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
      if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgPred for details\n")}

  
      ##rename values according to unmarked to extract from object
      ##lambda
      if(identical(parm.type, "lambda")) {
          parm.type1 <- "lambda"
          ##check mixture type for mixture models
          mixture.type <- sapply(X = cand.set, FUN = function(i) i@mixture)
          unique.mixture <- unique(mixture.type)
          if(length(unique.mixture) > 1) {
              if(any(unique.mixture == "ZIP")) stop("\nThis function does not yet support mixing ZIP with other distributions\n")
          } else {
              mixture.id <- unique(mixture.type)
              if(identical(unique.mixture, "ZIP")) {
                  if(identical(type, "link")) stop("\nLink scale not yet supported for ZIP mixtures\n")
                  if(identical(type, "response")) warning("\nModel-averaging linear predictor from a ZIP model not yet implemented\n")
              }
          }

      }



      ##gamma
      if(identical(parm.type, "gamma")) {
          parm.type1 <- "gamma"

      }
  
      ##omega
      if(identical(parm.type, "omega")) {
          parm.type1 <- "omega"

      }

      
      ##iota (for immigration = TRUE with dynamics = "autoreg", "trend", "ricker", or "gompertz")
      if(identical(parm.type, "iota")) {
          parm.type1 <- "iota"
          ##check that parameter appears in all models
          parfreq <- sum(sapply(cand.set, FUN = function(i) any(names(i@estimates@estimates) == parm.type1)))
          if(!identical(length(cand.set), parfreq)) {
              stop("\nParameter \'", parm.type1, "\' (parm.type = \"", parm.type, "\") does not appear in all models:",
                   "\ncannot compute model-averaged predictions across all models\n")
          }
      
      }
      

      ##detect
      if(identical(parm.type, "detect")) {
          parm.type1 <- "det"

         
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
 
      ##extract inverse link function
      if(identical(select.link, "logistic")) {
          link.inv <- plogis
      }
      if(identical(select.link, "exp")) {
          link.inv <- exp
      }

     
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  
  ##determine number of observations in new data set
  nobserv <- dim(newdata)[1]
  
  ##determine number of columns in new data set
  ncolumns <- dim(newdata)[2]

  ##if only 1 column, add an additional column to avoid problems in computation with predictSE.mer( )
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
 
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames,
                    second.ord = second.ord, nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = nobserv, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.est", "Uncond.SE")

  ##required for CI
  if(identical(type, "response")) {
    Mod.avg.out.link <- matrix(NA, nrow = nobserv, ncol = 2)
    colnames(Mod.avg.out.link) <- c("Mod.avg.est", "Uncond.SE")
  }
  

  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      AICctmp <- AICctab
        
      if(identical(type, "response")) {
        
        ##extract fitted value for observation obs
        if(identical(parm.type, "lambda") && identical(mixture.id, "ZIP")) {
          fit <- unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE,
                                               newdata = newdata[obs, ])$fit))
          SE <- unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE,
                                              newdata = newdata[obs, ])$se.fit))
          AICctmp$fit.link <- NA
          AICctmp$SE.link <- NA
          
        } else {
          
          fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                               type = parm.type1)$Predicted))
          
          ##extract SE for fitted value for observation obs
          SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                              type = parm.type1)$SE))
          
          
          if(length(unique.link) == 1) {
            AICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                              type = parm.type1, backTransform = FALSE)$Predicted))
            AICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                             type = parm.type1, backTransform = FALSE)$SE)) * sqrt(c.hat)
          } else {
            warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
            AICctmp$fit.link <- NA
            AICctmp$SE.link <- NA
          } 
        }
      }

      
      if(identical(type, "link")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = parm.type1, backTransform = FALSE)$Predicted))

      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = parm.type1, backTransform = FALSE)$SE))
    }
      
      ##create temporary data.frame to store fitted values and SE 
      AICctmp$fit <- fit
      AICctmp$SE <- SE
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }  
      }
    }
  }


    ##begin loop - QAICc
    if(second.ord == TRUE && c.hat > 1){
      for (obs in 1:nobserv) {

        ##create temporary data.frame to store fitted values and SE
        QAICctmp <- AICctab

        if(identical(type, "response")) {
          
          ##extract fitted value for observation obs
          if(identical(parm.type, "lambda") && identical(mixture.id, "ZIP")) {
            fit <- unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE,
                                                 newdata = newdata[obs, ])$fit))
            SE <- unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE,
                                                newdata = newdata[obs, ])$se.fit))
            QAICctmp$fit.link <- NA
            QAICctmp$SE.link <- NA

          } else {
        
            ##extract fitted value for observation obs
            fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                 type = parm.type1)$Predicted))
            
            ##extract SE for fitted value for observation obs
            SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                type = parm.type1)$SE))


            if(length(unique.link) == 1) {
              QAICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                 type = parm.type1, backTransform = FALSE)$Predicted))
              QAICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                type = parm.type1, backTransform = FALSE)$SE)) * sqrt(c.hat)
            } else {
              warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
              QAICctmp$fit.link <- NA
              QAICctmp$SE.link <- NA
            }
          }
        }       
              
        
        if(identical(type, "link")) {
          ##extract fitted value for observation obs
          fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                               type = parm.type1, backTransform = FALSE)$Predicted))

          ##extract SE for fitted value for observation obs
          SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                              type = parm.type1, backTransform = FALSE)$SE))
        }
      
        ##create temporary data.frame to store fitted values and SE 
        QAICctmp$fit <- fit
        QAICctmp$SE <- SE * sqrt(c.hat)
      
        ##compute model averaged prediction and store in output matrix
        Mod.avg.out[obs, 1] <- sum(QAICctmp$QAICcWt * QAICctmp$fit)
        ##compute unconditional SE and store in output matrix

        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
        if(identical(uncond.se, "old")) {
          Mod.avg.out[obs, 2] <- sum(QAICctmp$QAICcWt * sqrt(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2))
          if(identical(type, "response")) {
            Mod.avg.out.link[obs, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$fit.link)
            Mod.avg.out.link[obs, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.link^2 + (QAICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
          }
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
          Mod.avg.out[obs, 2] <- sqrt(sum(QAICctmp$QAICcWt *(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2)))
          if(identical(type, "response")) {
            Mod.avg.out.link[obs, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$fit.link)
            Mod.avg.out.link[obs, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.link^2 + (QAICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
          }
        }
      }
    }
  
  
  
    ##create temporary data.frame to store fitted values and SE - AIC
    if(second.ord == FALSE && c.hat == 1) {
      for (obs in 1:nobserv) {

        ##create temporary data.frame to store fitted values and SE 
        AICtmp <- AICctab
        
        if(identical(type, "response")) {
        
          if(identical(parm.type, "lambda") && identical(mixture.id, "ZIP")) {
            fit <- unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE,
                                                 newdata = newdata[obs, ])$fit))
            SE <- unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE,
                                                newdata = newdata[obs, ])$se.fit))
            AICtmp$fit.link <- NA
            AICtmp$SE.link <- NA
            
          } else {
            fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                 type = parm.type1)$Predicted))
            
            ##extract SE for fitted value for observation obs
            SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                type = parm.type1)$SE))

            if(length(unique.link) == 1) {
              AICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                               type = parm.type1, backTransform = FALSE)$Predicted))
              AICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                              type = parm.type1, backTransform = FALSE)$SE))
            } else {
              warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
              AICtmp$fit.link <- NA
              AICtmp$SE.link <- NA
            }
          }
        }


        if(identical(type, "link")) {
          ##extract fitted value for observation obs
          fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                               type = parm.type1, backTransform = FALSE)$Predicted))

          ##extract SE for fitted value for observation obs
          SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                              type = parm.type1, backTransform = FALSE)$SE))
        }
        

        AICtmp$fit <- fit
        AICtmp$SE <- SE

        ##compute model averaged prediction and store in output matrix
        Mod.avg.out[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit)
        
        ##compute unconditional SE and store in output matrix
        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
        if(identical(uncond.se, "old")) {
          Mod.avg.out[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2))
          if(identical(type, "response")) {
            Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
            Mod.avg.out.link[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
          }
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
          Mod.avg.out[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2)))
          if(identical(type, "response")) {
            Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
            Mod.avg.out.link[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
          }        
        }  
      }
    }

  
    ##begin loop - QAIC
    if(second.ord == FALSE && c.hat > 1){
      for (obs in 1:nobserv) {

        ##create temporary data.frame to store fitted values and SE 
        QAICtmp <- AICctab
        
        if(identical(type, "response")) {
          
          ##extract fitted value for observation obs
          if(identical(parm.type, "lambda") && identical(mixture.id, "ZIP")) {
            fit <- unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE,
                                                 newdata = newdata[obs, ])$fit))
            SE <- unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE,
                                                newdata = newdata[obs, ])$se.fit))
            QAICtmp$fit.link <- NA
            QAICtmp$SE.link <- NA
            
          } else {
            
            ##extract fitted value for observation obs
            fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                 type = parm.type1)$Predicted))

            ##extract SE for fitted value for observation obs
            SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                type = parm.type1)$SE))

            if(length(unique.link) == 1) {
              QAICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                type = parm.type1, backTransform = FALSE)$Predicted))
              QAICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                               type = parm.type1, backTransform = FALSE)$SE)) * sqrt(c.hat)
            } else {
              warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
              QAICtmp$fit.link <- NA
              QAICtmp$SE.link <- NA
            }
          }
        }
          
        if(identical(type, "link")) {
          ##extract fitted value for observation obs
          fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                               type = parm.type1, backTransform = FALSE)$Predicted))
          
          ##extract SE for fitted value for observation obs
          SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                              type = parm.type1, backTransform = FALSE)$SE))
        }
      
        ##create temporary data.frame to store fitted values and SE 
        QAICtmp$fit <- fit
        QAICtmp$SE <- SE * sqrt(c.hat)
        
        ##compute model averaged prediction and store in output matrix
        Mod.avg.out[obs, 1] <- sum(QAICtmp$QAICWt * QAICtmp$fit)
        ##compute unconditional SE and store in output matrix

        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
        if(identical(uncond.se, "old")) {
          Mod.avg.out[obs, 2] <- sum(QAICtmp$QAICWt * sqrt(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2))
          if(identical(type, "response")) {
            Mod.avg.out.link[obs, 1] <- sum(QAICtmp$QAICWt*QAICtmp$fit.link)
            Mod.avg.out.link[obs, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.link^2 + (QAICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
          }
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
          Mod.avg.out[obs, 2] <- sqrt(sum(QAICtmp$QAICWt *(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2)))
          if(identical(type, "response")) {
            Mod.avg.out.link[obs, 1] <- sum(QAICtmp$QAICWt*QAICtmp$fit.link)
            Mod.avg.out.link[obs, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.link^2 + (QAICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
          }        
        }
      }
    }


##############changes start here
    mod.avg.pred <- Mod.avg.out[, 1]
    uncond.se <- Mod.avg.out[, 2]
  
    ##compute confidence interval
    zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)

    if(identical(type, "link")) {
      lower.CL <- mod.avg.pred - zcrit * uncond.se
      upper.CL <- mod.avg.pred + zcrit * uncond.se
    } else {
      lower.CL <- link.inv(Mod.avg.out.link[, 1] - zcrit * Mod.avg.out.link[, 2])
      upper.CL <- link.inv(Mod.avg.out.link[, 1] + zcrit * Mod.avg.out.link[, 2])
    }
  
    ##create matrix
    matrix.output <- matrix(data = c(mod.avg.pred, uncond.se, lower.CL, upper.CL),
                            nrow = nrow(newdata), ncol = 4)
    colnames(matrix.output) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
  
    ##organize as list
    Mod.pred.list <- list("type" = type, "mod.avg.pred" = mod.avg.pred, "uncond.se" = uncond.se,
                          "conf.level" = conf.level, "lower.CL" = lower.CL, "upper.CL" = upper.CL,
                          "matrix.output" = matrix.output)
    class(Mod.pred.list) <- c("modavgPred", "list")
    return(Mod.pred.list)
  }



##distsamp
modavgPred.AICunmarkedFitDS <-
  function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           type = "response", c.hat = 1, parm.type = NULL, ...) {

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
  if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgPred for details\n")}

  
  ##rename values according to unmarked to extract from object
  ##lambda
  if(identical(parm.type, "lambda")) {
      parm.type1 <- "state"

  }
  
      ##detect
      if(identical(parm.type, "detect")){
          parm.type1 <- "det"
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

      ##extract inverse link function
      if(identical(select.link, "logistic")) {
          link.inv <- plogis
      }
      if(identical(select.link, "exp")) {
          link.inv <- exp
      }
  
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  
  ##determine number of observations in new data set
  nobserv <- dim(newdata)[1]
  
  ##determine number of columns in new data set
  ncolumns <- dim(newdata)[2]

  ##if only 1 column, add an additional column to avoid problems in computation with predictSE.mer( )
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
 
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames,
                    second.ord = second.ord, nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = nobserv, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.est", "Uncond.SE")

  ##required for CI
  if(identical(type, "response")) {
    Mod.avg.out.link <- matrix(NA, nrow = nobserv, ncol = 2)
    colnames(Mod.avg.out.link) <- c("Mod.avg.est", "Uncond.SE")
  }
  

  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
    
    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      AICctmp <- AICctab
      
      if(identical(type, "response")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1)$Predicted))
        
        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                            type = parm.type1)$SE))
        if(length(unique.link) == 1) {
          AICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))
          AICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                           type = parm.type1, backTransform = FALSE)$SE))
        } else {
          warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
          AICctmp$fit.link <- NA
          AICctmp$SE.link <- NA
        }
      }
    
      
      if(identical(type, "link")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))

        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          type = parm.type1, backTransform = FALSE)$SE))
      }
      

      AICctmp$fit <- fit
      AICctmp$SE <- SE
    
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }        
      }
    }
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
    
    for (obs in 1:nobserv) {
      
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
  
      if(identical(type, "response")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                           type = parm.type1)$Predicted))
      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          type = parm.type1)$SE))
      if(length(unique.link) == 1) {
        QAICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                           type = parm.type1, backTransform = FALSE)$Predicted))
        QAICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                          type = parm.type1, backTransform = FALSE)$SE)) * sqrt(c.hat)
      } else {
        warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
        QAICctmp$fit.link <- NA
        QAICctmp$SE.link <- NA
      }
    }
      
      if(identical(type, "link")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                           type = parm.type1, backTransform = FALSE)$Predicted))
      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          type = parm.type1, backTransform = FALSE)$SE))
    }

      
      QAICctmp$fit <- fit
      QAICctmp$SE <- SE * sqrt(c.hat)
      
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(QAICctmp$QAICcWt * QAICctmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(QAICctmp$QAICcWt * sqrt(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.link^2 + (QAICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(QAICctmp$QAICcWt *(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.link^2 + (QAICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }
    }
  }
  
  
  ##begin loop - AIC
  if(second.ord == FALSE && c.hat == 1) {

    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      AICtmp <- AICctab
      
      if(identical(type, "response")) {
          fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                               type = parm.type1)$Predicted))

          ##extract SE for fitted value for observation obs
          SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                              type = parm.type1)$SE))
          if(length(unique.link) == 1) {
            AICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                             type = parm.type1, backTransform = FALSE)$Predicted))
            AICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                            type = parm.type1, backTransform = FALSE)$SE))
            
          } else {
            warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
            AICtmp$fit.link <- NA
            AICtmp$SE.link <- NA
          }
        }
      

      if(identical(type, "link")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))
        
        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                            type = parm.type1, backTransform = FALSE)$SE))
      }
        
      
      AICtmp$fit <- fit
      AICtmp$SE <- SE

      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }  
    }
  }

  
  ##begin loop - QAIC
  if(second.ord == FALSE && c.hat > 1){
     
    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      QAICtmp <- AICctab
      
      if(identical(type, "response")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = parm.type1)$Predicted))

      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = parm.type1)$SE))
      
      if(length(unique.link) == 1) {
        QAICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                          type = parm.type1, backTransform = FALSE)$Predicted))
        QAICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                         type = parm.type1, backTransform = FALSE)$SE)) * sqrt(c.hat)
      } else {
        warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
        QAICtmp$fit.link <- NA
        QAICtmp$SE.link <- NA
      }
    }
      
      if(identical(type, "link")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = parm.type1, backTransform = FALSE)$Predicted))

      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = parm.type1, backTransform = FALSE)$SE))
    }
      
      
      QAICtmp$fit <- fit
      QAICtmp$SE <- SE * sqrt(c.hat)
      
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(QAICtmp$QAICWt * QAICtmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(QAICtmp$QAICWt * sqrt(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICtmp$QAICWt*QAICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.link^2 + (QAICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(QAICtmp$QAICWt *(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICtmp$QAICWt*QAICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.link^2 + (QAICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }
    }
  }

##############changes start here
  mod.avg.pred <- Mod.avg.out[, 1]
  uncond.se <- Mod.avg.out[, 2]
  
  ##compute confidence interval
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)

  if(identical(type, "link")) {
    lower.CL <- mod.avg.pred - zcrit * uncond.se
    upper.CL <- mod.avg.pred + zcrit * uncond.se
  } else {
    lower.CL <- link.inv(Mod.avg.out.link[, 1] - zcrit * Mod.avg.out.link[, 2])
    upper.CL <- link.inv(Mod.avg.out.link[, 1] + zcrit * Mod.avg.out.link[, 2])
  }
  
  ##create matrix
  matrix.output <- matrix(data = c(mod.avg.pred, uncond.se, lower.CL, upper.CL),
                          nrow = nrow(newdata), ncol = 4)
  colnames(matrix.output) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
  
  ##organize as list
  Mod.pred.list <- list("type" = type, "mod.avg.pred" = mod.avg.pred, "uncond.se" = uncond.se,
                        "conf.level" = conf.level, "lower.CL" = lower.CL, "upper.CL" = upper.CL,
                        "matrix.output" = matrix.output)
  class(Mod.pred.list) <- c("modavgPred", "list")
  return(Mod.pred.list)
}



##gdistsamp
modavgPred.AICunmarkedFitGDS <-
  function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           type = "response", c.hat = 1, parm.type = NULL, ...) {

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
  if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgPred for details\n")}

  
  ##rename values according to unmarked to extract from object

  ##lambda
  if(identical(parm.type, "lambda")) {
      parm.type1 <- "lambda"

  }
  
  ##detect
  if(identical(parm.type, "detect")) {
      parm.type1 <- "det"
      ##check for key function used
      keyid <- unique(sapply(cand.set, FUN = function(i) i@keyfun))
      if(any(keyid == "uniform")) stop("\nDetection parameter not found in some models\n")

  }

      
      ##availability
      if(identical(parm.type, "phi")) {
          parm.type1 <- "phi"
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

      ##extract inverse link function
      if(identical(select.link, "logistic")) {
          link.inv <- plogis
      }
      if(identical(select.link, "exp")) {
          link.inv <- exp
      }
  
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  
  ##determine number of observations in new data set
  nobserv <- dim(newdata)[1]
  
  ##determine number of columns in new data set
  ncolumns <- dim(newdata)[2]

  ##if only 1 column, add an additional column to avoid problems in computation with predictSE.mer( )
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
 
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames,
                    second.ord = second.ord, nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = nobserv, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.est", "Uncond.SE")
  

  ##required for CI
  if(identical(type, "response")) {
    Mod.avg.out.link <- matrix(NA, nrow = nobserv, ncol = 2)
    colnames(Mod.avg.out.link) <- c("Mod.avg.est", "Uncond.SE")
  }


  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
    
    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      AICctmp <- AICctab
      
      if(identical(type, "response")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1)$Predicted))
        
        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                            type = parm.type1)$SE))
        if(length(unique.link) == 1) {
          AICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))
          AICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                           type = parm.type1, backTransform = FALSE)$SE))
        } else {
          warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
          AICctmp$fit.link <- NA
          AICctmp$SE.link <- NA
        }
      }
    
      
      if(identical(type, "link")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))

        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          type = parm.type1, backTransform = FALSE)$SE))
      }
      

      AICctmp$fit <- fit
      AICctmp$SE <- SE
    
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }        
      }
    }
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
    
    for (obs in 1:nobserv) {
      
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
  
      if(identical(type, "response")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                           type = parm.type1)$Predicted))
      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          type = parm.type1)$SE))
      if(length(unique.link) == 1) {
        QAICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                           type = parm.type1, backTransform = FALSE)$Predicted))
        QAICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                          type = parm.type1, backTransform = FALSE)$SE)) * sqrt(c.hat)
      } else {
        warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
        QAICctmp$fit.link <- NA
        QAICctmp$SE.link <- NA
      }
    }
      
      if(identical(type, "link")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                           type = parm.type1, backTransform = FALSE)$Predicted))
      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          type = parm.type1, backTransform = FALSE)$SE))
    }

      
      QAICctmp$fit <- fit
      QAICctmp$SE <- SE * sqrt(c.hat)
      
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(QAICctmp$QAICcWt * QAICctmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(QAICctmp$QAICcWt * sqrt(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.link^2 + (QAICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(QAICctmp$QAICcWt *(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.link^2 + (QAICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }
    }
  }
  
  
  
  ##begin loop - AIC
  if(second.ord == FALSE && c.hat == 1) {

    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      AICtmp <- AICctab
      
      if(identical(type, "response")) {
          fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                               type = parm.type1)$Predicted))

          ##extract SE for fitted value for observation obs
          SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                              type = parm.type1)$SE))
          if(length(unique.link) == 1) {
            AICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                             type = parm.type1, backTransform = FALSE)$Predicted))
            AICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                            type = parm.type1, backTransform = FALSE)$SE))
            
          } else {
            warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
            AICtmp$fit.link <- NA
            AICtmp$SE.link <- NA
          }
        }
      

      if(identical(type, "link")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))
        
        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                            type = parm.type1, backTransform = FALSE)$SE))
      }
        
      
      AICtmp$fit <- fit
      AICtmp$SE <- SE

      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }  
    }
  }

  
  ##begin loop - QAIC
  if(second.ord == FALSE && c.hat > 1){
     
    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      QAICtmp <- AICctab
      
      if(identical(type, "response")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = parm.type1)$Predicted))

      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = parm.type1)$SE))
      
      if(length(unique.link) == 1) {
        QAICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                          type = parm.type1, backTransform = FALSE)$Predicted))
        QAICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                         type = parm.type1, backTransform = FALSE)$SE)) * sqrt(c.hat)
      } else {
        warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
        QAICtmp$fit.link <- NA
        QAICtmp$SE.link <- NA
      }
    }
      
      if(identical(type, "link")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = parm.type1, backTransform = FALSE)$Predicted))

      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = parm.type1, backTransform = FALSE)$SE))
    }
      
      
      QAICtmp$fit <- fit
      QAICtmp$SE <- SE * sqrt(c.hat)
      
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(QAICtmp$QAICWt * QAICtmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(QAICtmp$QAICWt * sqrt(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICtmp$QAICWt*QAICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.link^2 + (QAICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(QAICtmp$QAICWt *(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICtmp$QAICWt*QAICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.link^2 + (QAICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }
    }
  }

##############changes start here
  mod.avg.pred <- Mod.avg.out[, 1]
  uncond.se <- Mod.avg.out[, 2]
  
  ##compute confidence interval
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)

  if(identical(type, "link")) {
    lower.CL <- mod.avg.pred - zcrit * uncond.se
    upper.CL <- mod.avg.pred + zcrit * uncond.se
  } else {
    lower.CL <- link.inv(Mod.avg.out.link[, 1] - zcrit * Mod.avg.out.link[, 2])
    upper.CL <- link.inv(Mod.avg.out.link[, 1] + zcrit * Mod.avg.out.link[, 2])
  }
  
  ##create matrix
  matrix.output <- matrix(data = c(mod.avg.pred, uncond.se, lower.CL, upper.CL),
                          nrow = nrow(newdata), ncol = 4)
  colnames(matrix.output) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
  
  ##organize as list
  Mod.pred.list <- list("type" = type, "mod.avg.pred" = mod.avg.pred, "uncond.se" = uncond.se,
                        "conf.level" = conf.level, "lower.CL" = lower.CL, "upper.CL" = upper.CL,
                        "matrix.output" = matrix.output)
  class(Mod.pred.list) <- c("modavgPred", "list")
  return(Mod.pred.list)
}



##occuFP
modavgPred.AICunmarkedFitOccuFP <-
  function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           type = "response", c.hat = 1, parm.type = NULL, ...) {

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
  if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgPred for details\n")}

  
  ##rename values according to unmarked to extract from object
  ##psi
  if(identical(parm.type, "psi")) {
      parm.type1 <- "state"

  }

  
  ##detect
      if(identical(parm.type, "detect")) {
          parm.type1 <- "det"

      }
  
  ##false positives
      if(identical(parm.type, "falsepos") || identical(parm.type, "fp")) {
          parm.type1 <- "fp"

      }

  ##certain detections
  if(identical(parm.type, "certain")) {
      parm.type1 <- "b"
      ##check that parameter appears in all models
      parfreq <- sum(sapply(cand.set, FUN = function(i) any(names(i@estimates@estimates) == parm.type1)))
      if(!identical(length(cand.set), parfreq)) {
          stop("\nParameter \'", parm.type1, "\' (parm.type = \"", parm.type, "\") does not appear in all models:",
               "\ncannot compute model-averaged predictions across all models\n")
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
      
      ##extract inverse link function
      if(identical(select.link, "logistic")) {
          link.inv <- plogis
      }
  
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  
  ##determine number of observations in new data set
  nobserv <- dim(newdata)[1]
  
  ##determine number of columns in new data set
  ncolumns <- dim(newdata)[2]

  ##if only 1 column, add an additional column to avoid problems in computation with predictSE.mer( )
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
 
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames,
                    second.ord = second.ord, nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = nobserv, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.est", "Uncond.SE")

  
  ##required for CI
  if(identical(type, "response")) {
    Mod.avg.out.link <- matrix(NA, nrow = nobserv, ncol = 2)
    colnames(Mod.avg.out.link) <- c("Mod.avg.est", "Uncond.SE")
  }
  
  
  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
    
    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      AICctmp <- AICctab
      
      if(identical(type, "response")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1)$Predicted))
        
        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                            type = parm.type1)$SE))
        if(length(unique.link) == 1) {
          AICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))
          AICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                           type = parm.type1, backTransform = FALSE)$SE))
        } else {
          warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
          AICctmp$fit.link <- NA
          AICctmp$SE.link <- NA
        }
      }
    
      
      if(identical(type, "link")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))

        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          type = parm.type1, backTransform = FALSE)$SE))
      }
      

      AICctmp$fit <- fit
      AICctmp$SE <- SE
    
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }        
      }
    }
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
    
    for (obs in 1:nobserv) {
      
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
  
      if(identical(type, "response")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                           type = parm.type1)$Predicted))
      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          type = parm.type1)$SE))
      if(length(unique.link) == 1) {
        QAICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                           type = parm.type1, backTransform = FALSE)$Predicted))
        QAICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                          type = parm.type1, backTransform = FALSE)$SE)) * sqrt(c.hat)
      } else {
        warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
        QAICctmp$fit.link <- NA
        QAICctmp$SE.link <- NA
      }
    }
      
      if(identical(type, "link")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                           type = parm.type1, backTransform = FALSE)$Predicted))
      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          type = parm.type1, backTransform = FALSE)$SE))
    }

      
      QAICctmp$fit <- fit
      QAICctmp$SE <- SE * sqrt(c.hat)
      
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(QAICctmp$QAICcWt * QAICctmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(QAICctmp$QAICcWt * sqrt(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.link^2 + (QAICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(QAICctmp$QAICcWt *(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.link^2 + (QAICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }
    }
  }
  
  
  
  ##begin loop - AIC
  if(second.ord == FALSE && c.hat == 1) {

    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      AICtmp <- AICctab
      
      if(identical(type, "response")) {
          fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                               type = parm.type1)$Predicted))

          ##extract SE for fitted value for observation obs
          SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                              type = parm.type1)$SE))
          if(length(unique.link) == 1) {
            AICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                             type = parm.type1, backTransform = FALSE)$Predicted))
            AICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                            type = parm.type1, backTransform = FALSE)$SE))
            
          } else {
            warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
            AICtmp$fit.link <- NA
            AICtmp$SE.link <- NA
          }
        }
      

      if(identical(type, "link")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))
        
        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                            type = parm.type1, backTransform = FALSE)$SE))
      }
        
      
      AICtmp$fit <- fit
      AICtmp$SE <- SE

      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }  
    }
  }

  
  ##begin loop - QAIC
  if(second.ord == FALSE && c.hat > 1){
     
    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      QAICtmp <- AICctab
      
      if(identical(type, "response")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = parm.type1)$Predicted))

      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = parm.type1)$SE))
      
      if(length(unique.link) == 1) {
        QAICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                          type = parm.type1, backTransform = FALSE)$Predicted))
        QAICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                         type = parm.type1, backTransform = FALSE)$SE)) * sqrt(c.hat)
      } else {
        warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
        QAICtmp$fit.link <- NA
        QAICtmp$SE.link <- NA
      }
    }
      
      if(identical(type, "link")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = parm.type1, backTransform = FALSE)$Predicted))

      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = parm.type1, backTransform = FALSE)$SE))
    }
      
      
      QAICtmp$fit <- fit
      QAICtmp$SE <- SE * sqrt(c.hat)
      
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(QAICtmp$QAICWt * QAICtmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(QAICtmp$QAICWt * sqrt(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICtmp$QAICWt*QAICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.link^2 + (QAICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(QAICtmp$QAICWt *(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICtmp$QAICWt*QAICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.link^2 + (QAICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }
    }
  }

##############changes start here
  mod.avg.pred <- Mod.avg.out[, 1]
  uncond.se <- Mod.avg.out[, 2]
  
  ##compute confidence interval
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)

  if(identical(type, "link")) {
    lower.CL <- mod.avg.pred - zcrit * uncond.se
    upper.CL <- mod.avg.pred + zcrit * uncond.se
  } else {
    lower.CL <- link.inv(Mod.avg.out.link[, 1] - zcrit * Mod.avg.out.link[, 2])
    upper.CL <- link.inv(Mod.avg.out.link[, 1] + zcrit * Mod.avg.out.link[, 2])
  }
  
  ##create matrix
  matrix.output <- matrix(data = c(mod.avg.pred, uncond.se, lower.CL, upper.CL),
                          nrow = nrow(newdata), ncol = 4)
  colnames(matrix.output) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
  
  ##organize as list
  Mod.pred.list <- list("type" = type, "mod.avg.pred" = mod.avg.pred, "uncond.se" = uncond.se,
                        "conf.level" = conf.level, "lower.CL" = lower.CL, "upper.CL" = upper.CL,
                        "matrix.output" = matrix.output)
  class(Mod.pred.list) <- c("modavgPred", "list")
  return(Mod.pred.list)
}



##multinomPois
modavgPred.AICunmarkedFitMPois <-
  function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           type = "response", c.hat = 1, parm.type = NULL, ...) {

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
  if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgPred for details\n")}
  

  ##rename values according to unmarked to extract from object
  ##lambda
  if(identical(parm.type, "lambda")) {
    parm.type1 <- "state"

  }
  
  
      ##detect
      if(identical(parm.type, "detect")) {
          parm.type1 <- "det"

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

  
  ##extract inverse link function
  if(identical(select.link, "logistic")) {
    link.inv <- plogis
  }
  if(identical(select.link, "exp")) {
    link.inv <- exp
  }
####changes#############

  
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  
  ##determine number of observations in new data set
  nobserv <- dim(newdata)[1]
  
  ##determine number of columns in new data set
  ncolumns <- dim(newdata)[2]

  ##if only 1 column, add an additional column to avoid problems in computation with predictSE.mer( )
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
 
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames,
                    second.ord = second.ord, nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = nobserv, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.est", "Uncond.SE")
  
  ##required for CI
  if(identical(type, "response")) {
    Mod.avg.out.link <- matrix(NA, nrow = nobserv, ncol = 2)
    colnames(Mod.avg.out.link) <- c("Mod.avg.est", "Uncond.SE")
  }

  
  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
    
    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      AICctmp <- AICctab
      
      if(identical(type, "response")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1)$Predicted))
        
        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                            type = parm.type1)$SE))
        if(length(unique.link) == 1) {
          AICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))
          AICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                           type = parm.type1, backTransform = FALSE)$SE))
        } else {
          warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
          AICctmp$fit.link <- NA
          AICctmp$SE.link <- NA
        }
      }
    
      
      if(identical(type, "link")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))

        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          type = parm.type1, backTransform = FALSE)$SE))
      }
      

      AICctmp$fit <- fit
      AICctmp$SE <- SE
    
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }        
      }
    }
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
    
    for (obs in 1:nobserv) {
      
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
  
      if(identical(type, "response")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                           type = parm.type1)$Predicted))
      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          type = parm.type1)$SE))
      if(length(unique.link) == 1) {
        QAICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                           type = parm.type1, backTransform = FALSE)$Predicted))
        QAICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                          type = parm.type1, backTransform = FALSE)$SE)) * sqrt(c.hat)
      } else {
        warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
        QAICctmp$fit.link <- NA
        QAICctmp$SE.link <- NA
      }
    }
      
      if(identical(type, "link")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                           type = parm.type1, backTransform = FALSE)$Predicted))
      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          type = parm.type1, backTransform = FALSE)$SE))
    }

      
      QAICctmp$fit <- fit
      QAICctmp$SE <- SE * sqrt(c.hat)
      
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(QAICctmp$QAICcWt * QAICctmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(QAICctmp$QAICcWt * sqrt(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.link^2 + (QAICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(QAICctmp$QAICcWt *(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.link^2 + (QAICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }
    }
  }
  
  
  ##begin loop - AIC
  if(second.ord == FALSE && c.hat == 1) {

    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      AICtmp <- AICctab
      
      if(identical(type, "response")) {
          fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                               type = parm.type1)$Predicted))

          ##extract SE for fitted value for observation obs
          SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                              type = parm.type1)$SE))
          if(length(unique.link) == 1) {
            AICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                             type = parm.type1, backTransform = FALSE)$Predicted))
            AICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                            type = parm.type1, backTransform = FALSE)$SE))
            
          } else {
            warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
            AICtmp$fit.link <- NA
            AICtmp$SE.link <- NA
          }
        }
      

      if(identical(type, "link")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))
        
        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                            type = parm.type1, backTransform = FALSE)$SE))
      }
        
      
      AICtmp$fit <- fit
      AICtmp$SE <- SE

      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }  
    }
  }

  
  ##begin loop - QAIC
  if(second.ord == FALSE && c.hat > 1){
     
    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      QAICtmp <- AICctab
      
      if(identical(type, "response")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = parm.type1)$Predicted))

      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = parm.type1)$SE))
      
      if(length(unique.link) == 1) {
        QAICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                          type = parm.type1, backTransform = FALSE)$Predicted))
        QAICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                         type = parm.type1, backTransform = FALSE)$SE)) * sqrt(c.hat)
      } else {
        warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
        QAICtmp$fit.link <- NA
        QAICtmp$SE.link <- NA
      }
    }
      
      if(identical(type, "link")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = parm.type1, backTransform = FALSE)$Predicted))

      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = parm.type1, backTransform = FALSE)$SE))
    }
      
      
      QAICtmp$fit <- fit
      QAICtmp$SE <- SE * sqrt(c.hat)
      
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(QAICtmp$QAICWt * QAICtmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(QAICtmp$QAICWt * sqrt(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICtmp$QAICWt*QAICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.link^2 + (QAICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(QAICtmp$QAICWt *(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICtmp$QAICWt*QAICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.link^2 + (QAICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }
    }
  }

##############changes start here
  mod.avg.pred <- Mod.avg.out[, 1]
  uncond.se <- Mod.avg.out[, 2]
  
  ##compute confidence interval
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)

  if(identical(type, "link")) {
    lower.CL <- mod.avg.pred - zcrit * uncond.se
    upper.CL <- mod.avg.pred + zcrit * uncond.se
  } else {
    lower.CL <- link.inv(Mod.avg.out.link[, 1] - zcrit * Mod.avg.out.link[, 2])
    upper.CL <- link.inv(Mod.avg.out.link[, 1] + zcrit * Mod.avg.out.link[, 2])
  }
  
  ##create matrix
  matrix.output <- matrix(data = c(mod.avg.pred, uncond.se, lower.CL, upper.CL),
                          nrow = nrow(newdata), ncol = 4)
  colnames(matrix.output) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
  
  ##organize as list
  Mod.pred.list <- list("type" = type, "mod.avg.pred" = mod.avg.pred, "uncond.se" = uncond.se,
                        "conf.level" = conf.level, "lower.CL" = lower.CL, "upper.CL" = upper.CL,
                        "matrix.output" = matrix.output)
  class(Mod.pred.list) <- c("modavgPred", "list")
  return(Mod.pred.list)
}



##gmultmix
modavgPred.AICunmarkedFitGMM <-
  function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           type = "response", c.hat = 1, parm.type = NULL, ...) {

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
  if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgPred for details\n")}
  


  ##rename values according to unmarked to extract from object
  ##lambda
  if(identical(parm.type, "lambda")) {
    parm.type1 <- "lambda"
  }
  
      ##detect
      if(identical(parm.type, "detect")) {
          parm.type1 <- "det"
      }

      ##availability
      if(identical(parm.type, "phi")) {
          parm.type1 <- "phi"
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


  ##extract inverse link function
  if(identical(select.link, "logistic")) {
    link.inv <- plogis
  }
  if(identical(select.link, "exp")) {
    link.inv <- exp
  }
####changes#############

  
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  
  ##determine number of observations in new data set
  nobserv <- dim(newdata)[1]
  
  ##determine number of columns in new data set
  ncolumns <- dim(newdata)[2]

  ##if only 1 column, add an additional column to avoid problems in computation with predictSE.mer( )
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
 
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames,
                    second.ord = second.ord, nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = nobserv, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.est", "Uncond.SE")
  
  ##required for CI
  if(identical(type, "response")) {
    Mod.avg.out.link <- matrix(NA, nrow = nobserv, ncol = 2)
    colnames(Mod.avg.out.link) <- c("Mod.avg.est", "Uncond.SE")
  }


  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
    
    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      AICctmp <- AICctab
      
      if(identical(type, "response")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1)$Predicted))
        
        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                            type = parm.type1)$SE))
        if(length(unique.link) == 1) {
          AICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))
          AICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                           type = parm.type1, backTransform = FALSE)$SE))
        } else {
          warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
          AICctmp$fit.link <- NA
          AICctmp$SE.link <- NA
        }
      }
    
      
      if(identical(type, "link")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))

        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          type = parm.type1, backTransform = FALSE)$SE))
      }
      

      AICctmp$fit <- fit
      AICctmp$SE <- SE
    
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }        
      }
    }
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
    
    for (obs in 1:nobserv) {
      
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
  
      if(identical(type, "response")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                           type = parm.type1)$Predicted))
      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          type = parm.type1)$SE))
      if(length(unique.link) == 1) {
        QAICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                           type = parm.type1, backTransform = FALSE)$Predicted))
        QAICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                          type = parm.type1, backTransform = FALSE)$SE)) * sqrt(c.hat)
      } else {
        warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
        QAICctmp$fit.link <- NA
        QAICctmp$SE.link <- NA
      }
    }
      
      if(identical(type, "link")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                           type = parm.type1, backTransform = FALSE)$Predicted))
      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          type = parm.type1, backTransform = FALSE)$SE))
    }

      
      QAICctmp$fit <- fit
      QAICctmp$SE <- SE * sqrt(c.hat)
      
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(QAICctmp$QAICcWt * QAICctmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(QAICctmp$QAICcWt * sqrt(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.link^2 + (QAICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(QAICctmp$QAICcWt *(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.link^2 + (QAICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }
    }
  }
  
  
  ##begin loop - AIC
  if(second.ord == FALSE && c.hat == 1) {

    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      AICtmp <- AICctab
      
      if(identical(type, "response")) {
          fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                               type = parm.type1)$Predicted))

          ##extract SE for fitted value for observation obs
          SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                              type = parm.type1)$SE))
          if(length(unique.link) == 1) {
            AICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                             type = parm.type1, backTransform = FALSE)$Predicted))
            AICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                            type = parm.type1, backTransform = FALSE)$SE))
            
          } else {
            warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
            AICtmp$fit.link <- NA
            AICtmp$SE.link <- NA
          }
        }
      

      if(identical(type, "link")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))
        
        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                            type = parm.type1, backTransform = FALSE)$SE))
      }
        
      
      AICtmp$fit <- fit
      AICtmp$SE <- SE

      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }  
    }
  }

  
  ##begin loop - QAIC
  if(second.ord == FALSE && c.hat > 1){
     
    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      QAICtmp <- AICctab
      
      if(identical(type, "response")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = parm.type1)$Predicted))

      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = parm.type1)$SE))
      
      if(length(unique.link) == 1) {
        QAICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                          type = parm.type1, backTransform = FALSE)$Predicted))
        QAICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                         type = parm.type1, backTransform = FALSE)$SE)) * sqrt(c.hat)
      } else {
        warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
        QAICtmp$fit.link <- NA
        QAICtmp$SE.link <- NA
      }
    }
      
      if(identical(type, "link")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = parm.type1, backTransform = FALSE)$Predicted))

      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = parm.type1, backTransform = FALSE)$SE))
    }
      
      
      QAICtmp$fit <- fit
      QAICtmp$SE <- SE * sqrt(c.hat)
      
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(QAICtmp$QAICWt * QAICtmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(QAICtmp$QAICWt * sqrt(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICtmp$QAICWt*QAICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.link^2 + (QAICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(QAICtmp$QAICWt *(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICtmp$QAICWt*QAICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.link^2 + (QAICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }
    }
  }

##############changes start here
  mod.avg.pred <- Mod.avg.out[, 1]
  uncond.se <- Mod.avg.out[, 2]
  
  ##compute confidence interval
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)

  if(identical(type, "link")) {
    lower.CL <- mod.avg.pred - zcrit * uncond.se
    upper.CL <- mod.avg.pred + zcrit * uncond.se
  } else {
    lower.CL <- link.inv(Mod.avg.out.link[, 1] - zcrit * Mod.avg.out.link[, 2])
    upper.CL <- link.inv(Mod.avg.out.link[, 1] + zcrit * Mod.avg.out.link[, 2])
  }
  
  ##create matrix
  matrix.output <- matrix(data = c(mod.avg.pred, uncond.se, lower.CL, upper.CL),
                          nrow = nrow(newdata), ncol = 4)
  colnames(matrix.output) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
  
  ##organize as list
  Mod.pred.list <- list("type" = type, "mod.avg.pred" = mod.avg.pred, "uncond.se" = uncond.se,
                        "conf.level" = conf.level, "lower.CL" = lower.CL, "upper.CL" = upper.CL,
                        "matrix.output" = matrix.output)
  class(Mod.pred.list) <- c("modavgPred", "list")
  return(Mod.pred.list)
}



##gpcount
modavgPred.AICunmarkedFitGPC <-
  function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           type = "response", c.hat = 1, parm.type = NULL, ...) {

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
  if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgPred for details\n")}


  ##rename values according to unmarked to extract from object
  ##lambda
  if(identical(parm.type, "lambda")) {
    parm.type1 <- "lambda"
  }
  
      ##detect
      if(identical(parm.type, "detect")) {
          parm.type1 <- "det"
      }

      ##availability
      if(identical(parm.type, "phi")) {
          parm.type1 <- "phi"
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
      
  ##extract inverse link function
  if(identical(select.link, "logistic")) {
    link.inv <- plogis
  }
  if(identical(select.link, "exp")) {
    link.inv <- exp
  }
####changes#############

  
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  
  ##determine number of observations in new data set
  nobserv <- dim(newdata)[1]
  
  ##determine number of columns in new data set
  ncolumns <- dim(newdata)[2]

  ##if only 1 column, add an additional column to avoid problems in computation with predictSE.mer( )
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
 
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames,
                    second.ord = second.ord, nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = nobserv, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.est", "Uncond.SE")
  
  ##required for CI
  if(identical(type, "response")) {
    Mod.avg.out.link <- matrix(NA, nrow = nobserv, ncol = 2)
    colnames(Mod.avg.out.link) <- c("Mod.avg.est", "Uncond.SE")
  }

  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
    
    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      AICctmp <- AICctab
      
      if(identical(type, "response")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1)$Predicted))
        
        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                            type = parm.type1)$SE))
        if(length(unique.link) == 1) {
          AICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))
          AICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                           type = parm.type1, backTransform = FALSE)$SE))
        } else {
          warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
          AICctmp$fit.link <- NA
          AICctmp$SE.link <- NA
        }
      }
    
      
      if(identical(type, "link")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))

        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          type = parm.type1, backTransform = FALSE)$SE))
      }
      

      AICctmp$fit <- fit
      AICctmp$SE <- SE
    
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }        
      }
    }
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
    
    for (obs in 1:nobserv) {
      
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
  
      if(identical(type, "response")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                           type = parm.type1)$Predicted))
      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          type = parm.type1)$SE))
      if(length(unique.link) == 1) {
        QAICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                           type = parm.type1, backTransform = FALSE)$Predicted))
        QAICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                          type = parm.type1, backTransform = FALSE)$SE)) * sqrt(c.hat)
      } else {
        warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
        QAICctmp$fit.link <- NA
        QAICctmp$SE.link <- NA
      }
    }
      
      if(identical(type, "link")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                           type = parm.type1, backTransform = FALSE)$Predicted))
      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          type = parm.type1, backTransform = FALSE)$SE))
    }

      
      QAICctmp$fit <- fit
      QAICctmp$SE <- SE * sqrt(c.hat)
      
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(QAICctmp$QAICcWt * QAICctmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(QAICctmp$QAICcWt * sqrt(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.link^2 + (QAICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(QAICctmp$QAICcWt *(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.link^2 + (QAICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }
    }
  }
  
  
  ##begin loop - AIC
  if(second.ord == FALSE && c.hat == 1) {

    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      AICtmp <- AICctab
      
      if(identical(type, "response")) {
          fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                               type = parm.type1)$Predicted))

          ##extract SE for fitted value for observation obs
          SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                              type = parm.type1)$SE))
          if(length(unique.link) == 1) {
            AICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                             type = parm.type1, backTransform = FALSE)$Predicted))
            AICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                            type = parm.type1, backTransform = FALSE)$SE))
            
          } else {
            warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
            AICtmp$fit.link <- NA
            AICtmp$SE.link <- NA
          }
        }
      

      if(identical(type, "link")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))
        
        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                            type = parm.type1, backTransform = FALSE)$SE))
      }
        
      
      AICtmp$fit <- fit
      AICtmp$SE <- SE

      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }  
    }
  }

  
  ##begin loop - QAIC
  if(second.ord == FALSE && c.hat > 1){
     
    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      QAICtmp <- AICctab
      
      if(identical(type, "response")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = parm.type1)$Predicted))

      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = parm.type1)$SE))
      
      if(length(unique.link) == 1) {
        QAICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                          type = parm.type1, backTransform = FALSE)$Predicted))
        QAICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                         type = parm.type1, backTransform = FALSE)$SE)) * sqrt(c.hat)
      } else {
        warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
        QAICtmp$fit.link <- NA
        QAICtmp$SE.link <- NA
      }
    }
      
      if(identical(type, "link")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = parm.type1, backTransform = FALSE)$Predicted))

      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = parm.type1, backTransform = FALSE)$SE))
    }
      
      
      QAICtmp$fit <- fit
      QAICtmp$SE <- SE * sqrt(c.hat)
      
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(QAICtmp$QAICWt * QAICtmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(QAICtmp$QAICWt * sqrt(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICtmp$QAICWt*QAICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.link^2 + (QAICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(QAICtmp$QAICWt *(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICtmp$QAICWt*QAICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.link^2 + (QAICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }
    }
  }

##############changes start here
  mod.avg.pred <- Mod.avg.out[, 1]
  uncond.se <- Mod.avg.out[, 2]
  
  ##compute confidence interval
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)

  if(identical(type, "link")) {
    lower.CL <- mod.avg.pred - zcrit * uncond.se
    upper.CL <- mod.avg.pred + zcrit * uncond.se
  } else {
    lower.CL <- link.inv(Mod.avg.out.link[, 1] - zcrit * Mod.avg.out.link[, 2])
    upper.CL <- link.inv(Mod.avg.out.link[, 1] + zcrit * Mod.avg.out.link[, 2])
  }
  
  ##create matrix
  matrix.output <- matrix(data = c(mod.avg.pred, uncond.se, lower.CL, upper.CL),
                          nrow = nrow(newdata), ncol = 4)
  colnames(matrix.output) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
  
  ##organize as list
  Mod.pred.list <- list("type" = type, "mod.avg.pred" = mod.avg.pred, "uncond.se" = uncond.se,
                        "conf.level" = conf.level, "lower.CL" = lower.CL, "upper.CL" = upper.CL,
                        "matrix.output" = matrix.output)
  class(Mod.pred.list) <- c("modavgPred", "list")
  return(Mod.pred.list)
}  



##multmixOpen
modavgPred.AICunmarkedFitMMO <-
  function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           type = "response", c.hat = 1, parm.type = NULL, ...) {

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
      if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgPred for details\n")}

  
      ##rename values according to unmarked to extract from object
      ##lambda
      if(identical(parm.type, "lambda")) {
          parm.type1 <- "lambda"
          ##check mixture type for mixture models
          mixture.type <- sapply(X = cand.set, FUN = function(i) i@mixture)
          unique.mixture <- unique(mixture.type)
          if(length(unique.mixture) > 1) {
              if(any(unique.mixture == "ZIP")) stop("\nThis function does not yet support mixing ZIP with other distributions\n")
          } else {
              mixture.id <- unique(mixture.type)
              if(identical(unique.mixture, "ZIP")) {
                  if(identical(type, "link")) stop("\nLink scale not yet supported for ZIP mixtures\n")
                  if(identical(type, "response")) warning("\nModel-averaging linear predictor from a ZIP model not yet implemented\n")
              }
          }

      }



      ##gamma
      if(identical(parm.type, "gamma")) {
          parm.type1 <- "gamma"

      }
  
      ##omega
      if(identical(parm.type, "omega")) {
          parm.type1 <- "omega"

      }

      ##iota (for immigration = TRUE with dynamics = "autoreg", "trend", "ricker", or "gompertz")
      if(identical(parm.type, "iota")) {
          parm.type1 <- "iota"
          ##check that parameter appears in all models
          parfreq <- sum(sapply(cand.set, FUN = function(i) any(names(i@estimates@estimates) == parm.type1)))
          if(!identical(length(cand.set), parfreq)) {
              stop("\nParameter \'", parm.type1, "\' (parm.type = \"", parm.type, "\") does not appear in all models:",
                   "\ncannot compute model-averaged predictions across all models\n")
          }

      }

      ##detect
      if(identical(parm.type, "detect")) {

          parm.type1 <- "det"

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
      
      ##extract inverse link function
      if(identical(select.link, "logistic")) {
          link.inv <- plogis
      }
      if(identical(select.link, "exp")) {
          link.inv <- exp
      }
     
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  
  ##determine number of observations in new data set
  nobserv <- dim(newdata)[1]
  
  ##determine number of columns in new data set
  ncolumns <- dim(newdata)[2]

  ##if only 1 column, add an additional column to avoid problems in computation with predictSE.mer( )
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
 
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames,
                    second.ord = second.ord, nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = nobserv, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.est", "Uncond.SE")

  ##required for CI
  if(identical(type, "response")) {
    Mod.avg.out.link <- matrix(NA, nrow = nobserv, ncol = 2)
    colnames(Mod.avg.out.link) <- c("Mod.avg.est", "Uncond.SE")
  }
  

  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      AICctmp <- AICctab
        
      if(identical(type, "response")) {
        
        ##extract fitted value for observation obs
        if(identical(parm.type, "lambda") && identical(mixture.id, "ZIP")) {
          fit <- unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE,
                                               newdata = newdata[obs, ])$fit))
          SE <- unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE,
                                              newdata = newdata[obs, ])$se.fit))
          AICctmp$fit.link <- NA
          AICctmp$SE.link <- NA
          
        } else {
          
          fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                               type = parm.type1)$Predicted))
          
          ##extract SE for fitted value for observation obs
          SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                              type = parm.type1)$SE))
          
          
          if(length(unique.link) == 1) {
            AICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                              type = parm.type1, backTransform = FALSE)$Predicted))
            AICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                             type = parm.type1, backTransform = FALSE)$SE)) * sqrt(c.hat)
          } else {
            warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
            AICctmp$fit.link <- NA
            AICctmp$SE.link <- NA
          } 
        }
      }

      
      if(identical(type, "link")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = parm.type1, backTransform = FALSE)$Predicted))

      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = parm.type1, backTransform = FALSE)$SE))
    }
      
      ##create temporary data.frame to store fitted values and SE 
      AICctmp$fit <- fit
      AICctmp$SE <- SE
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }  
      }
    }
  }


    ##begin loop - QAICc
    if(second.ord == TRUE && c.hat > 1){
      for (obs in 1:nobserv) {

        ##create temporary data.frame to store fitted values and SE
        QAICctmp <- AICctab

        if(identical(type, "response")) {
          
          ##extract fitted value for observation obs
          if(identical(parm.type, "lambda") && identical(mixture.id, "ZIP")) {
            fit <- unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE,
                                                 newdata = newdata[obs, ])$fit))
            SE <- unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE,
                                                newdata = newdata[obs, ])$se.fit))
            QAICctmp$fit.link <- NA
            QAICctmp$SE.link <- NA

          } else {
        
            ##extract fitted value for observation obs
            fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                 type = parm.type1)$Predicted))
            
            ##extract SE for fitted value for observation obs
            SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                type = parm.type1)$SE))


            if(length(unique.link) == 1) {
              QAICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                 type = parm.type1, backTransform = FALSE)$Predicted))
              QAICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                type = parm.type1, backTransform = FALSE)$SE)) * sqrt(c.hat)
            } else {
              warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
              QAICctmp$fit.link <- NA
              QAICctmp$SE.link <- NA
            }
          }
        }       
              
        
        if(identical(type, "link")) {
          ##extract fitted value for observation obs
          fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                               type = parm.type1, backTransform = FALSE)$Predicted))

          ##extract SE for fitted value for observation obs
          SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                              type = parm.type1, backTransform = FALSE)$SE))
        }
      
        ##create temporary data.frame to store fitted values and SE 
        QAICctmp$fit <- fit
        QAICctmp$SE <- SE * sqrt(c.hat)
      
        ##compute model averaged prediction and store in output matrix
        Mod.avg.out[obs, 1] <- sum(QAICctmp$QAICcWt * QAICctmp$fit)
        ##compute unconditional SE and store in output matrix

        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
        if(identical(uncond.se, "old")) {
          Mod.avg.out[obs, 2] <- sum(QAICctmp$QAICcWt * sqrt(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2))
          if(identical(type, "response")) {
            Mod.avg.out.link[obs, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$fit.link)
            Mod.avg.out.link[obs, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.link^2 + (QAICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
          }
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
          Mod.avg.out[obs, 2] <- sqrt(sum(QAICctmp$QAICcWt *(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2)))
          if(identical(type, "response")) {
            Mod.avg.out.link[obs, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$fit.link)
            Mod.avg.out.link[obs, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.link^2 + (QAICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
          }
        }
      }
    }
  
  
  
    ##create temporary data.frame to store fitted values and SE - AIC
    if(second.ord == FALSE && c.hat == 1) {
      for (obs in 1:nobserv) {

        ##create temporary data.frame to store fitted values and SE 
        AICtmp <- AICctab
        
        if(identical(type, "response")) {
        
          if(identical(parm.type, "lambda") && identical(mixture.id, "ZIP")) {
            fit <- unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE,
                                                 newdata = newdata[obs, ])$fit))
            SE <- unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE,
                                                newdata = newdata[obs, ])$se.fit))
            AICtmp$fit.link <- NA
            AICtmp$SE.link <- NA
            
          } else {
            fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                 type = parm.type1)$Predicted))
            
            ##extract SE for fitted value for observation obs
            SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                type = parm.type1)$SE))

            if(length(unique.link) == 1) {
              AICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                               type = parm.type1, backTransform = FALSE)$Predicted))
              AICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                              type = parm.type1, backTransform = FALSE)$SE))
            } else {
              warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
              AICtmp$fit.link <- NA
              AICtmp$SE.link <- NA
            }
          }
        }


        if(identical(type, "link")) {
          ##extract fitted value for observation obs
          fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                               type = parm.type1, backTransform = FALSE)$Predicted))

          ##extract SE for fitted value for observation obs
          SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                              type = parm.type1, backTransform = FALSE)$SE))
        }
        

        AICtmp$fit <- fit
        AICtmp$SE <- SE

        ##compute model averaged prediction and store in output matrix
        Mod.avg.out[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit)
        
        ##compute unconditional SE and store in output matrix
        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
        if(identical(uncond.se, "old")) {
          Mod.avg.out[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2))
          if(identical(type, "response")) {
            Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
            Mod.avg.out.link[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
          }
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
          Mod.avg.out[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2)))
          if(identical(type, "response")) {
            Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
            Mod.avg.out.link[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
          }        
        }  
      }
    }

  
    ##begin loop - QAIC
    if(second.ord == FALSE && c.hat > 1){
      for (obs in 1:nobserv) {

        ##create temporary data.frame to store fitted values and SE 
        QAICtmp <- AICctab
        
        if(identical(type, "response")) {
          
          ##extract fitted value for observation obs
          if(identical(parm.type, "lambda") && identical(mixture.id, "ZIP")) {
            fit <- unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE,
                                                 newdata = newdata[obs, ])$fit))
            SE <- unlist(lapply(X = cand.set, FUN = function(i)predictSE(i, se.fit = TRUE,
                                                newdata = newdata[obs, ])$se.fit))
            QAICtmp$fit.link <- NA
            QAICtmp$SE.link <- NA
            
          } else {
            
            ##extract fitted value for observation obs
            fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                 type = parm.type1)$Predicted))

            ##extract SE for fitted value for observation obs
            SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                type = parm.type1)$SE))

            if(length(unique.link) == 1) {
              QAICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                                type = parm.type1, backTransform = FALSE)$Predicted))
              QAICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                               type = parm.type1, backTransform = FALSE)$SE)) * sqrt(c.hat)
            } else {
              warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
              QAICtmp$fit.link <- NA
              QAICtmp$SE.link <- NA
            }
          }
        }
          
        if(identical(type, "link")) {
          ##extract fitted value for observation obs
          fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                               type = parm.type1, backTransform = FALSE)$Predicted))
          
          ##extract SE for fitted value for observation obs
          SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                              type = parm.type1, backTransform = FALSE)$SE))
        }
      
        ##create temporary data.frame to store fitted values and SE 
        QAICtmp$fit <- fit
        QAICtmp$SE <- SE * sqrt(c.hat)
        
        ##compute model averaged prediction and store in output matrix
        Mod.avg.out[obs, 1] <- sum(QAICtmp$QAICWt * QAICtmp$fit)
        ##compute unconditional SE and store in output matrix

        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
        if(identical(uncond.se, "old")) {
          Mod.avg.out[obs, 2] <- sum(QAICtmp$QAICWt * sqrt(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2))
          if(identical(type, "response")) {
            Mod.avg.out.link[obs, 1] <- sum(QAICtmp$QAICWt*QAICtmp$fit.link)
            Mod.avg.out.link[obs, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.link^2 + (QAICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
          }
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
          Mod.avg.out[obs, 2] <- sqrt(sum(QAICtmp$QAICWt *(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2)))
          if(identical(type, "response")) {
            Mod.avg.out.link[obs, 1] <- sum(QAICtmp$QAICWt*QAICtmp$fit.link)
            Mod.avg.out.link[obs, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.link^2 + (QAICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
          }        
        }
      }
    }


##############changes start here
    mod.avg.pred <- Mod.avg.out[, 1]
    uncond.se <- Mod.avg.out[, 2]
  
    ##compute confidence interval
    zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)

    if(identical(type, "link")) {
      lower.CL <- mod.avg.pred - zcrit * uncond.se
      upper.CL <- mod.avg.pred + zcrit * uncond.se
    } else {
      lower.CL <- link.inv(Mod.avg.out.link[, 1] - zcrit * Mod.avg.out.link[, 2])
      upper.CL <- link.inv(Mod.avg.out.link[, 1] + zcrit * Mod.avg.out.link[, 2])
    }
  
    ##create matrix
    matrix.output <- matrix(data = c(mod.avg.pred, uncond.se, lower.CL, upper.CL),
                            nrow = nrow(newdata), ncol = 4)
    colnames(matrix.output) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
  
    ##organize as list
    Mod.pred.list <- list("type" = type, "mod.avg.pred" = mod.avg.pred, "uncond.se" = uncond.se,
                          "conf.level" = conf.level, "lower.CL" = lower.CL, "upper.CL" = upper.CL,
                          "matrix.output" = matrix.output)
    class(Mod.pred.list) <- c("modavgPred", "list")
    return(Mod.pred.list)
  }



##distsampOpen
modavgPred.AICunmarkedFitDSO <-
  function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           type = "response", c.hat = 1, parm.type = NULL, ...) {

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
  if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgPred for details\n")}

  
      ##rename values according to unmarked to extract from object
      ##lambda
      if(identical(parm.type, "lambda")) {
          parm.type1 <- "lambda"

          ##check mixture type for mixture models
          mixture.type <- sapply(X = cand.set, FUN = function(i) i@mixture)
          unique.mixture <- unique(mixture.type)
          if(length(unique.mixture) > 1) {
              if(any(unique.mixture == "ZIP")) stop("\nThis function does not yet support mixing ZIP with other distributions\n")
          } else {
              mixture.id <- unique(mixture.type)
              if(identical(unique.mixture, "ZIP")) {
                  if(identical(type, "link")) stop("\nLink scale not yet supported for ZIP mixtures\n")
                  if(identical(type, "response")) warning("\nModel-averaging linear predictor from a ZIP model not yet implemented\n")
              }
          }

      }
  
      ##gamma - recruitment
      if(identical(parm.type, "gamma")) {
          parm.type1 <- "gamma"

      }

      ##omega
      if(identical(parm.type, "omega")) {
          parm.type1 <- "omega"

      }

      ##iota (for immigration = TRUE with dynamics = "autoreg", "trend", "ricker", or "gompertz")
      if(identical(parm.type, "iota")) {
          parm.type1 <- "iota"
          ##check that parameter appears in all models
          parfreq <- sum(sapply(cand.set, FUN = function(i) any(names(i@estimates@estimates) == parm.type1)))
          if(!identical(length(cand.set), parfreq)) {
              stop("\nParameter \'", parm.type1, "\' (parm.type = \"", parm.type, "\") does not appear in all models:",
                   "\ncannot compute model-averaged predictions across all models\n")
          }

      }
      
      ##detect
      if(identical(parm.type, "detect")){
      parm.type1 <- "det"
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

      ##extract inverse link function
      if(identical(select.link, "logistic")) {
          link.inv <- plogis
      }
      if(identical(select.link, "exp")) {
          link.inv <- exp
      }
  
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  
  ##determine number of observations in new data set
  nobserv <- dim(newdata)[1]
  
  ##determine number of columns in new data set
  ncolumns <- dim(newdata)[2]

  ##if only 1 column, add an additional column to avoid problems in computation with predictSE.mer( )
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
 
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames,
                    second.ord = second.ord, nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = nobserv, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.est", "Uncond.SE")

  ##required for CI
  if(identical(type, "response")) {
    Mod.avg.out.link <- matrix(NA, nrow = nobserv, ncol = 2)
    colnames(Mod.avg.out.link) <- c("Mod.avg.est", "Uncond.SE")
  }
  

  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
    
    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      AICctmp <- AICctab
      
      if(identical(type, "response")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1)$Predicted))
        
        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                            type = parm.type1)$SE))
        if(length(unique.link) == 1) {
          AICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))
          AICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                           type = parm.type1, backTransform = FALSE)$SE))
        } else {
          warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
          AICctmp$fit.link <- NA
          AICctmp$SE.link <- NA
        }
      }
    
      
      if(identical(type, "link")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))

        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          type = parm.type1, backTransform = FALSE)$SE))
      }
      

      AICctmp$fit <- fit
      AICctmp$SE <- SE
    
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }        
      }
    }
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
    
    for (obs in 1:nobserv) {
      
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
  
      if(identical(type, "response")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                           type = parm.type1)$Predicted))
      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          type = parm.type1)$SE))
      if(length(unique.link) == 1) {
        QAICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                           type = parm.type1, backTransform = FALSE)$Predicted))
        QAICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                          type = parm.type1, backTransform = FALSE)$SE)) * sqrt(c.hat)
      } else {
        warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
        QAICctmp$fit.link <- NA
        QAICctmp$SE.link <- NA
      }
    }
      
      if(identical(type, "link")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                           type = parm.type1, backTransform = FALSE)$Predicted))
      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          type = parm.type1, backTransform = FALSE)$SE))
    }

      
      QAICctmp$fit <- fit
      QAICctmp$SE <- SE * sqrt(c.hat)
      
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(QAICctmp$QAICcWt * QAICctmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(QAICctmp$QAICcWt * sqrt(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.link^2 + (QAICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(QAICctmp$QAICcWt *(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.link^2 + (QAICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }
    }
  }
  
  
  ##begin loop - AIC
  if(second.ord == FALSE && c.hat == 1) {

    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      AICtmp <- AICctab
      
      if(identical(type, "response")) {
          fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                               type = parm.type1)$Predicted))

          ##extract SE for fitted value for observation obs
          SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                              type = parm.type1)$SE))
          if(length(unique.link) == 1) {
            AICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                             type = parm.type1, backTransform = FALSE)$Predicted))
            AICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                            type = parm.type1, backTransform = FALSE)$SE))
            
          } else {
            warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
            AICtmp$fit.link <- NA
            AICtmp$SE.link <- NA
          }
        }
      

      if(identical(type, "link")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))
        
        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                            type = parm.type1, backTransform = FALSE)$SE))
      }
        
      
      AICtmp$fit <- fit
      AICtmp$SE <- SE

      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }  
    }
  }

  
  ##begin loop - QAIC
  if(second.ord == FALSE && c.hat > 1){
     
    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      QAICtmp <- AICctab
      
      if(identical(type, "response")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = parm.type1)$Predicted))

      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = parm.type1)$SE))
      
      if(length(unique.link) == 1) {
        QAICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                          type = parm.type1, backTransform = FALSE)$Predicted))
        QAICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                         type = parm.type1, backTransform = FALSE)$SE)) * sqrt(c.hat)
      } else {
        warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
        QAICtmp$fit.link <- NA
        QAICtmp$SE.link <- NA
      }
    }
      
      if(identical(type, "link")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = parm.type1, backTransform = FALSE)$Predicted))

      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = parm.type1, backTransform = FALSE)$SE))
    }
      
      
      QAICtmp$fit <- fit
      QAICtmp$SE <- SE * sqrt(c.hat)
      
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(QAICtmp$QAICWt * QAICtmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(QAICtmp$QAICWt * sqrt(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICtmp$QAICWt*QAICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.link^2 + (QAICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(QAICtmp$QAICWt *(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICtmp$QAICWt*QAICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.link^2 + (QAICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }
    }
  }

##############changes start here
  mod.avg.pred <- Mod.avg.out[, 1]
  uncond.se <- Mod.avg.out[, 2]
  
  ##compute confidence interval
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)

  if(identical(type, "link")) {
    lower.CL <- mod.avg.pred - zcrit * uncond.se
    upper.CL <- mod.avg.pred + zcrit * uncond.se
  } else {
    lower.CL <- link.inv(Mod.avg.out.link[, 1] - zcrit * Mod.avg.out.link[, 2])
    upper.CL <- link.inv(Mod.avg.out.link[, 1] + zcrit * Mod.avg.out.link[, 2])
  }
  
  ##create matrix
  matrix.output <- matrix(data = c(mod.avg.pred, uncond.se, lower.CL, upper.CL),
                          nrow = nrow(newdata), ncol = 4)
  colnames(matrix.output) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
  
  ##organize as list
  Mod.pred.list <- list("type" = type, "mod.avg.pred" = mod.avg.pred, "uncond.se" = uncond.se,
                        "conf.level" = conf.level, "lower.CL" = lower.CL, "upper.CL" = upper.CL,
                        "matrix.output" = matrix.output)
  class(Mod.pred.list) <- c("modavgPred", "list")
  return(Mod.pred.list)
}



##occuTTD
modavgPred.AICunmarkedFitOccuTTD <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                             nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                             type = "response", c.hat = 1, parm.type = NULL, ...) {

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
  if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgPred for details\n")}

  
  ##rename values according to unmarked to extract from object
  ##psi
  if(identical(parm.type, "psi")) {
      parm.type1 <- "psi"

  }

    ##gamma - colonization
    if(identical(parm.type, "gamma")) {
        nseasons <- unique(sapply(cand.set, FUN = function(i) i@data@numPrimary))
        if(nseasons == 1) {
            stop("\nParameter \'gamma\' does not appear in single-season models\n")
        }
        parm.type1 <- "col"
        
    }

    ##epsilon - extinction
    if(identical(parm.type, "epsilon")) {
        nseasons <- unique(sapply(cand.set, FUN = function(i) i@data@numPrimary))
        if(nseasons == 1) {
            stop("\nParameter \'epsilon\' does not appear in single-season models\n")
        }
        parm.type1 <- "ext"
    }
        
          
    ##detect
    if(identical(parm.type, "detect")) {
        parm.type1 <- "det"

    } #lambda is rate parameter of species not detected at time t to be detected at next time step

    
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

    ##extract inverse link function
    if(identical(select.link, "logistic")) {
        link.inv <- plogis
    }

    ##extract inverse link function
    if(identical(select.link, "exp")) {
        link.inv <- exp
    }

    ##extract inverse link function
    if(identical(select.link, "cloglog")) {
        link.inv <- function(x) 1 - exp(-exp(x))
    }
    
####changes#############
  
  
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  
  ##determine number of observations in new data set
  nobserv <- dim(newdata)[1]
  
  ##determine number of columns in new data set
  ncolumns <- dim(newdata)[2]

  ##if only 1 column, add an additional column to avoid problems in computation with predict( )
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
 
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames,
                    second.ord = second.ord, nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = nobserv, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.est", "Uncond.SE")
  

  ##required for CI
  if(identical(type, "response")) {
    Mod.avg.out.link <- matrix(NA, nrow = nobserv, ncol = 2)
    colnames(Mod.avg.out.link) <- c("Mod.avg.est", "Uncond.SE")
  }
  
  
  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
    
    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      AICctmp <- AICctab
      
      if(identical(type, "response")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1)$Predicted))
        
        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                            type = parm.type1)$SE))
        if(length(unique.link) == 1) {
          AICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))
          AICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                             type = parm.type1, backTransform = FALSE)$SE))
        } else {
          warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
          AICctmp$fit.link <- NA
          AICctmp$SE.link <- NA
        }
      }
    
      
      if(identical(type, "link")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))

        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          type = parm.type1, backTransform = FALSE)$SE))
      }
      

      AICctmp$fit <- fit
      AICctmp$SE <- SE
    
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }        
      }
    }
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
    
    for (obs in 1:nobserv) {
      
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
  
      if(identical(type, "response")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                           type = parm.type1)$Predicted))
      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          type = parm.type1)$SE))
      if(length(unique.link) == 1) {
        QAICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                           type = parm.type1, backTransform = FALSE)$Predicted))
        QAICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                          type = parm.type1, backTransform = FALSE)$SE)) * sqrt(c.hat)
      } else {
        warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
        QAICctmp$fit.link <- NA
        QAICctmp$SE.link <- NA
      }
    }
      
      if(identical(type, "link")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                           type = parm.type1, backTransform = FALSE)$Predicted))
      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          type = parm.type1, backTransform = FALSE)$SE))
    }

      
      QAICctmp$fit <- fit
      QAICctmp$SE <- SE * sqrt(c.hat)
      
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(QAICctmp$QAICcWt * QAICctmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(QAICctmp$QAICcWt * sqrt(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.link^2 + (QAICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(QAICctmp$QAICcWt *(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.link^2 + (QAICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }
    }
  }
  
  
  
  ##begin loop - AIC
  if(second.ord == FALSE && c.hat == 1) {

    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      AICtmp <- AICctab
      
      if(identical(type, "response")) {
          fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                               type = parm.type1)$Predicted))

          ##extract SE for fitted value for observation obs
          SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                              type = parm.type1)$SE))
          if(length(unique.link) == 1) {
            AICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                             type = parm.type1, backTransform = FALSE)$Predicted))
            AICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                            type = parm.type1, backTransform = FALSE)$SE))
            
          } else {
            warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
            AICtmp$fit.link <- NA
            AICtmp$SE.link <- NA
          }
        }
      

      if(identical(type, "link")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))
        
        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                            type = parm.type1, backTransform = FALSE)$SE))
      }
        
      
      AICtmp$fit <- fit
      AICtmp$SE <- SE

      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }  
    }
  }

  
  ##begin loop - QAICc
  if(second.ord == FALSE && c.hat > 1){
     
    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      QAICtmp <- AICctab
      
      if(identical(type, "response")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = parm.type1)$Predicted))

      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = parm.type1)$SE))
      
      if(length(unique.link) == 1) {
        QAICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                          type = parm.type1, backTransform = FALSE)$Predicted))
        QAICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                         type = parm.type1, backTransform = FALSE)$SE)) * sqrt(c.hat)
      } else {
        warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
        QAICtmp$fit.link <- NA
        QAICtmp$SE.link <- NA
      }
    }
      
      if(identical(type, "link")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = parm.type1, backTransform = FALSE)$Predicted))

      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = parm.type1, backTransform = FALSE)$SE))
    }
      
      
      QAICtmp$fit <- fit
      QAICtmp$SE <- SE * sqrt(c.hat)
      
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(QAICtmp$QAICWt * QAICtmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(QAICtmp$QAICWt * sqrt(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICtmp$QAICWt*QAICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.link^2 + (QAICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(QAICtmp$QAICWt *(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICtmp$QAICWt*QAICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.link^2 + (QAICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }
    }
  }

##############changes start here
  mod.avg.pred <- Mod.avg.out[, 1]
  uncond.se <- Mod.avg.out[, 2]
  
  ##compute confidence interval
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)

  if(identical(type, "link")) {
    lower.CL <- mod.avg.pred - zcrit * uncond.se
    upper.CL <- mod.avg.pred + zcrit * uncond.se
  } else {
    lower.CL <- link.inv(Mod.avg.out.link[, 1] - zcrit * Mod.avg.out.link[, 2])
    upper.CL <- link.inv(Mod.avg.out.link[, 1] + zcrit * Mod.avg.out.link[, 2])
  }
  
  ##create matrix
  matrix.output <- matrix(data = c(mod.avg.pred, uncond.se, lower.CL, upper.CL),
                          nrow = nrow(newdata), ncol = 4)
  colnames(matrix.output) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
  
  ##organize as list
  Mod.pred.list <- list("type" = type, "mod.avg.pred" = mod.avg.pred, "uncond.se" = uncond.se,
                        "conf.level" = conf.level, "lower.CL" = lower.CL, "upper.CL" = upper.CL,
                        "matrix.output" = matrix.output)
  class(Mod.pred.list) <- c("modavgPred", "list")
  return(Mod.pred.list)
}



##occuMS
modavgPred.AICunmarkedFitOccuMS <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                            nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                            type = "response", c.hat = 1, parm.type = NULL, ...) {

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
  if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgPred for details\n")}


    ##check if type = "link"
    if(identical(type, "link")) stop("\nLink scale predictions not yet supported for this model type\n")
    
    ##rename values according to unmarked to extract from object
    ##psi
    if(identical(parm.type, "psi")) {
        parm.type1 <- "state"
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
        parm.type1 <- "transition"
        ##because the same elements have different labels in different parts of the results of this object type
        parm.type.alt <- parm.type
    }
          
          
    ##detect
    if(identical(parm.type, "detect")) {
        parm.type1 <- "det"
        ##because the same elements have different labels in different parts of the results of this object type
        parm.type.alt <- parm.type1
    }


####changes#############
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
    
    ##extract inverse link function
    #if(identical(select.link, "multinomial")) {
    #    link.inv <- TODO
    #}
    
    ##extract inverse link function
    #if(identical(select.link, "logistic")) {
    #    link.inv <- plogis
    #}
####changes#############
  
  
    ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  
    ##determine number of observations in new data set
    nobserv <- dim(newdata)[1]
  
    ##determine number of columns in new data set
    ncolumns <- dim(newdata)[2]

    ##if only 1 column, add an additional column to avoid problems in computation with predict( )
    if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
 
    ##store AICc table
    AICctab <- aictab(cand.set = cand.set, modnames = modnames,
                      second.ord = second.ord, nobs = nobs, sort = FALSE, c.hat = c.hat)


#################################################
#################################################
###CHANGES

    ##create object to hold Model-averaged estimates and unconditional SE's
    ##Mod.avg.out <- matrix(NA, nrow = nobserv, ncol = 2)
    ##colnames(Mod.avg.out) <- c("Mod.avg.est", "Uncond.SE")

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
    
    ##prepare list for model-averaged predictions and SE's
    predsOut <-  array(data = NA, dim = c(nobserv, 4, nparms),
                       dimnames = list(1:nobserv,
                                       c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL"),
                                       parmNames))

    ##begin loop - AICc
    if(second.ord == TRUE && c.hat == 1){
    
        for (obs in 1:nobserv) {

            ##create temporary data.frame to store fitted values and SE 
            AICctmp <- AICctab

            ##compute unconditional SE and store in output matrix
            ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    
            if(identical(uncond.se, "old")) {
                for(j in 1:nparms) { 
                    predsOut[obs, 1, j] <- sum(AICctmp$AICcWt * predsEst[obs, , j])
                    predsOut[obs, 2, j] <- sum(AICctmp$AICcWt * sqrt(predsSE[obs, , j]^2 + (predsEst[obs, , j] - predsOut[obs, 1, j])^2))
                }
            }

            ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
            if(identical(uncond.se, "revised")) {
                for(j in 1:nparms) { 
                    predsOut[obs, 1, j] <- sum(AICctmp$AICcWt * predsEst[obs, , j])
                    predsOut[obs, 2, j] <- sqrt(sum(AICctmp$AICcWt * (predsSE[obs, , j]^2 + (predsEst[obs, , j] - predsOut[obs, 1, j])^2)))
                }
            }
        }
    }

    


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){

    for (obs in 1:nobserv) {

            ##create temporary data.frame to store fitted values and SE 
            QAICctmp <- AICctab

            ##compute unconditional SE and store in output matrix
            ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    
            if(identical(uncond.se, "old")) {
                for(j in 1:nparms) { 
                    predsOut[obs, 1, j] <- sum(QAICctmp$QAICcWt * predsEst[obs, , j])
                    predsOut[obs, 2, j] <- sum(QAICctmp$QAICcWt * sqrt(predsSE[obs, , j]^2 + (predsEst[obs, , j] - predsOut[obs, 1, j])^2))
                }
            }

            ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
            if(identical(uncond.se, "revised")) {
                for(j in 1:nparms) { 
                    predsOut[obs, 1, j] <- sum(QAICctmp$QAICcWt * predsEst[obs, , j])
                    predsOut[obs, 2, j] <- sqrt(sum(QAICctmp$QAICcWt * (predsSE[obs, , j]^2 + (predsEst[obs, , j] - predsOut[obs, 1, j])^2)))
                }
            }
        }
    }

    
    ##begin loop - AIC
    if(second.ord == FALSE && c.hat == 1) {

    for (obs in 1:nobserv) {

            ##create temporary data.frame to store fitted values and SE 
            AICtmp <- AICctab

            ##compute unconditional SE and store in output matrix
            ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    
            if(identical(uncond.se, "old")) {
                for(j in 1:nparms) { 
                    predsOut[obs, 1, j] <- sum(AICtmp$AICWt * predsEst[obs, , j])
                    predsOut[obs, 2, j] <- sum(AICtmp$AICWt * sqrt(predsSE[obs, , j]^2 + (predsEst[obs, , j] - predsOut[obs, 1, j])^2))
                }
            }

            ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
            if(identical(uncond.se, "revised")) {
                for(j in 1:nparms) { 
                    predsOut[obs, 1, j] <- sum(AICtmp$AICWt * predsEst[obs, , j])
                    predsOut[obs, 2, j] <- sqrt(sum(AICtmp$AICWt * (predsSE[obs, , j]^2 + (predsEst[obs, , j] - predsOut[obs, 1, j])^2)))
                }
            }
        }
    }

    
    ##begin loop - QAICc
    if(second.ord == FALSE && c.hat > 1){

        for (obs in 1:nobserv) {

            ##create temporary data.frame to store fitted values and SE 
            QAICtmp <- AICctab

            ##compute unconditional SE and store in output matrix
            ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    
            if(identical(uncond.se, "old")) {
                for(j in 1:nparms) { 
                    predsOut[obs, 1, j] <- sum(QAICtmp$QAICWt * predsEst[obs, , j])
                    predsOut[obs, 2, j] <- sum(QAICtmp$QAICWt * sqrt(predsSE[obs, , j]^2 + (predsEst[obs, , j] - predsOut[obs, 1, j])^2))
                }
            }

            ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
            if(identical(uncond.se, "revised")) {
                for(j in 1:nparms) { 
                    predsOut[obs, 1, j] <- sum(QAICtmp$QAICWt * predsEst[obs, , j])
                    predsOut[obs, 2, j] <- sqrt(sum(QAICtmp$QAICWt * (predsSE[obs, , j]^2 + (predsEst[obs, , j] - predsOut[obs, 1, j])^2)))
                }
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
    arrayToDF <- expand.grid(dimnames(predsOut)[c(1, 3)])
    ##rename columns
    names(arrayToDF)[1:2] <- c("Observation", "Parameter")
    rawFrame <- data.frame(arrayToDF, arrayToMat)
    predsOut.mat <- as.matrix(rawFrame[, c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")])
    #rawFrame <- expand.grid(Obs = 1:matRows, Parms = names(predsOut))
    ##create label for rows
    rowLab <- paste(rawFrame$Parameter, rawFrame$Observation, sep = "-")
    rownames(predsOut.mat) <- rowLab

    ##convert array to list
    #predsOutList <- lapply(seq(dim(predsOut)[3]), function(i) predsOut[ , , i])
    #names(predsOutList) <- parmNames
    
    ##organize as list
    Mod.pred.list <- list("type" = type, "mod.avg.pred" = predsOut[, 1, ],
                          "uncond.se" = predsOut[, 2, ],
                          "conf.level" = conf.level,
                          "lower.CL" = predsOut[, 3, ],
                          "upper.CL" = predsOut[, 4, ],
                          "matrix.output" = predsOut.mat)
    class(Mod.pred.list) <- c("modavgPred", "list")
    return(Mod.pred.list)
}



##occuMulti
modavgPred.AICunmarkedFitOccuMulti <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                               nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                               type = "response", c.hat = 1, parm.type = NULL, ...) {

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
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgPred for details\n")}

    ##check if type = "link"
    if(identical(type, "link")) stop("\nLink scale predictions not yet supported for this model type\n")
    
    ##rename values according to unmarked to extract from object
    ##psi
    if(identical(parm.type, "psi")) {
        parm.type1 <- "state"
    }

    ##detect
    if(identical(parm.type, "detect")) {
        parm.type1 <- "det"
    }


####changes#############
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
    
    ##extract inverse link function
    #if(identical(select.link, "multinomial")) {
    #    link.inv <- TODO
    #}
    
    ##extract inverse link function
    #if(identical(select.link, "logistic")) {
    #    link.inv <- plogis
    #}
####changes#############
  
  
    ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  
    ##determine number of observations in new data set
    nobserv <- dim(newdata)[1]
  
    ##determine number of columns in new data set
    ncolumns <- dim(newdata)[2]

    ##if only 1 column, add an additional column to avoid problems in computation with predict( )
    if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
 
    ##store AICc table
    AICctab <- aictab(cand.set = cand.set, modnames = modnames,
                      second.ord = second.ord, nobs = nobs, sort = FALSE, c.hat = c.hat)


#################################################
#################################################
###CHANGES

    ##create object to hold Model-averaged estimates and unconditional SE's
    ##Mod.avg.out <- matrix(NA, nrow = nobserv, ncol = 2)
    ##colnames(Mod.avg.out) <- c("Mod.avg.est", "Uncond.SE")

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
    
    ##prepare list for model-averaged predictions and SE's
    predsOut <-  array(data = NA, dim = c(nobserv, 4, nparms),
                       dimnames = list(1:nobserv,
                                       c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL"),
                                       parmNames))
    
    ##begin loop - AICc
    if(second.ord == TRUE && c.hat == 1){
    
        for (obs in 1:nobserv) {

            ##create temporary data.frame to store fitted values and SE 
            AICctmp <- AICctab

            ##compute unconditional SE and store in output matrix
            ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    
            if(identical(uncond.se, "old")) {
                for(j in 1:nparms) { 
                    predsOut[obs, 1, j] <- sum(AICctmp$AICcWt * predsEst[obs, , j])
                    predsOut[obs, 2, j] <- sum(AICctmp$AICcWt * sqrt(predsSE[obs, , j]^2 + (predsEst[obs, , j] - predsOut[obs, 1, j])^2))
                }
            }

            ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
            if(identical(uncond.se, "revised")) {
                for(j in 1:nparms) { 
                    predsOut[obs, 1, j] <- sum(AICctmp$AICcWt * predsEst[obs, , j])
                    predsOut[obs, 2, j] <- sqrt(sum(AICctmp$AICcWt * (predsSE[obs, , j]^2 + (predsEst[obs, , j] - predsOut[obs, 1, j])^2)))
                }
            }
        }
    }

    


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){

    for (obs in 1:nobserv) {

            ##create temporary data.frame to store fitted values and SE 
            QAICctmp <- AICctab

            ##compute unconditional SE and store in output matrix
            ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    
            if(identical(uncond.se, "old")) {
                for(j in 1:nparms) { 
                    predsOut[obs, 1, j] <- sum(QAICctmp$QAICcWt * predsEst[obs, , j])
                    predsOut[obs, 2, j] <- sum(QAICctmp$QAICcWt * sqrt(predsSE[obs, , j]^2 + (predsEst[obs, , j] - predsOut[obs, 1, j])^2))
                }
            }

            ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
            if(identical(uncond.se, "revised")) {
                for(j in 1:nparms) { 
                    predsOut[obs, 1, j] <- sum(QAICctmp$QAICcWt * predsEst[obs, , j])
                    predsOut[obs, 2, j] <- sqrt(sum(QAICctmp$QAICcWt * (predsSE[obs, , j]^2 + (predsEst[obs, , j] - predsOut[obs, 1, j])^2)))
                }
            }
        }
    }

    
    ##begin loop - AIC
    if(second.ord == FALSE && c.hat == 1) {

    for (obs in 1:nobserv) {

            ##create temporary data.frame to store fitted values and SE 
            AICtmp <- AICctab

            ##compute unconditional SE and store in output matrix
            ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    
            if(identical(uncond.se, "old")) {
                for(j in 1:nparms) { 
                    predsOut[obs, 1, j] <- sum(AICtmp$AICWt * predsEst[obs, , j])
                    predsOut[obs, 2, j] <- sum(AICtmp$AICWt * sqrt(predsSE[obs, , j]^2 + (predsEst[obs, , j] - predsOut[obs, 1, j])^2))
                }
            }

            ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
            if(identical(uncond.se, "revised")) {
                for(j in 1:nparms) { 
                    predsOut[obs, 1, j] <- sum(AICtmp$AICWt * predsEst[obs, , j])
                    predsOut[obs, 2, j] <- sqrt(sum(AICtmp$AICWt * (predsSE[obs, , j]^2 + (predsEst[obs, , j] - predsOut[obs, 1, j])^2)))
                }
            }
        }
    }

    
    ##begin loop - QAICc
    if(second.ord == FALSE && c.hat > 1){

        for (obs in 1:nobserv) {

            ##create temporary data.frame to store fitted values and SE 
            QAICtmp <- AICctab

            ##compute unconditional SE and store in output matrix
            ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    
            if(identical(uncond.se, "old")) {
                for(j in 1:nparms) { 
                    predsOut[obs, 1, j] <- sum(QAICtmp$QAICWt * predsEst[obs, , j])
                    predsOut[obs, 2, j] <- sum(QAICtmp$QAICWt * sqrt(predsSE[obs, , j]^2 + (predsEst[obs, , j] - predsOut[obs, 1, j])^2))
                }
            }

            ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
            if(identical(uncond.se, "revised")) {
                for(j in 1:nparms) { 
                    predsOut[obs, 1, j] <- sum(QAICtmp$QAICWt * predsEst[obs, , j])
                    predsOut[obs, 2, j] <- sqrt(sum(QAICtmp$QAICWt * (predsSE[obs, , j]^2 + (predsEst[obs, , j] - predsOut[obs, 1, j])^2)))
                }
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
    arrayToDF <- expand.grid(dimnames(predsOut)[c(1, 3)])
    ##rename columns
    names(arrayToDF)[1:2] <- c("Observation", "Parameter")
    rawFrame <- data.frame(arrayToDF, arrayToMat)
    predsOut.mat <- as.matrix(rawFrame[, c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")])
    #rawFrame <- expand.grid(Obs = 1:matRows, Parms = names(predsOut))
    ##create label for rows
    rowLab <- paste(rawFrame$Parameter, rawFrame$Observation, sep = "-")
    rownames(predsOut.mat) <- rowLab

    ##convert array to list
    #predsOutList <- lapply(seq(dim(predsOut)[3]), function(i) predsOut[ , , i])
    #names(predsOutList) <- parmNames
    
    ##organize as list
    Mod.pred.list <- list("type" = type, "mod.avg.pred" = predsOut[, 1, ],
                          "uncond.se" = predsOut[, 2, ],
                          "conf.level" = conf.level,
                          "lower.CL" = predsOut[, 3, ],
                          "upper.CL" = predsOut[, 4, ],
                          "matrix.output" = predsOut.mat)
        
    class(Mod.pred.list) <- c("modavgPred", "list")
    return(Mod.pred.list)
}



##goccu
modavgPred.AICunmarkedFitGOccu <-
  function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           type = "response", c.hat = 1, parm.type = NULL, ...) {

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
  if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgPred for details\n")}


  ##rename values according to unmarked to extract from object
  ##psi
  if(identical(parm.type, "psi")) {
    parm.type1 <- "psi"
  }
  
      ##detect
      if(identical(parm.type, "detect")) {
          parm.type1 <- "det"
      }

      ##availability
      if(identical(parm.type, "phi")) {
          parm.type1 <- "phi"
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
      
  ##extract inverse link function
  if(identical(select.link, "logistic")) {
    link.inv <- plogis
  }
  ##extract inverse link function
  if(identical(select.link, "cloglog")) {
    link.inv <- function(x) 1 - exp(-exp(x))
  }
####changes#############

  
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  
  ##determine number of observations in new data set
  nobserv <- dim(newdata)[1]
  
  ##determine number of columns in new data set
  ncolumns <- dim(newdata)[2]

  ##if only 1 column, add an additional column to avoid problems in computation with predictSE.mer( )
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
 
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames,
                    second.ord = second.ord, nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = nobserv, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.est", "Uncond.SE")
  
  ##required for CI
  if(identical(type, "response")) {
    Mod.avg.out.link <- matrix(NA, nrow = nobserv, ncol = 2)
    colnames(Mod.avg.out.link) <- c("Mod.avg.est", "Uncond.SE")
  }

  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
    
    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      AICctmp <- AICctab
      
      if(identical(type, "response")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1)$Predicted))
        
        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                            type = parm.type1)$SE))
        if(length(unique.link) == 1) {
          AICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))
          AICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                           type = parm.type1, backTransform = FALSE)$SE))
        } else {
          warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
          AICctmp$fit.link <- NA
          AICctmp$SE.link <- NA
        }
      }
    
      
      if(identical(type, "link")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))

        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          type = parm.type1, backTransform = FALSE)$SE))
      }
      

      AICctmp$fit <- fit
      AICctmp$SE <- SE
    
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }        
      }
    }
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
    
    for (obs in 1:nobserv) {
      
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
  
      if(identical(type, "response")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                           type = parm.type1)$Predicted))
      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          type = parm.type1)$SE))
      if(length(unique.link) == 1) {
        QAICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                           type = parm.type1, backTransform = FALSE)$Predicted))
        QAICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                          type = parm.type1, backTransform = FALSE)$SE)) * sqrt(c.hat)
      } else {
        warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
        QAICctmp$fit.link <- NA
        QAICctmp$SE.link <- NA
      }
    }
      
      if(identical(type, "link")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                           type = parm.type1, backTransform = FALSE)$Predicted))
      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          type = parm.type1, backTransform = FALSE)$SE))
    }

      
      QAICctmp$fit <- fit
      QAICctmp$SE <- SE * sqrt(c.hat)
      
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(QAICctmp$QAICcWt * QAICctmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(QAICctmp$QAICcWt * sqrt(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.link^2 + (QAICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(QAICctmp$QAICcWt *(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.link^2 + (QAICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }
    }
  }
  
  
  ##begin loop - AIC
  if(second.ord == FALSE && c.hat == 1) {

    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      AICtmp <- AICctab
      
      if(identical(type, "response")) {
          fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                               type = parm.type1)$Predicted))

          ##extract SE for fitted value for observation obs
          SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                              type = parm.type1)$SE))
          if(length(unique.link) == 1) {
            AICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                             type = parm.type1, backTransform = FALSE)$Predicted))
            AICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                            type = parm.type1, backTransform = FALSE)$SE))
            
          } else {
            warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
            AICtmp$fit.link <- NA
            AICtmp$SE.link <- NA
          }
        }
      

      if(identical(type, "link")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))
        
        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                            type = parm.type1, backTransform = FALSE)$SE))
      }
        
      
      AICtmp$fit <- fit
      AICtmp$SE <- SE

      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }  
    }
  }

  
  ##begin loop - QAIC
  if(second.ord == FALSE && c.hat > 1){
     
    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      QAICtmp <- AICctab
      
      if(identical(type, "response")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = parm.type1)$Predicted))

      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = parm.type1)$SE))
      
      if(length(unique.link) == 1) {
        QAICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                          type = parm.type1, backTransform = FALSE)$Predicted))
        QAICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                         type = parm.type1, backTransform = FALSE)$SE)) * sqrt(c.hat)
      } else {
        warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
        QAICtmp$fit.link <- NA
        QAICtmp$SE.link <- NA
      }
    }
      
      if(identical(type, "link")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = parm.type1, backTransform = FALSE)$Predicted))

      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = parm.type1, backTransform = FALSE)$SE))
    }
      
      
      QAICtmp$fit <- fit
      QAICtmp$SE <- SE * sqrt(c.hat)
      
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(QAICtmp$QAICWt * QAICtmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(QAICtmp$QAICWt * sqrt(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICtmp$QAICWt*QAICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.link^2 + (QAICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(QAICtmp$QAICWt *(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICtmp$QAICWt*QAICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.link^2 + (QAICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }
    }
  }

##############changes start here
  mod.avg.pred <- Mod.avg.out[, 1]
  uncond.se <- Mod.avg.out[, 2]
  
  ##compute confidence interval
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)

  if(identical(type, "link")) {
    lower.CL <- mod.avg.pred - zcrit * uncond.se
    upper.CL <- mod.avg.pred + zcrit * uncond.se
  } else {
    lower.CL <- link.inv(Mod.avg.out.link[, 1] - zcrit * Mod.avg.out.link[, 2])
    upper.CL <- link.inv(Mod.avg.out.link[, 1] + zcrit * Mod.avg.out.link[, 2])
  }
  
  ##create matrix
  matrix.output <- matrix(data = c(mod.avg.pred, uncond.se, lower.CL, upper.CL),
                          nrow = nrow(newdata), ncol = 4)
  colnames(matrix.output) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
  
  ##organize as list
  Mod.pred.list <- list("type" = type, "mod.avg.pred" = mod.avg.pred, "uncond.se" = uncond.se,
                        "conf.level" = conf.level, "lower.CL" = lower.CL, "upper.CL" = upper.CL,
                        "matrix.output" = matrix.output)
  class(Mod.pred.list) <- c("modavgPred", "list")
  return(Mod.pred.list)
}  



##occuComm
modavgPred.AICunmarkedFitOccuComm <- function(cand.set, modnames = NULL, newdata, second.ord = TRUE,
                                              nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                                              type = "response", c.hat = 1, parm.type = NULL, ...) {

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
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?modavgPred for details\n")}

  
    ##rename values according to unmarked to extract from object
    ##psi
    if(identical(parm.type, "psi")) {
    parm.type1 <- "state"
    
    }

  
    ##detect
    if(identical(parm.type, "detect")) {
        parm.type1 <- "det"
        
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
    
    ##extract inverse link function
    if(identical(select.link, "logistic")) {
        link.inv <- plogis
    }

    ##extract inverse link function
    if(identical(select.link, "cloglog")) {
        link.inv <- function(x) 1 - exp(-exp(x))
    }

    
    ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  
    ##determine number of observations in new data set
    nobserv <- dim(newdata)[1]
  
  ##determine number of columns in new data set
  ncolumns <- dim(newdata)[2]

  ##if only 1 column, add an additional column to avoid problems in computation with predict( )
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
 
  ##store AICc table
  AICctab <- aictab(cand.set = cand.set, modnames = modnames,
                    second.ord = second.ord, nobs = nobs, sort = FALSE, c.hat = c.hat)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = nobserv, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.est", "Uncond.SE")
  

  ##required for CI
  if(identical(type, "response")) {
    Mod.avg.out.link <- matrix(NA, nrow = nobserv, ncol = 2)
    colnames(Mod.avg.out.link) <- c("Mod.avg.est", "Uncond.SE")
  }
  
  
  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
    
    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      AICctmp <- AICctab
      
      if(identical(type, "response")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1)$Predicted))
        
        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                            type = parm.type1)$SE))
        if(length(unique.link) == 1) {
          AICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))
          AICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                             type = parm.type1, backTransform = FALSE)$SE))
        } else {
          warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
          AICctmp$fit.link <- NA
          AICctmp$SE.link <- NA
        }
      }
    
      
      if(identical(type, "link")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))

        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          type = parm.type1, backTransform = FALSE)$SE))
      }
      

      AICctmp$fit <- fit
      AICctmp$SE <- SE
    
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE.link^2 + (AICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }        
      }
    }
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
    
    for (obs in 1:nobserv) {
      
      ##create temporary data.frame to store fitted values and SE 
      QAICctmp <- AICctab
  
      if(identical(type, "response")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                           type = parm.type1)$Predicted))
      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          type = parm.type1)$SE))
      if(length(unique.link) == 1) {
        QAICctmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                           type = parm.type1, backTransform = FALSE)$Predicted))
        QAICctmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                          type = parm.type1, backTransform = FALSE)$SE)) * sqrt(c.hat)
      } else {
        warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
        QAICctmp$fit.link <- NA
        QAICctmp$SE.link <- NA
      }
    }
      
      if(identical(type, "link")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                           type = parm.type1, backTransform = FALSE)$Predicted))
      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                          type = parm.type1, backTransform = FALSE)$SE))
    }

      
      QAICctmp$fit <- fit
      QAICctmp$SE <- SE * sqrt(c.hat)
      
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(QAICctmp$QAICcWt * QAICctmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(QAICctmp$QAICcWt * sqrt(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE.link^2 + (QAICctmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(QAICctmp$QAICcWt *(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICctmp$QAICcWt*QAICctmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE.link^2 + (QAICctmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }
    }
  }
  
  
  
  ##begin loop - AIC
  if(second.ord == FALSE && c.hat == 1) {

    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      AICtmp <- AICctab
      
      if(identical(type, "response")) {
          fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                               type = parm.type1)$Predicted))

          ##extract SE for fitted value for observation obs
          SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                              type = parm.type1)$SE))
          if(length(unique.link) == 1) {
            AICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                             type = parm.type1, backTransform = FALSE)$Predicted))
            AICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                            type = parm.type1, backTransform = FALSE)$SE))
            
          } else {
            warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
            AICtmp$fit.link <- NA
            AICtmp$SE.link <- NA
          }
        }
      

      if(identical(type, "link")) {
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                             type = parm.type1, backTransform = FALSE)$Predicted))
        
        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                            type = parm.type1, backTransform = FALSE)$SE))
      }
        
      
      AICtmp$fit <- fit
      AICtmp$SE <- SE

      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE.link^2 + (AICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }  
    }
  }

  
  ##begin loop - QAICc
  if(second.ord == FALSE && c.hat > 1){
     
    for (obs in 1:nobserv) {

      ##create temporary data.frame to store fitted values and SE 
      QAICtmp <- AICctab
      
      if(identical(type, "response")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = parm.type1)$Predicted))

      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = parm.type1)$SE))
      
      if(length(unique.link) == 1) {
        QAICtmp$fit.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                          type = parm.type1, backTransform = FALSE)$Predicted))
        QAICtmp$SE.link <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                                         type = parm.type1, backTransform = FALSE)$SE)) * sqrt(c.hat)
      } else {
        warning("\nIt is not appropriate to model-average linear predictors using different link functions\n")
        QAICtmp$fit.link <- NA
        QAICtmp$SE.link <- NA
      }
    }
      
      if(identical(type, "link")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = parm.type1, backTransform = FALSE)$Predicted))

      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = parm.type1, backTransform = FALSE)$SE))
    }
      
      
      QAICtmp$fit <- fit
      QAICtmp$SE <- SE * sqrt(c.hat)
      
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(QAICtmp$QAICWt * QAICtmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(QAICtmp$QAICWt * sqrt(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICtmp$QAICWt*QAICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE.link^2 + (QAICtmp$fit.link - Mod.avg.out.link[obs, 1])^2))
        }
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(QAICtmp$QAICWt *(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2)))
        if(identical(type, "response")) {
          Mod.avg.out.link[obs, 1] <- sum(QAICtmp$QAICWt*QAICtmp$fit.link)
          Mod.avg.out.link[obs, 2] <- sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE.link^2 + (QAICtmp$fit.link - Mod.avg.out.link[obs, 1])^2)))
        }
      }
    }
  }

##############changes start here
  mod.avg.pred <- Mod.avg.out[, 1]
  uncond.se <- Mod.avg.out[, 2]
  
  ##compute confidence interval
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)

  if(identical(type, "link")) {
    lower.CL <- mod.avg.pred - zcrit * uncond.se
    upper.CL <- mod.avg.pred + zcrit * uncond.se
  } else {
    lower.CL <- link.inv(Mod.avg.out.link[, 1] - zcrit * Mod.avg.out.link[, 2])
    upper.CL <- link.inv(Mod.avg.out.link[, 1] + zcrit * Mod.avg.out.link[, 2])
  }
  
  ##create matrix
  matrix.output <- matrix(data = c(mod.avg.pred, uncond.se, lower.CL, upper.CL),
                          nrow = nrow(newdata), ncol = 4)
  colnames(matrix.output) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
  
  ##organize as list
  Mod.pred.list <- list("type" = type, "mod.avg.pred" = mod.avg.pred, "uncond.se" = uncond.se,
                        "conf.level" = conf.level, "lower.CL" = lower.CL, "upper.CL" = upper.CL,
                        "matrix.output" = matrix.output)
  class(Mod.pred.list) <- c("modavgPred", "list")
  return(Mod.pred.list)
}



print.modavgPred <- function(x, digits = 3, ...) {

  if(any(names(x)=="type") ) {
    cat("\nModel-averaged predictions on the ", x$type, " scale\n",
        "based on entire model set and ",
        x$conf.level*100, "% confidence interval:\n\n", sep = "")
  } else {cat("\nModel-averaged predictions based on entire\n",
              "model set and ", x$conf.level*100,
              "% confidence interval:\n\n", sep = "")}

    nice.tab <- x$matrix.output
    colnames(nice.tab) <- c("mod.avg.pred", "uncond.se", "lower.CL", "upper.CL")
    ##add check to display output from occuMS or occuMulti
    ##check if rownames are present
    if(is.null(rownames(nice.tab))){
        nrows <- dim(nice.tab)[1]
        rownames(nice.tab) <- 1:nrows
    }
    print(round(nice.tab, digits = digits))
    cat("\n")
}

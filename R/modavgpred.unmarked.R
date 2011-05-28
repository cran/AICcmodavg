modavgpred.unmarked <- function(cand.set, modnames, newdata, second.ord = TRUE, type = "response",
                                c.hat = 1, nobs = NULL, uncond.se = "revised", parm.type = NULL) {

  ##check model types
  mod.class <- unlist(lapply(cand.set, FUN = function(i) class(i)[1]))
  mod.type <- unique(mod.class)

  if(length(mod.type) > 1) stop("This function is not appropriate to model-average parameters from different model types")

  ##check for supported mod.type
  supp.class <- c("unmarkedFitOccu", "unmarkedFitColExt", "unmarkedFitOccuRN", "unmarkedFitPCount", "unmarkedFitPCO")
                  
  if(!any(supp.class == mod.type)) {stop("\nFunction not yet defined for this object class\n")}
  
  ##rename values according to unmarked to extract from object
  ##psi
  if(identical(parm.type, "psi")) {
    if(identical(mod.type, "unmarkedFitOccu")) {parm.type1 <- "state"}
    if(identical(mod.type, "unmarkedFitColExt")) {parm.type1 <- "psi"}
  }

  ##gamma
  if(identical(parm.type, "gamma")) {
    if(identical(mod.type, "unmarkedFitColExt")) {parm.type1 <- "col"}
    if(identical(mod.type, "unmarkedFitPCO")) {parm.type1 <- "gamma"}
  }

  ##epsilon
  if(identical(parm.type, "epsilon")) {
    if(identical(mod.type, "unmarkedFitColExt")) {parm.type1 <- "ext"}
  }

  ##lambda
  if(identical(parm.type, "lambda")) {
    if(identical(mod.type, "unmarkedFitPCount")) {parm.type1 <- "state"}
    if(identical(mod.type, "unmarkedFitPCO")) {parm.type1 <- "lambda"}
    if(identical(mod.type, "unmarkedFitOccuRN")) {parm.type1 <- "state"}
  }

  ##omega
  if(identical(parm.type, "omega")) {
    if(identical(mod.type, "unmarkedFitPCO")) {parm.type1 <- "omega"}
  }

  ##detect
  if(identical(parm.type, "detect")) parm.type1 <- "det"

  
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  
  ##determine number of observations in new data set
  nobserv <- dim(newdata)[1]
  
  ##determine number of columns in new data set
  ncolumns <- dim(newdata)[2]

  ##if only 1 column, add an additional column to avoid problems in computation with predictSE.mer( )
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
 
  ##store AICc table
  AICctab <- aictab.unmarked(cand.set = cand.set, modnames = modnames, c.hat = c.hat, second.ord = second.ord,
                         nobs = nobs, sort = FALSE)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out <- matrix(NA, nrow = nobserv, ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.est", "Uncond.SE")
  

  ##begin loop - AICc
  if(second.ord == TRUE && c.hat == 1){
    for (obs in 1:nobserv) {

      if(identical(type, "response")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = parm.type1)$Predicted))

      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = parm.type1)$SE))
    }
      
      if(identical(type, "link")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = parm.type1, backTransform = FALSE)$Predicted))

      ##extract SE for fitted value for observation obs
      SE<-unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = parm.type1, backTransform = FALSE)$SE))
    }
      
      ##create temporary data.frame to store fitted values and SE 
      AICctmp <- AICctab
      AICctmp$fit <- fit
      AICctmp$SE <- SE
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICctmp$AICcWt*AICctmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICctmp$AICcWt*sqrt(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICctmp$AICcWt*(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2)))
      }
    }
  }


  ##begin loop - QAICc
  if(second.ord == TRUE && c.hat > 1){
    for (obs in 1:nobserv) {

      if(identical(type, "response")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = parm.type1)$Predicted))

      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = parm.type1)$SE))
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
      QAICctmp <- AICctab
      QAICctmp$fit <- fit
      QAICctmp$SE <- SE * sqrt(c.hat)
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(QAICctmp$QAICcWt * QAICctmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(QAICctmp$QAICcWt * sqrt(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(QAICctmp$QAICcWt *(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2)))
      }
    }
  }
  
  
  
  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord == FALSE && c.hat == 1) {
    for (obs in 1:nobserv) {

      if(identical(type, "response")) {
        
        ##extract fitted value for observation obs
        fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                         type = parm.type1)$Predicted))
        ##extract SE for fitted value for observation obs
        SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                        type = parm.type1)$SE))
      }

     if(identical(type, "link")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = parm.type1, backTransform = FALSE)$Predicted))

      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = parm.type1, backTransform = FALSE)$SE))
    }
        
      AICtmp <- AICctab
      AICtmp$fit <- fit
      AICtmp$SE <- SE

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(AICtmp$AICWt*AICtmp$fit)
      
      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(AICtmp$AICWt*sqrt(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(AICtmp$AICWt*(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2)))
      }  
      
    }
  }

  
  ##begin loop - QAICc
  if(second.ord == FALSE && c.hat > 1){
    for (obs in 1:nobserv) {

      if(identical(type, "response")) {
      ##extract fitted value for observation obs
      fit <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                       type = parm.type1)$Predicted))

      ##extract SE for fitted value for observation obs
      SE <- unlist(lapply(X = cand.set, FUN = function(i)predict(i, se.fit = TRUE, newdata = newdata[obs, ],
                                      type = parm.type1)$SE))
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
      QAICtmp <- AICctab
      QAICtmp$fit <- fit
      QAICtmp$SE <- SE * sqrt(c.hat)
      
      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1] <- sum(QAICtmp$QAICWt * QAICtmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2] <- sum(QAICtmp$QAICWt * sqrt(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2] <- sqrt(sum(QAICtmp$QAICWt *(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2)))
      }
    }
  }
  

  Mod.pred.list <- list("mod.avg.pred" = Mod.avg.out[,1], "uncond.se" = Mod.avg.out[,2])
  class(Mod.pred.list) <- c("modavgpred", "list")
  return(Mod.pred.list)
  
}

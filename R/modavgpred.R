modavgpred <- function(cand.set, modnames, newdata, type = "response", c.hat = 1,
                       gamdisp = NULL, second.ord = TRUE, nobs = NULL, uncond.se = "revised"){
  results <- NULL
  known <- rep(0, 3) #create an identifier of class type for lm, glm, lme, and mer
  ##extract classes
  mod.class <- unlist(lapply(X=cand.set, FUN=class))
  ##check if all are identical
  check.class <- unique(mod.class)

  ##determine if lm or glm
  if(identical(check.class, "lm") || identical(check.class, c("glm", "lm"))) {
    results <- modavgpred.glm(cand.set = cand.set, modnames = modnames, newdata = newdata,
                          type = type, c.hat = c.hat, gamdisp = gamdisp, second.ord = second.ord, 
                          nobs = nobs, uncond.se = uncond.se)
    known[1] <- 1
  }

  ##determine if lme
  if(identical(check.class, "lme"))  {
    results <- modavgpred.lme(cand.set = cand.set, modnames = modnames, newdata = newdata,
                          second.ord = second.ord, nobs = nobs, uncond.se = uncond.se)
    known[2] <- 1
  }

  ##determine if mer
  if(identical(check.class, "mer")) {
    results <- modavgpred.mer(cand.set = cand.set, modnames = modnames, newdata = newdata,
                              type = type, c.hat = c.hat, second.ord = second.ord, nobs = nobs,
                              uncond.se = uncond.se)
      known[3] <- 1
    }

    
#warn if class is neither lm, glm, nor lme
    if(sum(known) < 1) {stop("Function not yet defined for this object class")}

    return(results)
  }






modavgpred.glm <-
function(cand.set, modnames, newdata, type = "response", c.hat = 1, gamdisp = NULL, second.ord = TRUE,
         nobs = NULL, uncond.se = "revised") {
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  if(type=="terms") {stop("The terms argument is not defined for this function")}

  ##check family of glm to avoid problems when requesting predictions with argument 'dispersion'
  fam.type <- unlist(lapply(cand.set, FUN=function(i) family(i)$family))
  fam.unique <- unique(fam.type)
  if(identical(fam.unique, "gaussian")) {dispersion <- NULL} #set to NULL if gaussian is used
  ##poisson, binomial, and negative binomial defaults to 1 (no separate parameter for variance)
    
###################CHANGES####
##############################
  if(c.hat>1) {dispersion <- c.hat }
  if(!is.null(gamdisp)) {dispersion <- gamdisp}
  if(c.hat>1 && !is.null(gamdisp)) {stop("You cannot specify values for both \'c.hat\' and \'gamdisp\'")}
  ##dispersion is the dispersion parameter - this influences the SE's (to specify dispersion parameter for either overdispersed Poisson or Gamma glm)
  ##type enables to specify either "response" (original scale = point estimate) or "link" (linear predictor)
  
  ##check if object is of "lm" or "glm" class
  ##extract classes
  mod.class <- unlist(lapply(X=cand.set, FUN=class))
  ##check if all are identical
  check.class <- unique(mod.class)

  ##check that link function is the same for all models if linear predictor is used
  if(identical(type, "link")) {
    check.link <- unlist(lapply(X = cand.set, FUN=function(i) i$family$link))
    unique.link <- unique(x=check.link)
    if(length(unique.link) > 1) {stop(cat("It is not appropriate to compute a model averaged beta estimate\n",
                                          "with different link functions\n"))}
  }

 
  if(identical(check.class, "lm") || identical(check.class, c("glm", "lm")))  {


    ##check if model uses gamma distribution
    gam1<-unlist(lapply(cand.set, FUN=function(i) family(i)$family[1]=="Gamma")) #check for gamma regression models
    ##correct SE's for estimates of gamma regressions when gamdisp is specified
    if(any(gam1)==TRUE)  {
      ##check for specification of gamdisp argument
      if(is.null(gamdisp)) stop("You must specify a gamma dispersion parameter with gamma generalized linear models\n")
    }
  
    
    ##determine number of observations in new data set
    nobserv <- dim(newdata)[1]

    ##determine number of columns in new data set
    ncolumns <- dim(newdata)[2]

    ##if only 1 column, add an additional column to avoid problems in computation with predictSE.mer( )
    if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
 
    ##store AICc table
    AICctab<-aictab(cand.set=cand.set, modnames=modnames, c.hat=c.hat, second.ord=second.ord, nobs=nobs, sort=FALSE)

    ##create object to hold Model-averaged estimates and unconditional SE's
    Mod.avg.out<-matrix(NA, nrow=nobserv, ncol=2)
    colnames(Mod.avg.out)<-c("Mod.avg.est", "Uncond.SE")


    ##begin loop - AICc
    if(second.ord==TRUE && c.hat==1){
      for (obs in 1:nobserv) {

        ##extract fitted value for observation obs
        fit<-unlist(lapply(X=cand.set, FUN=function(i)predict(i, se.fit=TRUE, newdata=newdata[obs, ],
                                         type=type, dispersion=dispersion)$fit))

        ##extract SE for fitted value for observation obs
        SE<-unlist(lapply(X=cand.set, FUN=function(i)predict(i, se.fit=TRUE, newdata=newdata[obs, ],
                                        type=type, dispersion=dispersion)$se.fit))


        ##create temporary data.frame to store fitted values and SE 
        AICctmp<-AICctab
        AICctmp$fit<-fit
        AICctmp$SE<-SE

        ##compute model averaged prediction and store in output matrix
        Mod.avg.out[obs, 1]<-sum(AICctmp$AICcWt*AICctmp$fit)
        ##compute unconditional SE and store in output matrix

        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
        if(identical(uncond.se, "old")) {
          Mod.avg.out[obs, 2]<-sum(AICctmp$AICcWt*sqrt(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2))
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
          Mod.avg.out[obs, 2]<-sqrt(sum(AICctmp$AICcWt*(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2)))
        }
      }
    }




    ##create temporary data.frame to store fitted values and SE - QAICc
    if(second.ord==TRUE && c.hat > 1) {
      for (obs in 1:nobserv) {

        ##extract fitted value for observation obs
        fit<-unlist(lapply(X=cand.set, FUN=function(i)predict(i, se.fit=TRUE, newdata=newdata[obs, ],
                                         type=type, dispersion=dispersion)$fit))
        ##extract SE for fitted value for observation obs
        SE<-unlist(lapply(X=cand.set, FUN=function(i)predict(i, se.fit=TRUE, newdata=newdata[obs, ],
                                        type=type, dispersion=dispersion)$se.fit))

        QAICctmp<-AICctab
        QAICctmp$fit<-fit
        QAICctmp$SE<-SE

        ##compute model averaged prediction and store in output matrix
        Mod.avg.out[obs, 1]<-sum(QAICctmp$QAICcWt*QAICctmp$fit)

        ##compute unconditional SE and store in output matrix
        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
        if(identical(uncond.se, "old")) {
          Mod.avg.out[obs, 2]<-sum(QAICctmp$QAICcWt*sqrt(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2))
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
          Mod.avg.out[obs, 2]<-sqrt(sum(QAICctmp$QAICcWt*(QAICctmp$SE^2 + (QAICctmp$fit- Mod.avg.out[obs, 1])^2)))  
        }
      }
    }




    ##create temporary data.frame to store fitted values and SE - AIC
    if(second.ord==FALSE && c.hat==1) {
      for (obs in 1:nobserv) {

        ##extract fitted value for observation obs
        fit<-unlist(lapply(X=cand.set, FUN=function(i)predict(i, se.fit=TRUE, newdata=newdata[obs, ],
                                         type=type, dispersion=dispersion)$fit))
        ##extract SE for fitted value for observation obs
        SE<-unlist(lapply(X=cand.set, FUN=function(i)predict(i, se.fit=TRUE, newdata=newdata[obs, ], type=type,
                                        dispersion=dispersion)$se.fit))

        AICtmp<-AICctab
        AICtmp$fit<-fit
        AICtmp$SE<-SE

        ##compute model averaged prediction and store in output matrix
        Mod.avg.out[obs, 1]<-sum(AICtmp$AICWt*AICtmp$fit)

        ##compute unconditional SE and store in output matrix
        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
        if(identical(uncond.se, "old")) {
          Mod.avg.out[obs, 2]<-sum(AICtmp$AICWt*sqrt(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2))
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
          Mod.avg.out[obs, 2]<-sqrt(sum(AICtmp$AICWt*(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2)))
        }  
        
      }
    }




    ##create temporary data.frame to store fitted values and SE - QAIC
    if(second.ord==FALSE && c.hat > 1) {
      for (obs in 1:nobserv) {

        ##extract fitted value for observation obs
        fit<-unlist(lapply(X=cand.set, FUN=function(i)predict(i, se.fit=TRUE, newdata=newdata[obs, ],
                                         type=type, dispersion=dispersion)$fit))
        ##extract SE for fitted value for observation obs
        SE<-unlist(lapply(X=cand.set, FUN=function(i)predict(i, se.fit=TRUE, newdata=newdata[obs, ],
                                        type=type, dispersion=dispersion)$se.fit))

        QAICtmp<-AICctab
        QAICtmp$fit<-fit
        QAICtmp$SE<-SE

        ##compute model averaged prediction and store in output matrix
        Mod.avg.out[obs, 1]<-sum(QAICtmp$QAICWt*QAICtmp$fit)

        ##compute unconditional SE and store in output matrix
        if(identical(uncond.se, "old")) {
          Mod.avg.out[obs, 2]<-sum(QAICtmp$QAICWt*sqrt(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2))
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
          Mod.avg.out[obs, 2]<-sqrt(sum(QAICtmp$QAICWt*(QAICtmp$SE^2 + (QAICtmp$fit- Mod.avg.out[obs, 1])^2)))
        }  
        
      }
    }


    Mod.pred.list <- list("type" = type, "mod.avg.pred" = Mod.avg.out[,1], "uncond.se" = Mod.avg.out[,2])
    class(Mod.pred.list) <- c("modavgpred", "list")
    return(Mod.pred.list)

  } else {stop("This function is only appropriate with either \'lm\' or \'glm\' classes\n")}
}






modavgpred.lme <-
function(cand.set, modnames, newdata, second.ord = TRUE,
         nobs = NULL, uncond.se = "revised") {
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)

  ##check if object is of "lm" or "glm" class
  ##extract classes
  mod.class <- unlist(lapply(X=cand.set, FUN=class))
  ##check if all are identical
  check.class <- unique(mod.class)

  if(identical(check.class, "lme") )  {

    ##determine number of observations in new data set
    nobserv <- dim(newdata)[1]

    ##determine number of columns in new data set
    ncolumns <- dim(newdata)[2]

    ##if only 1 column, add an additional column to avoid problems in computation with predictSE.mer( )
    if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
 
    ##store AICc table
    AICctab<-aictab(cand.set=cand.set, modnames=modnames, second.ord=second.ord, nobs=nobs, sort=FALSE)

    ##create object to hold Model-averaged estimates and unconditional SE's
    Mod.avg.out<-matrix(NA, nrow=nobserv, ncol=2)
    colnames(Mod.avg.out)<-c("Mod.avg.est", "Uncond.SE")




    ##begin loop - AICc
    if(second.ord==TRUE){
      for (obs in 1:nobserv) {

        ##extract fitted value for observation obs
        fit<-unlist(lapply(X=cand.set, FUN=function(i)predictSE.lme(i, se.fit=TRUE, newdata=newdata[obs, ])$fit))

        ##extract SE for fitted value for observation obs
        SE<-unlist(lapply(X=cand.set, FUN=function(i)predictSE.lme(i, se.fit=TRUE, newdata=newdata[obs, ])$se.fit))


        ##create temporary data.frame to store fitted values and SE 
        AICctmp<-AICctab
        AICctmp$fit<-fit
        AICctmp$SE<-SE

        ##compute model averaged prediction and store in output matrix
        Mod.avg.out[obs, 1]<-sum(AICctmp$AICcWt*AICctmp$fit)
        ##compute unconditional SE and store in output matrix

        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
        if(identical(uncond.se, "old")) {
          Mod.avg.out[obs, 2]<-sum(AICctmp$AICcWt*sqrt(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2))
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
          Mod.avg.out[obs, 2]<-sqrt(sum(AICctmp$AICcWt*(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2)))
        }
      }
    }



    ##create temporary data.frame to store fitted values and SE - AIC
    if(second.ord==FALSE) {
      for (obs in 1:nobserv) {

        ##extract fitted value for observation obs
        fit<-unlist(lapply(X=cand.set, FUN=function(i)predict(i, se.fit=TRUE, newdata=newdata[obs, ])$fit))
        ##extract SE for fitted value for observation obs
        SE<-unlist(lapply(X=cand.set, FUN=function(i)predict(i, se.fit=TRUE, newdata=newdata[obs, ])$se.fit))

        AICtmp<-AICctab
        AICtmp$fit<-fit
        AICtmp$SE<-SE

        ##compute model averaged prediction and store in output matrix
        Mod.avg.out[obs, 1]<-sum(AICtmp$AICWt*AICtmp$fit)

        ##compute unconditional SE and store in output matrix
        ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
        if(identical(uncond.se, "old")) {
          Mod.avg.out[obs, 2]<-sum(AICtmp$AICWt*sqrt(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2))
        }

        ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
        if(identical(uncond.se, "revised")) {
          Mod.avg.out[obs, 2]<-sqrt(sum(AICtmp$AICWt*(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2)))
        }  
        
      }
    }



    Mod.pred.list <- list("mod.avg.pred" = Mod.avg.out[,1], "uncond.se" = Mod.avg.out[,2])
    class(Mod.pred.list) <- c("modavgpred", "list")
    return(Mod.pred.list)
    
  } else {stop("This function is only appropriate with \'lme\' class\n")}
}



modavgpred.mer <- function(cand.set, modnames, newdata, type = "response", c.hat = 1, second.ord = TRUE,
                           nobs = NULL, uncond.se = "revised") {
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  if(c.hat != 1) {warning("This function only allows \'c.hat = 1\' for \'mer\' class objects\n")}

  ##extract classes
  mod.class <- unlist(lapply(X=cand.set, FUN=class))
  ##check if all are identical
  check.class <- unique(mod.class)

  ##check that link function is the same for all models if linear predictor is used
  if(identical(type, "link")) {
    link.list <- unlist(lapply(X = cand.set, FUN = function(i) fam.link.mer(i)$link))
    check.link <- unique(link.list)
    if(length(check.link) > 1) stop(cat("\nIt is not appropriate to compute a model-averaged beta estimate\n",
                                        "from models using different link functions\n"))
  }

 
  if(!identical(check.class, "mer")) {stop("This function is only appropriate with \'mer\' classes\n")}

       
  ##determine number of observations in data set
  nobserv <- dim(newdata)[1]

  ##determine number of columns in data set
  ncolumns <- dim(newdata)[2]

  ##if only 1 column, add an additional column to avoid problems in computation with predictSE.mer( )
  if(ncolumns == 1) newdata$blank.fake.column.NAs <- NA
  
  ##store AICc table
  AICctab<-aictab(cand.set=cand.set, modnames=modnames, second.ord=second.ord, nobs=nobs, sort=FALSE)

  ##create object to hold Model-averaged estimates and unconditional SE's
  Mod.avg.out<-matrix(NA, nrow=nobserv, ncol=2)
  colnames(Mod.avg.out)<-c("Mod.avg.est", "Uncond.SE")


  ##begin loop - AICc
  if(second.ord==TRUE && c.hat==1){
    for (obs in 1:nobserv) {

      ##extract fitted value for observation obs
      fit<-unlist(lapply(X=cand.set, FUN=function(i)predictSE.mer(i, se.fit=TRUE, newdata=newdata[obs, ],
                                       type=type, level = 0)$fit))

      ##extract SE for fitted value for observation obs
      SE<-unlist(lapply(X=cand.set, FUN=function(i)predictSE.mer(i, se.fit=TRUE, newdata=newdata[obs, ],
                                      type=type, level = 0)$se.fit))


      ##create temporary data.frame to store fitted values and SE 
      AICctmp<-AICctab
      AICctmp$fit<-fit
      AICctmp$SE<-SE

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1]<-sum(AICctmp$AICcWt*AICctmp$fit)
      ##compute unconditional SE and store in output matrix

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2]<-sum(AICctmp$AICcWt*sqrt(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2]<-sqrt(sum(AICctmp$AICcWt*(AICctmp$SE^2 + (AICctmp$fit- Mod.avg.out[obs, 1])^2)))
      }
    }
  }




  ##create temporary data.frame to store fitted values and SE - AIC
  if(second.ord==FALSE && c.hat==1) {
    for (obs in 1:nobserv) {
      
      ##extract fitted value for observation obs
      fit<-unlist(lapply(X=cand.set, FUN=function(i)predictSE.mer(i, se.fit=TRUE, newdata=newdata[obs, ],
                                       type=type, level = 0)$fit))
      ##extract SE for fitted value for observation obs
      SE<-unlist(lapply(X=cand.set, FUN=function(i)predictSE.mer(i, se.fit=TRUE, newdata=newdata[obs, ],
                                      type=type, level = 0)$se.fit))

      AICtmp<-AICctab
      AICtmp$fit<-fit
      AICtmp$SE<-SE

      ##compute model averaged prediction and store in output matrix
      Mod.avg.out[obs, 1]<-sum(AICtmp$AICWt*AICtmp$fit)

      ##compute unconditional SE and store in output matrix
      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Mod.avg.out[obs, 2]<-sum(AICtmp$AICWt*sqrt(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Mod.avg.out[obs, 2]<-sqrt(sum(AICtmp$AICWt*(AICtmp$SE^2 + (AICtmp$fit- Mod.avg.out[obs, 1])^2)))
      }  
    }
  }


  Mod.pred.list <- list("type" = type, "mod.avg.pred" = Mod.avg.out[,1], "uncond.se" = Mod.avg.out[,2])
  class(Mod.pred.list) <- c("modavgpred", "list")
  return(Mod.pred.list)
  
}





print.modavgpred <- function(x, digits = 2, ...) {

  if(any(names(x)=="type") ) {
    cat("\nModel-averaged predictions on the", x$type, "scale based on entire model set:\n\n")
  } else {cat("\nModel-averaged predictions based on entire model set:\n\n")}

  nice.tab <- cbind(x$mod.avg.pred, x$uncond.se)
  colnames(nice.tab) <- c("mod.avg.pred", "uncond.se")
  nrows <- dim(nice.tab)[1]
  rownames(nice.tab) <- 1:nrows
  print(round(nice.tab, digits = digits))
  cat("\n")
}

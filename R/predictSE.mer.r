##########################################################################
##########################################################################
##FUNCTION STARTS HERE
##########################################################################
##########################################################################
##create function
predictSE.mer <- function(mod, newdata, se.fit = TRUE, type = "response", level = 0, print.matrix = FALSE){
  ##check for model type
  if(!identical(paste(class(mod)), "mer")) stop(cat("\nThis function is only appropriate with \'mer\' objects\n"))

  ##check for offset
  if(length(mod@offset) > 0) stop(cat("\nOffsets are not supported at the moment\n"))
  
##########################################################################
###determine characteristics of glmm
##########################################################################
  mod.details <- fam.link.mer(mod)
  fam.type <- mod.details$family
  link.type <- mod.details$link
  supp.link <- mod.details$supp
  
  if(identical(supp.link, "no")) stop("\nOnly canonical link is supported with current version of function\n")
    
  if(identical(link.type, "other")) stop("\nThis function is not yet defined for the specified link function\n")
##########################################################################      


  
  ##this part of code converts data.frame (including factors) into design matrix of model
  tt <- terms(mod)
  TT <- delete.response(tt)
  newdata <- as.data.frame(newdata)

#################################################################################################################
########################### This following clever piece of code is modified from predict.lme( ) from nlme package
#################################################################################################################  
  mfArgs <- list(formula = TT, data = newdata)
  dataMix <- do.call("model.frame", mfArgs)

  ## making sure factor levels are the same as in contrasts

###########this part creates a list to hold factors - changed from nlme
  orig.frame <- mod@frame

  ##matrix with info on factors
  fact.frame <- attr(attr(orig.frame, "terms"), "dataClasses")[-1]

  ##continue if factors
  if(any(fact.frame == "factor")) {
    id.factors <- which(fact.frame == "factor")
    fact.name <- names(fact.frame)[id.factors] #identify the rows for factors

    contr <- list( )
    for(j in fact.name) {
      contr[[j]] <- contrasts(orig.frame[, j])
    }
  }
##########end of code to create list changed from nlme

  
  for(i in names(dataMix)) {
    if (inherits(dataMix[,i], "factor") && !is.null(contr[[i]])) {
      levs <- levels(dataMix[,i])
      levsC <- rownames(contr[[i]])
      if (any(wch <- is.na(match(levs, levsC)))) {
        stop(paste("Levels", paste(levs[wch], collapse = ","),
                   "not allowed for", i))
      }
      attr(dataMix[,i], "contrasts") <- contr[[i]][levs, , drop = FALSE]
    }
  }
#################################################################################################################
########################### The previous clever piece of code is modified from predict.lme( ) from nlme package
#################################################################################################################

  m <- model.frame(TT, data=dataMix)
  des.matrix <- model.matrix(TT, m)
  newdata <- des.matrix  #we now have a design matrix 

  
##logical test for level
  if(!identical(level, 0)) stop(cat("This function does not support computation of predicted values\n",
                                    "or standard errors for higher levels of nesting\n"))


######START OF PREDICT FUNCTION
######
  fix.coef <- fixef(mod)
  ncoefs <- length(fix.coef)
  names.coef <- labels(fix.coef)
  nvals <- dim(newdata)[1]
  
  ##check for intercept fixed effect term in model
  int.yes <- any(names.coef == "(Intercept)")

  ##if no intercept term, return error
  if(!int.yes) stop("This function does not work with models excluding the intercept terms\n")
  
  formula <- character(length=ncoefs)

  nbetas <- ncoefs - 1
  
  if(int.yes & nbetas >= 1) {
    ##create loop to construct formula for derivative
    formula <- paste("Beta", 1:nbetas, sep="")
    formula <- c("Beta0", formula)
  } else {
    if(int.yes & nbetas == 0) {
      formula <- "Beta0"
    }
  }
  ##for models without intercept - formula <- paste("Beta", 1:ncoefs, sep="")
  

  ##a loop to assemble formula
  ##first, identify interaction terms
  inters <- rep(NA, ncoefs)
  for (m in 1:ncoefs) {
    inters[m] <- attr(regexpr(pattern = ":", text = names.coef[m]), "match.length")
  }

  ##change the name of the labels for flexibility
  names.cov <- paste("cov", 1:ncoefs-1, sep="")
  
  if(!int.yes) {names.cov <- paste("cov", 1:ncoefs, sep="")}
  
  id <- which(inters == 1)
  for (k in 1:length(id)) {
    names.cov[id[k]] <- paste("inter", k, sep="")
  }

  ##iterate and combine betas and covariates
  formula2 <- character(length = ncoefs)
  for(b in 1:ncoefs) {
    formula2[b] <- paste(formula[b], names.cov[b], sep="*")
  }
  ##replace with Beta0 if fixed intercept term present
  if(int.yes) {formula2[1] <- "Beta0"}
  
  ##collapse into a single equation and convert to expression
  equation <- parse(text  = as.expression(paste(formula2, collapse="+")),
                    srcfile = NULL)
  ##parse returns the unevaluated expression


  ##determine number of covariates excluding interaction terms
  ncovs <- ncoefs - length(id)

  ##assign values of covariates
  cov.values <- list( )

  ##if only intercept, then add column
  if(int.yes && ncovs == 1) {
    cov.values[[1]] <- 1
  }
  
  if(int.yes && ncovs > 1) {
    cov.values[[1]] <- rep(1, nvals)
    for (q in 2:ncoefs) {
      cov.values[[q]] <- newdata[, labels(fix.coef)[q]]
    }
  } else {
    for (q in 1:ncoefs) {
      cov.values[[q]] <- newdata[, labels(fix.coef)[q]]
    }
  }
    
  names(cov.values) <- names.cov
  cov.values.mat <- matrix(data = unlist(cov.values), nrow = nvals, ncol = ncoefs)

  

################################################################
####use the following code to compute predicted values and SE's
####on response scale if identity link is used OR link scale
  if((identical(type, "response") && identical(link.type, "identity")) || (identical(type, "link"))) {
    
    if(identical(se.fit, TRUE)) {
      ##determine number of partial derivatives to compute
      part.devs <- list( )
      for(j in 1:ncoefs) {
        part.devs[[j]] <- D(equation, formula[j])
      }
    }

    
    if(identical(se.fit, TRUE)) {
      ##substitute a given row for each covariate
      predicted.SE <- matrix(NA, nrow = nvals, ncol = 2)
      colnames(predicted.SE) <- c("Pred.value", "SE")
      rownames(predicted.SE) <- 1:nvals
      part.devs.eval <- list( )
      part.devs.eval[[1]] <- 1
      for (w in 1:nvals) {
        if(int.yes && ncovs > 1) {
          for (p in 2:ncoefs) {
            part.devs.eval[[p]] <-  cov.values[names.cov[p]][[1]][w]
          }
        } # else {  ##for cases without intercept
        ##for (p in 1:ncoefs) {
        ##  part.devs.eval[[p]] <-  cov.values[names.cov[p]][[1]][w]
        ##}
        ##}
    
        part.devs.solved <- unlist(part.devs.eval)
        
        ##extract vc matrix
        vcmat <- vcov(mod)
        
        mat_partialdevs<-as.matrix(part.devs.solved) #create matrix from vector of 2 rows by 1 column
        mat_tpartialdevs<-t(part.devs.solved)        #transpose of partial derivatives to have 2 columns by 1 row
        
        var_hat<-mat_tpartialdevs%*%vcmat%*%mat_partialdevs
        SE<-sqrt(var_hat)
        predicted.vals <- fix.coef%*%cov.values.mat[w,]
        predicted.SE[w, 1] <- predicted.vals
        predicted.SE[w, 2] <- SE@x #to extract only value computed
      }
      
      out.fit.SE <- list(fit = predicted.SE[,"Pred.value"], se.fit = predicted.SE[, "SE"])
      
    } else {
      predicted.SE <- matrix(NA, nrow = nvals, ncol = 1)
      colnames(predicted.SE) <- c("Pred.value")
      rownames(predicted.SE) <- 1:nvals
      for (w in 1:nvals) {
        predicted.vals <- fix.coef%*%cov.values.mat[w,]
        predicted.SE[w, 1] <- predicted.vals
      }
      
      out.fit.SE <- predicted.SE
      colnames(out.fit.SE) <- "fit"
      
    }
  }

###################################################################################
###################################################################################
####use the following code to compute predicted values and SE's
####on response scale if other than identity link is used
###################################################################################
###################################################################################
  
  if(identical(type, "response") && !identical(link.type, "identity")) {
    
    ##for binomial GLMM with logit link
    if(identical(link.type, "logit")) {
      ##build partial derivatives
      logit.eq <- parse(text  = as.expression(paste("exp(", equation, ")/(1 + exp(", equation, "))")),
                        srcfile = NULL)
      part.devs <- list( )
      for(j in 1:ncoefs) {
        part.devs[[j]] <- D(logit.eq, formula[j])
      }
    }

    ##for poisson, gaussian or Gamma GLMM with log link
    if(identical(link.type, "log")) {
      ##build partial derivatives
      log.eq <- parse(text  = as.expression(paste("exp(", equation, ")")),
                        srcfile = NULL)
      part.devs <- list( )
      for(j in 1:ncoefs) {
        part.devs[[j]] <- D(log.eq, formula[j])
      }
    }

    
    ##assign values of beta estimates to beta parameters
    beta.vals <- fix.coef
    names(beta.vals) <- formula

    ##neat way of assigning beta estimate values to objects using names in beta.vals
    for(d in 1:ncoefs) {  
      assign(names(beta.vals)[d], beta.vals[d])
    }

    if(identical(se.fit, TRUE)) {
      ##substitute a given row for each covariate
      predicted.SE <- matrix(NA, nrow = nvals, ncol = 2)
      colnames(predicted.SE) <- c("Pred.value", "SE")
      rownames(predicted.SE) <- 1:nvals
      part.devs.eval <- list( )
      for (w in 1:nvals) {
        if(int.yes && ncovs > 1) {
          for (p in 1:ncoefs) {
            cmds <- list( )
            for(r in 2:ncoefs) {
              ##create commands
              cmds[[r]] <- paste(names.cov[r], "=", "cov.values[[names.cov[", r, "]]][", w, "]")
            } 
            ##assemble commands
            cmd.arg <- paste(unlist(cmds), collapse = ", ")
            cmd.eval <- paste("eval(expr = part.devs[[", p, "]],", "envir = list(", cmd.arg, ")", ")")
            ##evaluate partial derivative
            part.devs.eval[[p]] <- eval(parse(text = cmd.eval))
          }
        }

        if(int.yes && ncovs == 1) {  #for cases with intercept only
          part.devs.eval[[1]] <- eval(part.devs[[1]])
        }
   
        ## else {  ##for cases without intercept
        ##for (p in 1:ncoefs) {
        ##  part.devs.eval[[p]] <-  cov.values[names.cov[p]][[1]][w]
        ##}
        ##}
    
        part.devs.solved <- unlist(part.devs.eval)
        
        ##extract vc matrix
        vcmat <- vcov(mod)
      
        mat_partialdevs<-as.matrix(part.devs.solved) #create matrix from vector of 2 rows by 1 column
        mat_tpartialdevs<-t(part.devs.solved)        #transpose of partial derivatives to have 2 columns by 1 row
      
        var_hat<-mat_tpartialdevs%*%vcmat%*%mat_partialdevs
        SE<-sqrt(var_hat)
        predicted.vals <- fix.coef%*%cov.values.mat[w,]
        if(identical(link.type, "logit")) {
          predicted.SE[w, 1] <- exp(predicted.vals)/(1 + exp(predicted.vals))
        } else {
          if(identical(link.type, "log")) {
            predicted.SE[w, 1] <- exp(predicted.vals)
          }
        }
        predicted.SE[w, 2] <- SE@x #to extract only value computed
      }
      out.fit.SE <- list(fit = predicted.SE[,"Pred.value"], se.fit = predicted.SE[, "SE"])
      
    } else {
      predicted.SE <- matrix(NA, nrow = nvals, ncol = 1)
      colnames(predicted.SE) <- c("Pred.value")
      rownames(predicted.SE) <- 1:nvals
      for (w in 1:nvals) {
        predicted.vals <- fix.coef%*%cov.values.mat[w,]
        if(identical(link.type, "logit")) {
          predicted.SE[w, 1] <- exp(predicted.vals)/(1 + exp(predicted.vals))
        } else {
          if(identical(link.type, "log")) {
            predicted.SE[w, 1] <- exp(predicted.vals)
          }
        }
      }

      out.fit.SE <- predicted.SE
      colnames(out.fit.SE) <- "fit"
     
    }
  }

  
###################################################################
    ##print as nice matrix, otherwise print as list
  if(identical(print.matrix, TRUE)) {
    out.fit.SE <- predicted.SE
    if(identical(se.fit, TRUE)) {
      colnames(out.fit.SE) <- c("fit", "se.fit")
    } else {
      colnames(out.fit.SE) <- c("fit")
    }
  }
  
  return(out.fit.SE)
  
}
###################################################################################
###################################################################################
##END OF PREDICTION FUNCTION
###################################################################################
###################################################################################

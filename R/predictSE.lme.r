predictSE.lme <- function(mod, newdata, se.fit = TRUE, level = 0, print.matrix = FALSE){
  ##check for model type
  if(!identical(class(mod), "lme")) stop(cat("This function is only appropriate with \'lme\' objects"))

  ##first part of code converts data.frame (including factors) into design matrix of model
  fixed <- mod$call$fixed[-2] #extract only fixed portion of model formula
  tt <- terms(mod)
  TT <- delete.response(tt)
  newdata <- as.data.frame(newdata)

#################################################################################################################
########################### This following clever piece of code is modified from predict.lme( ) from nlme package
#################################################################################################################  
  mfArgs <- list(formula = fixed, data = newdata)
  dataMix <- do.call("model.frame", mfArgs)

  ## making sure factor levels are the same as in contrasts
  contr <- mod$contrasts
  for(i in names(dataMix)) {
    if (inherits(dataMix[,i], "factor") && !is.null(contr[[i]])) {
      levs <- levels(dataMix[,i])
      levsC <- dimnames(contr[[i]])[[1]] ##could change to rownames(contr[[i]])
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


##

  if(identical(se.fit, TRUE)) {
  ##determine number of partial derivatives to compute
    part.devs <- list( )
    for(j in 1:ncoefs) {
      part.devs[[j]] <- D(equation, formula[j])
    }
    
  }

  ##determine number of covariates excluding interaction terms
  ncovs <- ncoefs - length(id)

  ##assign values of covariates
  cov.values <- list()

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
      
#####5)
      var_hat<-mat_tpartialdevs%*%vcmat%*%mat_partialdevs
      SE<-sqrt(var_hat)
      predicted.vals <- fix.coef%*%cov.values.mat[w,]
      predicted.SE[w, 1] <- predicted.vals
      predicted.SE[w, 2] <- SE
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

##
##END OF PREDICTION FUNCTION

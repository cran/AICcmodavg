modavg.mer <-
function(cand.set, parm, modnames, conf.level = 0.95, second.ord = TRUE, nobs = NULL, exclude = NULL,
         warn = TRUE, uncond.se = "revised"){

  ##check if class is appropriate
  ##extract classes
  mod.class <- unlist(lapply(X=cand.set, FUN=class))
  ##check if all are identical
  check.class <- unique(mod.class)

  if(!identical(check.class, "mer"))  {stop("\nThis function is only appropriate with the \'mer\' class\n")}

#####MODIFICATIONS BEGIN#######
  ##remove all leading and trailing white space and within parm
  parm <- gsub('[[:space:]]+', "", parm)
  
  ##reverse parm
  reversed.parm <- reverse.parm(parm)
  exclude <- reverse.exclude(exclude = exclude)
#####MODIFICATIONS END######

  

###################
  ##determine families of model
  fam.list <- unlist(lapply(X = cand.set, FUN = function(i) fam.link.mer(i)$family))
  check.fam <- unique(fam.list)
  if(length(check.fam) > 1) stop("\nIt is not appropriate to compute a model-averaged beta estimate\n",
                                       "from models using different families of distributions\n")
  ##determine link functions
  link.list <- unlist(lapply(X = cand.set, FUN = function(i) fam.link.mer(i)$link))
  check.link <- unique(link.list)
  if(length(check.link) > 1) stop("\nIt is not appropriate to compute a model-averaged beta estimate\n",
                                       "from models using different link functions\n")
###################       

  
  ##extract model formula for each model in cand.set
  mod_formula <- lapply(cand.set, FUN=function(i) labels(fixef(i)))

  nmods <- length(cand.set)
  
  ##setup matrix to indicate presence of parms in the model
  include <- matrix(NA, nrow=nmods, ncol=1)
  ##add a check for multiple instances of same variable in given model (i.e., interactions)
  include.check <- matrix(NA, nrow=nmods, ncol=1)

  ##iterate over each formula in mod_formula list
  for (i in 1:nmods) {
    idents <- NULL
    idents.check <- NULL
    form <- mod_formula[[i]]

######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
    ##iterate over each element of formula[[i]] in list
      if(is.null(reversed.parm)) {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)  
        }
      } else {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm, form[j],
                                                                  fixed=TRUE), "match.length")=="-1" , 0, 1)  
        }
      }
###MODIFICATIONS END
######################################################################################################
######################################################################################################

    include[i] <- ifelse(any(idents==1), 1, 0)
    include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
  }

#####################################################
###exclude == NULL; warn=TRUE:  warn that duplicates occur and stop
  if(is.null(exclude) && identical(warn, TRUE)) {
    ##check for duplicates in same model
    if(any(include.check == "duplicates")) {
      stop("\nSome models include more than one instance of the parameter of interest. \n",
           "This may be due to the presence of interaction/polynomial terms, or variables\n",
           "with similar names:\n",
           "\tsee \"?modavg\" for details on variable specification and \"exclude\" argument\n")
    }
    
  }
  
  ##exclude == NULL; warn=FALSE:  compute model-averaged beta estimate from models including variable of interest,
  ##assuming that the variable is not involved in interaction or higher order polynomial (x^2, x^3, etc...),
  ##warn that models were not excluded
  if(is.null(exclude) && identical(warn, FALSE)) {
    if(any(include.check == "duplicates")) {
      warning("\nMultiple instances of parameter of interest in given model is presumably\n",
              "not due to interaction or polynomial terms - these models will not be\n",
              "excluded from the computation of model-averaged estimate\n")
    }
    
  }
  
  ##warn if exclude is neither a list nor NULL
  if(!is.null(exclude)) {
    if(!is.list(exclude)) {stop("\nItems in \"exclude\" must be specified as a list\n")}
  }
  

  ##if exclude is list  
  if(is.list(exclude)) {
    
    ##determine number of elements in exclude
    nexcl <- length(exclude)

    ##check each formula for presence of exclude variable extracted with formula( )  
    not.include <- lapply(cand.set, FUN=formula)  #random effect portion is returned within parentheses
    ##because matching uses identical( ) to check fixed effects against formula( ),
    ##should not be problematic for variables included in random effects
    
    
    ##set up a new list with model formula
    forms <- list()
    for (i in 1:nmods) {
      form.tmp <- strsplit(as.character(not.include[i]), split="~")[[1]][-1]
      if(attr(regexpr("\\+", form.tmp), "match.length")==-1) {
        forms[i] <- form.tmp
      } else {forms[i] <- strsplit(form.tmp, split=" \\+ ")}
    }

    ######################################
    ##remove leading and trailing spaces as well as spaces within string
    ##forms <- lapply(forms.space, FUN = function(b) gsub('[[:space:]]+', "", b)) 
    ######################################
    
    ##additional check to see whether some variable names include "+"
    check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\+", i), "match.length")>0)[[1]]))
    if (any(check.forms==TRUE)) stop("\nPlease avoid \"+\" in variable names\n")
    
    ##additional check to determine if intercept was removed from models
    check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\- 1", i), "match.length")>0)[[1]]))
    if (any(check.forms==TRUE)) stop("\nModels without intercept are not supported in this version, please use alternative parameterization\n")

 
    ##search within formula for variables to exclude
    mod.exclude <- matrix(NA, nrow=nmods, ncol=nexcl)

    ##iterate over each element in exclude list
    for (var in 1:nexcl) {
      
      ##iterate over each formula in mod_formula list
      for (i in 1:nmods) {
        idents <- NULL
        form.excl <- forms[[i]]
                    
        ##iterate over each element of forms[[i]]
        for (j in 1:length(form.excl)) {
          idents[j] <- identical(exclude[var][[1]], gsub("(^ +)|( +$)", "", form.excl[j]))
        }
        mod.exclude[i,var] <- ifelse(any(idents==1), 1, 0)
      }
      
    }
  
    ##determine outcome across all variables to exclude
    to.exclude <- rowSums(mod.exclude)
    
    
    ##exclude models following models from model averaging  
    include[which(to.exclude>=1)] <- 0
    
    
  }


  
  ##add a check to determine if include always == 0
  if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }
  
  new.cand.set <- cand.set[which(include==1)] #select models including a given parameter
  new.mod.name <- modnames[which(include==1)]    #update model names
  
  new_table <- aictab.mer(cand.set=new.cand.set, modnames=new.mod.name, sort=FALSE,
                        second.ord=second.ord, nobs=nobs)  #recompute AIC table and associated measures
  new_table$Beta_est <- unlist(lapply(new.cand.set, FUN=function(i) fixef(i)[paste(parm)])) #extract beta estimate for parm
  new_table$SE <- unlist(lapply(new.cand.set, FUN=function(i) extractSE.mer(i)[paste(parm)]))
  
  ##compute model-averaged estimates, unconditional SE, and 95% CL
  if(second.ord==TRUE) {
    Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)
    
    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
    } 


  } else {
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
  
  zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
  Lower_CL<-Modavg_beta-zcrit*Uncond_SE
  Upper_CL<-Modavg_beta+zcrit*Uncond_SE
  out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                     "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                     "Upper.CL" = Upper_CL)
  
  class(out.modavg) <- c("modavg", "list")
  return(out.modavg)
}


##create SE extractor function that includes estimate labels
extractSE.mer <- function(mod){
  ##extract vcov matrix
  vcov.mat <- as.matrix(vcov(mod))
  se <- sqrt(diag(vcov.mat))
  fixed.labels <- names(fixef(mod))
  names(se) <- fixed.labels
  return(se)
}

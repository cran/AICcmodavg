mult.comp <- function(mod, factor.id, letter.labels = TRUE, sort = TRUE, c.hat = 1, second.ord = TRUE,
                      nobs = NULL) {
##2^(k - 1) combinations of group patterns

  ##check for model class
  ##works with glm, glmer, gls, lm, lme, lmer, rlm
  results <- NULL
  known <- 0
  ##extract classes
  check.class <- class(mod)

  ##determine if lm or glm
  if(identical(check.class, "lm") || identical(check.class, c("glm", "lm"))) {
    known <- 1
  }   

  ##determine if lme
  if(identical(check.class, "lme"))  {
    known <- 1
  }

  ##determine if gls
  if(identical(check.class, "gls"))  {
    known <- 1
  }

  ##determine if mer
  if(identical(check.class[1], "mer"))  {
    known <- 1
  }

  ##determine if merMod
  if(identical(check.class[1], "lmerMod") || identical(check.class[1], "glmerMod") ) {
    known <- 1
    }
  
  ##determine if rlm
  if(identical(check.class, c("rlm", "lm")))  {
    known <- 1
  }
  
  ##warn if class is not supported
  if(sum(known) < 1) {stop("\nFunction not yet defined for this object class\n")}

  ##extract data set from model
  ##unmarked models are not supported yet
  ##check if S4
  if(!isS4(mod)) {
    data.set <- eval(mod$call$data, environment(formula(mod)))
  } else {
    ##if mer or merMod objects
    if(identical(check.class[1], "mer") || identical(check.class[1], "lmerMod") || identical(check.class[1], "glmerMod")) {
      data.set <- eval(mod@call$data, environment(formula(mod)))
    }
#    ##if unmarked object
#    if(identical(attr(class(mod), "package"), "unmarked")) stop("\nModels of unmarked classes are not yet supported\n")
                                        #{data.set <- eval(mod@data, (environment(mod@formula)))}
  }

  ##identify column with factor.id
  parm.vals <- data.set[, factor.id]

  ##check that variable is factor
  if(!is.factor(parm.vals)) stop("\n'factor.id' must be a factor\n")

  ##identify response variable
  resp.id <- as.character(formula(mod))[2]
  ##check for response variable
  if(!any(names(data.set) == resp.id)) stop("\nThis function does not support transformations of the response variable in the formula\n")
  ##changed to allow when response variable is transformed
  ##in the formula to avoid error
  resp <- data.set[, resp.id]
  
  ##determine group identity
  groups.id <- levels(parm.vals)

  ##determine number of groups
  n.groups <- length(groups.id)
  
  ##order group means
  ord.means <- sort(tapply(X = resp, INDEX = parm.vals, FUN = mean))
  ##or sort(tapply(X = fitted(mod), INDEX = parm.vals, FUN = mean))
  
  ##order groups
  ord.groups <- names(ord.means)

  ##generic groups
  gen.groups <- paste(1:n.groups)

  ##create pattern matrix of group assignment
  if(n.groups == 2) {
    ##combinations for 2 groups
    ##{1, 2}
    ##{12}
    ##rows correspond to parameterization of different models
    pat.mat <- matrix(data =
                      c(1, 2,
                        1, 1),
                      byrow = TRUE,
                      ncol = 2)
  }

  if(n.groups == 3) {
    ##combinations for 3 groups
    ##{1, 2, 3}
    ##{12, 3}
    ##{1, 23}
    ##{123}
    pat.mat <- matrix(data =
                      c(1, 2, 3,
                        1, 1, 2,
                        1, 2, 2,
                        1, 1, 1),
                      byrow = TRUE,
                      ncol = 3)
  }

  if(n.groups == 4) {
    ##combinations for 4 groups
    ##{1, 2, 3, 4}
    ##{1, 2, 34}
    ##{12, 3, 4}
    ##{12, 34}
    ##{1, 23, 4}
    ##{1, 234}
    ##{123, 4}
    ##{1234}
    pat.mat <- matrix(data =
                      c(1, 2, 3, 4,
                        1, 2, 3, 3,
                        1, 1, 2, 3,
                        1, 1, 2, 2,
                        1, 2, 2, 3,
                        1, 2, 2, 2,
                        1, 1, 1, 2,
                        1, 1, 1, 1),
                      byrow = TRUE,
                      ncol = 4)
  }

  if(n.groups == 5) {
    ##combinations for 5 groups
    ##{1, 2, 3, 4, 5}
    ##{1, 2, 3, 45}
    ##{1, 2, 345}
    ##{1, 2, 34, 5}
    ##{12, 3, 4, 5}
    ##{12, 34, 5}
    ##{12, 3, 45}
    ##{12, 345}
    ##{1, 23, 4, 5}
    ##{1, 23, 45}
    ##{1, 234, 5}
    ##{123, 4, 5}
    ##{123, 45}
    ##{1234, 5}
    ##{1, 2345}
    ##{12345}
    pat.mat <- matrix(data =
                      c(1, 2, 3, 4, 5,
                        1, 2, 3, 4, 4,
                        1, 2, 3, 3, 3,
                        1, 2, 3, 3, 4,
                        1, 1, 2, 3, 4,
                        1, 1, 2, 2, 3,
                        1, 1, 2, 3, 3,
                        1, 1, 2, 2, 2,
                        1, 2, 2, 3, 4,
                        1, 2, 2, 3, 3,
                        1, 2, 2, 2, 3,
                        1, 1, 1, 2, 3,
                        1, 1, 1, 2, 2,
                        1, 1, 1, 1, 2,
                        1, 2, 2, 2, 2,
                        1, 1, 1, 1, 1),
                      byrow = TRUE,
                      ncol = 5)
  }
    
  ##combinations for 6 groups
  if(n.groups == 6) {
    pat.mat <- matrix(data =
                      c(1, 2, 2, 2, 2, 2,
                        1, 2, 2, 2, 2, 3,
                        1, 2, 2, 2, 3, 3,
                        1, 2, 2, 2, 3, 4,
                        1, 2, 2, 3, 3, 3,
                        1, 2, 3, 3, 3, 4,
                        1, 2, 3, 3, 3, 3,
                        1, 2, 3, 3, 4, 5,
                        1, 2, 3, 3, 4, 4,
                        1, 2, 3, 4, 4, 4,
                        1, 2, 3, 4, 4, 5,
                        1, 2, 3, 4, 5, 6,
                        1, 2, 3, 4, 5, 5,
                        1, 2, 2, 3, 4, 5,
                        1, 2, 2, 3, 4, 4,
                        1, 2, 2, 3, 3, 4,
                        1, 1, 2, 2, 2, 2,
                        1, 1, 2, 2, 2, 3,
                        1, 1, 2, 2, 3, 3,
                        1, 1, 2, 3, 3, 3,
                        1, 1, 1, 2, 2, 2,
                        1, 1, 1, 2, 2, 3,
                        1, 1, 1, 2, 3, 3,
                        1, 1, 1, 2, 3, 4,
                        1, 1, 2, 3, 4, 5,
                        1, 1, 1, 1, 2, 3,
                        1, 1, 1, 1, 2, 2,
                        1, 1, 1, 1, 1, 2,
                        1, 1, 1, 1, 1, 1,
                        1, 1, 2, 3, 4, 4,
                        1, 1, 2, 2, 3, 4,
                        1, 1, 2, 3, 3, 4),
                      byrow = TRUE,
                      ncol = 6)
  }

  if(n.groups > 6) stop("\nThis function supports a maximum of 6 groups\n")

  ##number of models
  n.mods <- nrow(pat.mat)
  data.iter <- data.set
  
  ##create list to store models
  mod.list <- vector("list", length = n.mods)
  
  if(n.groups == 2) {
  ##loop over matrix and change parameterization
  for(k in 1:nrow(pat.mat)) {
    orig.parm <- data.set[, factor.id]
    new.parm <- ifelse(orig.parm == ord.groups[1], pat.mat[k, 1], pat.mat[k, 2])   
  
    ##convert to factor only for cases with different groups
    if(length(unique(new.parm)) > 1) {
      data.iter[, factor.id] <- as.factor(new.parm)
      
      ##run model on updated data
      mod.list[[k]] <- update(mod, data = data.iter)
      
    } else {
      ##if constant, remove factor.id from model
      mod.list[[k]] <- update(mod, as.formula(paste(". ~ . -", factor.id)), data = data.iter)
    }
    
    
  }

}





if(n.groups == 3) {
  ##loop over matrix and change parameterization
  for(k in 1:nrow(pat.mat)) {
    orig.parm <- data.set[, factor.id]
    new.parm <- ifelse(orig.parm == ord.groups[1], pat.mat[k, 1],
                       ifelse(orig.parm == ord.groups[2], pat.mat[k, 2],
                              pat.mat[k, 3]))

    ##convert to factor only for cases with different groups
    if(length(unique(new.parm)) > 1) {
      data.iter[, factor.id] <- as.factor(new.parm)
      
      ##run model on updated data
      mod.list[[k]] <- update(mod, data = data.iter)
      
    } else {
      ##if constant, remove factor.id from model
      mod.list[[k]] <- update(mod, as.formula(paste(". ~ . -", factor.id)), data = data.iter)
    }
    
    
  }

}




if(n.groups == 4) {
  ##loop over matrix and change parameterization
  for(k in 1:nrow(pat.mat)) {
    orig.parm <- data.set[, factor.id]
    new.parm <- ifelse(orig.parm == ord.groups[1], pat.mat[k, 1],
                       ifelse(orig.parm == ord.groups[2], pat.mat[k, 2],
                              ifelse(orig.parm == ord.groups[3], pat.mat[k, 3],
                                     pat.mat[k, 4])))
    
                                        #    data.iter[, factor.id] <- new.parm
    
    ##convert to factor only for cases with different groups
    if(length(unique(new.parm)) > 1) {
      data.iter[, factor.id] <- as.factor(new.parm)
      
      ##run model on updated data
      mod.list[[k]] <- update(mod, data = data.iter)
      
    } else {
      ##if constant, remove factor.id from model
      mod.list[[k]] <- update(mod, as.formula(paste(". ~ . -", factor.id)), data = data.iter)
    }
    
    
  }

}





if(n.groups == 5) {
  ##loop over matrix and change parameterization
  for(k in 1:nrow(pat.mat)) {
    orig.parm <- data.set[, factor.id]
    new.parm <- ifelse(orig.parm == ord.groups[1], pat.mat[k, 1],
                       ifelse(orig.parm == ord.groups[2], pat.mat[k, 2],
                              ifelse(orig.parm == ord.groups[3], pat.mat[k, 3],
                                     ifelse(orig.parm == ord.groups[4], pat.mat[k, 4],
                                            pat.mat[k, 5]))))
    
    ##convert to factor only for cases with different groups
    if(length(unique(new.parm)) > 1) {
      data.iter[, factor.id] <- as.factor(new.parm)
      
      ##run model on updated data
      mod.list[[k]] <- update(mod, data = data.iter)
      
    } else {
      ##if constant, remove factor.id from model
      mod.list[[k]] <- update(mod, as.formula(paste(". ~ . -", factor.id)), data = data.iter)
    }
    
    
  }

}




if(n.groups == 6) {
  ##loop over matrix and change parameterization
  for(k in 1:nrow(pat.mat)) {
    orig.parm <- data.set[, factor.id]
    new.parm <- ifelse(orig.parm == ord.groups[1], pat.mat[k, 1],
                       ifelse(orig.parm == ord.groups[2], pat.mat[k, 2],
                              ifelse(orig.parm == ord.groups[3], pat.mat[k, 3],
                                     ifelse(orig.parm == ord.groups[4], pat.mat[k, 4],
                                            ifelse(orig.parm == ord.groups[5], pat.mat[k, 5],
                                                   pat.mat[k, 6])))))
       
    ##convert to factor only for cases with different groups
    if(length(unique(new.parm)) > 1) {
      data.iter[, factor.id] <- as.factor(new.parm)
      
      ##run model on updated data
      mod.list[[k]] <- update(mod, data = data.iter)
      
    } else {
      ##if constant, remove factor.id from model
      mod.list[[k]] <- update(mod, as.formula(paste(". ~ . -", factor.id)), data = data.iter)
    }
    
    
  }

}


  ##if group label should be letters
  if(letter.labels) {
    letter.mat <- matrix(data = letters[pat.mat], ncol = ncol(pat.mat))
    pat.mat <- letter.mat
  }
  
  ##create vector of names
  model.names <- apply(X = pat.mat, MARGIN = 1, FUN = function(i) paste(i, collapse = ""))##use contrasts in name
  #1111, 122, etc...
  model.names <- paste("m_", model.names, sep = "")
  
  ##compute AICc table
  out.table <- aictab(cand.set = mod.list, modnames = model.names,
                      sort = sort, c.hat = c.hat,
                      second.ord = second.ord,  nobs = nobs)

  ##arrange output
  results <- list(factor.id = factor.id, models = mod.list, model.names = model.names,
                  model.table = out.table, ordered.levels = ord.groups)
  class(results) <- c("mult.comp", "list")

  return(results)  
}



print.mult.comp <-
  function(x, digits = 2, LL = TRUE, ...) {
    ##extract model table
    parm <- x$factor.id
    ord.groups <- x$ordered.levels
    x <- x$model.table
    cat("\nModel selection for multiple comparisons of \"", parm, "\" based on", colnames(x)[3], ":\n")
    if (any(names(x) == "c_hat")) {cat("(c-hat estimate = ", x$c_hat[1], ")\n")}
    cat("\n")

    #check if Cum.Wt should be printed
    if(any(names(x) == "Cum.Wt")) {
      nice.tab <- cbind(x[, 2], x[, 3], x[, 4], x[, 6], x[, "Cum.Wt"], x[, 7])
      colnames(nice.tab) <- c(colnames(x)[c(2, 3, 4, 6)], "Cum.Wt", colnames(x)[7])
      rownames(nice.tab) <- x[, 1]
    } else {
      nice.tab <- cbind(x[, 2], x[, 3], x[, 4], x[, 6], x[, 7])
      colnames(nice.tab) <- c(colnames(x)[c(2, 3, 4, 6, 7)])
      rownames(nice.tab) <- x[, 1]
    }
    

    #if LL==FALSE
    if(identical(LL, FALSE)) {
      names.cols <- colnames(nice.tab)
      sel.LL <- which(attr(regexpr(pattern = "LL", text = names.cols), "match.length") > 1)
      nice.tab <- nice.tab[, -sel.LL]
    }
    
    print(round(nice.tab, digits = digits)) #select rounding off with digits argument
    cat("\n")
    cat("Order of \"", parm, "\" levels based on increasing means:", ord.groups, "\n")
    cat("\nLabels in model names denote grouping structure\n\n")
  }


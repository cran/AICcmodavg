mult.comp <- function(mod, factor.id, letter.labels = TRUE, sort = TRUE, c.hat = 1, second.ord = TRUE,
                          nobs = NULL, newdata = NULL, type = "response", uncond.se = "revised", gamdisp = NULL) {
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

      ##problem with merMod objects: factor becomes character in data set
      ##convert character to factor with default treatment contrasts
      #if(!is.factor(data.set[, factor.id])) {
      #  data.set[, factor.id] <- as.factor(data.set[, factor.id])
      #}
      ##if unmarked object
      ##    if(identical(attr(class(mod), "package"), "unmarked")) stop("\nModels of unmarked classes are not yet supported\n")
                                        #{data.set <- eval(mod@data, (environment(mod@formula)))}
    }
  }

  ##check for presence of interactions
  form.check <- formula(mod)
  form.mod <- strsplit(as.character(form.check), split="~")[[3]]
  if(attr(regexpr(paste(":", factor.id, sep = ""), form.mod), "match.length") != -1 || attr(regexpr(paste(factor.id, ":", sep = ""),
           form.mod), "match.length") != -1 ||  attr(regexpr(paste("\\*", factor.id), form.mod), "match.length") != -1 ||
     attr(regexpr(paste(factor.id, "\\*"), form.mod), "match.length") != -1 ) {
    stop("\nDo not involve the factor interest in interaction terms with this function,\n",
            "see \"?mult.comp\" for details on specifying interaction terms\n")
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

###################CHANGES####
##############################
  if(identical(class(mod), c("glm", "lm")))  {##check if model uses gamma distribution
    ##check family of glm to avoid problems when requesting predictions with argument 'dispersion'
    fam.type <- family(mod)$family
    if(identical(fam.type, "gaussian")) {
      dispersion <- NULL  #set to NULL if gaussian is used
    } else{dispersion <- c.hat}
    
    if(c.hat > 1) {dispersion <- c.hat }
    if(!is.null(gamdisp)) {dispersion <- gamdisp}
    if(c.hat > 1 && !is.null(gamdisp)) {stop("\nYou cannot specify values for both \'c.hat\' and \'gamdisp\'\n")}

    ##correct SE's for estimates of gamma regressions when gamdisp is specified
    if(identical(family(mod)$family[1], "Gamma"))  {
      ##check for specification of gamdisp argument
      if(is.null(gamdisp)) stop("\nYou must specify a gamma dispersion parameter with gamma generalized linear models\n")
    }
  }
  
  
  ##newdata is data frame with exact structure of the original data frame (same variable names and type)
  if(type == "terms") {stop("\nThe terms argument is not defined for this function\n")}

 
  ##check data frame for predictions
  if(is.null(newdata)) {
    ##if no data set is specified, use first observation in data set
    new.data <- data.set[1, ]
    preds.data <- new.data
    ##replicate first line n.groups time
    for(k in 1:(n.groups-1)){
      preds.data <- rbind(preds.data, new.data)
    }
    
    ##add ordered groups
    preds.data[, factor.id] <- as.factor(ord.groups)
    pred.parm <-  preds.data[, factor.id]
    
  } else {
    preds.data <- newdata

    ##check that no more rows are specified than group levels
    if(nrow(newdata) != n.groups) stop("\nnumber of rows in \'newdata\' must match number of factor levels\n")
    preds.data[, factor.id] <- as.factor(ord.groups)
    pred.parm <-  preds.data[, factor.id]
  }
    


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

  ##create list to store preds
  pred.list <- vector("list", length = n.mods)

  if(n.groups == 2) {
    ##loop over matrix and change parameterization
    for(k in 1:nrow(pat.mat)) {
      orig.parm <- data.set[, factor.id]
      ##for model
      new.parm <- ifelse(orig.parm == ord.groups[1], pat.mat[k, 1], pat.mat[k, 2])
      ##for predictions
      new.pred.parm <- ifelse(pred.parm == ord.groups[1], pat.mat[k, 1], pat.mat[k, 2])
      ##replace in preds.data
      preds.data[, factor.id] <- as.factor(new.pred.parm)
  
      ##convert to factor only for cases with different groups
      if(length(unique(new.parm)) > 1) {
        data.iter[, factor.id] <- as.factor(new.parm)
      
        ##run model on updated data
        mod.list[[k]] <- update(mod, data = data.iter)
 
        ##compute predictions
        ##determine if lm or rlm
        if(identical(class(mod), "lm") || identical(class(mod), c("rlm", "lm"))){
          pred.list[[k]] <- predict(mod.list[[k]], newdata = preds.data, se.fit = TRUE)
        }   

        ##determine if glm
        if(identical(class(mod), c("glm", "lm"))) {
          pred.list[[k]] <- predict(mod.list[[k]], type = type, newdata = preds.data, dispersion = dispersion,
                                    se.fit = TRUE)
        }   
      
        ##determine if lme
        if(identical(class(mod), "lme"))  {
          pred.list[[k]] <- predictSE.lme(mod.list[[k]], newdata = preds.data, se.fit = TRUE)
        }

        ##determine if gls
        if(identical(class(mod), "gls"))  {
          pred.list[[k]] <- predictSE.gls(mod.list[[k]], newdata = preds.data, se.fit = TRUE)
        }

        ##determine if mer or merMod
        if(identical(class(mod)[1], "mer") || identical(check.class[1], "lmerMod") || identical(check.class[1], "glmerMod")) {
          pred.list[[k]] <- predictSE.mer(mod.list[[k]], newdata = preds.data, se.fit = TRUE,
                                          type = type)
        }

       
      } else {
        ##if constant, remove factor.id from model
        mod.list[[k]] <- update(mod, as.formula(paste(". ~ . -", factor.id)), data = data.iter)
      
        ##compute predictions
        ##determine if lm or rlm
        if(identical(class(mod), "lm") || identical(class(mod), c("rlm", "lm"))){
          pred.list[[k]] <- predict(mod.list[[k]], newdata = preds.data, se.fit = TRUE)
        }   

        ##determine if glm
        if(identical(class(mod), c("glm", "lm"))) {
          pred.list[[k]] <- predict(mod.list[[k]], type = type, newdata = preds.data, dispersion = dispersion,
                                    se.fit = TRUE)
        }   
      
        ##determine if lme
        if(identical(class(mod), "lme"))  {
          pred.list[[k]] <- predictSE.lme(mod.list[[k]], newdata = preds.data, se.fit = TRUE)
        }

        ##determine if gls
        if(identical(class(mod), "gls"))  {
          pred.list[[k]] <- predictSE.gls(mod.list[[k]], newdata = preds.data, se.fit = TRUE)
        }
        
        ##determine if mer or merMod
        if(identical(class(mod)[1], "mer") || identical(check.class[1], "lmerMod") || identical(check.class[1], "glmerMod")) {
          pred.list[[k]] <- predictSE.mer(mod.list[[k]], newdata = preds.data, se.fit = TRUE,
                                          type = type)
        }
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
      
      ##for predictions
      new.pred.parm <- ifelse(pred.parm == ord.groups[1], pat.mat[k, 1],
                              ifelse(pred.parm == ord.groups[2], pat.mat[k, 2],
                                     pat.mat[k, 3]))
      ##replace in preds.data
      preds.data[, factor.id] <- as.factor(new.pred.parm)

      ##convert to factor only for cases with different groups
      if(length(unique(new.parm)) > 1) {
        data.iter[, factor.id] <- as.factor(new.parm)
        
        ##run model on updated data
        mod.list[[k]] <- update(mod, data = data.iter)
        ##compute predictions
        ##determine if lm or rlm
        if(identical(class(mod), "lm") || identical(class(mod), c("rlm", "lm"))){
          pred.list[[k]] <- predict(mod.list[[k]], newdata = preds.data, se.fit = TRUE)
        }   

        ##determine if glm
        if(identical(class(mod), c("glm", "lm"))) {
          pred.list[[k]] <- predict(mod.list[[k]], type = type, newdata = preds.data, dispersion = dispersion,
                                    se.fit = TRUE)
        }   

      
        ##determine if lme
        if(identical(class(mod), "lme"))  {
          pred.list[[k]] <- predictSE.lme(mod.list[[k]], newdata = preds.data, se.fit = TRUE)
        }

        ##determine if gls
        if(identical(class(mod), "gls"))  {
          pred.list[[k]] <- predictSE.gls(mod.list[[k]], newdata = preds.data, se.fit = TRUE)
        }

        ##determine if mer or merMod
        if(identical(class(mod)[1], "mer") || identical(check.class[1], "lmerMod") || identical(check.class[1], "glmerMod")) {
          pred.list[[k]] <- predictSE.mer(mod.list[[k]], newdata = preds.data, se.fit = TRUE,
                                          type = type)
        }
       
      } else {
        ##if constant, remove factor.id from model
        mod.list[[k]] <- update(mod, as.formula(paste(". ~ . -", factor.id)), data = data.iter)
        
        ##compute predictions
        ##determine if lm or rlm
        if(identical(class(mod), "lm") || identical(class(mod), c("rlm", "lm"))){
          pred.list[[k]] <- predict(mod.list[[k]], newdata = preds.data, se.fit = TRUE)
        }   
        
        ##determine if glm
        if(identical(class(mod), c("glm", "lm"))) {
          pred.list[[k]] <- predict(mod.list[[k]], type = type, newdata = preds.data, dispersion = dispersion,
                                    se.fit = TRUE)
        }   

        ##determine if lme
        if(identical(class(mod), "lme"))  {
          pred.list[[k]] <- predictSE.lme(mod.list[[k]], newdata = preds.data, se.fit = TRUE)
        }
        
        ##determine if gls
        if(identical(class(mod), "gls"))  {
          pred.list[[k]] <- predictSE.gls(mod.list[[k]], newdata = preds.data, se.fit = TRUE)
        }

        ##determine if mer or merMod
        if(identical(class(mod)[1], "mer") || identical(check.class[1], "lmerMod") || identical(check.class[1], "glmerMod")) {
          pred.list[[k]] <- predictSE.mer(mod.list[[k]], newdata = preds.data, se.fit = TRUE,
                                          type = type)
        }
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
      
      ##for predictions
      new.pred.parm <- ifelse(pred.parm == ord.groups[1], pat.mat[k, 1],
                              ifelse(pred.parm == ord.groups[2], pat.mat[k, 2],
                                     ifelse(pred.parm == ord.groups[3], pat.mat[k, 3],
                                            pat.mat[k, 4])))
      ##replace in preds.data
      preds.data[, factor.id] <- as.factor(new.pred.parm)
      ##    data.iter[, factor.id] <- new.parm
      
      ##convert to factor only for cases with different groups
      if(length(unique(new.parm)) > 1) {
        data.iter[, factor.id] <- as.factor(new.parm)
          
        ##run model on updated data
        mod.list[[k]] <- update(mod, data = data.iter)
        
        ##compute predictions
        ##determine if lm or rlm
        if(identical(class(mod), "lm") || identical(class(mod), c("rlm", "lm"))){
          pred.list[[k]] <- predict(mod.list[[k]], newdata = preds.data, se.fit = TRUE)
        }   

        ##determine if glm
        if(identical(class(mod), c("glm", "lm"))) {
          pred.list[[k]] <- predict(mod.list[[k]], type = type, newdata = preds.data, dispersion = dispersion,
                                    se.fit = TRUE)
        }   
   
        ##determine if lme
        if(identical(class(mod), "lme"))  {
          pred.list[[k]] <- predictSE.lme(mod.list[[k]], newdata = preds.data, se.fit = TRUE)
        }
        
        ##determine if gls
        if(identical(class(mod), "gls"))  {
          pred.list[[k]] <- predictSE.gls(mod.list[[k]], newdata = preds.data, se.fit = TRUE)
        }

        ##determine if mer or merMod
        if(identical(class(mod)[1], "mer") || identical(check.class[1], "lmerMod") || identical(check.class[1], "glmerMod")) {
          pred.list[[k]] <- predictSE.mer(mod.list[[k]], newdata = preds.data, se.fit = TRUE,
                                          type = type)
        }

      } else {
        ##if constant, remove factor.id from model
        mod.list[[k]] <- update(mod, as.formula(paste(". ~ . -", factor.id)), data = data.iter)
        ##compute predictions
        ##determine if lm or rlm
        if(identical(class(mod), "lm") || identical(class(mod), c("rlm", "lm"))){
          pred.list[[k]] <- predict(mod.list[[k]], newdata = preds.data, se.fit = TRUE)
        }   

        ##determine if glm
        if(identical(class(mod), c("glm", "lm"))) {
          pred.list[[k]] <- predict(mod.list[[k]], type = type, newdata = preds.data, dispersion = dispersion,
                                    se.fit = TRUE)
        }   

        ##determine if lme
        if(identical(class(mod), "lme"))  {
          pred.list[[k]] <- predictSE.lme(mod.list[[k]], newdata = preds.data, se.fit = TRUE)
        }

        ##determine if gls
        if(identical(class(mod), "gls"))  {
          pred.list[[k]] <- predictSE.gls(mod.list[[k]], newdata = preds.data, se.fit = TRUE)
        }

        ##determine if mer or merMod
        if(identical(class(mod)[1], "mer") || identical(check.class[1], "lmerMod") || identical(check.class[1], "glmerMod")) {
          pred.list[[k]] <- predictSE.mer(mod.list[[k]], newdata = preds.data, se.fit = TRUE,
                                          type = type)
        }
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

      ##for predictions
      new.pred.parm <- ifelse(pred.parm == ord.groups[1], pat.mat[k, 1],
                              ifelse(pred.parm == ord.groups[2], pat.mat[k, 2],
                                     ifelse(pred.parm == ord.groups[3], pat.mat[k, 3],
                                            ifelse(pred.parm == ord.groups[4], pat.mat[k, 4],
                                                   pat.mat[k, 5]))))
      ##replace in preds.data
      preds.data[, factor.id] <- as.factor(new.pred.parm)
      ##    data.iter[, factor.id] <- new.parm
    
    
      ##convert to factor only for cases with different groups
      if(length(unique(new.parm)) > 1) {
        data.iter[, factor.id] <- as.factor(new.parm)
        
        ##run model on updated data
        mod.list[[k]] <- update(mod, data = data.iter)

        ##compute predictions
        ##determine if lm or rlm
        if(identical(class(mod), "lm") || identical(class(mod), c("rlm", "lm"))){
          pred.list[[k]] <- predict(mod.list[[k]], newdata = preds.data, se.fit = TRUE)
        }   

        ##determine if glm
        if(identical(class(mod), c("glm", "lm"))) {
          pred.list[[k]] <- predict(mod.list[[k]], type = type, newdata = preds.data, dispersion = dispersion,
                                    se.fit = TRUE)
        }   

      
        ##determine if lme
        if(identical(class(mod), "lme"))  {
          pred.list[[k]] <- predictSE.lme(mod.list[[k]], newdata = preds.data, se.fit = TRUE)
        }

        ##determine if gls
        if(identical(class(mod), "gls"))  {
          pred.list[[k]] <- predictSE.gls(mod.list[[k]], newdata = preds.data, se.fit = TRUE)
        }

        ##determine if mer or merMod
        if(identical(class(mod)[1], "mer") || identical(check.class[1], "lmerMod") || identical(check.class[1], "glmerMod")) {
          pred.list[[k]] <- predictSE.mer(mod.list[[k]], newdata = preds.data, se.fit = TRUE,
                                          type = type)
        }

      } else {
        ##if constant, remove factor.id from model
        mod.list[[k]] <- update(mod, as.formula(paste(". ~ . -", factor.id)), data = data.iter)
        
        ##compute predictions
        ##determine if lm or rlm
        if(identical(class(mod), "lm") || identical(class(mod), c("rlm", "lm"))){
          pred.list[[k]] <- predict(mod.list[[k]], newdata = preds.data, se.fit = TRUE)
        }   

        ##determine if glm
        if(identical(class(mod), c("glm", "lm"))) {
          pred.list[[k]] <- predict(mod.list[[k]], type = type, newdata = preds.data, dispersion = dispersion,
                                    se.fit = TRUE)
        }

        ##determine if lme
        if(identical(class(mod), "lme"))  {
          pred.list[[k]] <- predictSE.lme(mod.list[[k]], newdata = preds.data, se.fit = TRUE)
        }

        ##determine if gls
        if(identical(class(mod), "gls"))  {
          pred.list[[k]] <- predictSE.gls(mod.list[[k]], newdata = preds.data, se.fit = TRUE)
        }

        ##determine if mer or merMod
        if(identical(class(mod)[1], "mer") || identical(check.class[1], "lmerMod") || identical(check.class[1], "glmerMod")) {
          pred.list[[k]] <- predictSE.mer(mod.list[[k]], newdata = preds.data, se.fit = TRUE,
                                          type = type)
        }
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
      
      ##for predictions
      new.pred.parm <- ifelse(pred.parm == ord.groups[1], pat.mat[k, 1],
                              ifelse(pred.parm == ord.groups[2], pat.mat[k, 2],
                                     ifelse(pred.parm == ord.groups[3], pat.mat[k, 3],
                                            ifelse(pred.parm == ord.groups[4], pat.mat[k, 4],
                                                   pat.mat[k, 5]))))
      ##replace in preds.data
      preds.data[, factor.id] <- as.factor(new.pred.parm)
      ##    data.iter[, factor.id] <- new.parm
    
      ##convert to factor only for cases with different groups
      if(length(unique(new.parm)) > 1) {
        data.iter[, factor.id] <- as.factor(new.parm)
      
        ##run model on updated data
        mod.list[[k]] <- update(mod, data = data.iter)
        
        ##compute predictions
        ##determine if lm or rlm
        if(identical(class(mod), "lm") || identical(class(mod), c("rlm", "lm"))){
          pred.list[[k]] <- predict(mod.list[[k]], newdata = preds.data, se.fit = TRUE)
        }   

        ##determine if glm
        if(identical(class(mod), c("glm", "lm"))) {
          pred.list[[k]] <- predict(mod.list[[k]], type = type, newdata = preds.data, dispersion = dispersion,
                                    se.fit = TRUE)
        }   

        ##determine if lme
        if(identical(class(mod), "lme"))  {
          pred.list[[k]] <- predictSE.lme(mod.list[[k]], newdata = preds.data, se.fit = TRUE)
        }

        ##determine if gls
        if(identical(class(mod), "gls"))  {
          pred.list[[k]] <- predictSE.gls(mod.list[[k]], newdata = preds.data, se.fit = TRUE)
        }

        ##determine if mer or merMod
        if(identical(class(mod)[1], "mer") || identical(check.class[1], "lmerMod") || identical(check.class[1], "glmerMod")) {
          pred.list[[k]] <- predictSE.mer(mod.list[[k]], newdata = preds.data, se.fit = TRUE,
                                          type = type)
        }

      } else {
        ##if constant, remove factor.id from model
        mod.list[[k]] <- update(mod, as.formula(paste(". ~ . -", factor.id)), data = data.iter)
        
        ##compute predictions
        ##determine if lm or rlm
        if(identical(class(mod), "lm") || identical(class(mod), c("rlm", "lm"))){
          pred.list[[k]] <- predict(mod.list[[k]], newdata = preds.data, se.fit = TRUE)
        }   

        ##determine if glm
        if(identical(class(mod), c("glm", "lm"))) {
          pred.list[[k]] <- predict(mod.list[[k]], type = type, newdata = preds.data, dispersion = dispersion,
                                    se.fit = TRUE)
        }   

        ##determine if lme
        if(identical(class(mod), "lme"))  {
          pred.list[[k]] <- predictSE.lme(mod.list[[k]], newdata = preds.data, se.fit = TRUE)
        }
        
        ##determine if gls
        if(identical(class(mod), "gls"))  {
          pred.list[[k]] <- predictSE.gls(mod.list[[k]], newdata = preds.data, se.fit = TRUE)
        }

        ##determine if mer or merMod
        if(identical(class(mod)[1], "mer") || identical(check.class[1], "lmerMod") || identical(check.class[1], "glmerMod")) {
          pred.list[[k]] <- predictSE.mer(mod.list[[k]], newdata = preds.data, se.fit = TRUE,
                                          type = type)
        }
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
  ##1111, 122, etc...
  model.names <- paste("m_", model.names, sep = "")

  
  ##compute AICc table for output
  out.table <- aictab(cand.set = mod.list, modnames = model.names,
                      sort = sort, c.hat = c.hat,
                      second.ord = second.ord,  nobs = nobs)

  out.table.avg <- aictab(cand.set = mod.list, modnames = model.names,
                          sort = FALSE, c.hat = c.hat,
                          second.ord = second.ord,  nobs = nobs)
  
  Mod.avg.out <- matrix(NA, nrow = nrow(preds.data), ncol = 2)
  colnames(Mod.avg.out) <- c("Mod.avg.est", "Uncond.SE")
  rownames(Mod.avg.out) <- ord.groups
  
  ##iterate over predictions
  for(m in 1:nrow(preds.data)){

    ##add fitted values in table
    out.table.avg$fit <- unlist(lapply(X = pred.list, FUN = function(i) i$fit[m]))
    out.table.avg$SE <- unlist(lapply(X = pred.list, FUN = function(i) i$se.fit[m]))

    ##model-averaged estimate
    mod.avg.est <- out.table.avg[, 6] %*% out.table.avg$fit
    Mod.avg.out[m, 1] <- mod.avg.est

    ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
    if(identical(uncond.se, "old")) {
      Mod.avg.out[m, 2] <- sum(out.table.avg[, 6] * sqrt(out.table.avg$SE^2 + (out.table.avg$fit - mod.avg.est)^2))
    }

    ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
    if(identical(uncond.se, "revised")) {
      Mod.avg.out[m, 2] <- sqrt(sum(out.table.avg[, 6] * (out.table.avg$SE^2 + (out.table.avg$fit - mod.avg.est)^2)))
    }
  }
  
#################################################
###NEW CODE
#################################################

##compute model-average estimates

#################################################
###NEW CODE
#################################################

  ##arrange output
  results <- list(factor.id = factor.id, models = mod.list, model.names = model.names,
                  model.table = out.table, ordered.levels = ord.groups, model.avg.est = Mod.avg.out)
  class(results) <- c("mult.comp", "list")
  
  return(results)  
}


print.mult.comp <-
  function(x, digits = 2, LL = TRUE, ...) {
    ##extract model table
    parm <- x$factor.id
    ord.groups <- x$ordered.levels
    mod.avg.out <- x$model.avg.est
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
    cat("\nLabels in model names denote grouping structure and\n")
    cat("are ordered based on increasing means:", ord.groups, "\n")
    cat("\nModel-averaged estimates of group means:", "\n")
    print(round(mod.avg.out, digits = digits))
    cat("\n")  
  }

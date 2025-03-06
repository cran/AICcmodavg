##extracting data characteristics for subsequent prediction typically with modavgPred or modavgEffect

##generic
extractX <- function(cand.set, ...){
    cand.set <- formatCands(cand.set)
    UseMethod("extractX", cand.set)
}



##default
extractX.default <- function(cand.set, ...){
  stop("\nFunction not yet defined for this object class\n")
}



##aov
extractX.AICaov.lm <- function(cand.set, ...) {
    
    ##extract predictors from list
    form.list <- as.character(lapply(cand.set, FUN = function(x) formula(x)[[3]]))
    ##extract based on "+"
    form.noplus <- unlist(sapply(form.list, FUN = function(i) strsplit(i, split = "\\+")))
    ##remove extra white space
    form.clean <- gsub("(^ +)|( +$)", "", form.noplus)
    unique.clean <- unique(form.clean)
    ##exclude empty strings and intercept
    unique.predictors <- unique.clean[nchar(unique.clean) != 0 & unique.clean != "1"]

    
    ##extract response
    resp <- unique(as.character(sapply(cand.set, FUN = function(x) formula(x)[[2]])))
    ##check if different response used
    if(length(resp) > 1) stop("\nThe response variable should be identical in all models\n")
    
    
    ##extract data from model objects
    dsets <- lapply(cand.set, FUN = function(i) i$model)

    ##remove model names from list
    names(dsets) <- NULL

    ##combine data sets
    combo <- do.call(what = "cbind", dsets)
    dframe <- combo[, unique(names(combo))]

    ##remove response from data frame
    dframe <- dframe[, names(dframe) != resp]

    
    ##assemble results
    result <- list("predictors" = unique.predictors,
                   "data" = dframe)
    class(result) <- "extractX"
    return(result)
}



##glm
extractX.AICglm.lm <- function(cand.set, ...) {
    
    ##extract predictors from list
    form.list <- as.character(lapply(cand.set, FUN = function(x) formula(x)[[3]]))
    ##extract based on "+"
    form.noplus <- unlist(sapply(form.list, FUN = function(i) strsplit(i, split = "\\+")))
    ##remove extra white space
    form.clean <- gsub("(^ +)|( +$)", "", form.noplus)
    unique.clean <- unique(form.clean)
    ##exclude empty strings and intercept
    unique.predictors <- unique.clean[nchar(unique.clean) != 0 & unique.clean != "1"]

    
    ##    ##extract response
    resp <- unique(as.character(sapply(cand.set, FUN = function(x) formula(x)[[2]])))
    ##    ##check if different response used
    if(length(resp) > 1) stop("\nThe response variable should be identical in all models\n")


    ##extract data from model objects
    dsets <- lapply(cand.set, FUN = function(i) i$model)

    ##remove model names from list
    names(dsets) <- NULL

    ##combine data sets
    combo <- do.call(what = "cbind", dsets)
    dframe <- combo[, unique(names(combo))]

    ##remove response from data frame
    dframe <- dframe[, names(dframe) != resp]

    
    ##assemble results
    result <- list("predictors" = unique.predictors,
                   "data" = dframe)
    class(result) <- "extractX"
    return(result)
}



##glmmTMB
extractX.AICglmmTMB <- function(cand.set, ...) {
    
    ##extract predictors from list
    form.list <- as.character(lapply(cand.set, FUN = function(x) formula(x)[[3]]))
    ##extract based on "+"
    form.noplus <- unlist(sapply(form.list, FUN = function(i) strsplit(i, split = "\\+")))
    ##remove extra white space
    form.clean <- gsub("(^ +)|( +$)", "", form.noplus)
    unique.clean <- unique(form.clean)
    ##exclude empty strings and intercept
    unique.predictors <- unique.clean[nchar(unique.clean) != 0 & unique.clean != "1"]

    ##extract response
    resp <- unique(as.character(sapply(cand.set, FUN = function(x) formula(x)[[2]])))
    ##check if different response used
    if(length(resp) > 1) stop("\nThe response variable should be identical in all models\n")
    
    
    ##check for | in variance terms
    pipe.id <- which(regexpr("\\|", unique.predictors) != -1)

    ##remove variance terms from string of predictors
    if(length(pipe.id) > 0) {unique.predictors <- unique.predictors[-pipe.id]}
    
    
    ##extract data from model objects
    dsets <- lapply(cand.set, FUN = function(i) (i$frame))

    ##remove model names from list
    names(dsets) <- NULL

    ##combine data sets
    combo <- do.call(what = "cbind", dsets)
    dframe <- combo[, unique(names(combo))]

    ##remove response from data frame
    dframe <- dframe[, names(dframe) != resp]


    ##check for interactions specified with *
    inter.star <- any(regexpr("\\*", unique.predictors) != -1)
    
    ##check for interaction terms
    inter.id <- any(regexpr("\\:", unique.predictors) != -1)
    
    ##inter.star and inter.id
    if(inter.star && inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.nostar.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(terms.nostar.clean, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }


    ##inter.star
    if(inter.star && !inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
    }
    
        
    ##inter.id
    if(!inter.star && inter.id) {
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }
        

    ##none
    if(!inter.star && !inter.id) {
        ##remove extra white space
        terms.clean <- unique.predictors
    }

        
    ##combine in single character vector
    final.predictors <- unique(terms.clean)
        
    ##check for I( ) custom variables in formula
    I.id <- which(regexpr("I\\(", final.predictors) != -1)
        
    ##if I( ) used
    if(length(I.id) > 0) {
        dframe <- dframe[, final.predictors[-I.id], drop = FALSE]
    } else {
        dframe <- dframe[, final.predictors, drop = FALSE]            
    }

    
    ##assemble results
    result <- list("predictors" = unique.predictors,
                   "data" = dframe)
    class(result) <- "extractX"
    return(result)
}



##gls
extractX.AICgls <- function(cand.set, ...) {
    
    ##extract predictors from list
    form.list <- as.character(lapply(cand.set, FUN = function(x) formula(x)[[3]]))
    ##extract based on "+"
    form.noplus <- unlist(sapply(form.list, FUN = function(i) strsplit(i, split = "\\+")))
    ##remove extra white space
    form.clean <- gsub("(^ +)|( +$)", "", form.noplus)
    unique.clean <- unique(form.clean)
    ##exclude empty strings and intercept
    unique.predictors <- unique.clean[nchar(unique.clean) != 0 & unique.clean != "1"]


    ##extract response
    resp <- unique(as.character(sapply(cand.set, FUN = function(x) formula(x)[[2]])))
    ##check if different response used
    if(length(resp) > 1) stop("\nThe response variable should be identical in all models\n")
    
    
    ##extract data from model objects
    dsets <- lapply(cand.set, FUN = function(i) getData(i))
    
    ##remove model names from list
    names(dsets) <- NULL
    
    ##combine data sets
    combo <- do.call(what = "cbind", dsets)
    dframe <- combo[, unique(names(combo))]

    ##remove response from data frame
    dframe <- dframe[, names(dframe) != resp]

    ##check for interactions specified with *
    inter.star <- any(regexpr("\\*", unique.predictors) != -1)
    
    ##check for interaction terms
    inter.id <- any(regexpr("\\:", unique.predictors) != -1)
    
    ##inter.star and inter.id
    if(inter.star && inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.nostar.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(terms.nostar.clean, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }


    ##inter.star
    if(inter.star && !inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
    }
    
        
    ##inter.id
    if(!inter.star && inter.id) {
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }
        

    ##none
    if(!inter.star && !inter.id) {
        ##remove extra white space
        terms.clean <- unique.predictors
    }

        
    ##combine in single character vector
    final.predictors <- unique(terms.clean)
        
    ##check for I( ) custom variables in formula
    I.id <- which(regexpr("I\\(", final.predictors) != -1)
        
    ##if I( ) used
    if(length(I.id) > 0) {
        dframe <- dframe[, final.predictors[-I.id], drop = FALSE]
    } else {
        dframe <- dframe[, final.predictors, drop = FALSE]            
    }

    
    ##assemble results
    result <- list("predictors" = unique.predictors,
                   "data" = dframe)
    class(result) <- "extractX"
    return(result)
}



##lm
extractX.AIClm <- function(cand.set, ...) {
    
    ##extract predictors from list
    form.list <- as.character(lapply(cand.set, FUN = function(x) formula(x)[[3]]))
    ##extract based on "+"
    form.noplus <- unlist(sapply(form.list, FUN = function(i) strsplit(i, split = "\\+")))
    ##remove extra white space
    form.clean <- gsub("(^ +)|( +$)", "", form.noplus)
    unique.clean <- unique(form.clean)
    ##exclude empty strings and intercept
    unique.predictors <- unique.clean[nchar(unique.clean) != 0 & unique.clean != "1"]
    

    ##extract response
    resp <- unique(as.character(sapply(cand.set, FUN = function(x) formula(x)[[2]])))
    ##check if different response used
    if(length(resp) > 1) stop("\nThe response variable should be identical in all models\n")


    ##extract data from model objects
    dsets <- lapply(cand.set, FUN = function(i) i$model)

    ##remove model names from list
    names(dsets) <- NULL

    ##combine data sets
    combo <- do.call(what = "cbind", dsets)
    dframe <- combo[, unique(names(combo))]

    ##remove response from data frame
    dframe <- dframe[, names(dframe) != resp]

    
    ##assemble results
    result <- list("predictors" = unique.predictors,
                   "data" = dframe)
    class(result) <- "extractX"
    return(result)
}



##lme
extractX.AIClme <- function(cand.set, ...) {
    
    ##extract predictors from list
    form.list <- as.character(lapply(cand.set, FUN = function(x) formula(x)[[3]]))
    ##extract based on "+"
    form.noplus <- unlist(sapply(form.list, FUN = function(i) strsplit(i, split = "\\+")))
    ##remove extra white space
    form.clean <- gsub("(^ +)|( +$)", "", form.noplus)
    unique.clean <- unique(form.clean)
    ##exclude empty strings and intercept
    unique.predictors <- unique.clean[nchar(unique.clean) != 0 & unique.clean != "1"]

    
    ##extract response
    resp <- unique(as.character(sapply(cand.set, FUN = function(x) formula(x)[[2]])))
    ##check if different response used
    if(length(resp) > 1) stop("\nThe response variable should be identical in all models\n")


    ##extract data from model objects
    dsets <- lapply(cand.set, FUN = function(i) getData(i))

    ##remove model names from list
    names(dsets) <- NULL

    ##combine data sets
    combo <- do.call(what = "cbind", dsets)
    dframe <- combo[, unique(names(combo))]

    ##remove response from data frame
    dframe <- dframe[, names(dframe) != resp]


    ##check for interactions specified with *
    inter.star <- any(regexpr("\\*", unique.predictors) != -1)
    
    ##check for interaction terms
    inter.id <- any(regexpr("\\:", unique.predictors) != -1)
    
    ##inter.star and inter.id
    if(inter.star && inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.nostar.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(terms.nostar.clean, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }


    ##inter.star
    if(inter.star && !inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
    }
    
        
    ##inter.id
    if(!inter.star && inter.id) {
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }
        

    ##none
    if(!inter.star && !inter.id) {
        ##remove extra white space
        terms.clean <- unique.predictors
    }

        
    ##combine in single character vector
    final.predictors <- unique(terms.clean)
        
    ##check for I( ) custom variables in formula
    I.id <- which(regexpr("I\\(", final.predictors) != -1)
        
    ##if I( ) used
    if(length(I.id) > 0) {
        dframe <- dframe[, final.predictors[-I.id], drop = FALSE]
    } else {
        dframe <- dframe[, final.predictors, drop = FALSE]            
    }

    
    ##assemble results
    result <- list("predictors" = unique.predictors,
                   "data" = dframe)
    class(result) <- "extractX"
    return(result)
}



##glmerMod
extractX.AICglmerMod <- function(cand.set, ...) {
    
    ##extract predictors from list
    form.list <- as.character(lapply(cand.set, FUN = function(x) formula(x)[[3]]))
    ##extract based on "+"
    form.noplus <- unlist(sapply(form.list, FUN = function(i) strsplit(i, split = "\\+")))
    ##remove extra white space
    form.clean <- gsub("(^ +)|( +$)", "", form.noplus)
    unique.clean <- unique(form.clean)
    ##exclude empty strings and intercept
    unique.predictors <- unique.clean[nchar(unique.clean) != 0 & unique.clean != "1"]

    ##extract response
    resp <- unique(as.character(sapply(cand.set, FUN = function(x) formula(x)[[2]])))
    ##check if different response used
    if(length(resp) > 1) stop("\nThe response variable should be identical in all models\n")
    
    
    ##check for | in variance terms
    pipe.id <- which(regexpr("\\|", unique.predictors) != -1)

    ##remove variance terms from string of predictors
    if(length(pipe.id) > 0) {unique.predictors <- unique.predictors[-pipe.id]}
    
    
    ##extract data from model objects
    dsets <- lapply(cand.set, FUN = function(i) (i@frame))

    ##remove model names from list
    names(dsets) <- NULL

    ##combine data sets
    combo <- do.call(what = "cbind", dsets)
    dframe <- combo[, unique(names(combo))]

    ##remove response from data frame
    dframe <- dframe[, names(dframe) != resp]


    ##check for interactions specified with *
    inter.star <- any(regexpr("\\*", unique.predictors) != -1)
    
    ##check for interaction terms
    inter.id <- any(regexpr("\\:", unique.predictors) != -1)
    
    ##inter.star and inter.id
    if(inter.star && inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.nostar.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(terms.nostar.clean, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }


    ##inter.star
    if(inter.star && !inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
    }
    
        
    ##inter.id
    if(!inter.star && inter.id) {
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }
        

    ##none
    if(!inter.star && !inter.id) {
        ##remove extra white space
        terms.clean <- unique.predictors
    }

        
    ##combine in single character vector
    final.predictors <- unique(terms.clean)
        
    ##check for I( ) custom variables in formula
    I.id <- which(regexpr("I\\(", final.predictors) != -1)
        
    ##if I( ) used
    if(length(I.id) > 0) {
        dframe <- dframe[, final.predictors[-I.id], drop = FALSE]
    } else {
        dframe <- dframe[, final.predictors, drop = FALSE]            
    }

    
    ##assemble results
    result <- list("predictors" = unique.predictors,
                   "data" = dframe)
    class(result) <- "extractX"
    return(result)
}



##lmerMod
extractX.AIClmerMod <- function(cand.set, ...) {
    
    ##extract predictors from list
    form.list <- as.character(lapply(cand.set, FUN = function(x) formula(x)[[3]]))
    ##extract based on "+"
    form.noplus <- unlist(sapply(form.list, FUN = function(i) strsplit(i, split = "\\+")))
    ##remove extra white space
    form.clean <- gsub("(^ +)|( +$)", "", form.noplus)
    unique.clean <- unique(form.clean)
    ##exclude empty strings and intercept
    unique.predictors <- unique.clean[nchar(unique.clean) != 0 & unique.clean != "1"]

    ##extract response
    resp <- unique(as.character(sapply(cand.set, FUN = function(x) formula(x)[[2]])))
    ##check if different response used
    if(length(resp) > 1) stop("\nThe response variable should be identical in all models\n")


    ##check for | in variance terms
    pipe.id <- which(regexpr("\\|", unique.predictors) != -1)

    ##remove variance terms from string of predictors
    if(length(pipe.id) > 0) {unique.predictors <- unique.predictors[-pipe.id]}
    
    
    ##extract data from model objects
    dsets <- lapply(cand.set, FUN = function(i) i@frame)

    ##remove model names from list
    names(dsets) <- NULL

    ##combine data sets
    combo <- do.call(what = "cbind", dsets)
    dframe <- combo[, unique(names(combo))]

    ##remove response from data frame
    dframe <- dframe[, names(dframe) != resp]


    ##check for interactions specified with *
    inter.star <- any(regexpr("\\*", unique.predictors) != -1)
    
    ##check for interaction terms
    inter.id <- any(regexpr("\\:", unique.predictors) != -1)
    
    ##inter.star and inter.id
    if(inter.star && inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.nostar.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(terms.nostar.clean, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }


    ##inter.star
    if(inter.star && !inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
    }
    
        
    ##inter.id
    if(!inter.star && inter.id) {
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }
        

    ##none
    if(!inter.star && !inter.id) {
        ##remove extra white space
        terms.clean <- unique.predictors
    }

        
    ##combine in single character vector
    final.predictors <- unique(terms.clean)
        
    ##check for I( ) custom variables in formula
    I.id <- which(regexpr("I\\(", final.predictors) != -1)
        
    ##if I( ) used
    if(length(I.id) > 0) {
        dframe <- dframe[, final.predictors[-I.id], drop = FALSE]
    } else {
        dframe <- dframe[, final.predictors, drop = FALSE]            
    }

    
    ##assemble results
    result <- list("predictors" = unique.predictors,
                   "data" = dframe)
    class(result) <- "extractX"
    return(result)
}



##lmerModLmerTest
extractX.AIClmerModLmerTest <- function(cand.set, ...) {
    
    ##extract predictors from list
    form.list <- as.character(lapply(cand.set, FUN = function(x) formula(x)[[3]]))
    ##extract based on "+"
    form.noplus <- unlist(sapply(form.list, FUN = function(i) strsplit(i, split = "\\+")))
    ##remove extra white space
    form.clean <- gsub("(^ +)|( +$)", "", form.noplus)
    unique.clean <- unique(form.clean)
    ##exclude empty strings and intercept
    unique.predictors <- unique.clean[nchar(unique.clean) != 0 & unique.clean != "1"]

    ##extract response
    resp <- unique(as.character(sapply(cand.set, FUN = function(x) formula(x)[[2]])))
    ##check if different response used
    if(length(resp) > 1) stop("\nThe response variable should be identical in all models\n")


    ##check for | in variance terms
    pipe.id <- which(regexpr("\\|", unique.predictors) != -1)

    ##remove variance terms from string of predictors
    if(length(pipe.id) > 0) {unique.predictors <- unique.predictors[-pipe.id]}
    
    
    ##extract data from model objects
    dsets <- lapply(cand.set, FUN = function(i) i@frame)

    ##remove model names from list
    names(dsets) <- NULL

    ##combine data sets
    combo <- do.call(what = "cbind", dsets)
    dframe <- combo[, unique(names(combo))]

    ##remove response from data frame
    dframe <- dframe[, names(dframe) != resp]


    ##check for interactions specified with *
    inter.star <- any(regexpr("\\*", unique.predictors) != -1)
    
    ##check for interaction terms
    inter.id <- any(regexpr("\\:", unique.predictors) != -1)
    
    ##inter.star and inter.id
    if(inter.star && inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.nostar.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(terms.nostar.clean, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }


    ##inter.star
    if(inter.star && !inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
    }
    
        
    ##inter.id
    if(!inter.star && inter.id) {
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }
        

    ##none
    if(!inter.star && !inter.id) {
        ##remove extra white space
        terms.clean <- unique.predictors
    }

        
    ##combine in single character vector
    final.predictors <- unique(terms.clean)
        
    ##check for I( ) custom variables in formula
    I.id <- which(regexpr("I\\(", final.predictors) != -1)
        
    ##if I( ) used
    if(length(I.id) > 0) {
        dframe <- dframe[, final.predictors[-I.id], drop = FALSE]
    } else {
        dframe <- dframe[, final.predictors, drop = FALSE]            
    }

    
    ##assemble results
    result <- list("predictors" = unique.predictors,
                   "data" = dframe)
    class(result) <- "extractX"
    return(result)
}



##rlm
extractX.AICrlm.lm <- function(cand.set, ...) {
    
    ##extract predictors from list
    form.list <- as.character(lapply(cand.set, FUN = function(x) formula(x)[[3]]))
    ##extract based on "+"
    form.noplus <- unlist(sapply(form.list, FUN = function(i) strsplit(i, split = "\\+")))
    ##remove extra white space
    form.clean <- gsub("(^ +)|( +$)", "", form.noplus)
    unique.clean <- unique(form.clean)
    ##exclude empty strings and intercept
    unique.predictors <- unique.clean[nchar(unique.clean) != 0 & unique.clean != "1"]
    

    ##extract response
    resp <- unique(as.character(sapply(cand.set, FUN = function(x) formula(x)[[2]])))
    ##check if different response used
    if(length(resp) > 1) stop("\nThe response variable should be identical in all models\n")


    ##extract data from model objects
    dsets <- lapply(cand.set, FUN = function(i) i$model)

    ##remove model names from list
    names(dsets) <- NULL

    ##combine data sets
    combo <- do.call(what = "cbind", dsets)
    dframe <- combo[, unique(names(combo))]

    ##remove response from data frame
    dframe <- dframe[, names(dframe) != resp]

    
    ##assemble results
    result <- list("predictors" = unique.predictors,
                   "data" = dframe)
    class(result) <- "extractX"
    return(result)
}



##survreg
extractX.AICsurvreg <- function(cand.set, ...) {
    
    ##extract predictors from list
    form.list <- as.character(lapply(cand.set, FUN = function(x) formula(x)[[3]]))
    ##extract based on "+"
    form.noplus <- unlist(sapply(form.list, FUN = function(i) strsplit(i, split = "\\+")))
    ##remove extra white space
    form.clean <- gsub("(^ +)|( +$)", "", form.noplus)
    unique.clean <- unique(form.clean)
    ##exclude empty strings and intercept
    unique.predictors <- unique.clean[nchar(unique.clean) != 0 & unique.clean != "1"]
    
    
    ##extract data from model object - identical for each because uses eval(data)
    dsets <- lapply(cand.set, FUN = function(i) eval(i$call$data))

    ##remove model names from list
    names(dsets) <- NULL

    ##combine data sets
    combo <- do.call(what = "cbind", dsets)
    dframe <- combo[, unique(names(combo))]

    ##extract response
    ##resp <- unique(as.character(sapply(cand.set, FUN = function(x) formula(x)[[2]])))
    ##check if different response used
    ##if(length(resp) > 1) stop("\nThe response variable should be identical in all models\n")

    ##check for interactions specified with *
    inter.star <- any(regexpr("\\*", unique.predictors) != -1)
        
    ##check for interaction terms
    inter.id <- any(regexpr("\\:", unique.predictors) != -1)

    ##inter.star and inter.id
    if(inter.star && inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.nostar.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(terms.nostar.clean, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }


    ##inter.star
    if(inter.star && !inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
    }
    
        
    ##inter.id
    if(!inter.star && inter.id) {
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }
        

    ##none
    if(!inter.star && !inter.id) {
        ##remove extra white space
        terms.clean <- unique.predictors
    }

        
    ##combine in single character vector
    final.predictors <- unique(terms.clean)
        

    ##check for I( ) custom variables in formula
    I.id <- which(regexpr("I\\(", final.predictors) != -1)
        
    ##if I( ) used
    if(length(I.id) > 0) {
        dframe <- dframe[, final.predictors[-I.id], drop = FALSE]
    } else {
        dframe <- dframe[, final.predictors, drop = FALSE]            
    }

    
    ##assemble results
    result <- list("predictors" = unique.predictors,
                   "data" = dframe)
    class(result) <- "extractX"
    return(result)
}



##unmarkedFitOccu
extractX.AICunmarkedFitOccu <- function(cand.set, parm.type = NULL, ...) {

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?extractX for details\n")}

    ##extract predictors from list
    ##psi
    if(identical(parm.type, "psi")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formula[[3]]))
    }

    ##detect
    if(identical(parm.type, "detect")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formula[[2]]))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }


    ##extract based on "+"
    form.noplus <- unlist(sapply(form.list, FUN = function(i) strsplit(i, split = "\\+")))
    ##remove extra white space
    form.clean <- gsub("(^ +)|( +$)", "", form.noplus)
    unique.clean <- unique(form.clean)
    ##exclude empty strings and intercept
    unique.predictors <- unique.clean[nchar(unique.clean) != 0 & unique.clean != "1"]
    

    ##extract data from model objects - identical for all models
    dsets <- lapply(cand.set, FUN = function(i) unmarked::getData(i))
    ##check that same data are used
    unique.dsets <- unique(dsets)
    if(length(unique.dsets) != 1) stop("\nData sets differ across models:\n check data carefully\n")
    unFrame <- unique.dsets[[1]]
    
    ##extract siteCovs
    siteVars <- siteCovs(unFrame)
    ##extract obsCovs
    obsVars <- obsCovs(unFrame)

    ##check for interactions specified with *
    inter.star <- any(regexpr("\\*", unique.predictors) != -1)
        
    ##check for interaction terms
    inter.id <- any(regexpr("\\:", unique.predictors) != -1)

    ##inter.star and inter.id
    if(inter.star && inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.nostar.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(terms.nostar.clean, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }


    ##inter.star
    if(inter.star && !inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
    }
    
        
    ##inter.id
    if(!inter.star && inter.id) {
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }
        

    ##none
    if(!inter.star && !inter.id) {
        ##remove extra white space
        terms.clean <- unique.predictors
    }

        
    ##combine in single character vector
    final.predictors <- unique(terms.clean)
        
    ##find where predictors occur
    if(!is.null(obsVars)) {
        obsID <- obsVars[, intersect(final.predictors, names(obsVars)), drop = FALSE]
        if(nrow(obsID) > 0) {
            obsID.info <- capture.output(str(obsID))[-1]
        } else {
            obsID.info <- NULL
        }
    } else {obsID.info <- NULL}
    
    if(!is.null(siteVars)) {
        siteID <- siteVars[, intersect(final.predictors, names(siteVars)), drop = FALSE]
        if(nrow(siteID) > 0) {
            siteID.info <- capture.output(str(siteID))[-1]
        } else {
            siteID.info <- NULL
        }
    } else {siteID.info <- NULL}

    ##store data sets
    data.out <- list( )
    if(is.null(obsVars)) {
        data.out$obsCovs <- NULL
    } else {data.out$obsCovs <- obsID}
    if(is.null(siteVars)) {
        data.out$siteCovs <- NULL
    } else {data.out$siteCovs <- siteID}

    
    ##assemble results
    result <- list("predictors" = unique.predictors,
                   "data" = data.out)
    class(result) <- "extractX"
    return(result)
}



##unmarkedFitColExt
extractX.AICunmarkedFitColExt <- function(cand.set, parm.type = NULL, ...) {

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?extractX for details\n")}

    ##extract predictors from list
    ##psi
    if(identical(parm.type, "psi")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@psiformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }

    ##gamma
    if(identical(parm.type, "gamma")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@gamformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }

    ##epsilon
    if(identical(parm.type, "epsilon")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@epsformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }

    ##detect
    if(identical(parm.type, "detect")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@detformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }


    ##extract based on "+"
    form.noplus <- unlist(sapply(form.list, FUN = function(i) strsplit(i, split = "\\+")))
    ##remove extra white space
    form.clean <- gsub("(^ +)|( +$)", "", form.noplus)
    unique.clean <- unique(form.clean)
    ##exclude empty strings and intercept
    unique.predictors <- unique.clean[nchar(unique.clean) != 0 & unique.clean != "1"]
    
    ##extract data from model objects - identical for all models
    dsets <- lapply(cand.set, FUN = function(i) unmarked::getData(i))
    ##check that same data are used
    unique.dsets <- unique(dsets)
    if(length(unique.dsets) != 1) stop("\nData sets differ across models:\n check data carefully\n")
    unFrame <- unique.dsets[[1]]
    
    ##extract siteCovs
    siteVars <- siteCovs(unFrame)
    ##extract obsCovs
    obsVars <- obsCovs(unFrame)
    ##extract yearlySiteCovs
    yearlyVars <- yearlySiteCovs(unFrame)
    
    ##check for interactions specified with *
    inter.star <- any(regexpr("\\*", unique.predictors) != -1)
        
    ##check for interaction terms
    inter.id <- any(regexpr("\\:", unique.predictors) != -1)

    ##inter.star and inter.id
    if(inter.star && inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.nostar.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(terms.nostar.clean, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }


    ##inter.star
    if(inter.star && !inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
    }
    
        
    ##inter.id
    if(!inter.star && inter.id) {
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }
        

    ##none
    if(!inter.star && !inter.id) {
        ##remove extra white space
        terms.clean <- unique.predictors
    }

        
    ##combine in single character vector
    final.predictors <- unique(terms.clean)
        
    ##find where predictors occur
    if(!is.null(obsVars)) {
        obsID <- obsVars[, intersect(final.predictors, names(obsVars)), drop = FALSE]
        if(nrow(obsID) > 0) {
            obsID.info <- capture.output(str(obsID))[-1]
        } else {
            obsID.info <- NULL
        }
    } else {obsID.info <- NULL}
    
    if(!is.null(siteVars)) {
        siteID <- siteVars[, intersect(final.predictors, names(siteVars)), drop = FALSE]
        if(nrow(siteID) > 0) {
            siteID.info <- capture.output(str(siteID))[-1]
        } else {
            siteID.info <- NULL
        }
    } else {siteID.info <- NULL}
    
    if(!is.null(yearlyVars)) {
        yearlyID <- yearlyVars[, intersect(final.predictors, names(yearlyVars)), drop = FALSE]
        if(nrow(yearlyID) > 0) {
            yearlyID.info <- capture.output(str(yearlyID))[-1]
        } else {
            yearlyID.info <- NULL
        }
    } else {yearlyID.info <- NULL}

    
    ##store data sets
    data.out <- list( )
    if(is.null(obsVars)) {
        data.out$obsCovs <- NULL
    } else {data.out$obsCovs <- obsID}
    if(is.null(siteVars)) {
        data.out$siteCovs <- NULL
    } else {data.out$siteCovs <- siteID}
    if(is.null(yearlyVars)) {
        data.out$yearlySiteCovs <- NULL
    } else {data.out$yearlySiteCovs <- yearlyID}

    
    ##assemble results
    result <- list("predictors" = unique.predictors,
                   "data" = data.out)
    class(result) <- "extractX"
    return(result)
}



##unmarkedFitOccurRN
extractX.AICunmarkedFitOccuRN <- function(cand.set, parm.type = NULL, ...) {

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?extractX for details\n")}

    ##extract predictors from list
    ##psi
    if(identical(parm.type, "psi")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formula[[3]]))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }

    ##detect
    if(identical(parm.type, "detect")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formula[[2]]))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }

    ##extract based on "+"
    form.noplus <- unlist(sapply(form.list, FUN = function(i) strsplit(i, split = "\\+")))
    ##remove extra white space
    form.clean <- gsub("(^ +)|( +$)", "", form.noplus)
    unique.clean <- unique(form.clean)
    ##exclude empty strings and intercept
    unique.predictors <- unique.clean[nchar(unique.clean) != 0 & unique.clean != "1"]
    
    ##extract data from model objects - identical for all models
    dsets <- lapply(cand.set, FUN = function(i) unmarked::getData(i))
    ##check that same data are used
    unique.dsets <- unique(dsets)
    if(length(unique.dsets) != 1) stop("\nData sets differ across models:\n check data carefully\n")
    unFrame <- unique.dsets[[1]]
    
    ##extract siteCovs
    siteVars <- siteCovs(unFrame)
    ##extract obsCovs
    obsVars <- obsCovs(unFrame)
        
    ##check for interactions specified with *
    inter.star <- any(regexpr("\\*", unique.predictors) != -1)
        
    ##check for interaction terms
    inter.id <- any(regexpr("\\:", unique.predictors) != -1)

    ##inter.star and inter.id
    if(inter.star && inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.nostar.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(terms.nostar.clean, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }


    ##inter.star
    if(inter.star && !inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
    }
    
        
    ##inter.id
    if(!inter.star && inter.id) {
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }
        

    ##none
    if(!inter.star && !inter.id) {
        ##remove extra white space
        terms.clean <- unique.predictors
    }

        
    ##combine in single character vector
    final.predictors <- unique(terms.clean)
        
    ##find where predictors occur
    if(!is.null(obsVars)) {
        obsID <- obsVars[, intersect(final.predictors, names(obsVars)), drop = FALSE]
        if(nrow(obsID) > 0) {
            obsID.info <- capture.output(str(obsID))[-1]
        } else {
            obsID.info <- NULL
        }
    } else {obsID.info <- NULL}
    
    if(!is.null(siteVars)) {
        siteID <- siteVars[, intersect(final.predictors, names(siteVars)), drop = FALSE]
        if(nrow(siteID) > 0) {
            siteID.info <- capture.output(str(siteID))[-1]
        } else {
            siteID.info <- NULL
        }
    } else {siteID.info <- NULL}
    

    ##store data sets
    data.out <- list( )
    if(is.null(obsVars)) {
        data.out$obsCovs <- NULL
    } else {data.out$obsCovs <- obsID}
    if(is.null(siteVars)) {
        data.out$siteCovs <- NULL
    } else {data.out$siteCovs <- siteID}
    

    ##assemble results
    result <- list("predictors" = unique.predictors,
                   "data" = data.out)
    class(result) <- "extractX"
    return(result)
}



##unmarkedFitPCO
extractX.AICunmarkedFitPCO <- function(cand.set, parm.type = NULL, ...) {

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?extractX for details\n")}

    ##extract predictors from list
    ##lambda
    if(identical(parm.type, "lambda")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formlist$lambdaformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }

    ##gamma
    if(identical(parm.type, "gamma")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formlist$gammaformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }

    ##omega
    if(identical(parm.type, "omega")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formlist$omegaformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }

    ##iota
    if(identical(parm.type, "iota")) {
        ##check that parameter appears in all models
        parfreq <- sum(sapply(cand.set, FUN = function(i) any(names(i@estimates@estimates) == parm.type)))
        if(!identical(length(cand.set), parfreq)) {
            stop("\nParameter \'iota\' does not appear in all models\n")
        }
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formlist$iotaformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }

    ##detect
    if(identical(parm.type, "detect")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formlist$pformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }


    ##extract based on "+"
    form.noplus <- unlist(sapply(form.list, FUN = function(i) strsplit(i, split = "\\+")))
    ##remove extra white space
    form.clean <- gsub("(^ +)|( +$)", "", form.noplus)
    unique.clean <- unique(form.clean)
    ##exclude empty strings and intercept
    unique.predictors <- unique.clean[nchar(unique.clean) != 0 & unique.clean != "1"]
    
    ##extract data from model objects - identical for all models
    dsets <- lapply(cand.set, FUN = function(i) unmarked::getData(i))
    ##check that same data are used
    unique.dsets <- unique(dsets)
    if(length(unique.dsets) != 1) stop("\nData sets differ across models:\n check data carefully\n")
    unFrame <- unique.dsets[[1]]
    
    ##extract siteCovs
    siteVars <- siteCovs(unFrame)
    ##extract obsCovs
    obsVars <- obsCovs(unFrame)
    ##extract yearlySiteCovs
    yearlyVars <- yearlySiteCovs(unFrame)
    
    ##check for interactions specified with *
    inter.star <- any(regexpr("\\*", unique.predictors) != -1)
        
    ##check for interaction terms
    inter.id <- any(regexpr("\\:", unique.predictors) != -1)

    ##inter.star and inter.id
    if(inter.star && inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.nostar.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(terms.nostar.clean, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }


    ##inter.star
    if(inter.star && !inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
    }
    
        
    ##inter.id
    if(!inter.star && inter.id) {
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }
        

    ##none
    if(!inter.star && !inter.id) {
        ##remove extra white space
        terms.clean <- unique.predictors
    }

        
    ##combine in single character vector
    final.predictors <- unique(terms.clean)
        
    ##find where predictors occur
    if(!is.null(obsVars)) {
        obsID <- obsVars[, intersect(final.predictors, names(obsVars)), drop = FALSE]
        if(nrow(obsID) > 0) {
            obsID.info <- capture.output(str(obsID))[-1]
        } else {
            obsID.info <- NULL
        }
    } else {obsID.info <- NULL}
    
    if(!is.null(siteVars)) {
        siteID <- siteVars[, intersect(final.predictors, names(siteVars)), drop = FALSE]
        if(nrow(siteID) > 0) {
            siteID.info <- capture.output(str(siteID))[-1]
        } else {
            siteID.info <- NULL
        }
    } else {siteID.info <- NULL}
    
    if(!is.null(yearlyVars)) {
        yearlyID <- yearlyVars[, intersect(final.predictors, names(yearlyVars)), drop = FALSE]
        if(nrow(yearlyID) > 0) {
            yearlyID.info <- capture.output(str(yearlyID))[-1]
        } else {
            yearlyID.info <- NULL
        }
    } else {yearlyID.info <- NULL}
    

    ##store data sets
    data.out <- list( )
    if(is.null(obsVars)) {
        data.out$obsCovs <- NULL
    } else {data.out$obsCovs <- obsID}
    if(is.null(siteVars)) {
        data.out$siteCovs <- NULL
    } else {data.out$siteCovs <- siteID}
    if(is.null(yearlyVars)) {
        data.out$yearlySiteCovs <- NULL
    } else {data.out$yearlySiteCovs <- yearlyID}

    
    ##assemble results
    result <- list("predictors" = unique.predictors,
                   "data" = data.out)
    class(result) <- "extractX"
    return(result)
}



##unmarkedFitPCount
extractX.AICunmarkedFitPCount <- function(cand.set, parm.type = NULL, ...) {

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?extractX for details\n")}

    ##extract predictors from list
    ##lambda
    if(identical(parm.type, "lambda")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formula[[3]]))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }

    ##detect
    if(identical(parm.type, "detect")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formula[[2]]))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }


    ##extract based on "+"
    form.noplus <- unlist(sapply(form.list, FUN = function(i) strsplit(i, split = "\\+")))
    ##remove extra white space
    form.clean <- gsub("(^ +)|( +$)", "", form.noplus)
    unique.clean <- unique(form.clean)
    ##exclude empty strings and intercept
    unique.predictors <- unique.clean[nchar(unique.clean) != 0 & unique.clean != "1"]
    
    ##extract data from model objects - identical for all models
    dsets <- lapply(cand.set, FUN = function(i) unmarked::getData(i))
    ##check that same data are used
    unique.dsets <- unique(dsets)
    if(length(unique.dsets) != 1) stop("\nData sets differ across models:\n check data carefully\n")
    unFrame <- unique.dsets[[1]]
    
    ##extract siteCovs
    siteVars <- siteCovs(unFrame)
    ##extract obsCovs
    obsVars <- obsCovs(unFrame)
    
    ##check for interactions specified with *
    inter.star <- any(regexpr("\\*", unique.predictors) != -1)
        
    ##check for interaction terms
    inter.id <- any(regexpr("\\:", unique.predictors) != -1)

    ##inter.star and inter.id
    if(inter.star && inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.nostar.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(terms.nostar.clean, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }


    ##inter.star
    if(inter.star && !inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
    }
    
        
    ##inter.id
    if(!inter.star && inter.id) {
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }
        

    ##none
    if(!inter.star && !inter.id) {
        ##remove extra white space
        terms.clean <- unique.predictors
    }

        
    ##combine in single character vector
    final.predictors <- unique(terms.clean)
        
    ##find where predictors occur
    if(!is.null(obsVars)) {
        obsID <- obsVars[, intersect(final.predictors, names(obsVars)), drop = FALSE]
        if(nrow(obsID) > 0) {
            obsID.info <- capture.output(str(obsID))[-1]
        } else {
            obsID.info <- NULL
        }
    } else {obsID.info <- NULL}
    
    if(!is.null(siteVars)) {
        siteID <- siteVars[, intersect(final.predictors, names(siteVars)), drop = FALSE]
        if(nrow(siteID) > 0) {
            siteID.info <- capture.output(str(siteID))[-1]
        } else {
            siteID.info <- NULL
        }
    } else {siteID.info <- NULL}
    

    ##store data sets
    data.out <- list( )
    if(is.null(obsVars)) {
        data.out$obsCovs <- NULL
    } else {data.out$obsCovs <- obsID}
    if(is.null(siteVars)) {
        data.out$siteCovs <- NULL
    } else {data.out$siteCovs <- siteID}
    

    ##assemble results
    result <- list("predictors" = unique.predictors,
                   "data" = data.out)
    class(result) <- "extractX"
    return(result)
}



##unmarkedFitDS
extractX.AICunmarkedFitDS <- function(cand.set, parm.type = NULL, ...) {

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?extractX for details\n")}

    ##extract predictors from list
    ##lambda
    if(identical(parm.type, "lambda")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formula[[3]]))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }

    ##detect
    if(identical(parm.type, "detect")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formula[[2]]))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }


    ##extract based on "+"
    form.noplus <- unlist(sapply(form.list, FUN = function(i) strsplit(i, split = "\\+")))
    ##remove extra white space
    form.clean <- gsub("(^ +)|( +$)", "", form.noplus)
    unique.clean <- unique(form.clean)
    ##exclude empty strings and intercept
    unique.predictors <- unique.clean[nchar(unique.clean) != 0 & unique.clean != "1"]
    
    ##extract data from model objects - identical for all models
    dsets <- lapply(cand.set, FUN = function(i) unmarked::getData(i))
    ##check that same data are used
    unique.dsets <- unique(dsets)
    if(length(unique.dsets) != 1) stop("\nData sets differ across models:\n check data carefully\n")
    unFrame <- unique.dsets[[1]]
    
    ##extract siteCovs
    siteVars <- siteCovs(unFrame)
    ##extract obsCovs
    obsVars <- obsCovs(unFrame)
        
    ##check for interactions specified with *
    inter.star <- any(regexpr("\\*", unique.predictors) != -1)
        
    ##check for interaction terms
    inter.id <- any(regexpr("\\:", unique.predictors) != -1)

    ##inter.star and inter.id
    if(inter.star && inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.nostar.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(terms.nostar.clean, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }


    ##inter.star
    if(inter.star && !inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
    }
    
        
    ##inter.id
    if(!inter.star && inter.id) {
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }
        

    ##none
    if(!inter.star && !inter.id) {
        ##remove extra white space
        terms.clean <- unique.predictors
    }

        
    ##combine in single character vector
    final.predictors <- unique(terms.clean)
        
    ##find where predictors occur
    if(!is.null(obsVars)) {
        obsID <- obsVars[, intersect(final.predictors, names(obsVars)), drop = FALSE]
        if(nrow(obsID) > 0) {
            obsID.info <- capture.output(str(obsID))[-1]
        } else {
            obsID.info <- NULL
        }
    } else {obsID.info <- NULL}
    
    if(!is.null(siteVars)) {
        siteID <- siteVars[, intersect(final.predictors, names(siteVars)), drop = FALSE]
        if(nrow(siteID) > 0) {
            siteID.info <- capture.output(str(siteID))[-1]
        } else {
            siteID.info <- NULL
        }
    } else {siteID.info <- NULL}
    

    ##store data sets
    data.out <- list( )
    if(is.null(obsVars)) {
        data.out$obsCovs <- NULL
    } else {data.out$obsCovs <- obsID}
    if(is.null(siteVars)) {
        data.out$siteCovs <- NULL
    } else {data.out$siteCovs <- siteID}
    

    ##assemble results
    result <- list("predictors" = unique.predictors,
                   "data" = data.out)
    class(result) <- "extractX"
    return(result)
}



##unmarkedFitGDS
extractX.AICunmarkedFitGDS <- function(cand.set, parm.type = NULL, ...) {

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?extractX for details\n")}

    ##extract predictors from list
    ##lambda
    if(identical(parm.type, "lambda")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formlist$lambdaformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }

    
    ##phi
    if(identical(parm.type, "phi")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formlist$phiformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }


    ##detect
    if(identical(parm.type, "detect")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formlist$pformula[[2]]))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }


    ##extract based on "+"
    form.noplus <- unlist(sapply(form.list, FUN = function(i) strsplit(i, split = "\\+")))
    ##remove extra white space
    form.clean <- gsub("(^ +)|( +$)", "", form.noplus)
    unique.clean <- unique(form.clean)
    ##exclude empty strings and intercept
    unique.predictors <- unique.clean[nchar(unique.clean) != 0 & unique.clean != "1"]
    
    ##extract data from model objects - identical for all models
    dsets <- lapply(cand.set, FUN = function(i) unmarked::getData(i))
    ##check that same data are used
    unique.dsets <- unique(dsets)
    if(length(unique.dsets) != 1) stop("\nData sets differ across models:\n check data carefully\n")
    unFrame <- unique.dsets[[1]]
    
    ##extract siteCovs
    siteVars <- siteCovs(unFrame)
    ##extract obsCovs
    obsVars <- obsCovs(unFrame)
    ##extract yearlySiteCovs
    yearlyVars <- yearlySiteCovs(unFrame)
    
    ##check for interactions specified with *
    inter.star <- any(regexpr("\\*", unique.predictors) != -1)
        
    ##check for interaction terms
    inter.id <- any(regexpr("\\:", unique.predictors) != -1)

    ##inter.star and inter.id
    if(inter.star && inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.nostar.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(terms.nostar.clean, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }


    ##inter.star
    if(inter.star && !inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
    }
    
        
    ##inter.id
    if(!inter.star && inter.id) {
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }
        

    ##none
    if(!inter.star && !inter.id) {
        ##remove extra white space
        terms.clean <- unique.predictors
    }

        
    ##combine in single character vector
    final.predictors <- unique(terms.clean)
        
    ##find where predictors occur
    if(!is.null(obsVars)) {
        obsID <- obsVars[, intersect(final.predictors, names(obsVars)), drop = FALSE]
        if(nrow(obsID) > 0) {
            obsID.info <- capture.output(str(obsID))[-1]
        } else {
            obsID.info <- NULL
        }
    } else {obsID.info <- NULL}
    
    if(!is.null(siteVars)) {
        siteID <- siteVars[, intersect(final.predictors, names(siteVars)), drop = FALSE]
        if(nrow(siteID) > 0) {
            siteID.info <- capture.output(str(siteID))[-1]
        } else {
            siteID.info <- NULL
        }
    } else {siteID.info <- NULL}
    
    if(!is.null(yearlyVars)) {
        yearlyID <- yearlyVars[, intersect(final.predictors, names(yearlyVars)), drop = FALSE]
        if(nrow(yearlyID) > 0) {
            yearlyID.info <- capture.output(str(yearlyID))[-1]
        } else {
            yearlyID.info <- NULL
        }
    } else {yearlyID.info <- NULL}
    

        ##store data sets
    data.out <- list( )
    if(is.null(obsVars)) {
        data.out$obsCovs <- NULL
    } else {data.out$obsCovs <- obsID}
    if(is.null(siteVars)) {
        data.out$siteCovs <- NULL
    } else {data.out$siteCovs <- siteID}
    if(is.null(yearlyVars)) {
        data.out$yearlySiteCovs <- NULL
    } else {data.out$yearlySiteCovs <- yearlyID}


    ##assemble results
    result <- list("predictors" = unique.predictors,
                   "data" = data.out)
    class(result) <- "extractX"
    return(result)
}



##unmarkedFitOccuFP
extractX.AICunmarkedFitOccuFP <- function(cand.set, parm.type = NULL, ...) {

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?extractX for details\n")}

    ##extract predictors from list
    ##psi
    if(identical(parm.type, "psi")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@stateformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }

    ##fp
    if(identical(parm.type, "falsepos") || identical(parm.type, "fp")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@FPformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }

    ##certain
    if(identical(parm.type, "certain")) {
        ##check that parameter appears in all models
        parfreq <- sum(sapply(cand.set, FUN = function(i) any(names(i@estimates@estimates) == parm.type)))
        if(!identical(length(cand.set), parfreq)) {
            stop("\nParameter \'b\' does not appear in all models\n")
        }
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@Bformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }
    
    ##detect
    if(identical(parm.type, "detect")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@detformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }

    ##extract based on "+"
    form.noplus <- unlist(sapply(form.list, FUN = function(i) strsplit(i, split = "\\+")))
    ##remove extra white space
    form.clean <- gsub("(^ +)|( +$)", "", form.noplus)
    unique.clean <- unique(form.clean)
    ##exclude empty strings and intercept
    unique.predictors <- unique.clean[nchar(unique.clean) != 0 & unique.clean != "1"]
    
    ##extract data from model objects - identical for all models
    dsets <- lapply(cand.set, FUN = function(i) unmarked::getData(i))
    ##check that same data are used
    unique.dsets <- unique(dsets)
    if(length(unique.dsets) != 1) stop("\nData sets differ across models:\n check data carefully\n")
    unFrame <- unique.dsets[[1]]
    
    ##extract siteCovs
    siteVars <- siteCovs(unFrame)
    ##extract obsCovs
    obsVars <- obsCovs(unFrame)
        
    ##check for interactions specified with *
    inter.star <- any(regexpr("\\*", unique.predictors) != -1)
        
    ##check for interaction terms
    inter.id <- any(regexpr("\\:", unique.predictors) != -1)

    ##inter.star and inter.id
    if(inter.star && inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.nostar.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(terms.nostar.clean, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }


    ##inter.star
    if(inter.star && !inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
    }
    
        
    ##inter.id
    if(!inter.star && inter.id) {
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }
        

    ##none
    if(!inter.star && !inter.id) {
        ##remove extra white space
        terms.clean <- unique.predictors
    }

        
    ##combine in single character vector
    final.predictors <- unique(terms.clean)
        
    ##find where predictors occur
    if(!is.null(obsVars)) {
        obsID <- obsVars[, intersect(final.predictors, names(obsVars)), drop = FALSE]
        if(nrow(obsID) > 0) {
            obsID.info <- capture.output(str(obsID))[-1]
        } else {
            obsID.info <- NULL
        }
    } else {obsID.info <- NULL}
    
    if(!is.null(siteVars)) {
        siteID <- siteVars[, intersect(final.predictors, names(siteVars)), drop = FALSE]
        if(nrow(siteID) > 0) {
            siteID.info <- capture.output(str(siteID))[-1]
        } else {
            siteID.info <- NULL
        }
    } else {siteID.info <- NULL}
    

    ##store data sets
    data.out <- list( )
    if(is.null(obsVars)) {
        data.out$obsCovs <- NULL
    } else {data.out$obsCovs <- obsID}
    if(is.null(siteVars)) {
        data.out$siteCovs <- NULL
    } else {data.out$siteCovs <- siteID}

    
    ##assemble results
    result <- list("predictors" = unique.predictors,
                   "data" = data.out)
    class(result) <- "extractX"
    return(result)
}



##unmarkedFitGMM
extractX.AICunmarkedFitMPois <- function(cand.set, parm.type = NULL, ...) {

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?extractX for details\n")}

    ##extract predictors from list
    ##lambda
    if(identical(parm.type, "lambda")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formula[[3]]))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }

    
    ##detect
    if(identical(parm.type, "detect")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formula[[2]]))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }


    ##extract based on "+"
    form.noplus <- unlist(sapply(form.list, FUN = function(i) strsplit(i, split = "\\+")))
    ##remove extra white space
    form.clean <- gsub("(^ +)|( +$)", "", form.noplus)
    unique.clean <- unique(form.clean)
    ##exclude empty strings and intercept
    unique.predictors <- unique.clean[nchar(unique.clean) != 0 & unique.clean != "1"]
    

    ##extract data from model objects - identical for all models
    dsets <- lapply(cand.set, FUN = function(i) unmarked::getData(i))
    ##check that same data are used
    unique.dsets <- unique(dsets)
    if(length(unique.dsets) != 1) stop("\nData sets differ across models:\n check data carefully\n")
    unFrame <- unique.dsets[[1]]
    
    ##extract siteCovs
    siteVars <- siteCovs(unFrame)
    ##extract obsCovs
    obsVars <- obsCovs(unFrame)
        
    ##check for interactions specified with *
    inter.star <- any(regexpr("\\*", unique.predictors) != -1)
        
    ##check for interaction terms
    inter.id <- any(regexpr("\\:", unique.predictors) != -1)

    ##inter.star and inter.id
    if(inter.star && inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.nostar.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(terms.nostar.clean, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }


    ##inter.star
    if(inter.star && !inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
    }
    
        
    ##inter.id
    if(!inter.star && inter.id) {
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }
        

    ##none
    if(!inter.star && !inter.id) {
        ##remove extra white space
        terms.clean <- unique.predictors
    }

        
    ##combine in single character vector
    final.predictors <- unique(terms.clean)
        
    ##find where predictors occur
    if(!is.null(obsVars)) {
        obsID <- obsVars[, intersect(final.predictors, names(obsVars)), drop = FALSE]
        if(nrow(obsID) > 0) {
            obsID.info <- capture.output(str(obsID))[-1]
        } else {
            obsID.info <- NULL
        }
    } else {obsID.info <- NULL}
    
    if(!is.null(siteVars)) {
        siteID <- siteVars[, intersect(final.predictors, names(siteVars)), drop = FALSE]
        if(nrow(siteID) > 0) {
            siteID.info <- capture.output(str(siteID))[-1]
        } else {
            siteID.info <- NULL
        }
    } else {siteID.info <- NULL}
    

    ##store data sets
    data.out <- list( )
    if(is.null(obsVars)) {
        data.out$obsCovs <- NULL
    } else {data.out$obsCovs <- obsID}
    if(is.null(siteVars)) {
        data.out$siteCovs <- NULL
    } else {data.out$siteCovs <- siteID}
    

    ##assemble results
    result <- list("predictors" = unique.predictors,
                   "data" = data.out)
    class(result) <- "extractX"
    return(result)
}



##unmarkedFitGMM
extractX.AICunmarkedFitGMM <- function(cand.set, parm.type = NULL, ...) {

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?extractX for details\n")}

    ##extract predictors from list
    ##lambda
    if(identical(parm.type, "lambda")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formlist$lambdaformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }


    ##phi
    if(identical(parm.type, "phi")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formlist$phiformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }

    
    ##detect
    if(identical(parm.type, "detect")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formlist$pformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }


    ##extract based on "+"
    form.noplus <- unlist(sapply(form.list, FUN = function(i) strsplit(i, split = "\\+")))
    ##remove extra white space
    form.clean <- gsub("(^ +)|( +$)", "", form.noplus)
    unique.clean <- unique(form.clean)
    ##exclude empty strings and intercept
    unique.predictors <- unique.clean[nchar(unique.clean) != 0 & unique.clean != "1"]
    

    ##extract data from model objects - identical for all models
    dsets <- lapply(cand.set, FUN = function(i) unmarked::getData(i))
    ##check that same data are used
    unique.dsets <- unique(dsets)
    if(length(unique.dsets) != 1) stop("\nData sets differ across models:\n check data carefully\n")
    unFrame <- unique.dsets[[1]]
    
    ##extract siteCovs
    siteVars <- siteCovs(unFrame)
    ##extract obsCovs
    obsVars <- obsCovs(unFrame)
    ##extract yearlySiteCovs
    yearlyVars <- yearlySiteCovs(unFrame)
    
    ##check for interactions specified with *
    inter.star <- any(regexpr("\\*", unique.predictors) != -1)
        
    ##check for interaction terms
    inter.id <- any(regexpr("\\:", unique.predictors) != -1)

    ##inter.star and inter.id
    if(inter.star && inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.nostar.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(terms.nostar.clean, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }


    ##inter.star
    if(inter.star && !inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
    }
    
        
    ##inter.id
    if(!inter.star && inter.id) {
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }
        

    ##none
    if(!inter.star && !inter.id) {
        ##remove extra white space
        terms.clean <- unique.predictors
    }

        
    ##combine in single character vector
    final.predictors <- unique(terms.clean)
        
    ##find where predictors occur
    if(!is.null(obsVars)) {
        obsID <- obsVars[, intersect(final.predictors, names(obsVars)), drop = FALSE]
        if(nrow(obsID) > 0) {
            obsID.info <- capture.output(str(obsID))[-1]
        } else {
            obsID.info <- NULL
        }
    } else {obsID.info <- NULL}
    
    if(!is.null(siteVars)) {
        siteID <- siteVars[, intersect(final.predictors, names(siteVars)), drop = FALSE]
        if(nrow(siteID) > 0) {
            siteID.info <- capture.output(str(siteID))[-1]
        } else {
            siteID.info <- NULL
        }
    } else {siteID.info <- NULL}
    
    if(!is.null(yearlyVars)) {
        yearlyID <- yearlyVars[, intersect(final.predictors, names(yearlyVars)), drop = FALSE]
        if(nrow(yearlyID) > 0) {
            yearlyID.info <- capture.output(str(yearlyID))[-1]
        } else {
            yearlyID.info <- NULL
        }
    } else {yearlyID.info <- NULL}
    

    ##store data sets
    data.out <- list( )
    if(is.null(obsVars)) {
        data.out$obsCovs <- NULL
    } else {data.out$obsCovs <- obsID}
    if(is.null(siteVars)) {
        data.out$siteCovs <- NULL
    } else {data.out$siteCovs <- siteID}
    if(is.null(yearlyVars)) {
        data.out$yearlySiteCovs <- NULL
    } else {data.out$yearlySiteCovs <- yearlyID}


    ##assemble results
    result <- list("predictors" = unique.predictors,
                   "data" = data.out)
    class(result) <- "extractX"
    return(result)
}



##unmarkedFitGPC
extractX.AICunmarkedFitGPC <- function(cand.set, parm.type = NULL, ...) {

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?extractX for details\n")}

    ##extract predictors from list
    ##lambda
    if(identical(parm.type, "lambda")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formlist$lambdaformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }


    ##phi
    if(identical(parm.type, "phi")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formlist$phiformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }

    
    ##detect
    if(identical(parm.type, "detect")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formlist$pformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }


    ##extract based on "+"
    form.noplus <- unlist(sapply(form.list, FUN = function(i) strsplit(i, split = "\\+")))
    ##remove extra white space
    form.clean <- gsub("(^ +)|( +$)", "", form.noplus)
    unique.clean <- unique(form.clean)
    ##exclude empty strings and intercept
    unique.predictors <- unique.clean[nchar(unique.clean) != 0 & unique.clean != "1"]
    

    ##extract data from model objects - identical for all models
    dsets <- lapply(cand.set, FUN = function(i) unmarked::getData(i))
    ##check that same data are used
    unique.dsets <- unique(dsets)
    if(length(unique.dsets) != 1) stop("\nData sets differ across models:\n check data carefully\n")
    unFrame <- unique.dsets[[1]]
    
    ##extract siteCovs
    siteVars <- siteCovs(unFrame)
    ##extract obsCovs
    obsVars <- obsCovs(unFrame)
    ##extract yearlySiteCovs
    yearlyVars <- yearlySiteCovs(unFrame)
    
    ##check for interactions specified with *
    inter.star <- any(regexpr("\\*", unique.predictors) != -1)
        
    ##check for interaction terms
    inter.id <- any(regexpr("\\:", unique.predictors) != -1)

    ##inter.star and inter.id
    if(inter.star && inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.nostar.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(terms.nostar.clean, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }


    ##inter.star
    if(inter.star && !inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
    }
    
        
    ##inter.id
    if(!inter.star && inter.id) {
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }
        

    ##none
    if(!inter.star && !inter.id) {
        ##remove extra white space
        terms.clean <- unique.predictors
    }

        
    ##combine in single character vector
    final.predictors <- unique(terms.clean)
        
    ##find where predictors occur
    if(!is.null(obsVars)) {
        obsID <- obsVars[, intersect(final.predictors, names(obsVars)), drop = FALSE]
        if(nrow(obsID) > 0) {
            obsID.info <- capture.output(str(obsID))[-1]
        } else {
            obsID.info <- NULL
        }
    } else {obsID.info <- NULL}
    
    if(!is.null(siteVars)) {
        siteID <- siteVars[, intersect(final.predictors, names(siteVars)), drop = FALSE]
        if(nrow(siteID) > 0) {
            siteID.info <- capture.output(str(siteID))[-1]
        } else {
            siteID.info <- NULL
        }
    } else {siteID.info <- NULL}
    
    if(!is.null(yearlyVars)) {
        yearlyID <- yearlyVars[, intersect(final.predictors, names(yearlyVars)), drop = FALSE]
        if(nrow(yearlyID) > 0) {
            yearlyID.info <- capture.output(str(yearlyID))[-1]
        } else {
            yearlyID.info <- NULL
        }
    } else {yearlyID.info <- NULL}
    

    ##store data sets
    data.out <- list( )
    if(is.null(obsVars)) {
        data.out$obsCovs <- NULL
    } else {data.out$obsCovs <- obsID}
    if(is.null(siteVars)) {
        data.out$siteCovs <- NULL
    } else {data.out$siteCovs <- siteID}
    if(is.null(yearlyVars)) {
        data.out$yearlySiteCovs <- NULL
    } else {data.out$yearlySiteCovs <- yearlyID}

    
    ##assemble results
    result <- list("predictors" = unique.predictors,
                   "data" = data.out)
    class(result) <- "extractX"
    return(result)
}



##unmarkedFitOccu
extractX.AICunmarkedFitOccuMulti <- function(cand.set, parm.type = NULL, ...) {

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?extractX for details\n")}

    ##extract predictors from list
    ##psi
    if(identical(parm.type, "psi")) {
        form.list <- lapply(cand.set, FUN = function(x) names(x@estimates@estimates$state@estimates))
    }

    ##detect
    if(identical(parm.type, "detect")) {
        form.list <- lapply(cand.set, FUN = function(x) names(x@estimates@estimates$det@estimates))
    }


    ##exclude empty strings and intercept
    formStrings <- unlist(form.list)
    notInclude <- grep(pattern = "(Intercept)", x = formStrings)
    formNoInt <- formStrings[-notInclude]

    ##extract only variable names
    formJustVars <- unlist(strsplit(formNoInt, split = "\\]"))
    formMat <- matrix(data = formJustVars, ncol = 2, byrow = TRUE)
    ##remove extra white space
    form.clean <- gsub("(^ +)|( +$)", "", formMat[, 2])
    unique.predictors <- unique(form.clean)    

    ##extract data from model objects - identical for all models
    dsets <- lapply(cand.set, FUN = function(i) unmarked::getData(i))
    ##check that same data are used
    unique.dsets <- unique(dsets)
    if(length(unique.dsets) != 1) stop("\nData sets differ across models:\n check data carefully\n")
    unFrame <- unique.dsets[[1]]
    
    ##extract siteCovs
    siteVars <- siteCovs(unFrame)
    ##extract obsCovs
    obsVars <- obsCovs(unFrame)

    ##check for interactions specified with *
    inter.star <- any(regexpr("\\*", unique.predictors) != -1)
        
    ##check for interaction terms
    inter.id <- any(regexpr("\\:", unique.predictors) != -1)

    ##inter.star and inter.id
    if(inter.star && inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.nostar.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(terms.nostar.clean, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }


    ##inter.star
    if(inter.star && !inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
    }
    
        
    ##inter.id
    if(!inter.star && inter.id) {
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }
        

    ##none
    if(!inter.star && !inter.id) {
        ##remove extra white space
        terms.clean <- unique.predictors
    }

        
    ##combine in single character vector
    final.predictors <- unique(terms.clean)
        
    ##find where predictors occur
    if(!is.null(obsVars)) {
        obsID <- obsVars[, intersect(final.predictors, names(obsVars)), drop = FALSE]
        if(nrow(obsID) > 0) {
            obsID.info <- capture.output(str(obsID))[-1]
        } else {
            obsID.info <- NULL
        }
    } else {obsID.info <- NULL}
    
    if(!is.null(siteVars)) {
        siteID <- siteVars[, intersect(final.predictors, names(siteVars)), drop = FALSE]
        if(nrow(siteID) > 0) {
            siteID.info <- capture.output(str(siteID))[-1]
        } else {
            siteID.info <- NULL
        }
    } else {siteID.info <- NULL}

    ##store data sets
    data.out <- list( )
    if(is.null(obsVars)) {
        data.out$obsCovs <- NULL
    } else {data.out$obsCovs <- obsID}
    if(is.null(siteVars)) {
        data.out$siteCovs <- NULL
    } else {data.out$siteCovs <- siteID}

    
    ##assemble results
    result <- list("predictors" = unique.predictors,
                   "data" = data.out)
    class(result) <- "extractX"
    return(result)
}



##unmarkedFitOccuMS
extractX.AICunmarkedFitOccuMS <- function(cand.set, parm.type = NULL, ...) {

  ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?extractX for details\n")}

    ##extract predictors from list
    ##psi
    if(identical(parm.type, "psi")) {
        form.list <- lapply(cand.set, FUN = function(x) x@psiformulas)
    }

    ##detect
    if(identical(parm.type, "detect")) {
        form.list <- lapply(cand.set, FUN = function(x) x@detformulas)
    }

    ##transition
    if(identical(parm.type, "phi")) {
        ##check that parameter appears in all models
        nseasons <- unique(sapply(cand.set, FUN = function(i) i@data@numPrimary))
        if(nseasons == 1) {
            stop("\nParameter \'phi\' does not appear in single-season models\n")
        }
        
        form.list <- lapply(cand.set, FUN = function(x) x@phiformulas)
    }

    ##exclude empty strings and intercept
    formStrings <- gsub(pattern = "~", replacement = "",
                        x = unlist(form.list))
    formNoInt <- formStrings[formStrings != "1"]
    ##remove extra white space
    form.clean <- gsub("(^ +)|( +$)", "", formNoInt)
    unique.predictors <- unique(form.clean)    

    ##extract data from model objects - identical for all models
    dsets <- lapply(cand.set, FUN = function(i) unmarked::getData(i))
    ##check that same data are used
    unique.dsets <- unique(dsets)
    if(length(unique.dsets) != 1) stop("\nData sets differ across models:\n check data carefully\n")
    unFrame <- unique.dsets[[1]]
    
    ##extract siteCovs
    siteVars <- siteCovs(unFrame)
    ##extract obsCovs
    obsVars <- obsCovs(unFrame)

    ##check for interactions specified with *
    inter.star <- any(regexpr("\\*", unique.predictors) != -1)
        
    ##check for interaction terms
    inter.id <- any(regexpr("\\:", unique.predictors) != -1)

    ##inter.star and inter.id
    if(inter.star && inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.nostar.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(terms.nostar.clean, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }


    ##inter.star
    if(inter.star && !inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
    }
    
        
    ##inter.id
    if(!inter.star && inter.id) {
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }
        

    ##none
    if(!inter.star && !inter.id) {
        ##remove extra white space
        terms.clean <- unique.predictors
    }

        
    ##combine in single character vector
    final.predictors <- unique(terms.clean)
        
    ##find where predictors occur
    if(!is.null(obsVars)) {
        obsID <- obsVars[, intersect(final.predictors, names(obsVars)), drop = FALSE]
        if(nrow(obsID) > 0) {
            obsID.info <- capture.output(str(obsID))[-1]
        } else {
            obsID.info <- NULL
        }
    } else {obsID.info <- NULL}
    
    if(!is.null(siteVars)) {
        siteID <- siteVars[, intersect(final.predictors, names(siteVars)), drop = FALSE]
        if(nrow(siteID) > 0) {
            siteID.info <- capture.output(str(siteID))[-1]
        } else {
            siteID.info <- NULL
        }
    } else {siteID.info <- NULL}

    ##store data sets
    data.out <- list( )
    if(is.null(obsVars)) {
        data.out$obsCovs <- NULL
    } else {data.out$obsCovs <- obsID}
    if(is.null(siteVars)) {
        data.out$siteCovs <- NULL
    } else {data.out$siteCovs <- siteID}

    
    ##assemble results
    result <- list("predictors" = unique.predictors,
                   "data" = data.out)
    class(result) <- "extractX"
    return(result)
}



##unmarkedFitOccuTTD
extractX.AICunmarkedFitOccuTTD <- function(cand.set, parm.type = NULL, ...) {

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?extractX for details\n")}

    ##extract predictors from list
    ##psi
    if(identical(parm.type, "psi")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@psiformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }

    ##gamma
    if(identical(parm.type, "gamma")) {
        nseasons <- unique(sapply(cand.set, FUN = function(i) i@data@numPrimary))
        if(nseasons == 1) {
            stop("\nParameter \'gamma\' does not appear in single-season models\n")
        }
        
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@gamformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }

    ##epsilon
    if(identical(parm.type, "epsilon")) {
        nseasons <- unique(sapply(cand.set, FUN = function(i) i@data@numPrimary))
        if(nseasons == 1) {
            stop("\nParameter \'epsilon\' does not appear in single-season models\n")
        }

        form.list <- as.character(lapply(cand.set, FUN = function(x) x@epsformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }

    ##detect
    if(identical(parm.type, "detect")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@detformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }


    ##extract based on "+"
    form.noplus <- unlist(sapply(form.list, FUN = function(i) strsplit(i, split = "\\+")))
    ##remove extra white space
    form.clean <- gsub("(^ +)|( +$)", "", form.noplus)
    unique.clean <- unique(form.clean)
    ##exclude empty strings and intercept
    unique.predictors <- unique.clean[nchar(unique.clean) != 0 & unique.clean != "1"]
    
    ##extract data from model objects - identical for all models
    dsets <- lapply(cand.set, FUN = function(i) unmarked::getData(i))
    ##check that same data are used
    unique.dsets <- unique(dsets)
    if(length(unique.dsets) != 1) stop("\nData sets differ across models:\n check data carefully\n")
    unFrame <- unique.dsets[[1]]
    
    ##extract siteCovs
    siteVars <- siteCovs(unFrame)
    ##extract obsCovs
    obsVars <- obsCovs(unFrame)
    ##extract yearlySiteCovs
    yearlyVars <- yearlySiteCovs(unFrame)
    
    ##check for interactions specified with *
    inter.star <- any(regexpr("\\*", unique.predictors) != -1)
        
    ##check for interaction terms
    inter.id <- any(regexpr("\\:", unique.predictors) != -1)

    ##inter.star and inter.id
    if(inter.star && inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.nostar.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(terms.nostar.clean, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }


    ##inter.star
    if(inter.star && !inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
    }
    
        
    ##inter.id
    if(!inter.star && inter.id) {
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }
        

    ##none
    if(!inter.star && !inter.id) {
        ##remove extra white space
        terms.clean <- unique.predictors
    }

        
    ##combine in single character vector
    final.predictors <- unique(terms.clean)
        
    ##find where predictors occur
    if(!is.null(obsVars)) {
        obsID <- obsVars[, intersect(final.predictors, names(obsVars)), drop = FALSE]
        if(nrow(obsID) > 0) {
            obsID.info <- capture.output(str(obsID))[-1]
        } else {
            obsID.info <- NULL
        }
    } else {obsID.info <- NULL}
    
    if(!is.null(siteVars)) {
        siteID <- siteVars[, intersect(final.predictors, names(siteVars)), drop = FALSE]
        if(nrow(siteID) > 0) {
            siteID.info <- capture.output(str(siteID))[-1]
        } else {
            siteID.info <- NULL
        }
    } else {siteID.info <- NULL}
    
    if(!is.null(yearlyVars)) {
        yearlyID <- yearlyVars[, intersect(final.predictors, names(yearlyVars)), drop = FALSE]
        if(nrow(yearlyID) > 0) {
            yearlyID.info <- capture.output(str(yearlyID))[-1]
        } else {
            yearlyID.info <- NULL
        }
    } else {yearlyID.info <- NULL}

    
    ##store data sets
    data.out <- list( )
    if(is.null(obsVars)) {
        data.out$obsCovs <- NULL
    } else {data.out$obsCovs <- obsID}
    if(is.null(siteVars)) {
        data.out$siteCovs <- NULL
    } else {data.out$siteCovs <- siteID}
    if(is.null(yearlyVars)) {
        data.out$yearlySiteCovs <- NULL
    } else {data.out$yearlySiteCovs <- yearlyID}

    
    ##assemble results
    result <- list("predictors" = unique.predictors,
                   "data" = data.out)
    class(result) <- "extractX"
    return(result)
}



##unmarkedFitMMO
extractX.AICunmarkedFitMMO <- function(cand.set, parm.type = NULL, ...) {

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?extractX for details\n")}

    ##extract predictors from list
    ##lambda
    if(identical(parm.type, "lambda")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formlist$lambdaformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }

    ##gamma
    if(identical(parm.type, "gamma")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formlist$gammaformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }

    ##omega
    if(identical(parm.type, "omega")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formlist$omegaformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }

    ##iota
    if(identical(parm.type, "iota")) {
        ##check that parameter appears in all models
        parfreq <- sum(sapply(cand.set, FUN = function(i) any(names(i@estimates@estimates) == parm.type)))
        if(!identical(length(cand.set), parfreq)) {
            stop("\nParameter \'iota\' does not appear in all models\n")
        }
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formlist$iotaformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }

    ##detect
    if(identical(parm.type, "detect")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formlist$pformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }


    ##extract based on "+"
    form.noplus <- unlist(sapply(form.list, FUN = function(i) strsplit(i, split = "\\+")))
    ##remove extra white space
    form.clean <- gsub("(^ +)|( +$)", "", form.noplus)
    unique.clean <- unique(form.clean)
    ##exclude empty strings and intercept
    unique.predictors <- unique.clean[nchar(unique.clean) != 0 & unique.clean != "1"]
    
    ##extract data from model objects - identical for all models
    dsets <- lapply(cand.set, FUN = function(i) unmarked::getData(i))
    ##check that same data are used
    unique.dsets <- unique(dsets)
    if(length(unique.dsets) != 1) stop("\nData sets differ across models:\n check data carefully\n")
    unFrame <- unique.dsets[[1]]
    
    ##extract siteCovs
    siteVars <- siteCovs(unFrame)
    ##extract obsCovs
    obsVars <- obsCovs(unFrame)
    ##extract yearlySiteCovs
    yearlyVars <- yearlySiteCovs(unFrame)
    
    ##check for interactions specified with *
    inter.star <- any(regexpr("\\*", unique.predictors) != -1)
        
    ##check for interaction terms
    inter.id <- any(regexpr("\\:", unique.predictors) != -1)

    ##inter.star and inter.id
    if(inter.star && inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.nostar.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(terms.nostar.clean, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }


    ##inter.star
    if(inter.star && !inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
    }
    
        
    ##inter.id
    if(!inter.star && inter.id) {
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }
        

    ##none
    if(!inter.star && !inter.id) {
        ##remove extra white space
        terms.clean <- unique.predictors
    }

        
    ##combine in single character vector
    final.predictors <- unique(terms.clean)
        
    ##find where predictors occur
    if(!is.null(obsVars)) {
        obsID <- obsVars[, intersect(final.predictors, names(obsVars)), drop = FALSE]
        if(nrow(obsID) > 0) {
            obsID.info <- capture.output(str(obsID))[-1]
        } else {
            obsID.info <- NULL
        }
    } else {obsID.info <- NULL}
    
    if(!is.null(siteVars)) {
        siteID <- siteVars[, intersect(final.predictors, names(siteVars)), drop = FALSE]
        if(nrow(siteID) > 0) {
            siteID.info <- capture.output(str(siteID))[-1]
        } else {
            siteID.info <- NULL
        }
    } else {siteID.info <- NULL}
    
    if(!is.null(yearlyVars)) {
        yearlyID <- yearlyVars[, intersect(final.predictors, names(yearlyVars)), drop = FALSE]
        if(nrow(yearlyID) > 0) {
            yearlyID.info <- capture.output(str(yearlyID))[-1]
        } else {
            yearlyID.info <- NULL
        }
    } else {yearlyID.info <- NULL}
    

    ##store data sets
    data.out <- list( )
    if(is.null(obsVars)) {
        data.out$obsCovs <- NULL
    } else {data.out$obsCovs <- obsID}
    if(is.null(siteVars)) {
        data.out$siteCovs <- NULL
    } else {data.out$siteCovs <- siteID}
    if(is.null(yearlyVars)) {
        data.out$yearlySiteCovs <- NULL
    } else {data.out$yearlySiteCovs <- yearlyID}

    
    ##assemble results
    result <- list("predictors" = unique.predictors,
                   "data" = data.out)
    class(result) <- "extractX"
    return(result)
}



##unmarkedFitDSO
extractX.AICunmarkedFitDSO <- function(cand.set, parm.type = NULL, ...) {

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?extractX for details\n")}

    ##extract predictors from list
    ##lambda
    if(identical(parm.type, "lambda")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formlist$lambdaformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }

    ##gamma
    if(identical(parm.type, "gamma")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formlist$gammaformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }

    ##omega
    if(identical(parm.type, "omega")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formlist$omegaformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }

    ##iota
    if(identical(parm.type, "iota")) {
        ##check that parameter appears in all models
        parfreq <- sum(sapply(cand.set, FUN = function(i) any(names(i@estimates@estimates) == parm.type)))
        if(!identical(length(cand.set), parfreq)) {
            stop("\nParameter \'iota\' does not appear in all models\n")
        }
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formlist$iotaformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }

    ##detect
    if(identical(parm.type, "detect")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formlist$pformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }


    ##extract based on "+"
    form.noplus <- unlist(sapply(form.list, FUN = function(i) strsplit(i, split = "\\+")))
    ##remove extra white space
    form.clean <- gsub("(^ +)|( +$)", "", form.noplus)
    unique.clean <- unique(form.clean)
    ##exclude empty strings and intercept
    unique.predictors <- unique.clean[nchar(unique.clean) != 0 & unique.clean != "1"]
    
    ##extract data from model objects - identical for all models
    dsets <- lapply(cand.set, FUN = function(i) unmarked::getData(i))
    ##check that same data are used
    unique.dsets <- unique(dsets)
    if(length(unique.dsets) != 1) stop("\nData sets differ across models:\n check data carefully\n")
    unFrame <- unique.dsets[[1]]
    
    ##extract siteCovs
    siteVars <- siteCovs(unFrame)
    ##extract obsCovs
    obsVars <- obsCovs(unFrame)
    ##extract yearlySiteCovs
    yearlyVars <- yearlySiteCovs(unFrame)
    
    ##check for interactions specified with *
    inter.star <- any(regexpr("\\*", unique.predictors) != -1)
        
    ##check for interaction terms
    inter.id <- any(regexpr("\\:", unique.predictors) != -1)

    ##inter.star and inter.id
    if(inter.star && inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.nostar.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(terms.nostar.clean, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }


    ##inter.star
    if(inter.star && !inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
    }
    
        
    ##inter.id
    if(!inter.star && inter.id) {
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }
        

    ##none
    if(!inter.star && !inter.id) {
        ##remove extra white space
        terms.clean <- unique.predictors
    }

        
    ##combine in single character vector
    final.predictors <- unique(terms.clean)
        
    ##find where predictors occur
    if(!is.null(obsVars)) {
        obsID <- obsVars[, intersect(final.predictors, names(obsVars)), drop = FALSE]
        if(nrow(obsID) > 0) {
            obsID.info <- capture.output(str(obsID))[-1]
        } else {
            obsID.info <- NULL
        }
    } else {obsID.info <- NULL}
    
    if(!is.null(siteVars)) {
        siteID <- siteVars[, intersect(final.predictors, names(siteVars)), drop = FALSE]
        if(nrow(siteID) > 0) {
            siteID.info <- capture.output(str(siteID))[-1]
        } else {
            siteID.info <- NULL
        }
    } else {siteID.info <- NULL}
    
    if(!is.null(yearlyVars)) {
        yearlyID <- yearlyVars[, intersect(final.predictors, names(yearlyVars)), drop = FALSE]
        if(nrow(yearlyID) > 0) {
            yearlyID.info <- capture.output(str(yearlyID))[-1]
        } else {
            yearlyID.info <- NULL
        }
    } else {yearlyID.info <- NULL}
    

    ##store data sets
    data.out <- list( )
    if(is.null(obsVars)) {
        data.out$obsCovs <- NULL
    } else {data.out$obsCovs <- obsID}
    if(is.null(siteVars)) {
        data.out$siteCovs <- NULL
    } else {data.out$siteCovs <- siteID}
    if(is.null(yearlyVars)) {
        data.out$yearlySiteCovs <- NULL
    } else {data.out$yearlySiteCovs <- yearlyID}

    
    ##assemble results
    result <- list("predictors" = unique.predictors,
                   "data" = data.out)
    class(result) <- "extractX"
    return(result)
}



##unmarkedFitGOccu
extractX.AICunmarkedFitGOccu <- function(cand.set, parm.type = NULL, ...) {

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?extractX for details\n")}

    ##extract predictors from list
    ##psi
    if(identical(parm.type, "psi")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formlist$psiformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }


    ##phi
    if(identical(parm.type, "phi")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formlist$phiformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }

    
    ##detect
    if(identical(parm.type, "detect")) {
        form.list <- as.character(lapply(cand.set, FUN = function(x) x@formlist$pformula))
        ##remove ~
        form.list <- gsub("~", replacement = "", x = form.list)
    }


    ##extract based on "+"
    form.noplus <- unlist(sapply(form.list, FUN = function(i) strsplit(i, split = "\\+")))
    ##remove extra white space
    form.clean <- gsub("(^ +)|( +$)", "", form.noplus)
    unique.clean <- unique(form.clean)
    ##exclude empty strings and intercept
    unique.predictors <- unique.clean[nchar(unique.clean) != 0 & unique.clean != "1"]
    

    ##extract data from model objects - identical for all models
    dsets <- lapply(cand.set, FUN = function(i) unmarked::getData(i))
    ##check that same data are used
    unique.dsets <- unique(dsets)
    if(length(unique.dsets) != 1) stop("\nData sets differ across models:\n check data carefully\n")
    unFrame <- unique.dsets[[1]]
    
    ##extract siteCovs
    siteVars <- siteCovs(unFrame)
    ##extract obsCovs
    obsVars <- obsCovs(unFrame)
    ##extract yearlySiteCovs
    yearlyVars <- yearlySiteCovs(unFrame)
    
    ##check for interactions specified with *
    inter.star <- any(regexpr("\\*", unique.predictors) != -1)
        
    ##check for interaction terms
    inter.id <- any(regexpr("\\:", unique.predictors) != -1)

    ##inter.star and inter.id
    if(inter.star && inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.nostar.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(terms.nostar.clean, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }


    ##inter.star
    if(inter.star && !inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
    }
    
        
    ##inter.id
    if(!inter.star && inter.id) {
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }
        

    ##none
    if(!inter.star && !inter.id) {
        ##remove extra white space
        terms.clean <- unique.predictors
    }

        
    ##combine in single character vector
    final.predictors <- unique(terms.clean)
        
    ##find where predictors occur
    if(!is.null(obsVars)) {
        obsID <- obsVars[, intersect(final.predictors, names(obsVars)), drop = FALSE]
        if(nrow(obsID) > 0) {
            obsID.info <- capture.output(str(obsID))[-1]
        } else {
            obsID.info <- NULL
        }
    } else {obsID.info <- NULL}
    
    if(!is.null(siteVars)) {
        siteID <- siteVars[, intersect(final.predictors, names(siteVars)), drop = FALSE]
        if(nrow(siteID) > 0) {
            siteID.info <- capture.output(str(siteID))[-1]
        } else {
            siteID.info <- NULL
        }
    } else {siteID.info <- NULL}
    
    if(!is.null(yearlyVars)) {
        yearlyID <- yearlyVars[, intersect(final.predictors, names(yearlyVars)), drop = FALSE]
        if(nrow(yearlyID) > 0) {
            yearlyID.info <- capture.output(str(yearlyID))[-1]
        } else {
            yearlyID.info <- NULL
        }
    } else {yearlyID.info <- NULL}
    

    ##store data sets
    data.out <- list( )
    if(is.null(obsVars)) {
        data.out$obsCovs <- NULL
    } else {data.out$obsCovs <- obsID}
    if(is.null(siteVars)) {
        data.out$siteCovs <- NULL
    } else {data.out$siteCovs <- siteID}
    if(is.null(yearlyVars)) {
        data.out$yearlySiteCovs <- NULL
    } else {data.out$yearlySiteCovs <- yearlyID}

    
    ##assemble results
    result <- list("predictors" = unique.predictors,
                   "data" = data.out)
    class(result) <- "extractX"
    return(result)
}



##unmarkedFitOccuComm
extractX.AICunmarkedFitOccuComm <- function(cand.set, parm.type = NULL, ...) {

    ##check for parm.type and stop if NULL
    if(is.null(parm.type)) {stop("\n'parm.type' must be specified for this model type, see ?extractX for details\n")}

    ##extract predictors from list
    ##get random intercepts and slopes
    allCoefs <- lapply(cand.set, FUN = function(i) randomTerms(i, addMean = TRUE))
    
    ##psi
    if(identical(parm.type, "psi")) {
        form.list <- lapply(allCoefs, FUN = function(x) x[x$Model == "psi", "Name"])
    }

    ##detect
    if(identical(parm.type, "detect")) {
        form.list <- lapply(allCoefs, FUN = function(x) x[x$Model == "p", "Name"])
    }


    ##exclude empty strings and intercept
    formStrings <- unlist(form.list)
    notInclude <- grep(pattern = "(Intercept)", x = formStrings)
    formNoInt <- formStrings[-notInclude]

    ##extract only variable names
    formJustVars <- unlist(strsplit(formNoInt, split = "\\]"))
    #formMat <- matrix(data = formJustVars, ncol = 2, byrow = TRUE)
    ##remove extra white space
    form.clean <- gsub("(^ +)|( +$)", "", formJustVars)
    unique.predictors <- unique(form.clean)    

    ##extract data from model objects - identical for all models
    dsets <- lapply(cand.set, FUN = function(i) unmarked::getData(i))
    ##check that same data are used
    unique.dsets <- unique(dsets)
    if(length(unique.dsets) != 1) stop("\nData sets differ across models:\n check data carefully\n")
    unFrame <- unique.dsets[[1]]
    
    ##extract siteCovs
    siteVars <- siteCovs(unFrame)
    ##extract obsCovs
    obsVars <- obsCovs(unFrame)

    ##check for interactions specified with *
    inter.star <- any(regexpr("\\*", unique.predictors) != -1)
        
    ##check for interaction terms
    inter.id <- any(regexpr("\\:", unique.predictors) != -1)

    ##inter.star and inter.id
    if(inter.star && inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.nostar.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(terms.nostar.clean, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }


    ##inter.star
    if(inter.star && !inter.id) {
        ##separate terms in interaction
        terms.nostar <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\*")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nostar)
    }
    
        
    ##inter.id
    if(!inter.star && inter.id) {
        ##separate terms in interaction
        terms.nointer <- unlist(sapply(unique.predictors, FUN = function(i) strsplit(i, split = "\\:")))
        ##remove extra white space
        terms.clean <- gsub("(^ +)|( +$)", "", terms.nointer)
    }
        

    ##none
    if(!inter.star && !inter.id) {
        ##remove extra white space
        terms.clean <- unique.predictors
    }

        
    ##combine in single character vector
    final.predictors <- unique(terms.clean)
        
    ##find where predictors occur
    if(!is.null(obsVars)) {
        obsID <- obsVars[, intersect(final.predictors, names(obsVars)), drop = FALSE]
        if(nrow(obsID) > 0) {
            obsID.info <- capture.output(str(obsID))[-1]
        } else {
            obsID.info <- NULL
        }
    } else {obsID.info <- NULL}
    
    if(!is.null(siteVars)) {
        siteID <- siteVars[, intersect(final.predictors, names(siteVars)), drop = FALSE]
        if(nrow(siteID) > 0) {
            siteID.info <- capture.output(str(siteID))[-1]
        } else {
            siteID.info <- NULL
        }
    } else {siteID.info <- NULL}

    ##store data sets
    data.out <- list( )
    if(is.null(obsVars)) {
        data.out$obsCovs <- NULL
    } else {data.out$obsCovs <- obsID}
    if(is.null(siteVars)) {
        data.out$siteCovs <- NULL
    } else {data.out$siteCovs <- siteID}

    
    ##assemble results
    result <- list("predictors" = unique.predictors,
                   "data" = data.out)
    class(result) <- "extractX"
    return(result)
}



##print method
print.extractX <- function(x, ...) {
    
    if(length(x$predictors) > 0) {
        cat("\nPredictors appearing in candidate models:\n")
        cat(x$predictors, sep = "    ")
        cat("\n")

        ##if unmarkedFit model
        if(!is.data.frame(x$data)) {
            ##determine number of elements
            nitems <- length(x$data)
            for(i in 1:nitems) {
                ##check if data frame contains data
                if(ncol(x$data[[i]]) > 0) {
                cat("\nStructure of predictors in ", names(x$data)[i], ":\n", sep = "")
                cat(capture.output(str(x$data[[i]]))[-1], sep = "\n")
                }
            }
            cat("\n")
        } else {
            
            cat("\nStructure of predictors:", "\n")
            cat(capture.output(str(x$data))[-1], sep = "\n")
            cat("\n")
            
        }
        
    } else {
        ##if only intercept is present
        cat("\nNo predictors appear in candidate models\n")
        cat("\n")
    }
}

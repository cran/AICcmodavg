##check SE's of parameters in model
##and identify SE's above a given threshold
##or with NaN

##mainly used with unmarkedFit objects, but also useful with classic GLM's

checkParms <- function(mod, se.max = 25, simplify = TRUE, ...) {
  UseMethod("checkParms", mod)
}



checkParms.default <- function(mod, se.max = 25, simplify = TRUE, ...) {
  stop("\nFunction not yet defined for this object class\n")
}



##unmarkedFit objects
checkParms.unmarkedFit <- function(mod, se.max = 25, simplify = TRUE, ...) {
    
    ##extract SE
    SEs <- sqrt(diag(vcov(mod)))

    ##extract names
    var.names <- names(SEs)

    ##extract parms, if unmarkedFrame, multinomial
    parm.names <- sapply(var.names, FUN = function(i) unlist(strsplit(i, split = "\\("))[1], simplify = TRUE)

    ##unique parms
    parm.id <- unique(parm.names)
    
    ##format to matrix
    ##models with several groups of parameters
    matSE <- data.frame(SEs = SEs, variable = var.names, parameter = parm.names)

    ##request simplified output
    if(identical(simplify, TRUE)) {

        ##create matrix to hold results for parm with highest SE
        out.result <- data.frame(variable = rep(NA, 1),
                                 max.se = rep(NA, 1),
                                 n.high.se = rep(NA, 1))

        ##identify maximum value of SE in model
        maxSE <- max(SEs)
        
        ##check if length = 0 (when NaN are present)
        if(is.nan(maxSE)) {
            nan.var <- as.character(matSE[is.nan(matSE$SEs), "variable"])
            if(length(nan.var) == 1) {
                nameSE <- nan.var
            } else {nameSE <- nan.var[1]}
            
            parmSE <- as.character(matSE[which(matSE$variable == nameSE), "parameter"])
            
        } else {        
            nameSE <- as.character(matSE[which(matSE$SEs == maxSE), "variable"])
            parmSE <- as.character(matSE[which(matSE$SEs == maxSE), "parameter"])
        }

        ##add to rowname
        rownames(out.result) <- parmSE
        
        out.result[, "variable"] <- nameSE
        out.result[, "max.se"] <- maxSE
        out.result[, "n.high.se"] <- length(which(matSE$SEs > se.max))
        
    }


    
    ##requesting long output
    if(identical(simplify, FALSE)) {
        
        ##create list to hold results for each parm type
        out.result <- data.frame(variable = rep(NA, length(parm.id)),
                                 max.se = rep(NA, length(parm.id)),
                                 n.high.se = rep(NA, length(parm.id)))
        rownames(out.result) <- parm.id
        
        ##for each parameter, identify maximum value of SE
        for(j in parm.id) {
            mat.parm <- matSE[matSE$parameter %in% j, ]
            maxSE <- max(mat.parm$SEs)
            out.result[j, "max.se"] <- maxSE
            nameSE <- as.character(mat.parm[which(mat.parm$SEs == maxSE), "variable"])
            ##check if NaN are present
            if(is.nan(maxSE)) {
                nan.var <- as.character(mat.parm[is.nan(mat.parm$SEs), "variable"])
                if(length(nan.var) == 1) {
                    nameSE <- nan.var
                } else {nameSE <- nan.var[1]}
            }

            out.result[j, "variable"] <- nameSE
            ##determine number of SE's > SE.limit
            out.result[j, "n.high.se"] <- length(which(mat.parm$SEs > se.max))
        }
    }

    out <- list(model.class = class(mod), se.max = se.max, result = out.result)
    class(out) <- "checkParms"
    return(out)
}



##betareg objects
checkParms.betareg <- function(mod, se.max = 25, ...) {
    
    ##extract SE
    SEs <- sqrt(diag(vcov(mod)))

    ##extract names
    var.names <- names(SEs)

    ##format as matrix
    matSE <- data.frame(SEs = SEs, variable = var.names)
    
    ##create matrix to hold results for parm with highest SE
    out.result <- data.frame(variable = rep(NA, 1),
                             max.se = rep(NA, 1),
                             n.high.se = rep(NA, 1))

    ##identify maximum value of SE in model
    maxSE <- max(SEs)
    
    ##check if length = 0 (when NaN are present)
    if(is.nan(maxSE)) {
        nan.var <- as.character(matSE[is.nan(matSE$SEs), "variable"])
        if(length(nan.var) == 1) {
            nameSE <- nan.var
        } else {nameSE <- nan.var[1]}
    } else {
        nameSE <- as.character(matSE[which(matSE$SEs == maxSE), "variable"])
    }
    
    rownames(out.result) <- "beta"
    
    out.result[, "variable"] <- nameSE
    out.result[, "max.se"] <- maxSE
    out.result[, "n.high.se"] <- length(which(matSE$SEs > se.max))

    out <- list(model.class = class(mod), se.max = se.max, result = out.result)
    class(out) <- "checkParms"
    return(out)
}



##clm objects
checkParms.clm <- function(mod, se.max = 25, ...) {
    
    ##extract SE
    SEs <- sqrt(diag(vcov(mod)))

    ##extract names
    var.names <- names(SEs)

    ##format as matrix
    matSE <- data.frame(SEs = SEs, variable = var.names)
    
    ##create matrix to hold results for parm with highest SE
    out.result <- data.frame(variable = rep(NA, 1),
                             max.se = rep(NA, 1),
                             n.high.se = rep(NA, 1))

    ##identify maximum value of SE in model
    maxSE <- max(SEs)
    
    ##check if length = 0 (when NaN are present)
    if(is.nan(maxSE)) {
        nan.var <- as.character(matSE[is.nan(matSE$SEs), "variable"])
        if(length(nan.var) == 1) {
            nameSE <- nan.var
        } else {nameSE <- nan.var[1]}
    } else {
        nameSE <- as.character(matSE[which(matSE$SEs == maxSE), "variable"])
    }
    
    rownames(out.result) <- "beta"
    
    out.result[, "variable"] <- nameSE
    out.result[, "max.se"] <- maxSE
    out.result[, "n.high.se"] <- length(which(matSE$SEs > se.max))

    out <- list(model.class = class(mod), se.max = se.max, result = out.result)
    class(out) <- "checkParms"
    return(out)
}



##clmm objects
checkParms.clmm <- function(mod, se.max = 25, ...) {
    
    ##extract SE
    SEs <- sqrt(diag(vcov(mod)))

    ##extract names
    var.names <- names(SEs)

    ##format as matrix
    matSE <- data.frame(SEs = SEs, variable = var.names)
    
    ##create matrix to hold results for parm with highest SE
    out.result <- data.frame(variable = rep(NA, 1),
                             max.se = rep(NA, 1),
                             n.high.se = rep(NA, 1))

    ##identify maximum value of SE in model
    maxSE <- max(SEs)
    
    ##check if length = 0 (when NaN are present)
    if(is.nan(maxSE)) {
        nan.var <- as.character(matSE[is.nan(matSE$SEs), "variable"])
        if(length(nan.var) == 1) {
            nameSE <- nan.var
        } else {nameSE <- nan.var[1]}
    } else {
        nameSE <- as.character(matSE[which(matSE$SEs == maxSE), "variable"])
    }
    
    rownames(out.result) <- "beta"
    
    out.result[, "variable"] <- nameSE
    out.result[, "max.se"] <- maxSE
    out.result[, "n.high.se"] <- length(which(matSE$SEs > se.max))

    out <- list(model.class = class(mod), se.max = se.max, result = out.result)
    class(out) <- "checkParms"
    return(out)
}



##coxme objects
checkParms.coxme <- function(mod, se.max = 25, ...) {
    
    ##extract SE
    SEs <- extractSE(mod)

    ##extract names
    var.names <- names(SEs)

    ##format as matrix
    matSE <- data.frame(SEs = SEs, variable = var.names)
    
    ##create matrix to hold results for parm with highest SE
    out.result <- data.frame(variable = rep(NA, 1),
                             max.se = rep(NA, 1),
                             n.high.se = rep(NA, 1))

    ##identify maximum value of SE in model
    maxSE <- max(SEs)
    
    ##check if length = 0 (when NaN are present)
    if(is.nan(maxSE)) {
        nan.var <- as.character(matSE[is.nan(matSE$SEs), "variable"])
        if(length(nan.var) == 1) {
            nameSE <- nan.var
        } else {nameSE <- nan.var[1]}
    } else {
        nameSE <- as.character(matSE[which(matSE$SEs == maxSE), "variable"])
    }
    
    rownames(out.result) <- "beta"
    
    out.result[, "variable"] <- nameSE
    out.result[, "max.se"] <- maxSE
    out.result[, "n.high.se"] <- length(which(matSE$SEs > se.max))

    out <- list(model.class = class(mod), se.max = se.max, result = out.result)
    class(out) <- "checkParms"
    return(out)
}



##coxph objects
checkParms.coxph <- function(mod, se.max = 25, ...) {
    
    ##extract SE
    SEs <- sqrt(diag(vcov(mod)))

    ##extract names
    var.names <- names(SEs)

    ##format as matrix
    matSE <- data.frame(SEs = SEs, variable = var.names)
    
    ##create matrix to hold results for parm with highest SE
    out.result <- data.frame(variable = rep(NA, 1),
                             max.se = rep(NA, 1),
                             n.high.se = rep(NA, 1))

    ##identify maximum value of SE in model
    maxSE <- max(SEs)
    
    ##check if length = 0 (when NaN are present)
    if(is.nan(maxSE)) {
        nan.var <- as.character(matSE[is.nan(matSE$SEs), "variable"])
        if(length(nan.var) == 1) {
            nameSE <- nan.var
        } else {nameSE <- nan.var[1]}
    } else {
        nameSE <- as.character(matSE[which(matSE$SEs == maxSE), "variable"])
    }
    
    rownames(out.result) <- "beta"
    
    out.result[, "variable"] <- nameSE
    out.result[, "max.se"] <- maxSE
    out.result[, "n.high.se"] <- length(which(matSE$SEs > se.max))

    out <- list(model.class = class(mod), se.max = se.max, result = out.result)
    class(out) <- "checkParms"
    return(out)
}



##GLM's
checkParms.glm <- function(mod, se.max = 25, ...) {
    
    ##extract SE
    SEs <- sqrt(diag(vcov(mod)))

    ##extract names
    var.names <- names(SEs)

    ##format as matrix
    matSE <- data.frame(SEs = SEs, variable = var.names)
    
    ##create matrix to hold results for parm with highest SE
    out.result <- data.frame(variable = rep(NA, 1),
                             max.se = rep(NA, 1),
                             n.high.se = rep(NA, 1))

    ##identify maximum value of SE in model
    maxSE <- max(SEs)
    
    ##check if length = 0 (when NaN are present)
    if(is.nan(maxSE)) {
        nan.var <- as.character(matSE[is.nan(matSE$SEs), "variable"])
        if(length(nan.var) == 1) {
            nameSE <- nan.var
        } else {nameSE <- nan.var[1]}
    } else {
        nameSE <- as.character(matSE[which(matSE$SEs == maxSE), "variable"])
    }
    
    rownames(out.result) <- "beta"
    
    out.result[, "variable"] <- nameSE
    out.result[, "max.se"] <- maxSE
    out.result[, "n.high.se"] <- length(which(matSE$SEs > se.max))

    out <- list(model.class = class(mod), se.max = se.max, result = out.result)
    class(out) <- "checkParms"
    return(out)
}



##GLM's
checkParms.glmmTMB <- function(mod, se.max = 25, ...) {
    
    ##extract SE
    SEs <- sqrt(diag(vcov(mod)$cond))

    ##extract names
    var.names <- names(SEs)

    ##format as matrix
    matSE <- data.frame(SEs = SEs, variable = var.names)
    
    ##create matrix to hold results for parm with highest SE
    out.result <- data.frame(variable = rep(NA, 1),
                             max.se = rep(NA, 1),
                             n.high.se = rep(NA, 1))

    ##identify maximum value of SE in model
    maxSE <- max(SEs)
    
    ##check if length = 0 (when NaN are present)
    if(is.nan(maxSE)) {
        nan.var <- as.character(matSE[is.nan(matSE$SEs), "variable"])
        if(length(nan.var) == 1) {
            nameSE <- nan.var
        } else {nameSE <- nan.var[1]}
    } else {
        nameSE <- as.character(matSE[which(matSE$SEs == maxSE), "variable"])
    }
    
    rownames(out.result) <- "beta"
    
    out.result[, "variable"] <- nameSE
    out.result[, "max.se"] <- maxSE
    out.result[, "n.high.se"] <- length(which(matSE$SEs > se.max))

    out <- list(model.class = class(mod), se.max = se.max, result = out.result)
    class(out) <- "checkParms"
    return(out)
}



##gls objects
checkParms.gls <- function(mod, se.max = 25, ...) {
    
    ##extract SE
    SEs <- sqrt(diag(vcov(mod)))

    ##extract names
    var.names <- names(SEs)

    ##format as matrix
    matSE <- data.frame(SEs = SEs, variable = var.names)
    
    ##create matrix to hold results for parm with highest SE
    out.result <- data.frame(variable = rep(NA, 1),
                             max.se = rep(NA, 1),
                             n.high.se = rep(NA, 1))

    ##identify maximum value of SE in model
    maxSE <- max(SEs)
    
    ##check if length = 0 (when NaN are present)
    if(is.nan(maxSE)) {
        nan.var <- as.character(matSE[is.nan(matSE$SEs), "variable"])
        if(length(nan.var) == 1) {
            nameSE <- nan.var
        } else {nameSE <- nan.var[1]}
    } else {
        nameSE <- as.character(matSE[which(matSE$SEs == maxSE), "variable"])
    }
    
    rownames(out.result) <- "beta"
    
    out.result[, "variable"] <- nameSE
    out.result[, "max.se"] <- maxSE
    out.result[, "n.high.se"] <- length(which(matSE$SEs > se.max))

    out <- list(model.class = class(mod), se.max = se.max, result = out.result)
    class(out) <- "checkParms"
    return(out)
}



##gnls objects
checkParms.gnls <- function(mod, se.max = 25, ...) {
    
    ##extract SE
    SEs <- sqrt(diag(vcov(mod)))

    ##extract names
    var.names <- names(SEs)

    ##format as matrix
    matSE <- data.frame(SEs = SEs, variable = var.names)
    
    ##create matrix to hold results for parm with highest SE
    out.result <- data.frame(variable = rep(NA, 1),
                             max.se = rep(NA, 1),
                             n.high.se = rep(NA, 1))

    ##identify maximum value of SE in model
    maxSE <- max(SEs)
    
    ##check if length = 0 (when NaN are present)
    if(is.nan(maxSE)) {
        nan.var <- as.character(matSE[is.nan(matSE$SEs), "variable"])
        if(length(nan.var) == 1) {
            nameSE <- nan.var
        } else {nameSE <- nan.var[1]}
    } else {
        nameSE <- as.character(matSE[which(matSE$SEs == maxSE), "variable"])
    }
    
    rownames(out.result) <- "beta"
    
    out.result[, "variable"] <- nameSE
    out.result[, "max.se"] <- maxSE
    out.result[, "n.high.se"] <- length(which(matSE$SEs > se.max))

    out <- list(model.class = class(mod), se.max = se.max, result = out.result)
    class(out) <- "checkParms"
    return(out)
}



##hurdle objects
checkParms.hurdle <- function(mod, se.max = 25, ...) {
    
    ##extract SE
    SEs <- sqrt(diag(vcov(mod)))

    ##extract names
    var.names <- names(SEs)

    ##format as matrix
    matSE <- data.frame(SEs = SEs, variable = var.names)
    
    ##create matrix to hold results for parm with highest SE
    out.result <- data.frame(variable = rep(NA, 1),
                             max.se = rep(NA, 1),
                             n.high.se = rep(NA, 1))

    ##identify maximum value of SE in model
    maxSE <- max(SEs)
    
    ##check if length = 0 (when NaN are present)
    if(is.nan(maxSE)) {
        nan.var <- as.character(matSE[is.nan(matSE$SEs), "variable"])
        if(length(nan.var) == 1) {
            nameSE <- nan.var
        } else {nameSE <- nan.var[1]}
    } else {
        nameSE <- as.character(matSE[which(matSE$SEs == maxSE), "variable"])
    }
    
    rownames(out.result) <- "beta"
    
    out.result[, "variable"] <- nameSE
    out.result[, "max.se"] <- maxSE
    out.result[, "n.high.se"] <- length(which(matSE$SEs > se.max))

    out <- list(model.class = class(mod), se.max = se.max, result = out.result)
    class(out) <- "checkParms"
    return(out)
}



##classic linear regression (lm)
checkParms.lm <- function(mod, se.max = 25, ...) {
    
    ##extract SE
    SEs <- sqrt(diag(vcov(mod)))

    ##extract names
    var.names <- names(SEs)

    ##format as matrix
    matSE <- data.frame(SEs = SEs, variable = var.names)
    
    ##create matrix to hold results for parm with highest SE
    out.result <- data.frame(variable = rep(NA, 1),
                             max.se = rep(NA, 1),
                             n.high.se = rep(NA, 1))

    ##identify maximum value of SE in model
    maxSE <- max(SEs)
    
    ##check if length = 0 (when NaN are present)
    if(is.nan(maxSE)) {
        nan.var <- as.character(matSE[is.nan(matSE$SEs), "variable"])
        if(length(nan.var) == 1) {
            nameSE <- nan.var
        } else {nameSE <- nan.var[1]}
    } else {
        nameSE <- as.character(matSE[which(matSE$SEs == maxSE), "variable"])
    }
    
    rownames(out.result) <- "beta"
    
    out.result[, "variable"] <- nameSE
    out.result[, "max.se"] <- maxSE
    out.result[, "n.high.se"] <- length(which(matSE$SEs > se.max))

    out <- list(model.class = class(mod), se.max = se.max, result = out.result)
    class(out) <- "checkParms"
    return(out)
}



##lme objects
checkParms.lme <- function(mod, se.max = 25, ...) {
    
    ##extract SE
    SEs <- sqrt(diag(vcov(mod)))

    ##extract names
    var.names <- names(SEs)

    ##format as matrix
    matSE <- data.frame(SEs = SEs, variable = var.names)
    
    ##create matrix to hold results for parm with highest SE
    out.result <- data.frame(variable = rep(NA, 1),
                             max.se = rep(NA, 1),
                             n.high.se = rep(NA, 1))

    ##identify maximum value of SE in model
    maxSE <- max(SEs)
    
    ##check if length = 0 (when NaN are present)
    if(is.nan(maxSE)) {
        nan.var <- as.character(matSE[is.nan(matSE$SEs), "variable"])
        if(length(nan.var) == 1) {
            nameSE <- nan.var
        } else {nameSE <- nan.var[1]}
    } else {
        nameSE <- as.character(matSE[which(matSE$SEs == maxSE), "variable"])
    }
    
    rownames(out.result) <- "beta"
    
    out.result[, "variable"] <- nameSE
    out.result[, "max.se"] <- maxSE
    out.result[, "n.high.se"] <- length(which(matSE$SEs > se.max))

    out <- list(model.class = class(mod), se.max = se.max, result = out.result)
    class(out) <- "checkParms"
    return(out)
}



##lmekin objects
checkParms.lmekin <- function(mod, se.max = 25, ...) {
    
    ##extract SE
    SEs <- extractSE(mod)

    ##extract names
    var.names <- names(SEs)
    
    ##format as matrix
    matSE <- data.frame(SEs = SEs, variable = var.names)
    
    ##create matrix to hold results for parm with highest SE
    out.result <- data.frame(variable = rep(NA, 1),
                             max.se = rep(NA, 1),
                             n.high.se = rep(NA, 1))

    ##identify maximum value of SE in model
    maxSE <- max(SEs)
    
    ##check if length = 0 (when NaN are present)
    if(is.nan(maxSE)) {
        nan.var <- as.character(matSE[is.nan(matSE$SEs), "variable"])
        if(length(nan.var) == 1) {
            nameSE <- nan.var
        } else {nameSE <- nan.var[1]}
    } else {
        nameSE <- as.character(matSE[which(matSE$SEs == maxSE), "variable"])
    }
    
    rownames(out.result) <- "beta"
    
    out.result[, "variable"] <- nameSE
    out.result[, "max.se"] <- maxSE
    out.result[, "n.high.se"] <- length(which(matSE$SEs > se.max))

    out <- list(model.class = class(mod), se.max = se.max, result = out.result)
    class(out) <- "checkParms"
    return(out)
}



##maxlike objects
checkParms.maxlikeFit <- function(mod, se.max = 25, ...) {
    
    ##extract SE
    SEs <- sqrt(diag(vcov(mod)))

    ##extract names
    var.names <- names(SEs)

    ##format as matrix
    matSE <- data.frame(SEs = SEs, variable = var.names)
    
    ##create matrix to hold results for parm with highest SE
    out.result <- data.frame(variable = rep(NA, 1),
                             max.se = rep(NA, 1),
                             n.high.se = rep(NA, 1))

    ##identify maximum value of SE in model
    maxSE <- max(SEs)
    
    ##check if length = 0 (when NaN are present)
    if(is.nan(maxSE)) {
        nan.var <- as.character(matSE[is.nan(matSE$SEs), "variable"])
        if(length(nan.var) == 1) {
            nameSE <- nan.var
        } else {nameSE <- nan.var[1]}
    } else {
        nameSE <- as.character(matSE[which(matSE$SEs == maxSE), "variable"])
    }
    
    rownames(out.result) <- "beta"
    
    out.result[, "variable"] <- nameSE
    out.result[, "max.se"] <- maxSE
    out.result[, "n.high.se"] <- length(which(matSE$SEs > se.max))

    out <- list(model.class = class(mod), se.max = se.max, result = out.result)
    class(out) <- "checkParms"
    return(out)
}



##mer objects - old versions of lme4
#checkParms.mer <- function(mod, se.max = 25, ...) {
#    
#    ##extract SE
#    SEs <- sqrt(diag(vcov(mod)))
#
#    ##extract names
#    var.names <- names(SEs)
#
#    ##format as matrix
#    matSE <- data.frame(SEs = SEs, variable = var.names)
#    
#    ##create matrix to hold results for parm with highest SE
#    out.result <- data.frame(variable = rep(NA, 1),
#                             max.se = rep(NA, 1),
#                             n.high.se = rep(NA, 1))
#
#    ##identify maximum value of SE in model
#    maxSE <- max(SEs)
#    
#    ##check if length = 0 (when NaN are present)
#    if(is.nan(maxSE)) {
#        nan.var <- as.character(matSE[is.nan(matSE$SEs), "variable"])
#        if(length(nan.var) == 1) {
#            nameSE <- nan.var
#        } else {nameSE <- nan.var[1]}
#    } else {
#        nameSE <- as.character(matSE[which(matSE$SEs == maxSE), "variable"])
#    }
#    
#    rownames(out.result) <- "beta"
#    
#    out.result[, "variable"] <- nameSE
#    out.result[, "max.se"] <- maxSE
#    out.result[, "n.high.se"] <- length(which(matSE$SEs > se.max))
#
#    out <- list(model.class = class(mod), se.max = se.max, result = out.result)
#    class(out) <- "checkParms"
#    return(out)
#}



##merMod objects
checkParms.merMod <- function(mod, se.max = 25, ...) {
    
    ##extract SE
    SEs <- extractSE(mod)

    ##extract names
    var.names <- names(SEs)

    ##format as matrix
    matSE <- data.frame(SEs = SEs, variable = var.names)
    
    ##create matrix to hold results for parm with highest SE
    out.result <- data.frame(variable = rep(NA, 1),
                             max.se = rep(NA, 1),
                             n.high.se = rep(NA, 1))

    ##identify maximum value of SE in model
    maxSE <- max(SEs)
    
    ##check if length = 0 (when NaN are present)
    if(is.nan(maxSE)) {
        nan.var <- as.character(matSE[is.nan(matSE$SEs), "variable"])
        if(length(nan.var) == 1) {
            nameSE <- nan.var
        } else {nameSE <- nan.var[1]}
    } else {
        nameSE <- as.character(matSE[which(matSE$SEs == maxSE), "variable"])
    }
    
    rownames(out.result) <- "beta"
    
    out.result[, "variable"] <- nameSE
    out.result[, "max.se"] <- maxSE
    out.result[, "n.high.se"] <- length(which(matSE$SEs > se.max))

    out <- list(model.class = class(mod), se.max = se.max, result = out.result)
    class(out) <- "checkParms"
    return(out)
}



##lmerModLmerTest objects
checkParms.lmerModLmerTest <- function(mod, se.max = 25, ...) {
    
    ##extract SE
    SEs <- extractSE(mod)

    ##extract names
    var.names <- names(SEs)

    ##format as matrix
    matSE <- data.frame(SEs = SEs, variable = var.names)
    
    ##create matrix to hold results for parm with highest SE
    out.result <- data.frame(variable = rep(NA, 1),
                             max.se = rep(NA, 1),
                             n.high.se = rep(NA, 1))

    ##identify maximum value of SE in model
    maxSE <- max(SEs)
    
    ##check if length = 0 (when NaN are present)
    if(is.nan(maxSE)) {
        nan.var <- as.character(matSE[is.nan(matSE$SEs), "variable"])
        if(length(nan.var) == 1) {
            nameSE <- nan.var
        } else {nameSE <- nan.var[1]}
    } else {
        nameSE <- as.character(matSE[which(matSE$SEs == maxSE), "variable"])
    }
    
    rownames(out.result) <- "beta"
    
    out.result[, "variable"] <- nameSE
    out.result[, "max.se"] <- maxSE
    out.result[, "n.high.se"] <- length(which(matSE$SEs > se.max))

    out <- list(model.class = class(mod), se.max = se.max, result = out.result)
    class(out) <- "checkParms"
    return(out)
}



##multinom objects
checkParms.multinom <- function(mod, se.max = 25, ...) {
    
    ##extract SE
    SEs <- sqrt(diag(vcov(mod)))

    ##extract names
    var.names <- names(SEs)

    ##format as matrix
    matSE <- data.frame(SEs = SEs, variable = var.names)
    
    ##create matrix to hold results for parm with highest SE
    out.result <- data.frame(variable = rep(NA, 1),
                             max.se = rep(NA, 1),
                             n.high.se = rep(NA, 1))

    ##identify maximum value of SE in model
    maxSE <- max(SEs)
    
    ##check if length = 0 (when NaN are present)
    if(is.nan(maxSE)) {
        nan.var <- as.character(matSE[is.nan(matSE$SEs), "variable"])
        if(length(nan.var) == 1) {
            nameSE <- nan.var
        } else {nameSE <- nan.var[1]}
    } else {
        nameSE <- as.character(matSE[which(matSE$SEs == maxSE), "variable"])
    }
    
    rownames(out.result) <- "beta"
    
    out.result[, "variable"] <- nameSE
    out.result[, "max.se"] <- maxSE
    out.result[, "n.high.se"] <- length(which(matSE$SEs > se.max))

    out <- list(model.class = class(mod), se.max = se.max, result = out.result)
    class(out) <- "checkParms"
    return(out)
}



##nlme objects
checkParms.nlme <- function(mod, se.max = 25, ...) {
    
    ##extract SE
    SEs <- sqrt(diag(vcov(mod)))

    ##extract names
    var.names <- names(SEs)

    ##format as matrix
    matSE <- data.frame(SEs = SEs, variable = var.names)
    
    ##create matrix to hold results for parm with highest SE
    out.result <- data.frame(variable = rep(NA, 1),
                             max.se = rep(NA, 1),
                             n.high.se = rep(NA, 1))

    ##identify maximum value of SE in model
    maxSE <- max(SEs)
    
    ##check if length = 0 (when NaN are present)
    if(is.nan(maxSE)) {
        nan.var <- as.character(matSE[is.nan(matSE$SEs), "variable"])
        if(length(nan.var) == 1) {
            nameSE <- nan.var
        } else {nameSE <- nan.var[1]}
    } else {
        nameSE <- as.character(matSE[which(matSE$SEs == maxSE), "variable"])
    }
    
    rownames(out.result) <- "beta"
    
    out.result[, "variable"] <- nameSE
    out.result[, "max.se"] <- maxSE
    out.result[, "n.high.se"] <- length(which(matSE$SEs > se.max))

    out <- list(model.class = class(mod), se.max = se.max, result = out.result)
    class(out) <- "checkParms"
    return(out)
}



##nls objects
checkParms.nls <- function(mod, se.max = 25, ...) {
    
    ##extract SE
    SEs <- sqrt(diag(vcov(mod)))

    ##extract names
    var.names <- names(SEs)

    ##format as matrix
    matSE <- data.frame(SEs = SEs, variable = var.names)
    
    ##create matrix to hold results for parm with highest SE
    out.result <- data.frame(variable = rep(NA, 1),
                             max.se = rep(NA, 1),
                             n.high.se = rep(NA, 1))

    ##identify maximum value of SE in model
    maxSE <- max(SEs)
    
    ##check if length = 0 (when NaN are present)
    if(is.nan(maxSE)) {
        nan.var <- as.character(matSE[is.nan(matSE$SEs), "variable"])
        if(length(nan.var) == 1) {
            nameSE <- nan.var
        } else {nameSE <- nan.var[1]}
    } else {
        nameSE <- as.character(matSE[which(matSE$SEs == maxSE), "variable"])
    }
    
    rownames(out.result) <- "beta"
    
    out.result[, "variable"] <- nameSE
    out.result[, "max.se"] <- maxSE
    out.result[, "n.high.se"] <- length(which(matSE$SEs > se.max))

    out <- list(model.class = class(mod), se.max = se.max, result = out.result)
    class(out) <- "checkParms"
    return(out)
}



##polr objects
checkParms.polr <- function(mod, se.max = 25, ...) {
    
    ##extract SE
    SEs <- sqrt(diag(vcov(mod)))

    ##extract names
    var.names <- names(SEs)

    ##format as matrix
    matSE <- data.frame(SEs = SEs, variable = var.names)
    
    ##create matrix to hold results for parm with highest SE
    out.result <- data.frame(variable = rep(NA, 1),
                             max.se = rep(NA, 1),
                             n.high.se = rep(NA, 1))

    ##identify maximum value of SE in model
    maxSE <- max(SEs)
    
    ##check if length = 0 (when NaN are present)
    if(is.nan(maxSE)) {
        nan.var <- as.character(matSE[is.nan(matSE$SEs), "variable"])
        if(length(nan.var) == 1) {
            nameSE <- nan.var
        } else {nameSE <- nan.var[1]}
    } else {
        nameSE <- as.character(matSE[which(matSE$SEs == maxSE), "variable"])
    }
    
    rownames(out.result) <- "beta"
    
    out.result[, "variable"] <- nameSE
    out.result[, "max.se"] <- maxSE
    out.result[, "n.high.se"] <- length(which(matSE$SEs > se.max))

    out <- list(model.class = class(mod), se.max = se.max, result = out.result)
    class(out) <- "checkParms"
    return(out)
}



##rlm objects
checkParms.rlm <- function(mod, se.max = 25, ...) {
    
    ##extract SE
    SEs <- sqrt(diag(vcov(mod)))

    ##extract names
    var.names <- names(SEs)

    ##format as matrix
    matSE <- data.frame(SEs = SEs, variable = var.names)
    
    ##create matrix to hold results for parm with highest SE
    out.result <- data.frame(variable = rep(NA, 1),
                             max.se = rep(NA, 1),
                             n.high.se = rep(NA, 1))

    ##identify maximum value of SE in model
    maxSE <- max(SEs)
    
    ##check if length = 0 (when NaN are present)
    if(is.nan(maxSE)) {
        nan.var <- as.character(matSE[is.nan(matSE$SEs), "variable"])
        if(length(nan.var) == 1) {
            nameSE <- nan.var
        } else {nameSE <- nan.var[1]}
    } else {
        nameSE <- as.character(matSE[which(matSE$SEs == maxSE), "variable"])
    }
    
    rownames(out.result) <- "beta"
    
    out.result[, "variable"] <- nameSE
    out.result[, "max.se"] <- maxSE
    out.result[, "n.high.se"] <- length(which(matSE$SEs > se.max))

    out <- list(model.class = class(mod), se.max = se.max, result = out.result)
    class(out) <- "checkParms"
    return(out)
}



##survreg objects
checkParms.survreg <- function(mod, se.max = 25, ...) {
    
    ##extract SE
    SEs <- sqrt(diag(vcov(mod)))

    ##extract names
    var.names <- names(SEs)

    ##format as matrix
    matSE <- data.frame(SEs = SEs, variable = var.names)
    
    ##create matrix to hold results for parm with highest SE
    out.result <- data.frame(variable = rep(NA, 1),
                             max.se = rep(NA, 1),
                             n.high.se = rep(NA, 1))

    ##identify maximum value of SE in model
    maxSE <- max(SEs)
    
    ##check if length = 0 (when NaN are present)
    if(is.nan(maxSE)) {
        nan.var <- as.character(matSE[is.nan(matSE$SEs), "variable"])
        if(length(nan.var) == 1) {
            nameSE <- nan.var
        } else {nameSE <- nan.var[1]}
    } else {
        nameSE <- as.character(matSE[which(matSE$SEs == maxSE), "variable"])
    }
    
    rownames(out.result) <- "beta"
    
    out.result[, "variable"] <- nameSE
    out.result[, "max.se"] <- maxSE
    out.result[, "n.high.se"] <- length(which(matSE$SEs > se.max))

    out <- list(model.class = class(mod), se.max = se.max, result = out.result)
    class(out) <- "checkParms"
    return(out)
}



##vglm objects
checkParms.vglm <- function(mod, se.max = 25, ...) {
    
    ##extract SE
    SEs <- sqrt(diag(vcov(mod)))

    ##extract names
    var.names <- names(SEs)

    ##format as matrix
    matSE <- data.frame(SEs = SEs, variable = var.names)
    
    ##create matrix to hold results for parm with highest SE
    out.result <- data.frame(variable = rep(NA, 1),
                             max.se = rep(NA, 1),
                             n.high.se = rep(NA, 1))

    ##identify maximum value of SE in model
    maxSE <- max(SEs)
    
    ##check if length = 0 (when NaN are present)
    if(is.nan(maxSE)) {
        nan.var <- as.character(matSE[is.nan(matSE$SEs), "variable"])
        if(length(nan.var) == 1) {
            nameSE <- nan.var
        } else {nameSE <- nan.var[1]}
    } else {
        nameSE <- as.character(matSE[which(matSE$SEs == maxSE), "variable"])
    }
    
    rownames(out.result) <- "beta"
    
    out.result[, "variable"] <- nameSE
    out.result[, "max.se"] <- maxSE
    out.result[, "n.high.se"] <- length(which(matSE$SEs > se.max))

    out <- list(model.class = class(mod), se.max = se.max, result = out.result)
    class(out) <- "checkParms"
    return(out)
}



##zeroinfl objects
checkParms.zeroinfl <- function(mod, se.max = 25, ...) {
    
    ##extract SE
    SEs <- sqrt(diag(vcov(mod)))

    ##extract names
    var.names <- names(SEs)

    ##format as matrix
    matSE <- data.frame(SEs = SEs, variable = var.names)
    
    ##create matrix to hold results for parm with highest SE
    out.result <- data.frame(variable = rep(NA, 1),
                             max.se = rep(NA, 1),
                             n.high.se = rep(NA, 1))

    ##identify maximum value of SE in model
    maxSE <- max(SEs)
    
    ##check if length = 0 (when NaN are present)
    if(is.nan(maxSE)) {
        nan.var <- as.character(matSE[is.nan(matSE$SEs), "variable"])
        if(length(nan.var) == 1) {
            nameSE <- nan.var
        } else {nameSE <- nan.var[1]}
    } else {
        nameSE <- as.character(matSE[which(matSE$SEs == maxSE), "variable"])
    }
    
    rownames(out.result) <- "beta"
    
    out.result[, "variable"] <- nameSE
    out.result[, "max.se"] <- maxSE
    out.result[, "n.high.se"] <- length(which(matSE$SEs > se.max))

    out <- list(model.class = class(mod), se.max = se.max, result = out.result)
    class(out) <- "checkParms"
    return(out)
}



##print method
print.checkParms <- function(x, digits = 2, ...) {
    ##result data frame
    out.frame <- x$result
    
    ##round numeric variables
    out.frame$max.se <- round(out.frame$max.se, digits = digits)
    print(out.frame)
}

##compute predicted values and SE
predictSE.zip <- function(mod, newdata, se.fit = TRUE, parm.type = "lambda", type = "response", c.hat = 1, print.matrix = FALSE) {

  ##only response scale is supported for ZIP models
  if(identical(type, "link")) stop("\nLink scale not supported for predictions of this model type\n")
  
  ##extract data from model object
  if(!is.data.frame(newdata)) stop("\n'newdata' must be a data frame\n")
  new.data.set <- newdata

  ##nobs
  nvals <- nrow(new.data.set)

  ##extract variables on lambda for pcount( ) model
  if(identical(class(mod)[1], "unmarkedFitPCount")) {
    lam.est <- coef(mod@estimates@estimates$state)
  }
  ##extract variables on lambda for pcountOpen( ) model  
  if(identical(class(mod)[1], "unmarkedFitPCO")) {
    lam.est <- coef(mod@estimates@estimates$lambda)
  } 
  lam.est.noint <- lam.est[-1]

  ##total parameters lambda + psi
  ncoefs <- length(lam.est) + 1
  ##number of parameters on lambda
  n.est.lam <- length(lam.est)
  
  ##extract variables on psi
  psi.est <- coef(mod@estimates@estimates$psi)

  ##full model labels
  mod.lab <- labels(coef(mod))

  ##extract labels
  lam.lab <- labels(lam.est)
  lam.lab.noint <- lam.lab[-1]
  psi.lab <- labels(psi.est)

  ##extract formula from model
  formula <- mod@formula

  ##if lambda
  if(identical(parm.type, "lambda")) {
    if(identical(class(mod)[1], "unmarkedFitPCount")) {
      form <- as.formula(paste("~", formula[3], sep="")) #state
    }
    if(identical(class(mod)[1], "unmarkedFitPCO")) {
      form <- mod@formlist$lambdaformula
    }
  } else {
    stop("\nThis function only supports predictions on lamba\n")
  }

  
  ##extract model frame matrix
  Mat <- model.frame(formula = form, data = new.data.set)      
  des.mat <- model.matrix(form, Mat)

##########################################    
##########################################
  ##check for offset
  X.offset <- model.offset(Mat)
  if(is.null(X.offset)) {
    X.offset <- rep(0, nrow(Mat))
  }
    
  ##check for intercept
  if(identical(parm.type, "lambda")) {int.yes <- any(lam.lab == "lam(Int)")}
  
  ##if no intercept term, return error
  if(!int.yes) stop("\nThis function does not work with models excluding the intercept terms: change model parameterization\n")
  
  ##number of estimates (not counting intercept)
  n.est <- n.est.lam - 1

  if(n.est.lam > 1) {
    ##create a list holding each cov
    covs <- list( )
    for(i in 1:n.est) {
      covs[[i]] <- paste("cov", i, sep = "")
    }

    ##covariate labels
    cov.labels <- unlist(covs)
    
    ##change names of columns in design matrix
    colnames(des.mat) <- c("(Int)", unlist(covs))

  } else {colnames(des.mat) <- "(Int)"}

  ##names of columns in design matrix
  design.names <- colnames(des.mat)

  ##extract values from new.data.set
  cov.values <- list( )
  for(i in 1:n.est.lam) {
    cov.values[[i]] <- des.mat[, design.names[i]]
  }
  
  names(cov.values) <- design.names

  ##build equation

  ##iterate over betas except first
  if(n.est.lam > 1) {
    ##betas
    betas <- paste("beta", 0:(n.est), sep = "")
    betas.noint <- betas[-1]
    temp.eq <- list( )
    for(i in 1:length(betas.noint)){
      temp.eq[i] <- paste(betas.noint[i], "*", covs[i], sep = " ")
    }

    ##linear predictor log scale
    lam.eq.log <- paste(c("beta0", unlist(temp.eq)), collapse = " + ")
  } else {
    betas <- "beta0"
    lam.eq.log <- betas
  }

  ##linear predictor log scale
  lam.eq.resp <- paste("exp(", lam.eq.log, "+ Val.offset", ")")
  
  ##logit scale for psi0 (zero-inflation intercept)
  psi.eq <- paste("(1 - (exp(psi0)/(1 + exp(psi0))))")
  
  ##combine both parameters to get abundance
  final.eq <- paste(lam.eq.resp, "*", psi.eq)

  ##total estimates
  tot.est.names <- c(betas, "psi0")
  if(n.est.lam > 1) {
    tot.est <- c(tot.est.names, cov.labels)
  } else {
    tot.est <- c(tot.est.names)
  }
  
  ##extract vcov matrix
  ##multiply by c.hat
  vcmat <- vcov(mod)[c(lam.lab, psi.lab), c(lam.lab, psi.lab)] * c.hat


##################################
  ##start modifications
##################################

  eq.space <- parse(text = as.expression(paste(final.eq, collapse = "+")), 
                    srcfile = NULL)
  no.space <- gsub("[[:space:]]+", "", as.character(eq.space))
  equation <- parse(text = as.expression(no.space))
  if (identical(se.fit, TRUE)) {
    part.devs <- list( )
    for (j in 1:ncoefs) {
      part.devs[[j]] <- D(equation, tot.est.names[j])
    }
  }

  ##assign values of betas and psi
  for(i in 1:n.est.lam) {  
    assign(betas[i], lam.est[i])
  }
  psi0 <- psi.est
  
  cov.values.mat <- matrix(data = unlist(cov.values), nrow = nvals, 
                           ncol = n.est.lam)
  if (identical(se.fit, TRUE)) {
    predicted.SE <- matrix(NA, nrow = nvals, ncol = 2)
    colnames(predicted.SE) <- c("Pred.value", "SE")
    rownames(predicted.SE) <- 1:nvals
    pred.eq <- list()
    ##extract columns
    for (w in 1:nvals) {
      if (int.yes) {
        for (p in 1:n.est.lam) {
          pred.eq[[p]] <- des.mat[w, design.names[p]]
        }
      }
      ##values from design matrix
      design.vals <- unlist(pred.eq)

      ##add value of offset for w
      Val.offset <- X.offset[w]
      
      ##compute values for betas
      exp.beta.pred <- exp(lam.est %*% design.vals + Val.offset)
 
      ##compute predictions including psi
      predicted.vals <- exp.beta.pred * (1 - (exp(psi.est)/(1 + exp(psi.est))))
      
      ##assign values for covariates

      ##columns for covariates only - exclude intercept
      if(n.est.lam > 1) {
        design.covs <- design.vals[-1]
      
        for (p in 1:length(cov.labels)) {
          assign(cov.labels[p], design.covs[p])
        }
      }
      
      ##evaluate partial derivative
      part.devs.solved <- list( )
      for (j in 1:ncoefs) {
        part.devs.solved[[j]] <- eval(part.devs[[j]]) 
      }

      mat_partialdevs <- as.matrix(unlist(part.devs.solved))
      mat_tpartialdevs <- t(mat_partialdevs)
      var_hat <- mat_tpartialdevs %*% vcmat %*% mat_partialdevs
      SE <- sqrt(var_hat)
      predicted.SE[w, 1] <- predicted.vals
      predicted.SE[w, 2] <- SE
    }
    out.fit.SE <- list(fit = predicted.SE[, "Pred.value"], 
                       se.fit = predicted.SE[, "SE"])
  } else {
    predicted.SE <- matrix(NA, nrow = nvals, ncol = 1)
    colnames(predicted.SE) <- c("Pred.value")
    rownames(predicted.SE) <- 1:nvals
    pred.eq <- list( )
    ##extract columns
    for (w in 1:nvals) {
      if (int.yes) {
        for (p in 1:n.est.lam) {
          pred.eq[[p]] <- des.mat[w, design.names[p]]
        }
      }
      ##values from design matrix
      design.vals <- unlist(pred.eq)
      
      ##compute values for betas
      exp.beta.pred <- exp(lam.est %*% design.vals)
      
      ##compute predictions including psi
      predicted.vals <- exp.beta.pred * (1 - (exp(psi.est)/(1 + exp(psi.est))))

      predicted.SE[w, 1] <- predicted.vals
    }
    out.fit.SE <- predicted.SE
    colnames(out.fit.SE) <- "fit"
  }

  if (identical(print.matrix, TRUE)) {
    out.fit.SE <- predicted.SE
    if (identical(se.fit, TRUE)) {
      colnames(out.fit.SE) <- c("fit", "se.fit")
    }
    else {
      colnames(out.fit.SE) <- c("fit")
    }
  }
  return(out.fit.SE)
}

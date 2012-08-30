aictab <-
  function(cand.set, modnames, sort = TRUE, c.hat = 1, second.ord = TRUE, nobs = NULL) {
    results <- NULL
    known <- rep(0, 13) #create an identifier of class type other than lm, glm, multinom, polr, lme, gls, mer, unmarked, nls, or coxph
    ##extract classes
    mod.class <- unlist(lapply(X = cand.set, FUN = class))
    ##check if all are identical
    check.class <- unique(mod.class)

    ##determine if lm or glm
    if(identical(check.class, "lm") || identical(check.class, c("glm", "lm"))) {
      results <- aictab.glm(cand.set = cand.set, modnames = modnames, sort = sort, c.hat = c.hat,
                            second.ord = second.ord, nobs = nobs)
      known[1] <- 1
    }   

    ##determine if multinom
    if(identical(sort(check.class), c("multinom", "nnet"))) {
      results <- aictab.mult(cand.set = cand.set, modnames = modnames, sort = sort, c.hat = c.hat,
                             second.ord = second.ord, nobs = nobs)
      known[2] <- 1
    }   

    ##determine if polr
    if(identical(check.class, "polr")) {
      results <- aictab.polr(cand.set = cand.set, modnames = modnames, sort = sort,
                             second.ord = second.ord, nobs = nobs)
      known[3] <- 1
    }   
      

    ##determine if lme
    if(identical(check.class, "lme"))  {
      results <- aictab.lme(cand.set = cand.set, modnames = modnames, sort = sort,
                            second.ord = second.ord, nobs = nobs)
      known[4] <- 1
    }


    ##determine if gls
    if(identical(check.class, "gls"))  {
      results <- aictab.gls(cand.set = cand.set, modnames = modnames, sort = sort,
                            second.ord = second.ord, nobs = nobs)
      known[5] <- 1
    }

    
    ##determine if mer
    if(identical(check.class, "mer"))  {
      results <- aictab.mer(cand.set = cand.set, modnames = modnames, sort = sort,
                            second.ord = second.ord, nobs = nobs)
      known[6] <- 1
    }


    ##determine if unmarked
    unmarked.class <- c("unmarkedFitOccu", "unmarkedFitColExt", "unmarkedFitOccuRN", "unmarkedFitPCount", "unmarkedFitPCO",
                        "unmarkedFitDS", "unmarkedFitGDS")
    if(any(sapply(unmarked.class, FUN = function(i) identical(i, check.class)))) {
      results <- aictab.unmarked(cand.set = cand.set, modnames = modnames, sort = sort,
                                 c.hat = c.hat, second.ord = second.ord, nobs = nobs)
      known[7] <- 1
    }


    ##determine if nls
    if(identical(check.class, "nls"))  {
      results <- aictab.nls(cand.set = cand.set, modnames = modnames, sort = sort,
                            second.ord = second.ord, nobs = nobs)
      known[8] <- 1
    }


    ##determine if coxph
    if(identical(check.class, "coxph") || identical(check.class, c("coxph.null", "coxph"))) {
      results <- aictab.coxph(cand.set = cand.set, modnames = modnames, sort = sort,
                              second.ord = second.ord, nobs = nobs)
      known[9] <- 1
    }

    
    ##determine if rlm
    if(identical(check.class, c("rlm", "lm")))  {
      results <- aictab.rlm(cand.set = cand.set, modnames = modnames, sort = sort,
                            second.ord = second.ord, nobs = nobs)
      known[10] <- 1
    }

    ##determine if rlm
    if(identical(check.class, c("sclm", "clm")))  {
      results <- aictab.clm(cand.set = cand.set, modnames = modnames, sort = sort,
                            second.ord = second.ord, nobs = nobs)
      known[11] <- 1
    }

    ##determine if clmm
    if(identical(check.class, "clmm"))  {
      results <- aictab.clmm(cand.set = cand.set, modnames = modnames, sort = sort,
                             second.ord = second.ord, nobs = nobs)
      known[12] <- 1
    }

        
    ##warn if models are from a mixture of lm and lme model classes
    if(identical(sort(check.class), c("lm", "lme"))) {
      stop("\nFunction not appropriate for mixture of object classes:", "\n",
               "avoid mixing objects of classes \'lm\' and \'lme\'\n")
      known[13] <- 1
    }


    
    ##warn if class is neither lm, glm, multinom, polr, lme, gls, nls, mer, nor unmarkedFit
    if(sum(known) < 1) {stop("\nFunction not yet defined for this object class\n")}

    return(results)
  }



print.aictab <-
  function(x, digits = 2, LL = TRUE, ...) {
    cat("\nModel selection based on", colnames(x)[3], ":\n")
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
  }


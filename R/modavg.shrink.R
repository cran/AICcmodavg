modavg.shrink <-
  function(cand.set, parm, modnames, c.hat = 1, gamdisp = NULL, conf.level = 0.95, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", parm.type = NULL){

    mod.avg.shrink <- NULL
    known <- rep(0, 12) #create an identifier of class type other than lm, glm, multinom, polr, lme, gls, mer, unmarked, coxph
    ##extract classes
    mod.class <- unlist(lapply(X = cand.set, FUN = class))
    ##check if all are identical
    check.class <- unique(mod.class)

    ##determine if lm or glm  
    if(identical(check.class, "lm") || identical(check.class, c("glm", "lm"))) {
      mod.avg.shrink <- modavg.shrink.glm(cand.set = cand.set, parm = parm, modnames = modnames, c.hat = c.hat, gamdisp = gamdisp,
                                   conf.level = conf.level, second.ord = second.ord, nobs = nobs, uncond.se = uncond.se)
      known[1] <- 1
    }   


    ##determine if multinom
    if(identical(check.class, c("multinom", "nnet"))) {
      mod.avg.shrink <- modavg.shrink.mult(cand.set = cand.set, parm = parm, modnames = modnames, c.hat = c.hat,
                                    conf.level = conf.level, second.ord = second.ord, nobs = nobs, uncond.se = uncond.se)
      known[2] <- 1
    } 

    
    ##determine if polr
    if(identical(check.class, "polr")) {
      mod.avg.shrink <- modavg.shrink.polr(cand.set = cand.set, parm = parm, modnames = modnames, conf.level = conf.level,
                             second.ord = second.ord, nobs = nobs, uncond.se = uncond.se)
      known[3] <- 1
    }   
      

    
    ##determine if lme
    if(identical(check.class, "lme"))  {
      mod.avg.shrink <- modavg.shrink.lme(cand.set = cand.set, parm = parm, modnames = modnames, conf.level = conf.level,
                            second.ord = second.ord, nobs = nobs, uncond.se = uncond.se)
      known[4] <- 1
    }      

    ##determine if gls
    if(identical(check.class, "gls"))  {
      mod.avg.shrink <- modavg.shrink.gls(cand.set = cand.set, parm = parm, modnames = modnames, conf.level = conf.level,
                            second.ord = second.ord, nobs = nobs, uncond.se = uncond.se)
      known[5] <- 1
    }

    
    ##determine if mer
    if(identical(check.class, "mer"))  {
      mod.avg.shrink <- modavg.shrink.mer(cand.set = cand.set, parm = parm, modnames = modnames, conf.level = conf.level,
                            second.ord = second.ord, nobs = nobs, uncond.se = uncond.se)
      known[6] <- 1
    }

    
    ##determine if unmarked
    unmarked.class <- c("unmarkedFitOccu", "unmarkedFitColExt", "unmarkedFitOccuRN", "unmarkedFitPCount", "unmarkedFitPCO",
                        "unmarkedFitDS", "unmarkedFitGDS")
    if(any(sapply(unmarked.class, FUN = function(i) identical(i, check.class)))) {
      mod.avg.shrink <- modavg.shrink.unmarked(cand.set = cand.set, parm = parm, modnames = modnames, c.hat = c.hat, conf.level = conf.level,
                                        second.ord = second.ord, nobs = nobs, uncond.se = uncond.se, parm.type = parm.type)
      known[7] <- 1
    }
    

    ##determine if coxph
    if(identical(check.class, "coxph") || identical(check.class, c("coxph.null", "coxph"))) {
      mod.avg.shrink <- modavg.shrink.coxph(cand.set = cand.set, parm = parm, modnames = modnames, conf.level = conf.level,
                              second.ord = second.ord, nobs = nobs, uncond.se = uncond.se)
      known[8] <- 1
    }



    ##determine if rlm
    if(identical(check.class, c("rlm", "lm")))  {
      mod.avg.shrink <- modavg.shrink.rlm(cand.set = cand.set, parm = parm, modnames = modnames, conf.level = conf.level,
                                          second.ord = second.ord, nobs = nobs, uncond.se = uncond.se)
      known[9] <- 1
    }      

    ##determine if clm
    if(identical(check.class, c("sclm", "clm")))  {
      mod.avg.shrink <- modavg.shrink.clm(cand.set = cand.set, parm = parm, modnames = modnames, conf.level = conf.level,
                                          second.ord = second.ord, nobs = nobs, uncond.se = uncond.se)
      known[10] <- 1
    }

    ##determine if clmm
    if(identical(check.class, c("clmm")))  {
      mod.avg.shrink <- modavg.shrink.clmm(cand.set = cand.set, parm = parm, modnames = modnames, conf.level = conf.level,
                                           second.ord = second.ord, nobs = nobs, uncond.se = uncond.se)
      known[11] <- 1
    }      


    
    ##warn if models are from a mixture of model classes
    if(identical(sort(check.class), c("lm", "lme"))) {
      stop("\nFunction not appropriate for mixture of object classes:", "\n",
               "avoid mixing objects of classes lm and lme", "\n")
      known[12] <- 1
    }


    ##warn if class is neither lm, glm, multinom, polr, lme, gls nor mer
    if(sum(known) < 1) {stop("\nFunction not yet defined for this object class\n")}


      
    ##determine if object from nlm optimizer
    ##if(class(mod)[1]=="nlm") {mod.avg <- modavg.nlm(cand.set, parm, modnames, c.hat)}      
    return(mod.avg.shrink)


  }


print.modavg.shrink <-
  function(x, digits = 2, ...) {
    ic <- colnames(x$Mod.avg.table)[3]
    cat("\nMultimodel inference on \"", x$Parameter, "\" based on", ic, "\n")
    cat("\n", ic, "table used to obtain model-averaged estimate with shrinkage:\n")
    oldtab <- x$Mod.avg.table
    if (any(names(oldtab) == "c_hat")) {cat("\t(c-hat estimate = ", oldtab$c_hat[1], ")\n")}
    cat("\n")
    if (any(names(oldtab)=="c_hat")) {
      nice.tab <- cbind(oldtab[, 2], oldtab[, 3], oldtab[, 4], oldtab[, 6],
                        oldtab[, 9], oldtab[, 10])
    } else {nice.tab <- cbind(oldtab[, 2], oldtab[, 3], oldtab[, 4], oldtab[, 6],
                              oldtab[, 8], oldtab[, 9])
          }

##modify printing style if multinomial model is used  
    if(length(x$Mod.avg.beta) == 1) {
      colnames(nice.tab) <- c(colnames(oldtab)[c(2, 3, 4, 6)], "Estimate", "SE")
      rownames(nice.tab) <- oldtab[, 1]
      print(round(nice.tab, digits = digits))
      cat("\nModel-averaged estimate with shrinkage:", eval(round(x$Mod.avg.beta, digits = digits)), "\n")
      cat("Unconditional SE:", eval(round(x$Uncond.SE, digits = digits)), "\n")
      cat("",x$Conf.level*100, "% Unconditional confidence interval:", round(x$Lower.CL, digits = digits),
          ",", round(x$Upper.CL, digits = digits), "\n\n")
    } else {
      col.ns <- ncol(nice.tab)
      nice.tab <- nice.tab[,-c(col.ns - 1, col.ns)]
      colnames(nice.tab) <- c(colnames(oldtab)[c(2, 3, 4, 6)])
      rownames(nice.tab) <- oldtab[, 1]
      print(round(nice.tab, digits = digits))
      cat("\n\nModel-averaged estimates with shrinkage for different levels of response variable:", "\n\n")
      resp.labels <- labels(x$Mod.avg.beta)
      mult.out <- matrix(NA, nrow = length(resp.labels), ncol = 4)
      colnames(mult.out) <- c("Model-averaged estimate with shrinkage", "Uncond. SE", paste(x$Conf.level*100,"% lower CL"),
                              paste(x$Conf.level*100, "% upper CL"))
      rownames(mult.out) <- resp.labels
      mult.out[, 1] <- round(x$Mod.avg.beta, digits = digits)
      mult.out[, 2] <- round(x$Uncond.SE, digits = digits)
      mult.out[, 3] <- round(x$Lower.CL, digits = digits)
      mult.out[, 4] <- round(x$Upper.CL, digits = digits)
      print(mult.out)
      cat("\n")
    }
  }

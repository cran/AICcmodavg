modavg.effect <-
  function(cand.set, modnames, newdata, type = "response", c.hat = 1, gamdisp = NULL, 
           conf.level = 0.95, second.ord = TRUE, nobs = NULL, uncond.se = "revised", parm.type = NULL){

    mod.avg.eff <- NULL
    known <- rep(0, 7) #create an identifier of class type other than lm, glm, lme, gls, mer, unmarked
    ##extract classes
    mod.class <- unlist(lapply(X = cand.set, FUN = class))
    ##check if all are identical
    check.class <- unique(mod.class)

    ##determine if lm or glm  
    if(identical(check.class, "lm") || identical(check.class, c("glm", "lm"))) {
      mod.avg.eff <- modavg.effect.glm(cand.set = cand.set, modnames = modnames, newdata = newdata,
                                       type = type, c.hat = c.hat, gamdisp = gamdisp, conf.level = conf.level,
                                       second.ord = second.ord, nobs = nobs, uncond.se = uncond.se)
      known[1] <- 1
    }   


    ##determine if lme
    if(identical(check.class, "lme"))  {
      mod.avg.eff <- modavg.effect.lme(cand.set = cand.set, modnames = modnames, newdata = newdata,
                                       conf.level = conf.level, second.ord = second.ord, nobs = nobs,
                                       uncond.se = uncond.se)
      known[2] <- 1
    }      

    ##determine if gls
    if(identical(check.class, "gls"))  {
      mod.avg.eff <- modavg.effect.gls(cand.set = cand.set, modnames = modnames, newdata = newdata,
                                conf.level = conf.level, second.ord = second.ord,
                                nobs = nobs, uncond.se = uncond.se)
      known[3] <- 1
    }

    
    ##determine if mer
    if(identical(check.class, "mer"))  {
      mod.avg.eff <- modavg.effect.mer(cand.set = cand.set, modnames = modnames, newdata = newdata,
                                type = type, c.hat = c.hat, conf.level = conf.level,
                                second.ord = second.ord, nobs = nobs, uncond.se = uncond.se)
      known[4] <- 1
    }


    ##determine if unmarked
    unmarked.class <- c("unmarkedFitOccu", "unmarkedFitColExt", "unmarkedFitOccuRN", "unmarkedFitPCount", "unmarkedFitPCO",
                        "unmarkedFitDS", "unmarkedFitGDS")
    if(any(sapply(unmarked.class, FUN = function(i) identical(i, check.class)))) {
      mod.avg.eff <- modavg.effect.unmarked(cand.set = cand.set, modnames = modnames, newdata = newdata,
                                     type = type, c.hat = c.hat, conf.level = conf.level,
                                     second.ord = second.ord, nobs = nobs, uncond.se = uncond.se, parm.type = parm.type)
      known[5] <- 1
    }


    ##determine if rlm
    if(identical(check.class, c("rlm", "lm")))  {
      mod.avg.eff <- modavg.effect.rlm(cand.set = cand.set, modnames = modnames, newdata = newdata,
                                       conf.level = conf.level, second.ord = second.ord, nobs = nobs,
                                       uncond.se = uncond.se)
      known[6] <- 1
    }      


     
    ##warn if models are from a mixture of model classes
    if(identical(sort(check.class), c("lm", "lme"))) {
      stop("\nFunction not appropriate for mixture of object classes:\navoid mixing objects of classes lm and lme\n")
      known[7] <- 1
    }


    ##warn if class is neither lm, glm, multinom, polr, lme, gls nor mer
    if(sum(known) < 1) {stop("\nFunction not yet defined for this object class\n")}


    ##determine if object from nlm optimizer
    ##if(class(mod)[1]=="nlm") {mod.avg <- modavg.nlm(cand.set, parm, modnames, c.hat)}      
    return(mod.avg.eff)


  }


print.modavg.effect <- function(x, digits = 2, ...) {

  ##rework Group.variable labels
  old.type <- x$Group.variable
  stripped.type <- unlist(strsplit(old.type, split = "\\("))

  ic <- colnames(x$Mod.avg.table)[3]
  
  cat("\nModel-averaged effect size on the", x$Type, "scale based on entire model set:\n\n")
  
  ##extract elements
  if(length(stripped.type) == 1) {
    cat("\nMultimodel inference on \"", paste(x$Group.variable, x$Group1, sep = ""), "-",
        paste(x$Group.variable, x$Group2, sep = ""), "\" based on", ic, "\n")
    
    ##if unmarkedFit model, then print differently
  } else {
    ##extract parameter name
    parm.type <- gsub("(^ +)|( +$)", "", stripped.type[1])
    
    ##extract Group.variable name
    var.id <- gsub("(^ +)|( +$)", "", unlist(strsplit(stripped.type[2], "\\)"))[1])

    cat("\nMultimodel inference on \"", paste(parm.type, "(", var.id, x$Group1, ")", sep = ""), "-",
        paste(parm.type, "(", var.id, x$Group2, ")", sep = ""), "\" based on", ic, "\n")
  }

  
  cat("\n", ic, "table used to obtain model-averaged effect size:\n")
  oldtab <- x$Mod.avg.table
  if (any(names(oldtab)=="c_hat")) {cat("\t(c-hat estimate = ", oldtab$c_hat[1], ")\n")}
  cat("\n")
  if (any(names(oldtab)=="c_hat")) {
    nice.tab <- cbind(oldtab[,2], oldtab[,3], oldtab[,4], oldtab[,6],
                      oldtab[,9], oldtab[,10])
  } else {nice.tab <- cbind(oldtab[,2], oldtab[,3], oldtab[,4], oldtab[,6],
                            oldtab[,8], oldtab[,9])
        }

  colnames(nice.tab) <- c(colnames(oldtab)[c(2,3,4,6)], paste("Effect (", x$Group1, " - ", x$Group2, ")", sep = ""), "SE")
  rownames(nice.tab) <- oldtab[,1]
  print(round(nice.tab, digits=digits))
  cat("\nModel-averaged effect size:", eval(round(x$Mod.avg.eff, digits=digits)), "\n")
  cat("Unconditional SE:", eval(round(x$Uncond.se, digits=digits)), "\n")
  cat("",x$Conf.level * 100, "% Unconditional confidence interval:", round(x$Lower.CL, digits=digits),
      ",", round(x$Upper.CL, digits=digits), "\n\n")
}


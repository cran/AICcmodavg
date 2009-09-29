modavg <-
  function(cand.set, parm, modnames, c.hat=1, gamdisp=NULL, conf.level=0.95, second.ord=TRUE, nobs=NULL){

    mod.avg <- NULL
    known <- rep(0, 4) #create an identifier of class type other than lm, glm, or lme
#extract classes
    mod.class <- unlist(lapply(X=cand.set, FUN=class))
#check if all are identical
    check.class <- unique(mod.class)

#determine if lm or glm  
    if(identical(check.class, "lm") || identical(check.class, c("glm", "lm"))) {
      mod.avg <- modavg.glm(cand.set=cand.set, parm=parm, modnames=modnames, c.hat=c.hat, gamdisp=gamdisp,
                            conf.level=conf.level, second.ord=second.ord, nobs=nobs)
      known[1] <- 1
    }   


#determine if multinom
    if(identical(check.class, c("multinom", "nnet"))) {
      mod.avg <- modavg.mult(cand.set=cand.set, parm=parm, modnames=modnames, c.hat=c.hat,
                             second.ord=second.ord, nobs=nobs)
      known[2] <- 1
    } 


#determine if lme
    if(identical(check.class, "lme"))  {
      mod.avg <- modavg.lme(cand.set=cand.set, parm=parm, modnames=modnames,
                            conf.level=conf.level, second.ord=second.ord, nobs=nobs)
      known[3] <- 1
    }      


#warn if models are from a mixture of model classes
    if(identical(sort(check.class), c("lm", "lme"))) {
      stop(cat("Function not appropriate for mixture of object classes:", "\n",
               "avoid mixing objects of classes lm and lme", "\n"))
      known[4] <- 1
    }


#warn if class is neither lm, glm, or lme
    if(sum(known) < 1) {stop("Function not defined for this object class")}


#if(class(mod)[1]=="nlm") {mod.avg <- modavg.nlm(cand.set, parm, modnames, c.hat)}      #determine if object from nlm optimizer
    return(mod.avg)


  }


print.modavg <-
  function(x, digits=4, ...) {
    ic <- colnames(x$Mod.avg.table)[3]
    cat("\nMultimodel inference on \"", x$Parameter, "\" based on", ic, "\n")
    cat("\n", ic, "table used to obtain model-averaged estimate:\n")
    oldtab <- x$Mod.avg.table
    if (ncol(oldtab) > 9) {cat("\t(c-hat estimate = ", oldtab$c_hat[1], ")\n")}
    cat("\n")
    if (ncol(oldtab) > 9) {
      nice.tab <- cbind(oldtab[,2], oldtab[,3], oldtab[,4], oldtab[,6],
                        oldtab[,7], oldtab[,9], oldtab[,10])
    } else {nice.tab <- cbind(oldtab[,2], oldtab[,3], oldtab[,4], oldtab[,6],
                              oldtab[,7], oldtab[,8], oldtab[,9])
          }

#modify printing style if multinomial model is used  
    if(length(x$Mod.avg.beta)==1) {
      colnames(nice.tab) <- c(colnames(oldtab)[c(2,3,4,6,7)], "Estimate", "SE")
      rownames(nice.tab) <- oldtab[,1]
      print(round(nice.tab, digits=digits))
      cat("\nModel-averaged estimate:", eval(round(x$Mod.avg.beta, digits=digits)), "\n")
      cat("Unconditional SE:", eval(round(x$Uncond.SE, digits=digits)), "\n")
      cat("",x$Conf.level*100, "% Unconditional confidence interval:", round(x$Lower.CL, digits=digits),
          ",", round(x$Upper.CL, digits=digits), "\n\n")
    } else {
      nice.tab <- nice.tab[,-c(6,7)]
      colnames(nice.tab) <- c(colnames(oldtab)[c(2,3,4,6,7)])
      rownames(nice.tab) <- oldtab[,1]
      print(round(nice.tab, digits=digits))
      cat("\n\nModel-averaged estimates for different levels of response variable:", "\n\n")
      resp.labels <- labels(x$Mod.avg.beta)
      mult.out <- matrix(NA, nrow=length(resp.labels), ncol=4)
      colnames(mult.out) <- c("Model-averaged estimate", "Uncond. SE", paste(x$Conf.level*100,"% lower CL"),
                              paste(x$Conf.level*100, "% upper CL"))
      rownames(mult.out) <- resp.labels
      mult.out[,1] <- round(x$Mod.avg.beta, digits=digits)
      mult.out[,2] <- round(x$Uncond.SE, digits=digits)
      mult.out[,3] <- round(x$Lower.CL, digits=digits)
      mult.out[,4] <- round(x$Upper.CL, digits=digits)
      print(mult.out)
      cat("\n")
    }
  }

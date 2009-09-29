aictab <-
function(cand.set, modnames, sort=TRUE, c.hat=1, second.ord=TRUE, nobs=NULL) {
  results <- NULL
  known <- rep(0, 4) #create an identifier of class type other than lm, glm, multinom, or lme
#extract classes
  mod.class <- unlist(lapply(X=cand.set, FUN=class))
#check if all are identical
  check.class <- unique(mod.class)

#determine if lm or glm
  if(identical(check.class, "lm") || identical(check.class, c("glm", "lm"))) {
    results <- aictab.glm(cand.set=cand.set, modnames=modnames, sort=sort, c.hat=c.hat, second.ord=second.ord, nobs=nobs)
    known[1] <- 1
  }   


#determine if multinom
  if(identical(sort(check.class), c("multinom", "nnet"))) {
    results <- aictab.mult(cand.set=cand.set, modnames=modnames, sort=sort, c.hat=c.hat, second.ord=second.ord, nobs=nobs)
    known[2] <- 1
  }   
  

#determine if lme
  if(identical(check.class, "lme"))  {
    results <- aictab.lme(cand.set=cand.set, modnames=modnames, sort=sort, second.ord=second.ord, nobs=nobs)
    known[3] <- 1
  }

#warn if models are from a mixture of model classes
  if(identical(sort(check.class), c("lm", "lme"))) {
    stop(cat("Function not appropriate for mixture of object classes:", "\n",
             "avoid mixing objects of classes lm and lme", "\n"))
    known[4] <- 1
  }

#warn if class is neither lm, glm, multinom, nor lme
  if(sum(known) < 1) {stop("Function not defined for this object class")}

  return(results)
}

print.aictab <-
  function(x, digits=4, ...) {
    cat("\nModel selection based on", colnames(x)[3], ":\n")
    if (ncol(x) > 7) {cat("(c-hat estimate = ", x$c_hat[1], ")\n")}
    cat("\n")
    nice.tab <- cbind(x[,2], x[,3], x[,4], x[,6], x[,7])
    colnames(nice.tab) <- colnames(x)[c(2,3,4,6,7)]
    rownames(nice.tab) <- x[,1]
    print(round(nice.tab, digits=digits)) #select rounding off with digits argument
    cat("\n")
  }

